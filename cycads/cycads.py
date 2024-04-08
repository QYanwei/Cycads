#!/usr/bin/env python3
# coding: utf-8


"""
@name: Cycads long reads quality analyser
@description: All-in-one parsing fasta format file
@author: Yanwei Qi & Jiayuan Zhang
@email: qiyanwei1@genomics.cn, zhangjiayuan@genomics.cn

change log:
    2024/03/25  Coding completed.
    2024/02/18  Sketch finished.
    2024/01/18  Project created.
"""

import os,re,sys,time
import argparse
from warnings import warn

from . import fq_index, fq_datum, fq_figure, fq_filter, fq_align, bam_datum, bam_figure, all_report, helpers
from . import __version__, __description__

def check_binary_dependencies(args, dependencies):
    cycads_dir = os.path.abspath(os.path.dirname(__file__))
    for executable in dependencies:
        cli_path = args.get(executable + "_path")
        symlink_path = os.path.join(cycads_dir, 'tool', executable)  
        system_path = helpers.which(executable)   
        if cli_path and os.path.isfile(cli_path):
            args[executable] = cli_path
        elif os.path.isfile(symlink_path):
            args[executable] = symlink_path
        elif system_path:
            args[executable] = system_path
        else:
            raise RuntimeError(f"Unable to find {executable}")
        print(f"Using {executable} from {args[executable]}")


def parse_command_line_arguments():
    parser = argparse.ArgumentParser(description=__description__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    io_group = parser.add_argument_group('I/O', 'Input/output arguments.')
    io_group.add_argument("-f",   "--fastq", metavar="FASTQ_PATH", required=False, default=None, help="Input FASTQ file. Supported extensions include *.fastq and *.fastq.gz.")
    io_group.add_argument("-b",  "--bam", metavar="BAM_PATH", required=False, default=None, help="Input BAM file.")
    io_group.add_argument("-r",  "--reference", metavar="REFERENCE_PATH", required=False, default=None, help="Reference FASTA file.")
    io_group.add_argument("-o", "--output_dir", default='cycads_output', required=False, help="Output direcotry.")
    io_group.add_argument("-n", "--sample_name", default='sample', required=False, help="Sample name displayed in output reports.")
    io_group.add_argument('-p', "--platform", default='cyclone', required=False, help="Design for CycloneSEQ data, also adopt to ONT and PB data.")
    
    fastq_group = parser.add_argument_group('FASTQ', 'Arguments for FASTQ analyses. Only effective when FASTQ_PATH is supplied.')
    fastq_group.add_argument("-s", "--sample", metavar="N", type=int, default=10000, required=False, help="Only include a random sample of N reads from the input FASTQ file to accelerate evaluation.")
    fastq_group.add_argument("--seed", metavar="SEED", type=int, default=1, required=False, help="Random seed for sampling.")
    fastq_group.add_argument("-T", "--check_terminal_bases", metavar="N", type=int, default=200, required=False, help="Analyze N bases at both ends of each read.")

    filter_group = parser.add_argument_group('Filtering', 'Arguments for filtering the input FASTQ file. Only effective when FASTQ_PATH is supplied.')
    filter_group.add_argument("-F", "--filter", action='store_true', required=False, help="Output filtered FASTQ file. Analyses are always based on the input FASTQ file.")
    filter_group.add_argument("-e", "--extract", metavar="N", required=False, default=None, help="Randomly extract N reads from the input FASTQ file.")
    filter_group.add_argument("-Q", "--min_base_quality", metavar="MIN_BASE_QUALITY", type=float, default=7,   required=False, help="Remove reads with mean base quality less than MIN_BASE_QUALITY.")
    filter_group.add_argument("--min_length",  metavar="MIN_READ_LENGTH", type=int,   default=1, required=False, help="Remove reads shorter than MIN_READ_LENGTH.")
    filter_group.add_argument("--max_length",  metavar="MAX_READ_LENGTH", type=int,   default=1000000000, required=False, help="Remove reads longer than MAX_READ_LENGTH.")
    filter_group.add_argument("--trim_5_end", metavar="N", type=int,   default=0,  required=False, help="Trim N bases from the 5' end of each read.")
    filter_group.add_argument("--trim_3_end", metavar="N", type=int,   default=0,  required=False, help="Trim N bases from the 3' end of each read.")
    filter_group.add_argument("-d", "--target_depth", metavar="TARGET_DEPTH", type=float, default=None, required=False, help="Downsample FASTQ file to TARGET_DEPTH. Requires GENOME_SIZE to be supplied.")
    filter_group.add_argument("-g", "--genome_size", metavar="GENOME_SIZE", type=str, default=None, required=False, help="Genome size of sequenced sample. Required if TARGET_DEPTH is set.")

    hp_group = parser.add_argument_group('Homopolymers', 'Arguments related to homopolymer analyses. ')
    hp_group.add_argument("--min_homopolymer_size", metavar='MIN_HOMOPOLYMER_SIZE', type=int, default=2, required=False, help="Do not analyze homopolymers shorter than MIN_HOMOPOLYMER_SIZE.")
    hp_group.add_argument("--max_homopolymer_size", metavar='MAX_HOMOPOLYMER_SIZE', type=int, default=9, required=False, help="Do not analyze homopolymers longer than MAX_HOMOPOLYMER_SIZE.")
    hp_group.add_argument("--max_homopolymer_indel_size", metavar='MAX_HOMOPOLYMER_INDEL_SIZE', type=int, default=4, required=False, help="Analyze homopolymer expansion/contraction up to MAX_HOMOPOLYMER_INDEL_SIZE.")

    alignment_group = parser.add_argument_group('Alignment', 'Arguments for read alignment. Only effective when FASTQ_PATH and REFERENCE_PATH are supplied.')
    alignment_group.add_argument("--alignment_threads", metavar="THREADS", required=False, default=4, help="Number of threads used in read alignment.")
    alignment_group.add_argument("--sort_threads", metavar="THREADS", required=False, default=1, help="Number of threads used in sorting aligned segments.")
    alignment_group.add_argument("--minimap2_arguments", metavar="ARGUMENTS", required=False, type=str, default="-ax map-ont --secondary=no --MD --eqx -I 10G", help="Alignment arguments to be passed to minimap2.")

    denpendency_group = parser.add_argument_group('Dependencies', 'Arguments for custom paths to external binary dependencies. Cycads searches for binary dependencies in the following order: 1. arguments specified here; 2. the `dependencies` folder in Cycads installation path; 3. the system $PATH environmental variable.')
    denpendency_group.add_argument("--minimap2", required=False, default=None, help="Path to Minimap2.")
    denpendency_group.add_argument("--samtools", required=False, default=None, help="Path to samtools.")
    denpendency_group.add_argument("--pyfastx", required=False, default=None, help="Path to pyfastx.")

    args = vars(parser.parse_args())
    return args




banner = f"""                      === Cycads {__version__} ===
============================================================================
            Quality control & Data filtering & Error analysis
                            for Long-read sequencing
============================================================================
"""


def main():
    print(banner)
    args = parse_command_line_arguments()

    for argument in ("fastq", "bam", "reference"):
        path = args[argument]
        if path is not None and not os.path.isfile(path):
            raise IOError(f"Failed to find input file: {path}")

    output_dir = args["output_dir"]
    output_folder =  os.path.join(args["output_dir"], args["sample_name"])
    
    output_dir = args["output_dir"]
    sample_output_dir = os.path.join(args['output_dir'], args["sample_name"])
    args['sample_output_dir'] = sample_output_dir
    args['fastq_pickle_path'] = os.path.join(sample_output_dir, "fq.pickle")
    args['bam_pickle_path'] = os.path.join(sample_output_dir, "bam.pickle")
    args['fastq_summary_path'] = os.path.join(sample_output_dir, "fastq_summary.txt")
    args['filtered_fastq_path'] = os.path.join(sample_output_dir, "filtered_reads.fastq")
    args['output_bam_path'] = os.path.join(sample_output_dir, "aligned_reads.bam")

    report_dir = os.path.join(sample_output_dir, "HTML_report")
    args['report_dir'] = report_dir

    if os.path.isdir(output_dir):
        print(f"Output folder {output_dir} already exists.")
    else:
        os.makedirs(output_dir, exist_ok=True)
    try:
        os.makedirs(report_dir, exist_ok=True)    
    except Exception:
        raise IOError(f"Failed to create output folder: {report_dir}")

    check_binary_dependencies(args, ('pyfastx', 'minimap2', 'samtools'))

    if not args["fastq"] and not args["bam"]:
        raise ValueError("Invalid arguments. Please provide an input FASTQ or BAM file!")

    if args["fastq"] and not args["filter"] and not args["bam"] and not args["reference"] :
        if os.path.exists(args["fastq"]):
            fq_index.fq_index_action(args)           
            fq_datum.fq_datum_action(args)
            fq_figure.fq_figure_action(args)
            all_report.generate_html(args)
        else:
            print(args["fastq"] + " does not exist!")
    elif args["fastq"] and args["filter"] and not args["bam"] and not args["reference"]:
        if os.path.exists(args["fastq"]):
            fq_index.fq_index_action(args)
#            fq_datum.fq_datum_action(args)
#            fq_figure.fq_figure_action(args)
            fq_filter.fq_filter_action(args)
#            all_report.generate_html(args)
        else:
            print(args["fastq"] + " does not exist!")
    elif args["fastq"] and args["reference"] and not args["filter"] and not args["bam"]:
        if os.path.exists(args["fastq"]) and os.path.exists(args["reference"]):
            print(args['platform'])
            print(args['minimap2'])
            if args['platform'] in ['Cyclone', 'cyclone', 'CycloneSEQ'] or args['platform'] in ['ont', 'ONT']:
                pass
            elif args['platform'] in ['pb', 'PB', 'PB-CLR'] and 'minimap' in args['minimap2']:
                args['minimap2_arguments'] = '-ax map-pb --secondary=no --MD --eqx -I 10G'
            elif args['platform'] in ['hifi', 'pb-hifi', 'PB-HiFi'] and 'minimap2' in args['minimap2']:
                args['minimap2_arguments'] = '-ax map-hifi --secondary=no --MD --eqx -I 10G'
            else:
                platform = args['platform']
                raise IOError(f'Wrong sequencing platform input:{platform}, please check the platform parameter!')
            fq_index.fq_index_action(args)
            fq_datum.fq_datum_action(args)
            fq_figure.fq_figure_action(args)
            fq_align.fq_align_action(args)
            bam_datum.bam_datum_action(args)
            bam_figure.bam_figure_action(args)
            all_report.generate_html(args)
        elif not os.path.exists(args["fastq"]) and os.path.exists(args["reference"]):
            print(args["fastq"] + " does not exist!")  
        elif os.path.exists(args["fastq"]) and not os.path.exists(args["reference"]):
            print(args["reference"] + " does not exist!")
        else:
            print( "Both " + args["fastq"] + " " +args["reference"] + " are not exist!")
    elif args["bam"] and not args["fastq"] and not args["reference"] and not args["filter"]:
        bam_datum.bam_datum_action(args)
        bam_figure.bam_figure_action(args)
        all_report.generate_html(args)
        
    else:
        raise ValueError("Invalid arguments.")

if __name__ == '__main__':
    main()

