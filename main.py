#!/usr/bin/env python3
# coding: utf-8

# Copyright (C) 2024, CycloneSEQ.
# qiyanwei1@genomics.cn, or qiyanweii@icloud.com.

# Cycads is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# Cycads is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
@name: Cycads long reads quality analyser
@description: All-in-one parsing fasta format file
@author: Yanwei Qi & Jiayuan Zhang
@email: qiyanwei1@genomics.cn, zhangjiayuan@genomics.cn
@last modified by: qiyanwei

change log:
    2024/03/25  Coding completed.
    2024/02/18  Sketch finished.
    2024/01/18  Project created.
"""

import os,re,sys,time
import argparse

from cycad import fq_index
from cycad import fq_datum
from cycad import fq_figure
from cycad import fq_filter
from cycad import fq_align
from cycad import bam_datum
from cycad import bam_figure
from cycad import all_report
# basic stent function
## check file
def file_checkin(infile, func):
    if os.path.exists(infile):
        print(func + ': check...' + infile + ' file existed.')
        return 1
    else:
        print(func + ': check...' + infile + ' file not found.')
        return 0
## version control
def version():
    version_file = open('%sVERSION.txt' % './')
    return version_file.readline().strip()
## help document
def print_helpdoc():
    help_message = """
              ........:::=== Cycads V%s ===:::........
    ============================================================================
                Quality control & Data filtering & Error analysis
                             for Long-read sequencing
    ============================================================================
    """ % version()
    print(help_message)

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    ############################################## initialize subparsers ###############################################

    # initialize the options parser
    parser = argparse.ArgumentParser(description='Cycads:Long reads sequencing quality analyzer')
    parser.add_argument("-fq",   "--fastq",     required=False,  help="sequences.fq/fq.gz")
    parser.add_argument("-bam",  "--alignment", required=False, help="alignment.bam")
    parser.add_argument("-ref",  "--reference", required=False, help="reference.fasta")
    
    # for quality control.
    parser.add_argument("-P", "--platform", required=False, help="cyclone")
    parser.add_argument("-M", "--mode", type=str, default="overall", required=False, help="if you want fast, please set to sampling")
    parser.add_argument("-Hshift", "--head_shift_length", type=int, default=200, required=False, help="check head bases quality")
    parser.add_argument("-Tshift", "--tail_shift_length", type=int, default=200, required=False, help="check tail bases quality")
    parser.add_argument("-kmer", "--kmer_size_frequency",     type=int, default=5, required=False, help="observe kmer size specturm")
    parser.add_argument("-hpmin", "--homopolymer_min_length", type=int, default=2, required=False, help="observe minium homopolymer")
    parser.add_argument("-hpmax", "--homopolymer_max_length", type=int, default=9, required=False, help="observe maxium homopolymer")

    # for data filtering.
    parser.add_argument("-filter", "--filtering", required=False, help="get clean data")
    parser.add_argument("-Qmin", "--minium_quality",  type=float, default='10',   required=False, help="filter low quality reads")
    parser.add_argument("-Lmin", "--minium_length",   type=int,   default='1000', required=False, help="filter short reads")
    parser.add_argument("-Lmax", "--maxium_length",   type=int,   default='1000', required=False, help="filter large reads")
    parser.add_argument("-Hcut", "--cut_head_length", type=int,   default='200',  required=False, help="trim head sequences")
    parser.add_argument("-Tcut", "--cut_tail_length", type=int,   default='200', required=False,  help="trim head sequences")
    parser.add_argument("-Htrim", "--trim_head_homopolymer", type=int, default='200', required=False, help="trim head sequences homopolymer")
    parser.add_argument("-Ttrim", "--trim_tail_homopolymer", type=int, default='200', required=False, help="trim head sequences homopolymer")
    parser.add_argument("-Dx", "--downsampling_dx", type=str, default='10', required=False, help="depth")
    parser.add_argument("-Dg", "--downsampling_gs", type=str, default='4m', required=False, help="genome size")
    # third party tool
    # for error analysis.
    parser.add_argument("-t", "--thread", required=False, default=4, help="thread number")
    parser.add_argument("-hpmindelmax", "--homopolymer_indel_max", type=int, default=4, required=False, help="homopolymer max indel shift")
    # for storing result.
    parser.add_argument("-name", "--sample_name", default='cycads_report', required=False, help="prefix of output file name")
    ############################## parse provided arguments and run corresponding function #############################

    # get and check options
    #    args = None
    args = parser.parse_args()
    if (len(sys.argv) == 1) or (sys.argv[1] == '-h') or (sys.argv[1] == '-help') or (sys.argv[1] == '--help'):
        print_helpdoc()
        sys.exit(0)
    else:
        args = vars(parser.parse_args())
    if os.path.exists("./" + args["sample_name"]):
        print("output folder existed!")
    else:
        os.mkdir("./" + args["sample_name"])
    if os.path.exists("./" + args["sample_name"] + "/report_html"):
        print("report folder existed")
    else:
        os.mkdir("./" + args["sample_name"] + "/report_html")
    pwd_config_file = os.path.realpath(__file__)
    args["pyfastx"] = '/'.join(pwd_config_file.split('/')[:-1]) + '/tool/pyfastx'
    args["minimap2"] = '/'.join(pwd_config_file.split('/')[:-1]) + '/tool/minimap2'
    args["samtools"] = '/'.join(pwd_config_file.split('/')[:-1]) + '/tool/samtools'
    if not os.path.exists(args["pyfastx"]):
        print("pyfastx: not found in the ./tool/")
    else:
        pass
    if not os.path.exists(args["minimap2"]):
        print("minimap2: not found in the ./tool/")
    else:
        pass
    if not os.path.exists(args["samtools"]):
        print("samtools: not found in the ./tool/")
    else:
        pass
    if args["fastq"] and not args["filtering"] and not args["alignment"] and not args["reference"] :
        if os.path.exists(args["fastq"]):
            fq_index.fq_index_action(args)           
            fq_datum.fq_datum_action(args)
            fq_figure.fq_figure_action(args)
            all_report.generate_html(args)
        else:
            print(args["fastq"] + " does not exist!")
    elif args["fastq"] and args["filtering"] and not args["alignment"] and not args["reference"]:
        if os.path.exists(args["fastq"]):
            fq_index.fq_index_action(args)
            fq_datum.fq_datum_action(args)
            fq_figure.fq_figure_action(args)
            fq_filter.fq_filter_action(args)
            all_report.generate_html(args)
        else:
            print(args["fastq"] + " does not exist!")
    elif args["fastq"] and args["reference"] and not args["filtering"] and not args["alignment"]:
        if os.path.exists(args["fastq"]) and os.path.exists(args["reference"]):
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
    elif args["alignment"] and not args["fastq"] and not args["reference"] and not args["filtering"]:
        if os.path.exists(args["alignment"]):
            bam_datum.bam_datum_action(args)
            bam_figure.bam_figure_action(args)
            all_report.generate_html(args)
        else:
            print(args["alignment"] + " does not exist!")
    else:
        print("please input correct file: fq/fastq/fq.gz & fastq+referecence.fasta & alignment.bam")

