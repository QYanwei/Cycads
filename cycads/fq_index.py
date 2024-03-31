#!/usr/bin/env python3
# coding: utf-8
import os
import subprocess
import argparse


def fq_index_action(args):
    input_fastq_path = os.path.abspath(args['fastq'])
    output_dir = args['sample_output_dir']
    symlink_fastq_path = os.path.abspath(os.path.join(output_dir, "input.symlink.fastq"))
    if input_fastq_path.endswith(".gz"):
        symlink_fastq_path += '.gz'
    os.symlink(input_fastq_path, symlink_fastq_path)
    args['fastq'] = symlink_fastq_path

    pyfastx_path = args['pyfastx']
    print("Creating FASTQ index")
    command = [pyfastx_path, 'index', symlink_fastq_path]
    print(" ".join(command))
    subprocess.run(command)
    summary_path = args['fastq_summary_path']
    os.system(f'{pyfastx_path} stat {input_fastq_path} |awk \'BEGIN{{OFS="\t"}}{{NF=7}}1\' > {summary_path}')
    print("------ Sequencing data summary ------")
    with open(summary_path, 'rt') as f:
        print(f.read())

# if __name__ == '__main__':
#     parser = argparse.ArgumentParser()
#     parser.add_argument("-fq",  "--fastq",   required=True, help="sequences.fq/fq.gz")
#     parser.add_argument("-pfx", "--pyfastx", required=False, help="pyfastx index fastq")
#     parser.add_argument("-o", "--output_dir", default='./', required=False, help="Output direcotry")
#     parser.add_argument("-n", "--sample_name", default='cycads_report', required=False, help="prefix of output file name")
#     args = vars(parser.parse_args())
#     fq_index_action(args)

