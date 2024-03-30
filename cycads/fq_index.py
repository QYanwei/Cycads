#!/usr/bin/env python3
# coding: utf-8
import os
import subprocess
import argparse
from .helpers import which


def fq_index_action(args):
    sample_name = args['sample_name']
    infastq = args['fastq']
    pyfastx_path = args['pyfastx']
    subprocess.run([pyfastx_path, 'index', infastq])
    summary_path = args['fastq_summary_path']
    os.system(f'{pyfastx_path} stat {infastq} |awk \'BEGIN{{OFS="\t"}}{{NF=7}}1\' > {summary_path}')
    print("##sequencing summary overview")
    with open(summary_path, 'rt') as f:
        print(f.read())

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-fq",  "--fastq",   required=True, help="sequences.fq/fq.gz")
    parser.add_argument("-pfx", "--pyfastx", required=False, help="pyfastx index fastq")
    parser.add_argument("-o", "--output_dir", default='./', required=False, help="Output direcotry")
    parser.add_argument("-n", "--sample_name", default='cycads_report', required=False, help="prefix of output file name")
    args = vars(parser.parse_args())
    fq_index_action(args)

