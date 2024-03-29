#!/usr/bin/env python3
# coding: utf-8
import os
import subprocess
import argparse

def fq_index_action(args):
    output = args['sample_name']
    infastq = args['fastq']
    subprocess.run(['pyfastx', 'index', infastq])
    print("## Data overview")
    os.system('pyfastx stat ' +  infastq  + ' |awk \'BEGIN{OFS="\t"}{NF=7}1\' ' ) # TODO: use subprocess.run()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-fq",  "--fastq",   required=True, help="sequences.fq/fq.gz")
    parser.add_argument("-pfx", "--pyfastx", required=True, help="pyfastx index fastq")
    parser.add_argument("-name", "--sample_name", default='cycads_report', required=False, help="prefix of output file name")
    args = vars(parser.parse_args())
    fq_index_action(args)

