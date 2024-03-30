#!/usr/bin/env python3
# coding: utf-8

import os
import argparse

def fq_align_action(args):
    infastq = args['fastq']
    reference = args['reference']
    output = args['sample_name']
    minimap2 = args['minimap2']
    samtools = args['samtools']
    thread = args['thread']
    alignment_threads = args['alignment_threads']
    sort_threads = args['sort_threads']
    minimap2_arguments = args['minimap2_arguments']
    output_bam_path = args['output_bam_path']
    os.system(f'{minimap2} {minimap2_arguments} -t {alignment_threads} {reference} {infastq} | {samtools} sort -@ {sort_threads} -o {output_bam_path}')
    os.system(f'{samtools} index {output_bam_path}')

# if __name__ == '__main__':
#     parser = argparse.ArgumentParser()
#     parser.add_argument("-fq",  "--fastq",   required=True, help="sequences.fq/fq.gz")
#     parser.add_argument("-ref", "--reference", required=True, help="reference.fasta")
#     parser.add_argument("-mmp", "--minimap2", required=True, help="minimap2 tool")
#     parser.add_argument("-sam", "--samtools", required=True, help="samtools tool")
#     parser.add_argument("-t", "--thread", required=False, default=4, help="thread number")
#     parser.add_argument("-o", "--output_dir", default='./', required=False, help="Output direcotry")
#     parser.add_argument("-n", "--sample_name", default='cycads_report', required=False, help="prefix of output file name")
#     args = vars(parser.parse_args())
#     fq_align_action(args)

