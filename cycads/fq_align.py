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
    os.system(f'{minimap2} --secondary=no --MD --eqx -t {thread} -ax map-ont {reference} {infastq} |{samtools} sort -@ {thread} |{samtools} view -@ {thread} -bS > {output}/{output}.bam')
    os.system(f'{samtools} index {output}/{output}.bam')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-fq",  "--fastq",   required=True, help="sequences.fq/fq.gz")
    parser.add_argument("-ref", "--reference", required=True, help="reference.fasta")
    parser.add_argument("-mmp", "--minimap2", required=True, help="minimap2 tool")
    parser.add_argument("-sam", "--samtools", required=True, help="samtools tool")
    parser.add_argument("-t", "--thread", required=False, default=4, help="thread number")
    parser.add_argument("-O", "--output_dir", default='cycads_report', required=False, help="Output direcotry")
    parser.add_argument("-name", "--sample_name", default='cycads_report', required=False, help="prefix of output file name")
    args = vars(parser.parse_args())
    fq_align_action(args)

