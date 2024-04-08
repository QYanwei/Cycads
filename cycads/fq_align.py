#!/usr/bin/env python3
# coding: utf-8
import os
import argparse

def fq_align_action(args):
    input_fastq_path = args['fastq']
    reference = args['reference']
    minimap2 = args['minimap2']
    samtools = args['samtools']
    alignment_threads = args['alignment_threads']
    sort_threads = args['sort_threads']
    minimap2_arguments = args['minimap2_arguments']
    output_bam_path = args['output_bam_path']
    command = f'{minimap2} {minimap2_arguments} -t {alignment_threads} {reference} {input_fastq_path} | {samtools} sort -@ {sort_threads} -o {output_bam_path}'
    print(command)
    os.system(command)
    command = f'{samtools} index {output_bam_path}'
    print(command)
    os.system(command)
    if os.path.isfile(output_bam_path):
        args['bam'] = output_bam_path
    else:
        raise RuntimeError(f"Failed to align {input_fastq_path} to {reference}")
