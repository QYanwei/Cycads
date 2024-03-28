#!/usr/bin/env python3
# coding: utf-8
import os

input_fastq = '../test/ecoli.fq.gz'
pyfastx_tool = "../tool/pyfastx"

def pyfastx_stat_fq_and_idx(pyfastx_tool, input_fastq):
    os.system(pyfastx_tool + ' stat ' +  input_fastq  + ' > '+ '../test/test.pyfqx.txt')

pyfastx_stat_fq_and_idx(pyfastx_tool, input_fastq)