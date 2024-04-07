#!/usr/bin/env python3
# coding: utf-8
import re,os,sys,time
import numpy as np
import gzip
import pyfastx
import subprocess
from collections import Counter

def ends_cutting(raw_read, raw_qual, head_cut_length, tail_cut_length):
    new_read, new_qual = '', ''
    if head_cut_length > 0 and tail_cut_length == 0:
        new_read = raw_read[head_cut_length:]
        new_qual = raw_qual[head_cut_length:]
    elif head_cut_length == 0 and tail_cut_length > 0:
        new_read = raw_read[:-tail_cut_length]
        new_qual = raw_qual[:-tail_cut_length]
    elif head_cut_length > 0 and tail_cut_length > 0:
        new_read = raw_read[head_cut_length:-tail_cut_length]
        new_qual = raw_qual[head_cut_length:-tail_cut_length]
    else:
        pass
    return new_read, new_qual
def reads_buffer_io(queries, outfile):
    outfile.write(queries)

def readAvgQscore(quali, seq_len):
    value_list = Counter(list(quali))
    value_sum = sum([k * v for k, v in value_list.items()])
    seq_qscore = value_sum / seq_len
    return seq_qscore

def Downsampling(fq, args, out_file):
    read_number = 0
    base_number = 0
    fq = pyfastx.Fastq(fq, build_index=False)
    for name, raw_seq, raw_qual in fq:
        read_number += 1
        read_length = len(raw_seq)
        quali = [ord(i) - 33 for i in raw_qual]
        read_q = readAvgQscore(quali, read_length)
        if base_number <= (float(args['target_depth']) * float(args['genome_size'])):
            if read_length >= args['min_length'] and read_length <= args['max_length'] and read_q >= args['min_base_quality']:
                if args['trim_5_end'] or args['trim_3_end']:
                    if (args['trim_5_end'] + args['trim_3_end']) < read_length:
                        base_number += read_length - args['trim_5_end'] - args['trim_3_end']
                        new_read, new_qual = ends_cutting(raw_seq, raw_qual, args['trim_5_end'], args['trim_3_end'])
                        query = '%s\n%s\n%s\n%s\n' % ('@'+name, new_read, '+', new_qual)
                        reads_buffer_io(query, out_file)
                    else:
                        pass
                else:
                    base_number += read_length
                    query = '%s\n%s\n%s\n%s\n' % ('@'+name, raw_seq, '+', raw_qual)
                    reads_buffer_io(query, out_file)
            else:
                pass

def Filtering(fq, args, out_file):
    read_number = 0
    base_number = 0
    fq = pyfastx.Fastq(fq, build_index=False)
    for name, raw_seq, raw_qual in fq:
        read_number += 1
        read_length = len(raw_seq)
        quali = [ord(i) - 33 for i in raw_qual]
        read_q = readAvgQscore(quali, read_length)
        if read_length >= args['min_length'] and read_length <= args['max_length'] and read_q >= args['min_base_quality']:
            if (args['trim_5_end'] or args['trim_3_end']):
                if (args['trim_5_end'] + args['trim_3_end']) < read_length:
                    base_number += read_length
                    new_read, new_qual = ends_cutting(raw_seq, raw_qual, args['trim_5_end'], args['trim_3_end'])
                    query = '%s\n%s\n%s\n%s\n' % ('@' + name, new_read, '+', new_qual)
                    reads_buffer_io(query, out_file)
                else:
                    continue
            else:
                base_number += read_length
                query = '%s\n%s\n%s\n%s\n' % ('@'+name, raw_seq, '+', raw_qual)
                reads_buffer_io(query, out_file)

def Extracting(fastq, args, outfile):
    read_number = 0
    base_number = 0
    fq = pyfastx.Fastq(fastq)
    read_count = len(fq)
    seed = args['seed']
    sample_size = args['extract']
    selected_read_indices = np.random.default_rng(seed).integers(low=0, high=read_count, size=sample_size, dtype=np.uint64)
    for i in selected_read_indices:
        read = fq[i]
        name, raw_seq, raw_qual, raw_quali = read.name, read.seq, read.qual, read.quali
        read_number += 1
        read_length = len(raw_seq)
        read_q = readAvgQscore(quali, read_length)
        if read_length >= args['min_length'] and read_length <= args['max_length'] and read_q >= args['min_base_quality']:
            if (args['trim_5_end'] or args['trim_3_end']):
                if (args['trim_5_end'] + args['trim_3_end']) < read_length:
                    base_number += read_length
                    new_read, new_qual = ends_cutting(raw_seq, raw_qual, args['trim_5_end'], args['trim_3_end'])
                    query = '%s\n%s\n%s\n%s\n' % ('@' + name, new_read, '+', new_qual)
                    reads_buffer_io(query, outfile)
                else:
                    continue
            else:
                base_number += read_length
                query = '%s\n%s\n%s\n%s\n' % ('@'+name, raw_seq, '+', raw_qual)
                reads_buffer_io(query, outfile)

def filtered_fq_stat(args):
    input_fastq_path = os.path.abspath(args['filtered_fastq_path'])
    pyfastx_path = args['pyfastx']
    command = [pyfastx_path, 'index', input_fastq_path]
    subprocess.run(command)
    print("------ Sequencing data summary ------")
    os.system(f'{pyfastx_path} stat {input_fastq_path} |awk \'BEGIN{{OFS="\t"}}{{NF=7}}1\'')

def fq_filter_action(args):
    input_fastq_path = args['fastq']
    with open(args['filtered_fastq_path'], 'wt', buffering=1000) as f:
        if args['target_depth'] and args['genome_size']:
            Downsampling(input_fastq_path, args, f)
            filtered_fq_stat(args)
            f.close()
            sys.exit()
        elif args['extract']:
            Extracting(input_fastq_path, args, f)
            f.cloes()
            filtered_fq_stat(args)
            sys.exit()
        elif args['min_length'] or args['max_length']:
            Filtering(input_fastq_path, args, f)
            f.close()
            filtered_fq_stat(args)
            sys.exit()
        else:
            raise ValueError("Invalid arguments")
