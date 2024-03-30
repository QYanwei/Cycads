#!/usr/bin/env python3
# coding: utf-8
import re,os,sys,time
import numpy as np
import gzip
import pyfastx


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

def Downsampling(fq, args, out_file):
    read_number = 0
    base_number = 0
    buffer_iter = 0
    query = ''
    fq = pyfastx.Fastq(fq, build_index=False)
    for name, raw_seq, raw_qual in fq:
        read_number += 1
        read_length = len(raw_seq)
        if base_number <= args['target_depth'] * args['genome_size']:
            if read_length >= args['min_length']:
                if (args['trim_5_end'] or args['trim_3_end']):
                    if (args['trim_5_end'] + args['trim_3_end']) < read_length:
                        buffer_iter += 1
                        base_number += read_length
                        new_read, new_qual = ends_cutting(raw_seq, raw_qual, args['trim_5_end'], args['trim_3_end'])
                        query += '%s\n%s\n%s\n%s\n' % ('@'+name, new_read, '+', new_qual)
                        if buffer_iter == buffersize:
                            reads_buffer_io(query, out_file)
                            buffer_iter = 0
                            query = ''
                        else:
                            continue
                    else:
                        pass
                else:
                    buffer_iter += 1
                    base_number += read_length
                    query += '%s\n%s\n%s\n%s\n' % ('@'+name, raw_seq, '+', raw_qual)
                    if buffer_iter == buffersize:
                        reads_buffer_io(query, out_file)
                        buffer_iter = 0
                        query = ''
                    else:
                        continue
            else:
                pass
        else:
            outfile = gzip.open(out_file, 'wt+')
            outfile.write(query)
            outfile.close()
            break
def Filtering(fq, args, out_file):
    read_number = 0
    base_number = 0
    buffer_iter = 0
    query = ''
    fq = pyfastx.Fastq(fq, build_index=False)
    for name, raw_seq, raw_qual in fq:
        read_number += 1
        read_length = len(raw_seq)
        if read_length >= args['min_length'] and read_length <= args['max_length']:
            if (args['trim_5_end'] or args['trim_3_end']):
                if (args['trim_5_end'] + args['trim_3_end']) < read_length:
                    buffer_iter += 1
                    base_number += read_length
                    new_read, new_qual = ends_cutting(raw_seq, raw_qual, args['trim_5_end'], args['trim_3_end'])
                    query += '%s\n%s\n%s\n%s\n' % ('@' + name, new_read, '+', new_qual)
                    if buffer_iter == buffersize:
                        reads_buffer_io(query, out_file)
                        buffer_iter = 0
                        query = ''
                else:
                    continue
            else:
                buffer_iter += 1
                base_number += read_length
                query += '%s\n%s\n%s\n%s\n' % ('@'+name, raw_seq, '+', raw_qual)
                if buffer_iter == buffersize:
                    reads_buffer_io(query, out_file)
                    buffer_iter = 0
                    query = ''
def Extracting(fastq, args, outfile):
    read_number = 0
    base_number = 0
    buffer_iter = 0
    query = ''
    fq = pyfastx.Fastq(fastq)
    read_count = len(fq)
    seed = args['seed']
    sample_size = args['extract']
    selected_read_indices = np.random.default_rng(seed).randint(low=0, high=read_count, size=sample_size, dtype=np.uint64)
    for i in selected_read_indices:
        read = fq[i]
        name, raw_seq, raw_qual = read.name, read.seq, read.qual
        read_number += 1
        read_length = len(raw_seq)
        if read_length >= args['min_length'] and read_length <= args['max_length']:
            if (args['head_cut_length'] or args['tail_cut_length']):
                if (args['head_cut_length'] + args['tail_cut_length']) < read_length:
                    buffer_iter += 1
                    base_number += read_length
                    new_read, new_qual = ends_cutting(raw_seq, raw_qual, args['head_cut_length'], args['tail_cut_length'])
                    query += '%s\n%s\n%s\n%s\n' % ('@' + name, new_read, '+', new_qual)
                    if buffer_iter == buffersize:
                        reads_buffer_io(query, outfile)
                        buffer_iter = 0
                        query = ''
                else:
                    continue
            else:
                buffer_iter += 1
                base_number += read_length
                query += '%s\n%s\n%s\n%s\n' % ('@'+name, raw_seq, '+', raw_qual)
                if buffer_iter == buffersize:
                    reads_buffer_io(query, outfile)
                    buffer_iter = 0
                    query = ''


def fq_filter_action(args):
    input_fastq_path = args['fastq']
    with open(args['filtered_fastq_path'], 'wt') as f:
        if args['target_depth'] and args['genome_size']:
            Downsampling(input_fastq_path, args, f)
        elif args['extract']:
            Extracting(input_fastq_path, args, f)
        elif args['min_length'] or args['max_length']:
            Filtering(input_fastq_path, args, f)
        else:
            raise ValueError("Invalid arguments")


# if __name__ == "__main__":
#     args = {
#         'head_cut_length': 10,
#         'tail_cut_length': 0,
#         'random_read_size': 100,
#         'min_length': 1,
#         'max_length': 10000000,
#         'max_depth': 0,
#         'genome_size': 10000000,
#         'seed': 0,
#         'random_read': 0,
#     }
#     buffersize = 100
#     fq = '../test/ecoli.fq.gz'
#     out = '../test/ecoli.filter.gz'
#     outfile = gzip.open(out, 'wt')
#     if args['max_depth'] and args['genome_size']:
#         Downsampling(fq, args, outfile)
#         outfile.close()
#         sys.exit()
#     elif args['sample'] > 0:
#         Extracting(fq, args, outfile)
#         outfile.close()
#         sys.exit()
#     elif args['min_length'] or args['max_length']:
#         Filtering(fq, args, outfile)
#         outfile.close()
#         sys.exit()
#     else:
#         print('stop')
    