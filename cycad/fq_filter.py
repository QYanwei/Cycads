#!/usr/bin/env python3
# coding: utf-8
import re,os,sys,time
import numpy as np
import gzip
import pyfastx

def random_readnum(seed_num, read_size, sample_num):
    np.random.seed(seed_num)
    sample_list = np.random.randint(1, read_size, (1, sample_num))
    return sample_list
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
        if base_number <= args['max_depth'] * args['genome_size']:
            if read_length >= args['min_length']:
                if (args['head_cut_length'] or args['tail_cut_length']):
                    if (args['head_cut_length'] + args['tail_cut_length']) < read_length:
                        buffer_iter += 1
                        base_number += read_length
                        new_read, new_qual = ends_cutting(raw_seq, raw_qual, args['head_cut_length'], args['tail_cut_length'])
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
            if (args['head_cut_length'] or args['tail_cut_length']):
                if (args['head_cut_length'] + args['tail_cut_length']) < read_length:
                    buffer_iter += 1
                    base_number += read_length
                    new_read, new_qual = ends_cutting(raw_seq, raw_qual, args['head_cut_length'], args['tail_cut_length'])
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
    sample_list = random_readnum(args['seed'], len(fq), args['random_read'])[0]
    for i in sample_list:
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

if __name__ == "__main__":
    args = {
        'head_cut_length': 10,
        'tail_cut_length': 0,
        'random_read_size': 100,
        'min_length': 1,
        'max_length': 10000000,
        'max_depth': 0,
        'genome_size': 10000000,
        'seed': 0,
        'random_read': 0,
    }
    buffersize = 100
    fq = '../test/ecoli.fq.gz'
    out = '../test/ecoli.filter.gz'
    outfile = gzip.open(out, 'wt')
    if args['max_depth'] and args['genome_size']:
        Downsampling(fq, args, outfile)
        outfile.close()
        sys.exit()
    elif args['random_read']:
        Extracting(fq, args, outfile)
        outfile.close()
        sys.exit()
    elif args['min_length'] or args['max_length']:
        Filtering(fq, args, outfile)
        outfile.close()
        sys.exit()
    else:
        print('stop')
    