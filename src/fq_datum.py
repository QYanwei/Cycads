#!/usr/bin/env python3
# coding: utf-8
import re,os,sys,time
import pyfastx

import numpy as np
from threading import Thread
from collections import Counter

def readGCcontent(seq):
    C = seq.count("C")
    G = seq.count("G")
    GC = round( (C+G)/len(seq), 3)
    return GC

def readAvgQscore(quali):
    read_length = len(quali)
    value_list = Counter( list(quali) )
    value_sum = sum([ k * v for k, v in value_list.items() ])
    read_qscore = value_sum / read_length
    return read_qscore

def endBaseHeadParse(seq, shift_length, endBaseQual_dict):
    i = 0
    for b in seq[:shift_length]:
        endBaseQual_dict['HeadBaseContent_dict'][b][i] += 1
        i += 1
    return endBaseQual_dict
def endBaseTailParse(seq, shift_length, endBaseQual_dict):
    i = 0
    for b in seq[-shift_length:]:
        endBaseQual_dict['TailBaseContent_dict'][b][i] += 1
        i += 1
    return endBaseQual_dict
def endQualHeadParse(seq, quali, shift_length, endBaseQual_dict):
    i = 0
    for b in seq[:shift_length]:
        endBaseQual_dict['HeadQualContent_dict'][b][i] += quali[:shift_length][i]
        i += 1
    return endBaseQual_dict
def endQualTailParse(seq, quali, shift_length, endBaseQual_dict):
    i = 0
    for b in seq[-shift_length:]:
        endBaseQual_dict['TailQualContent_dict'][b][i] += quali[-shift_length:][i]
        i += 1
    return endBaseQual_dict
def read_quality_to_bin_score(base_qual_list: list, split_number_of_read_length: int):
    div, mod = divmod(len(base_qual_list), split_number_of_read_length)
    qs_in_percentage_pos = [sum(base_qual_list[i * div + min(i, mod):(i+1) * div + min(i+1, mod)])/len(base_qual_list[i * div + min(i, mod):(i+1) * div + min(i+1, mod)]) for i in range(split_number_of_read_length)]
    return qs_in_percentage_pos
def allBaseQualParse(seq, quali, split_part_num, allBaseQual_dict):
    if read_length >= split_part_num:
        i = 0
        for b in seq[:shift_length]:
            allBaseQual_dict['HeadQualContent_dict'][b][i] += quali[:shift_length][i]
            i += 1
    return endBaseQual_dict

def endBaseQualParse(seq, quali, shift_length, endBaseQual_dict):
    if len(seq) > shift_length * 2:
        threads = []
        threads.append(Thread(target=endBaseHeadParse, args=(seq, shift_length, endBaseQual_dict)))
        threads.append(Thread(target=endBaseTailParse, args=(seq, shift_length, endBaseQual_dict)))
        threads.append(Thread(target=endQualHeadParse, args=(seq, quali, shift_length, endBaseQual_dict)))
        threads.append(Thread(target=endQualTailParse, args=(seq, quali, shift_length, endBaseQual_dict)))
        for fun in threads:
            fun.start()
    return endBaseQual_dict
def kmerSpectrumParse(fq_path, kmer_size, output_dir):
    kmersize = kmer_size
    output = output_dir
    pwd_config_file = os.path.realpath(__file__)
    meryl = '/'.join(pwd_config_file.split('/')[:-1]) + '/tools/meryl'
    if fq_path.endswith('gz'):
        os.system('gunzip -c {}|awk \'NR %4 == 1 || NR %4 == 2 \'   > {}.fasta '.format(fq_path, output))
        os.system('{} count  k={} {}.fasta output {}.meryl'.format(meryl, str(kmersize), output, output))
        os.system('{} print {}.meryl |sort -k2nr > {}_kmer_{}_freq.txt'.format(meryl, output, output, str(kmersize)))
        os.system('rm -f {}.fasta'.format(output))
        os.system('rm -rf {}.meryl'.format(output))
    elif fq_path.endswith('fastq') or self.args.input.endswith('fq'):
        os.system('awk \'NR %4 == 1 || NR %4 == 2 \'   > {}.fasta '.format(fq_path, output))
        os.system('{} count  k={} {}.fasta output {}.meryl'.format(meryl, str(kmersize), output, output))
        os.system('{} print {}.meryl | sort -k2nr > {}_kmer_{}_freq.txt'.format(meryl, output, output, str(kmersize)))
        os.system('rm -f {}.fasta'.format(output))
        os.system('rm -rf {}.meryl'.format(output))
    else:
        print('please determining the input file suffix is fastq or fq or fq.gz!')
# kmerSpectrumParse('../test/ecoli.fq.gz', 5, '../test/')
def homopolymerParse(seq, homopolymer_size_min, homopolymer_dict):
    ATCG = ['A', 'G', 'C', 'T']
    for base in ATCG:
        patterns = re.compile(r"%s{%d,}" % (base, homopolymer_size_min))
        homostring_dict = dict(Counter(patterns.findall(seq)))
        if len(homostring_dict ):
            for homostring, homofreq in homostring_dict.items():
                homolen = len(homostring)
                if homolen in  homopolymer_dict[base].keys():
                    homopolymer_dict[base][homolen] += homofreq
                else:
                    homopolymer_dict[base][homolen] = homofreq
    return homopolymer_dict

def readParse(read, seqdict):
    seqdict['ID'].append(read.name)
    seqdict['GC'].append(readGCcontent(read.seq))
    seqdict['LEN'].append(len(read.seq))
    seqdict['QUAL1'].append(float(read.name.split('_')[-1]))
    seqdict['QUAL2'].append(readAvgQscore(read.quali))
    return seqdict
def random_readnum(seed_num, read_size, sample_num):
    np.random.seed(seed_num)
    sample_list = np.random.randint(1, read_size, (1, sample_num) )
    return sample_list
def sampling_analyser(fq, seed_num, read_num):
    seqdict = dict( {'ID': [], 'GC': [], 'LEN': [], 'QUAL1': [], 'QUAL2':[]} ) # QUAL1: read basecall Q, QUAL2: read average Q
    homopolymer_size_min = 5
    homopolymer_dict = {'A':{}, 'G':{}, 'C':{}, 'T':{}}
    shift_length = 200
    endBaseQual_dict = {
                    'HeadBaseContent_dict': {'A': [0]*shift_length, 'G': [0]*shift_length, 'C': [0]*shift_length, 'T': [0]*shift_length},
                    'TailBaseContent_dict': {'A': [0]*shift_length, 'G': [0]*shift_length, 'C': [0]*shift_length, 'T': [0]*shift_length},
                    'HeadQualContent_dict': {'A': [0]*shift_length, 'G': [0]*shift_length, 'C': [0]*shift_length, 'T': [0]*shift_length},
                    'TailQualContent_dict': {'A': [0]*shift_length, 'G': [0]*shift_length, 'C': [0]*shift_length, 'T': [0]*shift_length},
                    }
    sample_list = random_readnum(seed_num, len(fq), read_num)[0]
    for i in sample_list:
        read = fq[i]
        seqdict = readParse(read, seqdict)
        homopolymer_dict = homopolymerParse(read.seq, homopolymer_size_min, homopolymer_dict)
        endBaseQual_dict = endBaseQualParse(read.seq, read.quali, shift_length, endBaseQual_dict)
    return seqdict, homopolymer_dict, endBaseQual_dict

def overall_analyser(fq):
    seqdict = dict( {'ID': [], 'GC': [], 'LEN': [], 'QUAL1': [], 'QUAL2': []} )  # QUAL1: read basecall Q, QUAL2: read average Q
    homopolymer_size_min = 5
    homopolymer_dict = {'A':{}, 'G':{}, 'C':{}, 'T':{}}
    shift_length = 200
    endBaseQual_dict = {
                    'HeadBaseContent_dict': {'A': [0]*shift_length, 'G': [0]*shift_length, 'C': [0]*shift_length, 'T': [0]*shift_length},
                    'TailBaseContent_dict': {'A': [0]*shift_length, 'G': [0]*shift_length, 'C': [0]*shift_length, 'T': [0]*shift_length},
                    'HeadQualContent_dict': {'A': [0]*shift_length, 'G': [0]*shift_length, 'C': [0]*shift_length, 'T': [0]*shift_length},
                    'TailQualContent_dict': {'A': [0]*shift_length, 'G': [0]*shift_length, 'C': [0]*shift_length, 'T': [0]*shift_length},
                    }
    for read in fq:
        seqdict = readParse(read, seqdict)
        homopolymer_dict = homopolymerParse(read.seq, homopolymer_size_min, homopolymer_dict)
        endBaseQual_dict = endBaseQualParse(read.seq, read.quali, shift_length, endBaseQual_dict)
    return seqdict, homopolymer_dict, endBaseQual_dict

fq = pyfastx.Fastq('../test/ecoli.fq.gz')
seqdict1, homopolymer_dict1, endBaseQual_dict1  = sampling_analyser(fq, 1, 100)
seqdict2, homopolymer_dict2, endBaseQual_dict2 = overall_analyser(fq)

for k, v in endBaseQual_dict1.items():
    for i, j in v.items():
        print(k, i, j)

for k, v in endBaseQual_dict2.items():
    for i, j in v.items():
        print(k, i, j)
