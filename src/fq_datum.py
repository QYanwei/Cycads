#!/usr/bin/env python3
# coding: utf-8
import re,os,sys,time
import json
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
    read_qscore = round(value_sum / read_length, 3)
    return read_qscore

def endBaseHeadParse(seq, shift_length, endBaseQual_dict):
    i = 0
    for b in seq[:shift_length]:
        endBaseQual_dict['HeadBaseContent_dict'][b][i] += 1
        i += 1
    endBaseQual_dict['HeadBaseContent_dict']['S'] += 1
    return endBaseQual_dict
def endBaseTailParse(seq, shift_length, endBaseQual_dict):
    i = 0
    for b in seq[-shift_length:]:
        endBaseQual_dict['TailBaseContent_dict'][b][i] += 1
        i += 1
    endBaseQual_dict['TailBaseContent_dict']['S'] += 1
    return endBaseQual_dict
def endQualHeadParse(seq, quali, shift_length, endBaseQual_dict):
    i = 0
    for b in seq[:shift_length]:
        endBaseQual_dict['HeadQualContent_dict'][b][i] += quali[:shift_length][i]
        i += 1
    endBaseQual_dict['HeadQualContent_dict']['S'] += 1
    return endBaseQual_dict
def endQualTailParse(seq, quali, shift_length, endBaseQual_dict):
    i = 0
    for b in seq[-shift_length:]:
        endBaseQual_dict['TailQualContent_dict'][b][i] += quali[-shift_length:][i]
        i += 1
    endBaseQual_dict['TailQualContent_dict']['S'] += 1
    return endBaseQual_dict
def read_quality_to_bin_score(base_qual_list: list, split_number_of_read_length: int):
    div, mod = divmod(len(base_qual_list), split_number_of_read_length)
    qs_in_percentage_pos = [sum(base_qual_list[i * div + min(i, mod):(i+1) * div + min(i+1, mod)])/len(base_qual_list[i * div + min(i, mod):(i+1) * div + min(i+1, mod)]) for i in range(split_number_of_read_length)]
    return qs_in_percentage_pos
def allBaseQualParse(quali, split_part_num, allBaseQual_dict):
    if len(quali) >= split_part_num:
        qv_in_percentage_pos = read_quality_to_bin_score(quali, split_part_num)
        for j in range(split_part_num):
            allBaseQual_dict['PercentBaseQual_dict']['Q'][j] += qv_in_percentage_pos[j]
        allBaseQual_dict['PercentBaseQual_dict']['S'] += 1
    return allBaseQual_dict

def endBaseQualParse(seq, quali, shift_length, endBaseQual_dict):
    if len(seq) > shift_length * 2 and '@' not in seq:
        endBaseQual_dict = endBaseHeadParse(seq, shift_length, endBaseQual_dict)
        endBaseQual_dict = endBaseTailParse(seq, shift_length, endBaseQual_dict)
        endBaseQual_dict = endQualHeadParse(seq, quali, shift_length, endBaseQual_dict)
        endBaseQual_dict = endQualTailParse(seq, quali, shift_length, endBaseQual_dict)
    return endBaseQual_dict
def homopolymerParse(seq, homopolymer_size_min, homopolymer_dict):
    ATCG = ['A', 'G', 'C', 'T']
    for base in ATCG:
        patterns = re.compile(r"%s{%d,}" % (base, homopolymer_size_min))
        homostring_dict = dict(Counter(patterns.findall(seq)))
        if len(homostring_dict):
            for homostring, homofreq in homostring_dict.items():
                homolen = len(homostring)
                if homolen in  homopolymer_dict[base].keys():
                    homopolymer_dict[base][homolen] += homofreq
                else:
                    homopolymer_dict[base][homolen] = homofreq
    return homopolymer_dict
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


def readParse(read, seqdict):
    seqdict['ID'].append(read.id)
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
                    'HeadBaseContent_dict': {'A': [0]*shift_length, 'G': [0]*shift_length, 'C': [0]*shift_length, 'T': [0]*shift_length, 'S':0},
                    'TailBaseContent_dict': {'A': [0]*shift_length, 'G': [0]*shift_length, 'C': [0]*shift_length, 'T': [0]*shift_length, 'S':0},
                    'HeadQualContent_dict': {'A': [0]*shift_length, 'G': [0]*shift_length, 'C': [0]*shift_length, 'T': [0]*shift_length, 'S':0},
                    'TailQualContent_dict': {'A': [0]*shift_length, 'G': [0]*shift_length, 'C': [0]*shift_length, 'T': [0]*shift_length, 'S':0},
                    }
    split_part_num = 100
    allBaseQual_dict = {
                    'PercentBaseQual_dict': {'Q':[0]* split_part_num, 'S':0} # Q: average quality, S: base count
    }
    sample_list = random_readnum(seed_num, len(fq), read_num)[0]
    for i in sample_list:
        read = fq[i]
        seqdict = readParse(read, seqdict)
        homopolymer_dict = homopolymerParse(read.seq, homopolymer_size_min, homopolymer_dict)
        endBaseQual_dict = endBaseQualParse(read.seq, read.quali, shift_length, endBaseQual_dict)
        allBaseQual_dict = allBaseQualParse(read.quali, split_part_num, allBaseQual_dict)
    return seqdict, homopolymer_dict, endBaseQual_dict, allBaseQual_dict

def overall_analyser(fq):
    seqdict = dict( {'ID': [], 'GC': [], 'LEN': [], 'QUAL1': [], 'QUAL2': []} )  # QUAL1: read basecall Q, QUAL2: read average Q
    homopolymer_size_min = 5
    homopolymer_dict = {'A':{}, 'G':{}, 'C':{}, 'T':{}}
    shift_length = 200
    endBaseQual_dict = {
                    'HeadBaseContent_dict': {'A': [0]*shift_length, 'G': [0]*shift_length, 'C': [0]*shift_length, 'T': [0]*shift_length, 'S':0},
                    'TailBaseContent_dict': {'A': [0]*shift_length, 'G': [0]*shift_length, 'C': [0]*shift_length, 'T': [0]*shift_length, 'S':0},
                    'HeadQualContent_dict': {'A': [0]*shift_length, 'G': [0]*shift_length, 'C': [0]*shift_length, 'T': [0]*shift_length, 'S':0},
                    'TailQualContent_dict': {'A': [0]*shift_length, 'G': [0]*shift_length, 'C': [0]*shift_length, 'T': [0]*shift_length, 'S':0},
                    }
    split_part_num = 100
    allBaseQual_dict = {
                    'PercentBaseQual_dict': {'Q':[0]* split_part_num, 'S':0} # Q: average quality, S: base count
    }
    for read in fq:
        seqdict = readParse(read, seqdict)
        homopolymer_dict = homopolymerParse(read.seq, homopolymer_size_min, homopolymer_dict)
        endBaseQual_dict = endBaseQualParse(read.seq, read.quali, shift_length, endBaseQual_dict)
        allBaseQual_dict = allBaseQualParse(read.quali, split_part_num, allBaseQual_dict)
        # print(len(seqdict['GC']))
        # print(sys.getsizeof(allBaseQual_dict))
    return seqdict, homopolymer_dict, endBaseQual_dict, allBaseQual_dict

# fq = pyfastx.Fastq('../test/ecoli.fq.gz')
# seq_qual_dict1, homopolymer_dict1, endBaseQual_dict1, allBaseQual_dict1  = sampling_analyser(fq, 1, 100)
# seq_qual_dict2, homopolymer_dict2, endBaseQual_dict2, allBaseQual_dict2 = overall_analyser(fq)

def get_fq_datum(fastq, mode):
    if mode == 'sampling':
        fq = pyfastx.Fastq(fastq)
        seq_qual_dict, homopolymer_dict, endBaseQual_dict, allBaseQual_dict = sampling_analyser(fq, 1, 100)
        sampling_fq_datum_dict = {
            'seq_qual_dict': seq_qual_dict,
            'homopolymer_dict': homopolymer_dict,
            'endBaseQual_dict': endBaseQual_dict,
            'allBaseQual_dict': allBaseQual_dict
        }
        return sampling_fq_datum_dict
    elif mode == 'overall':
        fq = pyfastx.Fastq(fastq)
        seq_qual_dict, homopolymer_dict, endBaseQual_dict, allBaseQual_dict = overall_analyser(fq)
        overall_fq_datum_dict = {
            'seq_qual_dict': seq_qual_dict,
            'homopolymer_dict': homopolymer_dict,
            'endBaseQual_dict': endBaseQual_dict,
            'allBaseQual_dict': allBaseQual_dict
        }
        return overall_fq_datum_dict
    else:
        print('please check the mode of parsing target data!')
fastq = '../test/ecoli.fq.gz'
# mode = 'sampling'
mode = 'overall'
merged_fq_datum_dict = get_fq_datum(fastq, mode)

import pprint

with open( '../test/ecoli.seq.json', 'w') as jsonfile:
    filewidth = len(merged_fq_datum_dict['seq_qual_dict']['ID'])+ 60
    pprint.pprint(merged_fq_datum_dict, jsonfile, indent=4, width=filewidth, depth = 5, compact=True)