#!/usr/bin/env python3
# coding: utf-8

import re
import os
import pprint
import pyfastx
import argparse
import pickle
import numpy as np
from collections import Counter

def readGCcontent(seq, seq_len):
    C = seq.count("C")
    G = seq.count("G")
    return (C+G)/len(seq)

def readAvgQscore(quali, seq_len):
    value_list = Counter(list(quali))
    value_sum = sum([k * v for k, v in value_list.items()])
    seq_qscore = value_sum / seq_len
    return seq_qscore

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
        endBaseQual_dict['HeadQualContent_dict']['S'][b][i] += 1
        i += 1
    return endBaseQual_dict
def endQualTailParse(seq, quali, shift_length, endBaseQual_dict):
    i = 0
    for b in seq[-shift_length:]:
        endBaseQual_dict['TailQualContent_dict'][b][i] += quali[-shift_length:][i]
        endBaseQual_dict['TailQualContent_dict']['S'][b][i] += 1
        i += 1
    return endBaseQual_dict
def read_quality_to_bin_score(base_qual_list: list, split_number_of_read_length: int):
    div, mod = divmod(len(base_qual_list), split_number_of_read_length)
    qs_in_percentage_pos = [sum(base_qual_list[i * div + min(i, mod):(i + 1) * div + min(i + 1, mod)]) / len(
        base_qual_list[i * div + min(i, mod):(i + 1) * div + min(i + 1, mod)]) for i in
                            range(split_number_of_read_length)]
    return qs_in_percentage_pos
def allBaseQualParse(quali, split_part_num, allBaseQual_dict):
    if len(quali) >= split_part_num:
        qv_in_percentage_pos = read_quality_to_bin_score(quali, split_part_num)
        for j in range(split_part_num):
            allBaseQual_dict['PercentBaseQual_dict']['Q'][j] += qv_in_percentage_pos[j]
        allBaseQual_dict['PercentBaseQual_dict']['S'] += 1
    return allBaseQual_dict

def endBaseQualParse(seq, quali, args, endBaseQual_dict):
    if len(seq) > (args["check_terminal_bases"] + args["check_terminal_bases"]):
        endBaseQual_dict = endBaseHeadParse(seq, args["check_terminal_bases"], endBaseQual_dict)
        endBaseQual_dict = endBaseTailParse(seq, args["check_terminal_bases"], endBaseQual_dict)
        endBaseQual_dict = endQualHeadParse(seq, quali, args["check_terminal_bases"], endBaseQual_dict)
        endBaseQual_dict = endQualTailParse(seq, quali, args["check_terminal_bases"], endBaseQual_dict)
        endBaseQual_dict['HeadBaseContent_dict']['S'] += 1
        endBaseQual_dict['TailBaseContent_dict']['S'] += 1
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
    elif fq_path.endswith('fastq') or fq_path.endswith('fq'):
        os.system('cat {} | awk \'NR %4 == 1 || NR %4 == 2 \'   > {}.fasta '.format(fq_path, output))
        os.system('{} count  k={} {}.fasta output {}.meryl'.format(meryl, str(kmersize), output, output))
        os.system('{} print {}.meryl | sort -k2nr > {}_kmer_{}_freq.txt'.format(meryl, output, output, str(kmersize)))
        os.system('rm -f {}.fasta'.format(output))
        os.system('rm -rf {}.meryl'.format(output))
    else:
        print('please determining the input file suffix is fastq or fq or fq.gz!')
# kmerSpectrumParse('../test/ecoli.fq.gz', 5, '../test/')

def readParse(number, name, seq, quali, seqdict):
    seq_len = len(seq)
    seqdict['ID'].append(number)
    seqdict['GC'].append(readGCcontent(seq, seq_len))
    seqdict['LEN'].append(seq_len)
    seqdict['QUAL1'].append(float(name.split('_')[-1]))
    seqdict['QUAL2'].append(readAvgQscore(quali, seq_len))
    return seqdict
def sampling_analyser(args):
    seqdict = dict( {'ID': [], 'GC': [], 'LEN': [], 'QUAL1': [], 'QUAL2':[]} ) # QUAL1: read basecall Q, QUAL2: read average Q
    homopolymer_size_min = args["min_homopolymer_size"]
    homopolymer_dict = {'A':{}, 'G':{}, 'C':{}, 'T':{}}
    head_shift_length, tail_shift_length = args["check_terminal_bases"], args["check_terminal_bases"]
    endBaseQual_dict = {
                    'HeadBaseContent_dict': {'A': [0]*head_shift_length, 'G': [0]*head_shift_length, 'C': [0]*head_shift_length, 'T': [0]*head_shift_length, 'S':0},
                    'TailBaseContent_dict': {'A': [0]*tail_shift_length, 'G': [0]*tail_shift_length, 'C': [0]*tail_shift_length, 'T': [0]*tail_shift_length, 'S':0},
                    'HeadQualContent_dict': {'A': [0]*head_shift_length, 'G': [0]*head_shift_length, 'C': [0]*head_shift_length, 'T': [0]*head_shift_length,
                                             'S':{'A':[0]*head_shift_length, 'T':[0]*head_shift_length, 'G':[0]*head_shift_length, 'C': [0]*head_shift_length}},
                    'TailQualContent_dict': {'A': [0]*tail_shift_length, 'G': [0]*tail_shift_length, 'C': [0]*tail_shift_length, 'T': [0]*tail_shift_length,
                                             'S':{'A':[0]*tail_shift_length, 'T':[0]*tail_shift_length, 'G':[0]*tail_shift_length, 'C':[0]*tail_shift_length}},
                    }
    split_part_num = 100
    allBaseQual_dict = {'PercentBaseQual_dict': {'Q':[0]* split_part_num, 'S':0}} # Q: average quality, S: base count
    fq = pyfastx.Fastq(args["fastq"], build_index=False)
    seed = int(args["seed"])
    sample_size = int(args["sample"])
    read_count = len(fq)
    if sample_size<len(fq):
        selected_read_indices = np.random.default_rng(seed).integers(low=0, high=read_count, size=sample_size, dtype=np.uint64)
        print('random number')
    else:
        selected_read_indices = range(read_count)
        print('read number')

    for i in selected_read_indices:
        read = fq[i]
        seqdict = readParse(i, read.name, read.seq, read.quali, seqdict)
        homopolymer_dict = homopolymerParse(read.seq, homopolymer_size_min, homopolymer_dict)
        endBaseQual_dict = endBaseQualParse(read.seq, read.quali, args, endBaseQual_dict)
        allBaseQual_dict = allBaseQualParse(read.quali, split_part_num, allBaseQual_dict)
    return seqdict, homopolymer_dict, endBaseQual_dict, allBaseQual_dict

def overall_analyser(args):
    seqdict = dict({'ID': [], 'GC': [], 'LEN': [], 'QUAL1': [], 'QUAL2': []})  # QUAL1: read basecall Q, QUAL2: read average Q
    homopolymer_size_min = args["min_homopolymer_size"]
    homopolymer_dict = {'A':{}, 'G':{}, 'C':{}, 'T':{}}
    head_shift_length, tail_shift_length = args["check_terminal_bases"], args["check_terminal_bases"]
    endBaseQual_dict = {
                    'HeadBaseContent_dict': {'A': [0]*head_shift_length, 'G': [0]*head_shift_length, 'C': [0]*head_shift_length, 'T': [0]*head_shift_length, 'S':0},
                    'TailBaseContent_dict': {'A': [0]*tail_shift_length, 'G': [0]*tail_shift_length, 'C': [0]*tail_shift_length, 'T': [0]*tail_shift_length, 'S':0},
                    'HeadQualContent_dict': {'A': [0]*head_shift_length, 'G': [0]*head_shift_length, 'C': [0]*head_shift_length, 'T': [0]*head_shift_length,
                                             'S':{'A':[0]*head_shift_length, 'T':[0]*head_shift_length, 'G':[0]*head_shift_length, 'C': [0]*head_shift_length}},
                    'TailQualContent_dict': {'A': [0]*tail_shift_length, 'G': [0]*tail_shift_length, 'C': [0]*tail_shift_length, 'T': [0]*tail_shift_length,
                                             'S':{'A':[0]*tail_shift_length, 'T':[0]*tail_shift_length, 'G':[0]*tail_shift_length, 'C':[0]*tail_shift_length}},
                    }
    split_part_num = 100
    allBaseQual_dict = {'PercentBaseQual_dict': {'Q':[0]* split_part_num, 'S':0}} # Q: average quality, S: base count
    number = 0
    fq = pyfastx.Fastq(args["fastq"], build_index=False)
    for name, seq, qual in fq:
        number += 1
        quali = [ord(i) - 33 for i in qual]
        seqdict = readParse(number, name, seq, quali, seqdict)
        homopolymer_dict = homopolymerParse(seq, homopolymer_size_min, homopolymer_dict)
        endBaseQual_dict = endBaseQualParse(seq, quali, args, endBaseQual_dict)
        allBaseQual_dict = allBaseQualParse(quali, split_part_num, allBaseQual_dict)
    return seqdict, homopolymer_dict, endBaseQual_dict, allBaseQual_dict

def get_fq_datum(args):
    if args['sample'] > 0:
        seq_qual_dict, homopolymer_dict, endBaseQual_dict, allBaseQual_dict = sampling_analyser(args)
        sampling_fq_datum_dict = {
            'seq_qual_dict': seq_qual_dict,
            'homopolymer_dict': homopolymer_dict,
            'endBaseQual_dict': endBaseQual_dict,
            'allBaseQual_dict': allBaseQual_dict
        }
        return sampling_fq_datum_dict
    else:
        seq_qual_dict, homopolymer_dict, endBaseQual_dict, allBaseQual_dict = overall_analyser(args)
        overall_fq_datum_dict = {
            'seq_qual_dict': seq_qual_dict,
            'homopolymer_dict': homopolymer_dict,
            'endBaseQual_dict': endBaseQual_dict,
            'allBaseQual_dict': allBaseQual_dict
        }
        return overall_fq_datum_dict

def fq_datum_action(args):
    merged_fq_datum_dict = get_fq_datum(args)
    with open(args['fastq_pickle_path'], "wb") as f:
        filewidth = len(merged_fq_datum_dict['seq_qual_dict']['ID']) + 60
        pickle.dump(merged_fq_datum_dict, f)

