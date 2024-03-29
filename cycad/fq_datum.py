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
    gc_count = 0
    for base in seq:
        if base == 'G' or base == 'C':
            gc_count += 1
    return round(gc_count/seq_len, 3)

def readAvgQscore(quali, seq_len):
    value_list = Counter(list(quali))
    value_sum = sum([k * v for k, v in value_list.items()])
    seq_qscore = round(value_sum / seq_len, 3)
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
    if len(seq) > (args["head_shift_length"] + args["tail_shift_length"]):
        endBaseQual_dict = endBaseHeadParse(seq, args["head_shift_length"], endBaseQual_dict)
        endBaseQual_dict = endBaseTailParse(seq, args["tail_shift_length"], endBaseQual_dict)
        endBaseQual_dict = endQualHeadParse(seq, quali, args["head_shift_length"], endBaseQual_dict)
        endBaseQual_dict = endQualTailParse(seq, quali, args["tail_shift_length"], endBaseQual_dict)
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
def random_readnum(seed_num, read_size, sample_num):
    np.random.seed(seed_num)
    sample_list = np.random.randint(1, read_size, (1, sample_num) )
    return sample_list
def sampling_analyser(args):
    seqdict = dict( {'ID': [], 'GC': [], 'LEN': [], 'QUAL1': [], 'QUAL2':[]} ) # QUAL1: read basecall Q, QUAL2: read average Q
    homopolymer_size_min = args["homopolymer_min_length"]
    homopolymer_dict = {'A':{}, 'G':{}, 'C':{}, 'T':{}}
    head_shift_length, tail_shift_length = args["head_shift_length"], args["tail_shift_length"]
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
    seed_num = 1
    read_num = 10000
    print(len(fq))
    if read_num<len(fq):
        sample_list = random_readnum(seed_num, len(fq), read_num)[0]
        print('random number')
    else:
        sample_list = range(len(fq))
        print('read number')

    for i in sample_list:
        read = fq[i]
        seqdict = readParse(i, read.name, read.seq, read.quali, seqdict)
        homopolymer_dict = homopolymerParse(read.seq, homopolymer_size_min, homopolymer_dict)
        endBaseQual_dict = endBaseQualParse(read.seq, read.quali, args, endBaseQual_dict)
        allBaseQual_dict = allBaseQualParse(read.quali, split_part_num, allBaseQual_dict)
    return seqdict, homopolymer_dict, endBaseQual_dict, allBaseQual_dict

def overall_analyser(args):
    seqdict = dict({'ID': [], 'GC': [], 'LEN': [], 'QUAL1': [], 'QUAL2': []})  # QUAL1: read basecall Q, QUAL2: read average Q
    homopolymer_size_min = args["homopolymer_min_length"]
    homopolymer_dict = {'A':{}, 'G':{}, 'C':{}, 'T':{}}
    head_shift_length, tail_shift_length = args["head_shift_length"], args["tail_shift_length"]
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
    if args["mode"] == 'sampling':
        seq_qual_dict, homopolymer_dict, endBaseQual_dict, allBaseQual_dict = sampling_analyser(args)
        sampling_fq_datum_dict = {
            'seq_qual_dict': seq_qual_dict,
            'homopolymer_dict': homopolymer_dict,
            'endBaseQual_dict': endBaseQual_dict,
            'allBaseQual_dict': allBaseQual_dict
        }
        return sampling_fq_datum_dict
    elif args["mode"] == 'overall':
        seq_qual_dict, homopolymer_dict, endBaseQual_dict, allBaseQual_dict = overall_analyser(args)
        overall_fq_datum_dict = {
            'seq_qual_dict': seq_qual_dict,
            'homopolymer_dict': homopolymer_dict,
            'endBaseQual_dict': endBaseQual_dict,
            'allBaseQual_dict': allBaseQual_dict
        }
        return overall_fq_datum_dict
    else:
        print('please check the mode of parsing target data!')

def fq_datum_action(args):
    if os.path.exists(args["fastq"]):
        merged_fq_datum_dict = get_fq_datum(args)
        output_folder = args["sample_name"]
        output_path = os.path.join(output_folder, "fq.pickle")
        with open(output_path, "wb") as f:
            filewidth = len(merged_fq_datum_dict['seq_qual_dict']['ID']) + 60
            pickle.dump(merged_fq_datum_dict, f)

if __name__ == "__main__" :
    parser = argparse.ArgumentParser()
    parser.add_argument("-fq",    "--fastq",     required=True,  help="sequences.fq/fq.gz")
    parser.add_argument("-P", "--platform", required=False, help="cyclone")
    parser.add_argument("-M", "--mode", type=str, default="overall", required=False, help="if you want fast, please set to sampling")
    parser.add_argument("-Hshift", "--head_shift_length", type=int, default=200, required=False, help="check head bases quality")
    parser.add_argument("-Tshift", "--tail_shift_length", type=int, default=200, required=False, help="check tail bases quality")
    parser.add_argument("-kmer", "--kmer_size_frequency",     type=int, default=5, required=False, help="observe kmer size specturm")
    parser.add_argument("-hpmin", "--homopolymer_min_length", type=int, default=2, required=False, help="observe minium homopolymer")
    parser.add_argument("-hpmax", "--homopolymer_max_length", type=int, default=9, required=False, help="observe maxium homopolymer")
    parser.add_argument("-name", "--sample_name", default='cycads_report', required=False, help="prefix of output file name")
    args = vars(parser.parse_args())
    fq_datum_action(args)
