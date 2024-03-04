#!/usr/bin/env python3
# coding: utf-8
import re,os,sys,time
import pysam
import numpy as np

def mutation_summary(cigar_tuples, dict_var, dict_ins, dict_del):
    cigar='''M  BAM_CMATCH      0
       I        BAM_CINS        1
       D        BAM_CDEL        2
       N        BAM_CREF_SKIP   3
       S        BAM_CSOFT_CLIP  4
       H        BAM_CHARD_CLIP  5
       P        BAM_CPAD        6
       =        BAM_CEQUAL      7
       X        BAM_CDIFF       8
       B        BAM_CBACK       9'''
    total_items = 0
    total_match = 0
    read_insertion = 0
    read_deletion = 0
    for i in cigar_tuples:
        if i[0] == 0:
            dict_var['match'][0] += 1
            dict_var['match'][1] += i[1]
            total_match += i[1]
            total_items += i[1]
        elif i[0] == 8:
            dict_var['mismatch'][0] += 1
            dict_var['mismatch'][1] += i[1]
            total_items += i[1]
        elif i[0] == 1:
            dict_var['insertion'][0] += 1
            dict_var['insertion'][1] += i[1]
            total_items += i[1]
            if i[1] in dict_ins.keys():
                dict_ins[ i[1] ] += 1
            else:
                dict_ins[ i[1] ] = 1
        elif i[0] == 2:
            dict_var['deletion'][0] += 1
            dict_var['deletion'][1] += i[1]
            total_items += i[1]
            if i[1] in dict_del.keys():
                dict_del[ i[1] ] += 1
            else:
                dict_del[ i[1] ] = 1
        elif i[0] == 4:
            dict_var['softclip'][0] += 1
            dict_var['softclip'][1] += i[1]
        elif i[0] ==  5:
            dict_var['hardclip'][0] += 1
            dict_var['hardclip'][1] += i[1]
        else:
            continue
    read_identity_rate = total_match / total_items
    return dict_var, dict_ins, dict_del, read_identity_rate

def mapping_summary(bam):
    dict_var = {'match'     : [0, 0],
                'mismatch'  : [0, 0],
                'insertion' : [0, 0],
                'deletion'  : [0, 0],
                'softclip'  : [0, 0],
                'hardclip'  : [0, 0]}
    dict_ins, dict_del = {}, {}
    list_read_identity_rate = []
    total_count, mapped_count = 0, 0
    mapped, unmap = 0, 0
    samfile = pysam.AlignmentFile(bam, "rb")
    for read in samfile.fetch():
        name = read.query_name
        flag = read.flag
        mapq = read.mapping_quality
        total_count += 1
        if flag != 4 and flag != 256 and flag != 2048:
            mapped_count += 1
        if flag == 0:
            mapped += 1
        if flag == 4:
            unmap += 1
        if len(read.get_tags()) >= 6 and flag <= 16:
            dict_var, dict_ins, dict_del, read_identity_rate = mutation_summary(read.cigar, dict_var, dict_ins, dict_del)
            list_read_identity_rate.append(read_identity_rate)

    mapping_rate = mapped_count / total_count
    total_items = dict_var['match'][1] + dict_var['mismatch'][1] + dict_var['insertion'][1] +  dict_var['deletion'][1]
    identity_rate = dict_var['match'][1] / total_items
    mismatch_error_rate = dict_var['mismatch'][1] / total_items
    deletion_error_rate = dict_var['deletion'][1] / total_items
    insertion_error_rate = dict_var['insertion'][1] / total_items
    print("Mapping rate/ mapped_read /total_count: ", mapping_rate, mapped_count, total_count)
    print("Mapped/Unmapped read count:", mapped, unmap)
    print("Identity rate: ", identity_rate)
    print("Mismatch rate: ", mismatch_error_rate)
    print("Deletion rate: ", deletion_error_rate)
    print("Insertion rate: ", insertion_error_rate)
    # for i in np.unique(np.array(sorted(dict_ins.keys()) + sorted(dict_del.keys()))):
    #     if i in dict_ins.keys() and i in dict_del.keys():
    #         print(i, dict_ins[i], dict_del[i])
    #     elif i in dict_ins.keys() and i not in dict_del.keys():
    #         print(i, dict_ins[i], 0)
    #     elif i not in dict_ins.keys() and i in dict_del.keys():
    #         print(i, 0, dict_del[i])
    #     else:
    #         continue
mapping_file_report = {
"Sample name" : '',
"Identity rate": 0,
"Mismatch rate": 0,
"Deletion rate": 0,
"Insertion rate": 0,
"Identity rate list": [],
"Indel distribution dict": {},
"Mean mapping quality": 0,
"Mapping quality list": [],
"Mapping rate": 0,
"Coverage rate": 0,
"Coverage depth": 0,
}
bam = '../test/ecoli.sorted.bam'
mapping_summary(bam)

#!/usr/bin/env python3
# coding: utf-8
import re,os,sys,time
import json
import pyfastx

import numpy as np
from collections import Counter

def readGCcontent(seq):
    gc_count = 0
    for base in seq:
        if base == 'G' or base == 'C':
            gc_count += 1
    return round(gc_count / len(seq),3)

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
    endBaseQual_dict['HeadQualContent_dict']['S'][b] += 1
    return endBaseQual_dict
def endQualTailParse(seq, quali, shift_length, endBaseQual_dict):
    i = 0
    for b in seq[-shift_length:]:
        endBaseQual_dict['TailQualContent_dict'][b][i] += quali[-shift_length:][i]
        i += 1
    endBaseQual_dict['TailQualContent_dict']['S'][b] += 1
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
    if len(seq) > shift_length * 2:
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
                    'HeadQualContent_dict': {'A': [0]*shift_length, 'G': [0]*shift_length, 'C': [0]*shift_length, 'T': [0]*shift_length, 'S':{'A':0, 'T':0, 'G':0, 'C':0}},
                    'TailQualContent_dict': {'A': [0]*shift_length, 'G': [0]*shift_length, 'C': [0]*shift_length, 'T': [0]*shift_length, 'S':{'A':0, 'T':0, 'G':0, 'C':0}},
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
                    'HeadQualContent_dict': {'A': [0]*shift_length, 'G': [0]*shift_length, 'C': [0]*shift_length, 'T': [0]*shift_length, 'S':{'A':0, 'T':0, 'G':0, 'C':0}},
                    'TailQualContent_dict': {'A': [0]*shift_length, 'G': [0]*shift_length, 'C': [0]*shift_length, 'T': [0]*shift_length, 'S':{'A':0, 'T':0, 'G':0, 'C':0}},
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
    return seqdict, homopolymer_dict, endBaseQual_dict, allBaseQual_dict

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
    pprint.pprint(merged_fq_datum_dict, jsonfile, indent=4, width=filewidth, depth = 6, compact=True)