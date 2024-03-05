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
    insertion_frequency_dict, deletion_frequency_dict = {}, {}
    list_read_identity_rate = []
    total_count, uniq_mapped_count = 0, 0
    f0, f1, f2, f4, f8, f16, f32, f64, f128, f256, f512, f1024, f2048 = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    readprofile = pysam.AlignmentFile(bam, "rb")
    for read in readprofile.fetch():
        name = read.query_name
        flag = read.flag
        mapq = read.mapping_quality
        total_count += 1
        if flag == 0 or flag == 16:
            uniq_mapped_count += 1
            dict_var, insertion_frequency_dict, deletion_frequency_dict, read_identity_rate = mutation_summary(
                read.cigar, dict_var, insertion_frequency_dict, deletion_frequency_dict)
            list_read_identity_rate.append(read_identity_rate)
            uniq_mapped_count += 1
    mapping_rate = round(uniq_mapped_count * 100 / total_count, 3)
    total_items = dict_var['match'][1] + dict_var['mismatch'][1] + dict_var['insertion'][1] +  dict_var['deletion'][1]
    identity_rate = round((dict_var['match'][1] * 100) / total_items, 3)
    mismatch_error_rate = round((dict_var['mismatch'][1] * 100) / total_items, 3)
    deletion_error_rate = round((dict_var['deletion'][1] * 100) / total_items, 3)
    insertion_error_rate = round((dict_var['insertion'][1] * 100)/ total_items, 3)
    mapping_summary_dict = {
        'alignment_result_statistics_dict': {
            'mapping_rate': mapping_rate,
            'uniq_mapped_count': uniq_mapped_count,
            'total_count': total_count,
            'identity_rate': identity_rate,
            'mismatch_error_rate': mismatch_error_rate,
            'deletion_error_rate': deletion_error_rate,
            'insertion_error_rate': insertion_error_rate
        },
       'insertion_frequency_dict': [list(insertion_frequency_dict.keys()), list(insertion_frequency_dict.values())],
       'deletion_frequency_dict': [list(deletion_frequency_dict.keys()), list(deletion_frequency_dict.values())]
    }
    return mapping_summary_dict

bam = '../test/ecoli.sorted.bam'
mapping_summary_dict = mapping_summary(bam)
import pprint

with open('../test/ecoli.bam.json', 'w') as jsonfile:
    file_width = len(mapping_summary_dict['deletion_frequency_dict'][0]) + 1000
    pprint.pprint(mapping_summary_dict, jsonfile, indent=4, width=file_width, depth=5, compact=True)
