
import re,os,sys,time
import pysam
import numpy as np
import argparse

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
    read_deltion = 0
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
    for i in np.unique(np.array(sorted(dict_ins.keys()) + sorted(dict_del.keys()))):
        if i in dict_ins.keys() and i in dict_del.keys():
            print(i, dict_ins[i], dict_del[i])
        elif i in dict_ins.keys() and i not in dict_del.keys():
            print(i, dict_ins[i], 0)
        elif i not in dict_ins.keys() and i in dict_del.keys():
            print(i, 0, dict_del[i])
        else:
            continue
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
