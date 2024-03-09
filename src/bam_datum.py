#!/usr/bin/env python3
# coding: utf-8
import re
import pysam
import pyfastx

aln_event_sum_dict = {
    'mis_typ_dict' : { 'AT': 0, 'AC': 0, 'AG': 0, 'TC': 0, 'CG': 0, 'TG': 0, 'TA': 0, 'CA': 0, 'GA': 0, 'CT': 0, 'GC': 0, 'GT': 0},
    'del_len_dict' : dict(),
    'ins_len_dict' : dict(),
    'map_len_dict' : dict(),
    'mis_len_dict' : dict(),
    'hpm_len_dict' : dict(),
    'hpm_mis_len_dict' : dict(),
    'hpm_ins_len_dict' : dict(),
    'hpm_del_len_dict' : dict(),
    'query_idy_rate' : list(),
    'query_dif_rate' : list(),
    'query_ins_rate' : list(),
    'query_del_rate' : list()
}

def parsing_alignment_events(raw_ref, raw_seq, cigar_tuples, aln_event_sum_dict):
    new_ref, new_seq = '', ''
    new_map, new_mis, new_del, new_ins = 0, 0, 0, 0
    for i in cigar_tuples:
        if i[0] == 0: # matched
            new_seq += raw_seq[: i[1]]
            new_ref += raw_seq[: i[1]]
            raw_seq = raw_seq[i[1]:]
            raw_ref = raw_ref[i[1]:]
            new_map += 1
            if i[1] in aln_event_sum_dict['map_len_dict'].keys():
                aln_event_sum_dict['map_len_dict'][ i[1] ] += 1
            else:
                aln_event_sum_dict['map_len_dict'][ i[1] ] = 1
        elif i[0] == 1: # insertion
            new_seq += raw_seq[:i[1]]
            new_ref += i[1] * '-'
            raw_seq = raw_seq[i[1]:]
            raw_ref = raw_ref
            new_mis += 1
            if i[1] in aln_event_sum_dict['ins_len_dict'].keys():
                aln_event_sum_dict['ins_len_dict'][ i[1] ] += 1
            else:
                aln_event_sum_dict['ins_len_dict'][ i[1] ] = 1
        elif i[0] == 2: # deletion
            new_seq += i[1] * '-'
            new_ref += raw_ref[:i[1]]
            raw_seq = raw_seq
            raw_ref = raw_ref[i[1]:]
            new_del += 1
            if i[1] in aln_event_sum_dict['del_len_dict'].keys():
                aln_event_sum_dict['del_len_dict'][ i[1] ] += 1
            else:
                aln_event_sum_dict['del_len_dict'][ i[1] ] = 1
        elif i[0] == 8: # mismatch
            new_seq += raw_seq[:i[1]]
            new_ref += raw_ref[:i[1]]
            raw_seq = raw_seq[i[1]:]
            raw_ref = raw_ref[i[1]:]
            new_mis += 1
            aln_event_sum_dict['mis_typ_dict'][raw_ref[i[1]:] + raw_seq[i[1]:]] += 1
        elif i[0] == 4 or i[0] == 5:
            raw_seq = raw_seq[i[1]:]
        else:
            pass
    query_ins_rate = round(new_ins/(new_map + new_ins + new_del + new_mis), 6)
    query_del_rate = round(new_del/(new_map + new_ins + new_del + new_mis), 6)
    query_idy_rate = round(new_map/(new_map + new_ins + new_del + new_mis), 6)
    query_dif_rate = round(1 - query_idy_rate, 6)
    aln_event_sum_dict['query_ins_rate'].append(query_ins_rate)
    aln_event_sum_dict['query_del_rate'].append(query_del_rate)
    aln_event_sum_dict['query_idy_rate'].append(query_idy_rate)
    aln_event_sum_dict['query_dif_rate'].append(query_dif_rate)
    return new_ref, new_seq, aln_event_sum_dict

def short_softclip_hardclip_discard(cigar_tuples):
    soft_hard_length = 0
    for i in cigar_tuples:
        if i[0] == 4 or i[0] == 5:
            soft_hard_length += i[1]
    return soft_hard_length

def mapping_summary(bam, aln_event_sum_dict):
    readprofile = pysam.AlignmentFile(bam, "rb")
    for read in readprofile.fetch():
        flag = read.flag
        cigar_string = read.cigarstring
        cigar_tuples = read.cigartuples
        query_length = read.qlen
        if read.is_forward and flag == 0 or flag == 16:
            skipped_length = short_softclip_hardclip_discard(cigar_tuples)
            if skipped_length/query_length <= 0.5 and query_length > 100:
                raw_seq = read.seq
                raw_ref = read.get_reference_sequence()
                new_seq, new_ref, aln_event_sum_dict = parsing_alignment_events(raw_ref, raw_seq, cigar_tuples, aln_event_sum_dict)
                print('new_seq', new_seq[:10], new_seq[-10:], new_seq)
                print('new_ref', new_ref[:10], new_ref[-10:], new_ref)
            else:
                pass
        else:
            continue

bam = '../test/ecoli.sorted.bam'
fasta = '../ref/Reference_Ecoli.fasta'
fa = pyfastx.Fasta(fasta)
mapping_summary(bam, aln_event_sum_dict)
print(aln_event_sum_dict)
import pprint
#
with open('../test/ecoli.bam.json', 'w') as jsonfile:
    file_width = len(aln_event_sum_dict['query_idy_rate'] ) + 60
    pprint.pprint(aln_event_sum_dict, jsonfile, indent=4, width=file_width, depth=4, compact=True)
