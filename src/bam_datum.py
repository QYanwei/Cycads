#!/usr/bin/env python3
# coding: utf-8
import re
import pysam
import pyfastx

mismatch_event_dict = { 'AT': 0, 'AC': 0, 'AG': 0, 'TC': 0, 'CG': 0, 'TG': 0,
                        'TA': 0, 'CA': 0, 'GA': 0, 'CT': 0, 'GC': 0, 'GT': 0}
deletion_event_dict = { }
insertion_event_dict = { }


def alignment_ref_seq(raw_ref, raw_seq, cigar_tuples, ):
    new_ref = ''
    new_seq = ''
    for i in cigar_tuples:
        if i[0] == 0: # matched
            new_seq += raw_seq[: i[1]]
            new_ref += raw_seq[: i[1]]
            raw_seq = raw_seq[i[1]:]
            raw_ref = raw_ref[i[1]:]

        elif i[0] == 1: # insertion
            new_seq += raw_seq[:i[1]]
            new_ref += i[1] * '-'
            raw_seq = raw_seq[i[1]:]
            raw_ref = raw_ref
        elif i[0] == 2: # deletion
            new_seq += i[1] * '-'
            new_ref += raw_ref[:i[1]]
            raw_seq = raw_seq
            raw_ref = raw_ref[i[1]:]
        elif i[0] == 8: #mismatch
            new_seq += raw_seq[:i[1]]
            new_ref += raw_ref[:i[1]]
            raw_seq = raw_seq[i[1]:]
            raw_seq = raw_seq[i[1]:]
        elif i[0] == 4 or i[0] == 5:
            raw_seq = raw_seq[i[1]:]
        else:
            print(i)
    return new_ref, new_seq

def softclip_hardclip_remove(read_length, cigarstring):
    cigar_list = []
    soft_hard_length = 0
    pattern = r"([0-9]*[S=XMDIH]{1})"
    for match in re.finditer(pattern, cigarstring):
        cigar_list += [match.group()]
    # print(cigar_list)
    for cigar in cigar_list:
        if "S" in cigar or "H" in cigar:
            # print(cigar, int(cigar[:-1]))
            soft_hard_length += int(cigar[:-1])
    # print(soft_hard_length)
    return soft_hard_length/read_length

def mapping_summary(bam):
    readprofile = pysam.AlignmentFile(bam, "rb")
    for read in readprofile.fetch():
        flag = read.flag
        qlen = read.qlen
        alen = read.reference_length
        cstr = read.cigarstring
        haln = softclip_hardclip_remove(qlen, cstr)

        if read.is_forward and (flag == 0 or flag == 16) and haln <= 0.5 and qlen > 100:
            raw_seq = read.seq
            raw_ref = read.get_reference_sequence()
            new_seq, new_ref = alignment_ref_seq(raw_ref, raw_seq, read.cigartuples)
            print('new_seq', new_seq[:10], new_seq[-10:], new_seq)
            print('new_ref', new_ref[:10], new_ref[-10:], new_ref)
            print(read.cigarstring)
        else:
            continue


bam = '../test/ecoli.sorted.bam'
fasta = '../ref/Reference_Ecoli.fasta'
fa = pyfastx.Fasta(fasta)
# mapping_summary_dict = mapping_summary(bam)
# import pprint
#
# with open('../test/ecoli.bam.json', 'w') as jsonfile:
#     file_width = len(mapping_summary_dict['deletion_frequency_dict'][0]) + 60
#     pprint.pprint(mapping_summary_dict, jsonfile, indent=4, width=file_width, depth=4, compact=True)
