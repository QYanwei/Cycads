#!/usr/bin/env python3
# coding: utf-8
import re
import pysam
import pprint

def try_extend_reference_homopolymer(start_pos, end_pos, sequence, base_homopolymer):
    if sequence[start_pos] in ("-", base_homopolymer):
        while start_pos-1 >= 0 and sequence[start_pos-1] in ("-", base_homopolymer):
            start_pos -= 1
    if sequence[end_pos-1] in ("-", base_homopolymer):
        while end_pos < len(sequence) and sequence[end_pos] in ("-", base_homopolymer):
            end_pos += 1
    while sequence[start_pos] == "-":
        start_pos += 1
    while sequence[end_pos-1] == "-":
        end_pos -= 1
    return start_pos, end_pos
def try_extend_sequence_homopolymer(start_pos, end_pos, sequence, base_homopolymer):
    if sequence[start_pos: end_pos].replace("-", "") == "":
        return start_pos, end_pos
    if sequence[start_pos] == base_homopolymer:
        while start_pos-1 >= 0 and sequence[start_pos-1] == base_homopolymer:
            start_pos -= 1
    if sequence[end_pos-1] == base_homopolymer:
        while end_pos < len(sequence) and sequence[end_pos] == base_homopolymer:
            end_pos += 1
    return start_pos, end_pos
def parsing_homopolymer_error_event(new_ref, new_seq):
    start_init, end_init = -1, -1
    hpm_map, hpm_mis, hpm_ins, hpm_del = 0, 0, 0, 0
    # homopolymers regular expression pattern
    homopolymer_pattern = "A{2,}|C{2,}|G{2,}|T{2,}"
    for homo_pos in re.finditer(homopolymer_pattern, new_ref):
        start_loci, end_loci = homo_pos.span()
        if start_loci in range(start_init, end_init):
            continue
        base_homopolymer = homo_pos.group()[0]
        start_ref, end_ref = try_extend_reference_homopolymer(start_loci, end_loci, new_ref, base_homopolymer)
        homo_ref = new_ref[start_ref: end_ref]
        start_seq, end_seq = try_extend_sequence_homopolymer(start_loci, end_loci, new_seq, base_homopolymer)
        start_init = min(start_ref, start_seq)
        end_init = max(end_ref, end_seq)
        homopolymer_ref_segment = new_ref[start_init:end_init]
        homopolymer_seq_segment = new_seq[start_init:end_init]
        homo_ref_length = len(homo_ref.replace("-", ""))
        if homopolymer_ref_segment != homopolymer_seq_segment:
            if "-" in homopolymer_ref_segment:
                hpm_ins += 1
                ins_size = homopolymer_ref_segment.count("-")
                if homo_ref_length < homo_length_max:  # [···, shift_length+1=ins1, shift_length+2=ins2, ···]
                    if ins_size < shift_length:
                        homopolymer_aln_event_stat_dict[base_homopolymer][homo_ref_length][
                            shift_length + ins_size + 2] += 1
                        homopolymer_aln_event_stat_dict['S'][homo_ref_length][shift_length + ins_size + 2] += 1
                    else:
                        homopolymer_aln_event_stat_dict[base_homopolymer][homo_ref_length][shift_length * 2 + 1] += 1
                        homopolymer_aln_event_stat_dict['S'][homo_ref_length][shift_length * 2 + 1] += 1
                else:
                    if ins_size < shift_length:
                        homopolymer_aln_event_stat_dict[base_homopolymer][homo_length_max][
                            shift_length + ins_size + 2] += 1
                        homopolymer_aln_event_stat_dict['S'][homo_length_max][shift_length + ins_size + 2] += 1
                    else:
                        homopolymer_aln_event_stat_dict[base_homopolymer][homo_length_max][shift_length * 2 + 1] += 1
                        homopolymer_aln_event_stat_dict['S'][homo_length_max][shift_length * 2 + 1] += 1
            elif "-" in homopolymer_seq_segment:
                hpm_del += 1
                del_size = homopolymer_seq_segment.count("-")
                if homo_ref_length < homo_length_max:
                    if del_size < shift_length:
                        homopolymer_aln_event_stat_dict[base_homopolymer][homo_ref_length][shift_length - del_size] += 1
                        homopolymer_aln_event_stat_dict['S'][homo_ref_length][shift_length - del_size] += 1
                    else:
                        homopolymer_aln_event_stat_dict[base_homopolymer][homo_ref_length][
                            shift_length - shift_length] += 1
                        homopolymer_aln_event_stat_dict['S'][homo_ref_length][shift_length - shift_length] += 1
                else:
                    if del_size < shift_length:
                        homopolymer_aln_event_stat_dict[base_homopolymer][homo_length_max][shift_length - del_size] += 1
                        homopolymer_aln_event_stat_dict['S'][homo_length_max][shift_length - del_size] += 1
                    else:
                        homopolymer_aln_event_stat_dict[base_homopolymer][homo_length_max][
                            shift_length - shift_length] += 1
                        homopolymer_aln_event_stat_dict['S'][homo_length_max][shift_length - shift_length] += 1
            else:
                hpm_mis += 1
                if homo_ref_length < homo_length_max:
                    homopolymer_aln_event_stat_dict[base_homopolymer][homo_ref_length][shift_length] += 1
                    homopolymer_aln_event_stat_dict['S'][homo_ref_length][shift_length] += 1
                else:
                    homopolymer_aln_event_stat_dict[base_homopolymer][homo_length_max][shift_length] += 1
                    homopolymer_aln_event_stat_dict['S'][homo_length_max][shift_length] += 1
        else:
            hpm_map += 1
            if homo_ref_length < homo_length_max:
                homopolymer_aln_event_stat_dict[base_homopolymer][homo_ref_length][shift_length + 1] += 1
                homopolymer_aln_event_stat_dict['S'][homo_ref_length][shift_length + 1] += 1
            else:
                homopolymer_aln_event_stat_dict[base_homopolymer][homo_length_max][shift_length] += 1
                homopolymer_aln_event_stat_dict['S'][homo_length_max][shift_length] += 1
    return hpm_map, hpm_mis, hpm_ins, hpm_del
def parsing_alignment_events(raw_ref, raw_seq, cigar_tuples):
    new_ref, new_seq = '', ''
    qry_map, qry_mis, qry_del, qry_ins = 0, 0, 0, 0
    for i in cigar_tuples:
        if i[0] == 0:  # matched
            new_seq += raw_seq[: i[1]]
            new_ref += raw_seq[: i[1]]
            raw_seq = raw_seq[i[1]:]
            raw_ref = raw_ref[i[1]:]
            qry_map += i[1]
            if i[1] in overall_aln_event_stat_dict['all_map_len_dict'].keys():
                overall_aln_event_stat_dict['all_map_len_dict'][i[1]] += 1
            else:
                overall_aln_event_stat_dict['all_map_len_dict'][i[1]] = 1
        elif i[0] == 1:  # insertion
            new_seq += raw_seq[:i[1]]
            new_ref += i[1] * '-'
            raw_seq = raw_seq[i[1]:]
            raw_ref = raw_ref
            qry_ins += 1
            if i[1] in overall_aln_event_stat_dict['all_ins_len_dict'].keys():
                overall_aln_event_stat_dict['all_ins_len_dict'][i[1]] += 1
            else:
                overall_aln_event_stat_dict['all_ins_len_dict'][i[1]] = 1
        elif i[0] == 2:  # deletion
            new_seq += i[1] * '-'
            new_ref += raw_ref[:i[1]]
            raw_seq = raw_seq
            raw_ref = raw_ref[i[1]:]
            qry_del += 1
            if i[1] in overall_aln_event_stat_dict['all_del_len_dict'].keys():
                overall_aln_event_stat_dict['all_del_len_dict'][i[1]] += 1
            else:
                overall_aln_event_stat_dict['all_del_len_dict'][i[1]] = 1
        elif i[0] == 8:  # mismatch
            new_seq += raw_seq[:i[1]]
            new_ref += raw_ref[:i[1]]
            raw_seq = raw_seq[i[1]:]
            raw_ref = raw_ref[i[1]:]
            qry_mis += 1
            overall_aln_event_stat_dict['all_mis_typ_dict'][raw_ref[i[1]:] + raw_seq[i[1]:]] += 1
        elif i[0] == 4 or i[0] == 5:
            raw_seq = raw_seq[i[1]:]
        else:
            pass
    # overall error and identity events
    query_mis_rate = round(qry_mis/(qry_map + qry_ins + qry_del + qry_mis), 5)
    query_ins_rate = round(qry_ins/(qry_map + qry_ins + qry_del + qry_mis), 5)
    query_del_rate = round(qry_del/(qry_map + qry_ins + qry_del + qry_mis), 5)
    query_idy_rate = round(qry_map/(qry_map + qry_ins + qry_del + qry_mis), 5)
    query_dif_rate = round(1 - query_idy_rate, 5)
    query_aln_event_stat_dict['qry_mis_rate'].append(query_mis_rate)
    query_aln_event_stat_dict['qry_ins_rate'].append(query_ins_rate)
    query_aln_event_stat_dict['qry_del_rate'].append(query_del_rate)
    query_aln_event_stat_dict['qry_dif_rate'].append(query_dif_rate)
    query_aln_event_stat_dict['qry_idy_rate'].append(query_idy_rate)
    overall_aln_event_sum_dict['identity'] += qry_map
    overall_aln_event_sum_dict['substitution'] += qry_mis
    overall_aln_event_sum_dict['expansion'] += qry_ins
    overall_aln_event_sum_dict['contraction'] += qry_del
    # homopolymer events statistics
    hpm_map, hpm_mis, hpm_del, hpm_ins = parsing_homopolymer_error_event(new_ref, new_seq)
    qry_hpm_mis_rate = round(hpm_mis/(qry_map + qry_ins + qry_del + qry_mis), 5)
    qry_hpm_ins_rate = round(hpm_ins/(qry_map + qry_ins + qry_del + qry_mis), 5)
    qry_hpm_del_rate = round(hpm_del/(qry_map + qry_ins + qry_del + qry_mis), 5)
    qry_hpm_idy_rate = round(hpm_map/(hpm_map + hpm_mis + hpm_del + hpm_ins), 5)
    qry_hpm_dif_rate = round(1-qry_hpm_idy_rate, 5)
    query_aln_event_stat_dict['qry_hpm_mis_rate'].append(qry_hpm_mis_rate)
    query_aln_event_stat_dict['qry_hpm_ins_rate'].append(qry_hpm_ins_rate)
    query_aln_event_stat_dict['qry_hpm_del_rate'].append(qry_hpm_del_rate)
    query_aln_event_stat_dict['qry_hpm_idy_rate'].append(qry_hpm_idy_rate)
    query_aln_event_stat_dict['qry_hpm_dif_rate'].append(qry_hpm_dif_rate)
    overall_aln_event_sum_dict['hpm_substitution'] += hpm_mis
    overall_aln_event_sum_dict['hpm_expansion'] += hpm_ins
    overall_aln_event_sum_dict['hpm_contraction'] += hpm_del
    # non-homopolymer events statistics
    qry_non_hpm_mis_rate = round((qry_mis - hpm_mis)/(qry_map + qry_ins + qry_del + qry_mis), 5)
    qry_non_hpm_ins_rate = round((qry_ins - hpm_ins)/(qry_map + qry_ins + qry_del + qry_mis), 5)
    qry_non_hpm_del_rate = round((qry_del - hpm_del)/(qry_map + qry_ins + qry_del + qry_mis), 5)
    hpm_mid_number = hpm_mis - hpm_ins - hpm_del
    qry_non_hpm_idy_rate = round((qry_map - hpm_mid_number)/(qry_map + qry_ins + qry_del + qry_mis - hpm_mid_number), 5)
    qry_non_hpm_dif_rate = round(1 - qry_non_hpm_idy_rate, 5)
    overall_aln_event_sum_dict['non_hpm_substitution'] += qry_mis - hpm_mis
    overall_aln_event_sum_dict['non_hpm_expansion'] += qry_ins - hpm_ins
    overall_aln_event_sum_dict['non_hpm_contraction'] += qry_del - hpm_del
    query_aln_event_stat_dict['qry_non_hpm_mis_rate'].append(qry_non_hpm_mis_rate)
    query_aln_event_stat_dict['qry_non_hpm_ins_rate'].append(qry_non_hpm_ins_rate)
    query_aln_event_stat_dict['qry_non_hpm_del_rate'].append(qry_non_hpm_del_rate)
    query_aln_event_stat_dict['qry_non_hpm_dif_rate'].append(qry_non_hpm_dif_rate)
    query_aln_event_stat_dict['qry_non_hpm_idy_rate'].append(qry_non_hpm_idy_rate)

def short_softclip_hardclip_discard(cigar_tuples):
    soft_hard_length = 0
    for i in cigar_tuples:
        if i[0] == 4 or i[0] == 5:
            soft_hard_length += i[1]
    return soft_hard_length

def mapping_summary(in_bam):
    readprofile = pysam.AlignmentFile(in_bam, "rb")
    for read in readprofile.fetch():
        flag = read.flag
        cigar_tuples = read.cigartuples
        query_length = read.qlen
        if read.is_forward and flag == 0 or flag == 16:
            skipped_length = short_softclip_hardclip_discard(cigar_tuples)
            if skipped_length/query_length <= 0.5 and query_length > 100:
                raw_seq = read.seq
                raw_ref = read.get_reference_sequence()
                parsing_alignment_events(raw_ref, raw_seq, cigar_tuples)
            else:
                pass
        else:
            continue


if __name__ == "__main__":
    # parameter list
    bam = '../test/ecoli.sorted.bam'
    homo_length_max = 16
    shift_length = 4
    # initiated dict
    overall_aln_event_stat_dict = {
        'all_mis_typ_dict': {'AT': 0, 'AC': 0, 'AG': 0, 'TC': 0, 'CG': 0, 'TG': 0, 'TA': 0, 'CA': 0, 'GA': 0, 'CT': 0,
                             'GC': 0, 'GT': 0},
        'all_hpm_len_dict': dict(),
        'all_map_len_dict': dict(),
        'all_ins_len_dict': dict(),
        'all_del_len_dict': dict()
    }
    query_aln_event_stat_dict = {
        'qry_idy_rate': list(),
        'qry_dif_rate': list(),
        'qry_mis_rate': list(),
        'qry_ins_rate': list(),
        'qry_del_rate': list(),
        'qry_hpm_idy_rate': list(),
        'qry_hpm_dif_rate': list(),
        'qry_hpm_mis_rate': list(),
        'qry_hpm_ins_rate': list(),
        'qry_hpm_del_rate': list(),
        'qry_non_hpm_idy_rate': list(),
        'qry_non_hpm_dif_rate': list(),
        'qry_non_hpm_mis_rate': list(),
        'qry_non_hpm_ins_rate': list(),
        'qry_non_hpm_del_rate': list(),
    }
    homopolymer_aln_event_stat_dict = {
        'S': {i: [0] * (shift_length * 2 + 2) for i in range(2, homo_length_max + 1, 1)},
        'A': {i: [0] * (shift_length * 2 + 2) for i in range(2, homo_length_max + 1, 1)},
        'T': {i: [0] * (shift_length * 2 + 2) for i in range(2, homo_length_max + 1, 1)},
        'C': {i: [0] * (shift_length * 2 + 2) for i in range(2, homo_length_max + 1, 1)},
        'G': {i: [0] * (shift_length * 2 + 2) for i in range(2, homo_length_max + 1, 1)},
    }
    overall_aln_event_sum_dict = {
        'identity': 0,
        'substitution': 0,
        'expansion': 0,
        'contraction': 0,
        'hpm_substitution': 0,
        'hpm_expansion': 0,
        'hpm_contraction': 0,
        'non_hpm_substitution': 0,
        'non_hpm_expansion': 0,
        'non_hpm_contraction': 0,
    }
    # analyzing data
    mapping_summary(bam)
    # alignment dict
    merge_alignment_dict = {
        'overall_aln_event_sum_dict': overall_aln_event_sum_dict,
        'overall_aln_event_stat_dict': overall_aln_event_stat_dict,
        'query_aln_event_stat_dict': query_aln_event_stat_dict,
        'homopolymer_aln_event_stat_dict': homopolymer_aln_event_stat_dict
    }
    with open('../test/ecoli.bam.json', 'w') as jsonfile:
        file_width = len(merge_alignment_dict['query_aln_event_stat_dict']['qry_idy_rate']) * 6 + 60
        pprint.pprint(merge_alignment_dict, jsonfile, indent=4, width=file_width, depth=6, compact=True)
