#!/usr/bin/env python3
# coding: utf-8

import os,re,sys,time, pickle
import argparse
import pysam

def short_softclip_hardclip_discard(cigar_tuples):
    soft_hard_length = 0
    for i in cigar_tuples:
        if i[0] == 4 or i[0] == 5:
            soft_hard_length += i[1]
    return soft_hard_length
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
    if 0 <= end_pos-1 < len(sequence) and sequence[end_pos-1] == base_homopolymer:
        while end_pos < len(sequence) and sequence[end_pos] == base_homopolymer:
            end_pos += 1
    return start_pos, end_pos

def parsing_homopolymer_error_event(hpm_max_length, shift_length, new_ref, new_seq, homopolymer_aln_event_stat_dict):
#    print(new_ref)
#    print(new_seq)
#    print(len(new_ref) - len(new_seq))
    start_init, end_init = -1, -1
    hpm_map, hpm_mis, hpm_ins, hpm_del = 0, 0, 0, 0
    # homopolymers regular expression pattern
    homopolymer_pattern = "A{2,}|C{2,}|G{2,}|T{2,}" # TODO: use min_homopolymer_size
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
                if homo_ref_length < hpm_max_length:  # [···, shift_length+1=ins1, shift_length+2=ins2, ···]
                    if ins_size < shift_length:
                        homopolymer_aln_event_stat_dict[base_homopolymer][homo_ref_length][
                            shift_length + ins_size + 2] += 1
                        homopolymer_aln_event_stat_dict['S'][homo_ref_length][shift_length + ins_size + 2] += 1
                    else:
                        homopolymer_aln_event_stat_dict[base_homopolymer][homo_ref_length][shift_length * 2 + 1] += 1
                        homopolymer_aln_event_stat_dict['S'][homo_ref_length][shift_length * 2 + 1] += 1
                else:
                    if ins_size < shift_length:
                        homopolymer_aln_event_stat_dict[base_homopolymer][hpm_max_length][
                            shift_length + ins_size + 2] += 1
                        homopolymer_aln_event_stat_dict['S'][hpm_max_length][shift_length + ins_size + 2] += 1
                    else:
                        homopolymer_aln_event_stat_dict[base_homopolymer][hpm_max_length][shift_length * 2 + 1] += 1
                        homopolymer_aln_event_stat_dict['S'][hpm_max_length][shift_length * 2 + 1] += 1
            elif "-" in homopolymer_seq_segment:
                hpm_del += 1
                del_size = homopolymer_seq_segment.count("-")
                if homo_ref_length < hpm_max_length:
                    if del_size < shift_length:
                        homopolymer_aln_event_stat_dict[base_homopolymer][homo_ref_length][shift_length - del_size] += 1
                        homopolymer_aln_event_stat_dict['S'][homo_ref_length][shift_length - del_size] += 1
                    else:
                        homopolymer_aln_event_stat_dict[base_homopolymer][homo_ref_length][
                            shift_length - shift_length] += 1
                        homopolymer_aln_event_stat_dict['S'][homo_ref_length][shift_length - shift_length] += 1
                else:
                    if del_size < shift_length:
                        homopolymer_aln_event_stat_dict[base_homopolymer][hpm_max_length][shift_length - del_size] += 1
                        homopolymer_aln_event_stat_dict['S'][hpm_max_length][shift_length - del_size] += 1
                    else:
                        homopolymer_aln_event_stat_dict[base_homopolymer][hpm_max_length][
                            shift_length - shift_length] += 1
                        homopolymer_aln_event_stat_dict['S'][hpm_max_length][shift_length - shift_length] += 1
            else:
                hpm_mis += 1
                if homo_ref_length < hpm_max_length:
                    homopolymer_aln_event_stat_dict[base_homopolymer][homo_ref_length][shift_length] += 1
                    homopolymer_aln_event_stat_dict['S'][homo_ref_length][shift_length] += 1
                else:
                    homopolymer_aln_event_stat_dict[base_homopolymer][hpm_max_length][shift_length] += 1
                    homopolymer_aln_event_stat_dict['S'][hpm_max_length][shift_length] += 1
        else:
            hpm_map += 1
            if homo_ref_length < hpm_max_length:
                homopolymer_aln_event_stat_dict[base_homopolymer][homo_ref_length][shift_length + 1] += 1
                homopolymer_aln_event_stat_dict['S'][homo_ref_length][shift_length + 1] += 1
            else:
                homopolymer_aln_event_stat_dict[base_homopolymer][hpm_max_length][shift_length + 1] += 1
                homopolymer_aln_event_stat_dict['S'][hpm_max_length][shift_length + 1] += 1
    return hpm_map, hpm_mis, hpm_ins, hpm_del
def parsing_alignment_events(hpm_max_length, hpm_shift_length, raw_ref, raw_seq, cigar_tuples, overall_aln_event_sum_dict, overall_aln_event_stat_dict, query_aln_event_stat_dict, homopolymer_aln_event_stat_dict):
    new_ref, new_seq = '', ''
    qry_map, qry_mis, qry_del, qry_ins = 0, 0, 0, 0
    for i in cigar_tuples:
        if i[0] == 7:  # matched
            new_seq += raw_seq[: i[1]]
            new_ref += raw_seq[: i[1]]
            # print(raw_seq[: i[1]], raw_seq[: i[1]])
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
            qry_ins += i[1]
            if i[1] in overall_aln_event_stat_dict['all_ins_len_dict'].keys():
                overall_aln_event_stat_dict['all_ins_len_dict'][i[1]] += 1
            else:
                overall_aln_event_stat_dict['all_ins_len_dict'][i[1]] = 1
        elif i[0] == 2:  # deletion
            new_seq += i[1] * '-'
            new_ref += raw_ref[:i[1]]
            raw_seq = raw_seq
            raw_ref = raw_ref[i[1]:]
            qry_del += i[1]
            if i[1] in overall_aln_event_stat_dict['all_del_len_dict'].keys():
                overall_aln_event_stat_dict['all_del_len_dict'][i[1]] += 1
            else:
                overall_aln_event_stat_dict['all_del_len_dict'][i[1]] = 1
        elif i[0] == 8:  # mismatch
            new_seq += raw_seq[:i[1]]
            new_ref += raw_ref[:i[1]]
            if i[1] == 1 and (raw_ref[:i[1]] + raw_seq[:i[1]]) in overall_aln_event_stat_dict['all_mis_typ_dict'].keys():
                overall_aln_event_stat_dict['all_mis_typ_dict'][raw_ref[:i[1]] + raw_seq[:i[1]]] += 1
            raw_seq = raw_seq[i[1]:]
            raw_ref = raw_ref[i[1]:]
            qry_mis += i[1]
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
    hpm_map, hpm_mis, hpm_del, hpm_ins = parsing_homopolymer_error_event(hpm_max_length, hpm_shift_length, new_ref, new_seq, homopolymer_aln_event_stat_dict)
    qry_hpm_mis_rate = round(hpm_mis/(qry_map + qry_ins + qry_del + qry_mis), 5)
    qry_hpm_ins_rate = round(hpm_ins/(qry_map + qry_ins + qry_del + qry_mis), 5)
    qry_hpm_del_rate = round(hpm_del/(qry_map + qry_ins + qry_del + qry_mis), 5)
    qry_hpm_idy_rate = round(hpm_map/(qry_map + qry_ins + qry_del + qry_mis), 5)
    qry_hpm_dif_rate = round(1-qry_hpm_idy_rate, 5)
    query_aln_event_stat_dict['qry_hpm_mis_rate'].append(qry_hpm_mis_rate)
    query_aln_event_stat_dict['qry_hpm_ins_rate'].append(qry_hpm_ins_rate)
    query_aln_event_stat_dict['qry_hpm_del_rate'].append(qry_hpm_del_rate)
    query_aln_event_stat_dict['qry_hpm_idy_rate'].append(qry_hpm_idy_rate)
    query_aln_event_stat_dict['qry_hpm_dif_rate'].append(qry_hpm_dif_rate)
    overall_aln_event_sum_dict['hpm_identity'] += hpm_map
    overall_aln_event_sum_dict['hpm_substitution'] += hpm_mis
    overall_aln_event_sum_dict['hpm_expansion'] += hpm_ins
    overall_aln_event_sum_dict['hpm_contraction'] += hpm_del
    # non-homopolymer events statistics
#    qry_non_hpm_mis_rate = round((qry_mis - hpm_mis)/(qry_map + qry_ins + qry_del + qry_mis), 5)
#    qry_non_hpm_ins_rate = round((qry_ins - hpm_ins)/(qry_map + qry_ins + qry_del + qry_mis), 5)
#    qry_non_hpm_del_rate = round((qry_del - hpm_del)/(qry_map + qry_ins + qry_del + qry_mis), 5)
#    hpm_mid_number = hpm_mis + hpm_ins + hpm_del
#    qry_non_hpm_idy_rate = round((qry_map - hpm_mid_number)/(qry_map + qry_ins + qry_del + qry_mis), 5)
#    qry_non_hpm_dif_rate = round(1 - qry_non_hpm_idy_rate, 5)
#    overall_aln_event_sum_dict['non_hpm_substitution'] += qry_mis - hpm_mis
#    overall_aln_event_sum_dict['non_hpm_expansion'] += qry_ins - hpm_ins
#    overall_aln_event_sum_dict['non_hpm_contraction'] += qry_del - hpm_del
#    query_aln_event_stat_dict['qry_non_hpm_mis_rate'].append(qry_non_hpm_mis_rate)
#    query_aln_event_stat_dict['qry_non_hpm_ins_rate'].append(qry_non_hpm_ins_rate)
#    query_aln_event_stat_dict['qry_non_hpm_del_rate'].append(qry_non_hpm_del_rate)
#    query_aln_event_stat_dict['qry_non_hpm_dif_rate'].append(qry_non_hpm_dif_rate)
#    query_aln_event_stat_dict['qry_non_hpm_idy_rate'].append(qry_non_hpm_idy_rate)

def init_reverse_complement():
    TRANSLATION_TABLE = str.maketrans("ACTGactg", "TGACtgac")

    def reverse_complement(sequence):
        """
        >>> reverse_complement("AATC")
        'GATT'
        >>> reverse_complement("CCANT")
        'ANTGG'
        """
        sequence = str(sequence)
        return sequence.translate(TRANSLATION_TABLE)[::-1]

    return reverse_complement


reverse_complement = init_reverse_complement()


def bam_datum_action(args):
    hpm_max_length = args["max_homopolymer_size"]
    hpm_shift_length = args["max_homopolymer_indel_size"]
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
    }
    homopolymer_aln_event_stat_dict = {
        'S': {i: [0] * (hpm_shift_length * 2 + 2) for i in range(2, hpm_max_length + 1, 1)},
        'A': {i: [0] * (hpm_shift_length * 2 + 2) for i in range(2, hpm_max_length + 1, 1)},
        'T': {i: [0] * (hpm_shift_length * 2 + 2) for i in range(2, hpm_max_length + 1, 1)},
        'C': {i: [0] * (hpm_shift_length * 2 + 2) for i in range(2, hpm_max_length + 1, 1)},
        'G': {i: [0] * (hpm_shift_length * 2 + 2) for i in range(2, hpm_max_length + 1, 1)},
    }
    overall_aln_event_sum_dict = {
        'total_reads': 0,
        'mapped_reads': 0,
        'identity': 0,
        'substitution': 0,
        'expansion': 0,
        'contraction': 0,
        'hpm_identity': 0,
        'hpm_substitution': 0,
        'hpm_expansion': 0,
        'hpm_contraction': 0,
    }
    
    # analyzing data
    bam_path = args['bam']
    print(f"Analyzing {bam_path}")
    readprofile = pysam.AlignmentFile(bam_path, "rb")
    for read in readprofile:
        flag = read.flag
        cigar_tuples = read.cigartuples
        
        overall_aln_event_sum_dict['total_reads'] += 1
        if read.is_mapped:
            overall_aln_event_sum_dict['mapped_reads'] += 1
        else:
            continue
        if read.query_alignment_length < 100:
            continue
        # reverse HP compute
        if read.is_reverse:
            raw_seq = read.seq
            raw_ref = read.get_reference_sequence()
            raw_seq, raw_ref = raw_seq.upper(), raw_ref.upper()
            raw_seq, raw_ref = reverse_complement(raw_seq), reverse_complement(raw_ref)
            parsing_alignment_events(hpm_max_length, hpm_shift_length, raw_ref, raw_seq, tuple(reversed(cigar_tuples)), overall_aln_event_sum_dict, overall_aln_event_stat_dict, query_aln_event_stat_dict, homopolymer_aln_event_stat_dict)

        # forward HP compute
        if read.is_forward:
            raw_seq = read.seq
            raw_ref = read.get_reference_sequence()
            raw_seq, raw_ref = raw_seq.upper(), raw_ref.upper()
            parsing_alignment_events(hpm_max_length, hpm_shift_length, raw_ref, raw_seq, cigar_tuples, overall_aln_event_sum_dict, overall_aln_event_stat_dict, query_aln_event_stat_dict, homopolymer_aln_event_stat_dict)

    # structure write
    merge_alignment_dict = {
        'overall_aln_event_sum_dict': overall_aln_event_sum_dict,
        'overall_aln_event_stat_dict': overall_aln_event_stat_dict,
        'query_aln_event_stat_dict': query_aln_event_stat_dict,
        'homopolymer_aln_event_stat_dict': homopolymer_aln_event_stat_dict
    }

    output_path = args['bam_pickle_path']
    with open(output_path, 'wb') as picklefile:
        pickle.dump(merge_alignment_dict, picklefile)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-hpmax", "--homopolymer_max_length", type=int, default=9, required=False, help="observe maxium homopolymer")
    parser.add_argument("-hpmindelmax", "--homopolymer_indel_max", type=int, default=4, required=False, help="homopolymer max indel shift")

    parser.add_argument("-o", "--output_dir", default='./', required=False, help="Output direcotry")
    parser.add_argument("-n", "--sample_name", default='cycads_report', required=False, help="Prefix of output file name")

    # initiated dict
    args = vars(parser.parse_args())
    bam_datum_action(args)


