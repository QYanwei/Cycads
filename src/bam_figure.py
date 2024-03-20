#!/usr/bin/env python3
# coding: utf-8
import re, os, sys, time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import json
import warnings
warnings.filterwarnings("ignore", "is_categorical_dtype")
warnings.filterwarnings("ignore", "use_inf_as_na")

plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['figure.dpi'] = 300
palette = dict(A="tab:red", T="tab:green", C="tab:blue", G="tab:purple", S='tab:black')

def plot_query_event_rate_overlapping_densities(**bam_datum_dict):
    # identity rate
    qry_idy_rate = np.array(bam_datum_dict['qry_idy_rate'])
    qry_hpm_idy_rate = np.array(bam_datum_dict['qry_hpm_idy_rate'])
    qry_non_hpm_idy_rate = np.array(bam_datum_dict['qry_non_hpm_idy_rate'])
    qry_idy_df = pd.DataFrame([qry_idy_rate,qry_hpm_idy_rate, qry_non_hpm_idy_rate])
    qry_idy_df.index = ['overall identity rate', 'homopolymer identity rate', 'non-homopolymer identity rate']
    plt.clf()
    sns.displot(data=qry_idy_df.T, kind="kde")
    plt.savefig('../test/query_events_curve_idy' + '.dispplot.png')
    # error rate
    qry_dif_rate = np.array(bam_datum_dict['qry_dif_rate'])
    qry_hpm_dif_rate = np.array(bam_datum_dict['qry_hpm_dif_rate'])
    qry_non_hpm_dif_rate = np.array(bam_datum_dict['qry_non_hpm_dif_rate'])
    qry_dif_df = pd.DataFrame([qry_dif_rate,qry_hpm_dif_rate, qry_non_hpm_dif_rate])
    qry_dif_df.index = ['overall error rate', 'homopolymer error rate', 'non-homopolymer error rate']
    plt.clf()
    sns.displot(data=qry_dif_df.T, kind="kde")
    plt.savefig('../test/query_events_curve_dif' + '.dispplot.png')

def plot_insertion_deletion_frequency(**bam_datum_dict):
    indel_len_max = 10
    indel_range_dict = {'Insertion': {}, 'Deletion':{} }
    ins_sum = sum(bam_datum_dict['all_ins_len_dict'].values())
    del_sum = sum(bam_datum_dict['all_del_len_dict'].values())
    for i in bam_datum_dict['all_ins_len_dict'].keys():
        if int(i) < indel_len_max:
            indel_range_dict['Insertion'][int(i)] = bam_datum_dict['all_ins_len_dict'][i] * 100 / ins_sum
        elif int(i) >= indel_len_max:
            if indel_len_max in indel_range_dict['Insertion'].keys():
                indel_range_dict['Insertion'][indel_len_max] += bam_datum_dict['all_ins_len_dict'][i]
            else:
                indel_range_dict['Insertion'][indel_len_max] = bam_datum_dict['all_ins_len_dict'][i]
        else:
            pass
    for i in bam_datum_dict['all_del_len_dict'].keys():
        if int(i) < indel_len_max:
            indel_range_dict['Deletion'][int(i)] = bam_datum_dict['all_del_len_dict'][i] * 100 / del_sum
        elif int(i) >= indel_len_max:
            if indel_len_max in indel_range_dict['Deletion'].keys():
                indel_range_dict['Deletion'][indel_len_max] += bam_datum_dict['all_del_len_dict'][i]
            else:
                indel_range_dict['Deletion'][indel_len_max] = bam_datum_dict['all_del_len_dict'][i]
        else:
            pass
    indel_range_dict['Insertion'][indel_len_max] = indel_range_dict['Insertion'][indel_len_max] * 100 / ins_sum
    indel_range_dict['Deletion'][indel_len_max] = indel_range_dict['Deletion'][indel_len_max] * 100 / del_sum
    plt.clf()
    Ins_dataframe = pd.DataFrame([ indel_range_dict['Insertion'] ])
    Ins_dataframe.index = ['Insertion']
    sns.barplot(data=Ins_dataframe)
    plt.savefig('../test/query_insertion_frequency' + '.barplot.png')
    plt.clf()
    Del_dataframe = pd.DataFrame([ indel_range_dict['Deletion'] ])
    Del_dataframe.index = ['Deletion']
    sns.barplot(data=Del_dataframe)
    plt.savefig('../test/query_deletion_frequency' + '.barplot.png')

def plot_overall_homopolymer_length_event_frequency(homopolymer_aln_event_stat_dict):
    qry_hpm_event = np.array([homopolymer_aln_event_stat_dict['S'][i] for i in homopolymer_aln_event_stat_dict['S'].keys()])
    qry_hpm_event2= []
    for i in qry_hpm_event:
        lt = []
        p = 0
        for j in i /sum(i):
            p += j
            lt.append(p)
        qry_hpm_event2.append(lt)
    homopolymer_dataframe = pd.DataFrame(qry_hpm_event2)
    homopolymer_dataframe.index = range(2,len(qry_hpm_event2)+2,1)
    homopolymer_dataframe.columns = ["contraction 4bp+", "contraction 3bp", "contraction 2bp", "contraction 1bp",
                                     "correct segment",  "mismatch segment",
                                     "expansion 1bp", "expansion 2bp", "expansion 3bp", "expansion 4bp+"]
    plt.clf()
    sns.lineplot(data=homopolymer_dataframe)
    plt.savefig('../test/query_homopolymer_length_event' + '.lineplot.png')

def plot_overall_alignment_frequency(**overall_aln_event_sum_dict):
    all_event = overall_aln_event_sum_dict['identity'] + overall_aln_event_sum_dict['substitution'] + overall_aln_event_sum_dict['contraction'] + overall_aln_event_sum_dict['expansion']
    all_dif_event = overall_aln_event_sum_dict['substitution'] + overall_aln_event_sum_dict['contraction'] + overall_aln_event_sum_dict['expansion']
    all_dif = all_dif_event / all_event *100
    all_mis = overall_aln_event_sum_dict['substitution'] / all_event *100
    all_del = overall_aln_event_sum_dict['contraction'] / all_event *100
    all_ins = overall_aln_event_sum_dict['expansion'] / all_event *100
    non_hpm_dif_event = overall_aln_event_sum_dict['non_hpm_substitution'] + overall_aln_event_sum_dict['non_hpm_contraction'] + overall_aln_event_sum_dict['non_hpm_expansion']
    non_hpm_dif = non_hpm_dif_event / all_event *100
    non_hpm_mis = overall_aln_event_sum_dict['non_hpm_substitution'] / all_event *100
    non_hpm_del = overall_aln_event_sum_dict['non_hpm_contraction'] / all_event *100
    non_hpm_ins = overall_aln_event_sum_dict['non_hpm_expansion'] / all_event *100
    plt.clf()
    all_aln_rates = [all_dif, all_mis, all_del, all_ins]
    non_hpm_aln_rates = [non_hpm_dif, non_hpm_mis, non_hpm_del, non_hpm_ins]
    labels = ['Error rate', 'Mismatch', 'Deletion', 'Insertion']
    x = np.arange(len(labels))  # the label locations
    width = 0.35  # the width of the bars
    fig, ax = plt.subplots()
    ax.bar(x - width/2, all_aln_rates, width, label='Total error rate')
    ax.bar(x + width/2, non_hpm_aln_rates, width, label='Non-homopolymer error rate')
    ax.set_ylabel('Percentage(%)')
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.legend()
    plt.savefig('../test/query_all_error_item' + '.barplot.png')

bam_json_file_path = '../test/ecoli.bam.json'
with open(bam_json_file_path) as jsonfile:
    bam_datum_dict = json.loads(json.dumps(eval(jsonfile.read()))) # pprint: ' -> json: "
    plot_query_event_rate_overlapping_densities(**bam_datum_dict['query_aln_event_stat_dict'])
    plot_overall_homopolymer_length_event_frequency(bam_datum_dict['homopolymer_aln_event_stat_dict'])
    plot_insertion_deletion_frequency(**bam_datum_dict['overall_aln_event_stat_dict'])
    plot_overall_alignment_frequency(**bam_datum_dict['overall_aln_event_sum_dict'])
