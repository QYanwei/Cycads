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
    plt.style.use('ggplot')
    sns.distplot()
    sns.distplot(qry_idy_rate, hist=False, kde=False, fit=stats.norm,
                 fit_kws={'color': 'red', 'label': 'overall identity rate', 'linestyle': '-'})
    sns.distplot(qry_hpm_idy_rate, hist=False, kde=False, fit=stats.norm,
                 fit_kws={'color': 'purple', 'label': 'homopolymer identity rate', 'linestyle': '-'})
    sns.distplot(qry_non_hpm_idy_rate, hist=False, kde=False, fit=stats.norm,
                 fit_kws={'color': 'green', 'label': 'non-homopolymer identity rate', 'linestyle': '-'})
    plt.legend()
    plt.savefig('../test/query_events_curve_idy_dif.' + '.distpplot.png')
    # error rate
    qry_dif_rate = np.array(bam_datum_dict['qry_dif_rate'])
    qry_hpm_dif_rate = np.array(bam_datum_dict['qry_hpm_dif_rate'])
    qry_non_hpm_dif_rate = np.array(bam_datum_dict['qry_non_hpm_dif_rate'])
    plt.clf()
    sns.distplot(qry_dif_rate, hist=False, kde=False, fit=stats.norm,
                 fit_kws={'color': 'red', 'label': 'overall error rate', 'linestyle': ':'})
    sns.distplot(qry_hpm_dif_rate, hist=False, kde=False, fit=stats.norm,
                 fit_kws={'color': 'purple', 'label': 'homopolymer error rate', 'linestyle': ':'})
    sns.distplot(qry_non_hpm_dif_rate, hist=False, kde=False, fit=stats.norm,
                 fit_kws={'color': 'green', 'label': 'non-homopolymer error rate', 'linestyle': ':'})
    plt.legend()
    plt.savefig('../test/query_events_curve_dif.' + '.distpplot.png')
    

def plot_homopolymer_frequency(**bam_datum_dict):
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

def overall_length_event_frequency(bam_datum_dict):
    a = 0

bam_json_file_path = '../test/ecoli.bam.json'
with open(bam_json_file_path) as jsonfile:
    bam_datum_dict = json.loads(json.dumps(eval(jsonfile.read()))) # pprint: ' -> json: "
    all_ins_len_dict = bam_datum_dict['all_ins_len_dict']
    all_del_len_dict = bam_datum_dict['all_del_len_dict']
    plot_query_event_rate_overlapping_densities(**bam_datum_dict)
    plot_homopolymer_frequency(**bam_datum_dict)


