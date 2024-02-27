#!/usr/bin/env python3
# coding: utf-8
import re, os, sys, time
import numpy as np
import pandas as pd
import scipy
import matplotlib.pyplot as plt
import seaborn as sns
import json

import warnings
warnings.filterwarnings("ignore", "is_categorical_dtype")
warnings.filterwarnings("ignore", "use_inf_as_na")

plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['figure.dpi'] = 300
palette = dict(A="tab:red", T="tab:green", C="tab:blue", G="tab:purple", S='tab:black')

def plot_length_Nx_average_bar(**seq_qual_dict):
    def Nx_reads_length(list_read_length):
        # Calculate the total length of the sequences
        total_length = sum(list_read_length)
        list_read_length.sort(reverse=True)
        # Calculate the N10, N20, N30, N40, N50, N60, N70, N80 and N90 values
        n_values = []
        for i in range(10, 100, 10):
            n_value = total_length * i / 100
            n_sum = 0
            for length in list_read_length:
                n_sum += length
                if n_sum >= n_value:
                    n_values.append(length)
                    break
        return n_values
    average_length = round(sum(seq_qual_dict['LEN']) / len(seq_qual_dict['LEN']), 3)
    LengthFrame = pd.DataFrame( np.array( Nx_reads_length( seq_qual_dict['LEN'])  + [ average_length ] ).astype(int).T, columns=['read_length'])
    LengthFrame.loc[: , 'nx_and_average'] = [ 'N10', 'N20', 'N30', 'N40', 'N50', 'N60', 'N70', 'N80', 'N90', 'Average']
    plt.clf()
    g = sns.barplot(LengthFrame, x='nx_and_average', y='read_length')
    g.bar_label(g.containers[0], fontsize=10)
    plt.savefig('../test/read_length_biostat' + '.barplot.png')
def plot_gc_content_frequency_distribution(**seq_qual_dict):
    gcFrame = pd.DataFrame( np.array( seq_qual_dict['GC'] ).astype(float), columns=['GC'])
    plt.clf()
    sns.displot(gcFrame, color="steelblue", bins=50, stat="percent")
    plt.savefig('../test/read_gc_histplot' + '.barplot.png')
def plot_read_quality_frequency_distritution(**seq_qual_dict):
    QualityFrame = pd.DataFrame( np.array( seq_qual_dict['QUAL1'] ).astype(int).T, columns=['read_quality'])
    plt.clf()
    sns.displot(data=QualityFrame, x="read_quality", color="steelblue", stat="percent")
    plt.savefig('../test/read_quality_histplot' + '.barplot.png')
def plot_read_length_frequency_distribution(**seq_qual_dict):
    LengthFrame = pd.DataFrame( np.array( seq_qual_dict['LEN'] ).astype(int).T, columns=['read_length'])
    plt.clf()
    sns.histplot(data=LengthFrame, x="read_length", color="steelblue", stat="percent")
    # sns.histplot(data = LengthFrame, x="read_length", stat="percent", log_scale=True)
    plt.savefig('../test/read_length_histplot_nolog' + '.barplot.png')
def plot_read_length_cumulative_distribution(**seq_qual_dict):
    LengthFrame = pd.DataFrame(np.array(seq_qual_dict['LEN']).astype(int).T, columns=['read_length'])
    sns.displot(LengthFrame, x="read_length")
    sns.displot(LengthFrame, x="read_length", kind='ecdf')
    plt.savefig('../test/read_length_cumulative' + '.barplot.png')

def plot_length_quality_cross_scatter(**seq_qual_dict):
    LenQualFrame = pd.DataFrame( np.array([ seq_qual_dict['LEN'], seq_qual_dict['QUAL1'] ] ).T, columns = [ 'read_length', 'read_averageQ' ])
    plt.clf()
    sns.set_style("ticks")
    ax = sns.scatterplot(data=LenQualFrame, x="read_length", y="read_averageQ", linewidth=0, size=2, alpha=0.3, legend=None)
    bottom, top = ax.get_ylim()
    ax.set_ylim(bottom + 0.2, top - 0.2)
    # sns.scatterplot(data=LenQualFrame, x="read_length", y="read_averageQ", linewidth=0, size=2, alpha=0.3, legend=None)
    plt.savefig('../test/read_length_quality_cross' + '.scatterplot.png')

# def plot_kmer_spectrum_heatmap(kmer_list, kmer_count):
#
# def plot_homopolymer_frequency(homopolymer_dict):
#
def plot_ends_base_content_curve(end_base_content_dict):
    # plot head base content
    plt.clf()
    fig, ax = plt.subplots(figsize=(8, 6))
    atgc_head_df = pd.DataFrame(atgc_head_dt).astype(float)
    if int(lens) <= 500:
        ax_head = atgc_head_df.plot(linewidth=0.5)
    elif int(lens) <= 5000:
        ax_head = atgc_head_df.plot(linewidth=0.2)
    else:
        ax_head = atgc_head_df.plot(linewidth=0.1)
    #    plt.ylim(0, 1)
    plt.title(' '.join([tech, 'head', lens + 'bp']))
    plt.ylabel('Frequency of base content')
    plt.xlabel('Position in the head of read')
    plt.savefig('_'.join([tech, 'head', lens]) + '.png')
    
    # plot tail base content
    plt.clf()
    fig, ax = plt.subplots(figsize=(8, 6))
    atgc_tail_df = pd.DataFrame(atgc_tail_dt).astype(float)
    if int(lens) <= 500:
        ax_tail = atgc_tail_df.plot(linewidth=0.5)
    elif int(lens) <= 5000:
        ax_tail = atgc_tail_df.plot(linewidth=0.2)
    else:
        ax_tail = atgc_tail_df.plot(linewidth=0.1)
    #    plt.ylim(0, 1)
    plt.title(' '.join([tech, 'tail', lens + 'bp']))
    plt.ylabel('Frequency of base content')
    plt.xlabel('Position in the tail of read')
    plt.gca().invert_xaxis()
    plt.savefig('_'.join([tech, 'tail', lens]) + '.png')
def plot_ends_base_quality_curve(end_base_quality_dict):
    # plot head base content
    plt.clf()
    fig, ax = plt.subplots(figsize=(8, 6))
    atgc_head_df = pd.DataFrame(atgc_head_dt).astype(float)
    if int(lens) <= 500:
        ax_head = atgc_head_df.plot(linewidth=0.5)
    elif int(lens) <= 5000:
        ax_head = atgc_head_df.plot(linewidth=0.2)
    else:
        ax_head = atgc_head_df.plot(linewidth=0.1)
    #    plt.ylim(0, 1)
    plt.title(' '.join([tech, 'head', lens + 'bp']))
    plt.ylabel('Frequency of base content')
    plt.xlabel('Position in the head of read')
    plt.savefig('_'.join([tech, 'head', lens]) + '.png')
    
    # plot tail base content
    plt.clf()
    fig, ax = plt.subplots(figsize=(8, 6))
    atgc_tail_df = pd.DataFrame(atgc_tail_dt).astype(float)
    if int(lens) <= 500:
        ax_tail = atgc_tail_df.plot(linewidth=0.5)
    elif int(lens) <= 5000:
        ax_tail = atgc_tail_df.plot(linewidth=0.2)
    else:
        ax_tail = atgc_tail_df.plot(linewidth=0.1)
    #    plt.ylim(0, 1)
    plt.title(' '.join([tech, 'tail', lens + 'bp']))
    plt.ylabel('Frequency of base content')
    plt.xlabel('Position in the tail of read')
    plt.gca().invert_xaxis()
    plt.savefig('_'.join([tech, 'tail', lens]) + '.png')
# def plot_read_percent_qualtiy_curve(all_base_percent_average_quality_dict):


fastq_json_file_path = '../test/ecoli.seq.json'

with open(fastq_json_file_path) as jsonfile:
    fq_datum_dict = json.loads(json.dumps(eval(jsonfile.read()))) # pprint: ' -> json: "
    seq_qual_dict = fq_datum_dict['seq_qual_dict']
    homopolymer_dict = fq_datum_dict['homopolymer_dict']
    endBaseQual_dict = fq_datum_dict['endBaseQual_dict']
    allBaseQual_dict = fq_datum_dict['allBaseQual_dict']
    print(fq_datum_dict['seq_qual_dict'].keys())
    print(fq_datum_dict['homopolymer_dict'].keys())
    print(fq_datum_dict['endBaseQual_dict'].keys())
    print(fq_datum_dict['allBaseQual_dict'].keys())
    plot_length_Nx_average_bar(**seq_qual_dict)
    plot_gc_content_frequency_distribution(**seq_qual_dict)
    plot_read_quality_frequency_distritution(**seq_qual_dict)
    plot_read_length_frequency_distribution(**seq_qual_dict)
    plot_read_length_cumulative_distribution(**seq_qual_dict)
    plot_length_quality_cross_scatter(**seq_qual_dict)





