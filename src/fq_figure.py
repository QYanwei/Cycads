#!/usr/bin/env python3
# coding: utf-8
import re, os, sys, time
import numpy as np
import pandas as pd
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
    plt.clf()
    sns.displot(LengthFrame, x="read_length")
    sns.displot(LengthFrame, x="read_length", kind='ecdf')
    plt.savefig('../test/read_length_cumulative' + '.barplot.png')

def plot_length_quality_cross_scatter(**seq_qual_dict):
    LenQualFrame = pd.DataFrame( np.array([ seq_qual_dict['LEN'], seq_qual_dict['QUAL1'] ] ).T, columns = [ 'read_length', 'read_averageQ' ])
    plt.clf()
    sns.set_style("ticks")
    ax = sns.scatterplot(data=LenQualFrame, x="read_length", y="read_averageQ", linewidth=0, size=1, alpha=0.5, legend=None)
    bottom, top = ax.get_ylim()
    ax.set_ylim(bottom + 0.2, top - 0.2)
    # sns.scatterplot(data=LenQualFrame, x="read_length", y="read_averageQ", linewidth=0, size=2, alpha=0.3, legend=None)
    plt.savefig('../test/read_length_quality_cross' + '.scatterplot.png')

def plot_homopolymer_frequency(**homopolymer_dict):
    homopolymer_len_max = 20
    homopolymer_range_dict = {'A': {}, 'G': {}, 'C': {}, 'T': {} }
    for b in homopolymer_dict.keys():
        for i in homopolymer_dict[b].keys():
            if int(i) < homopolymer_len_max:
                homopolymer_range_dict[b][int(i)] = homopolymer_dict[b][i]
            elif int(i) >= homopolymer_len_max:
                homopolymer_range_dict[b][str(homopolymer_len_max)] += homopolymer_dict[b][i]
            else:
                pass
    plt.clf()
    A_dataframe = pd.DataFrame([homopolymer_range_dict['A']])
    A_dataframe.index = ['A']
    T_dataframe = pd.DataFrame([homopolymer_range_dict['T']])
    T_dataframe.index = ['T']
    G_dataframe = pd.DataFrame([homopolymer_range_dict['G']])
    G_dataframe.index = ['G']
    C_dataframe = pd.DataFrame([homopolymer_range_dict['C']])
    C_dataframe.index = ['C']
    homopolymerRangeFrame = pd.concat([A_dataframe, T_dataframe, G_dataframe, C_dataframe], axis=0)
    sns.lineplot(data=homopolymerRangeFrame.T, dashes=False)
    plt.savefig('../test/read_homopolymer_frequency' + '.lineplot.png')

# def plot_kmer_spectrum_heatmap(kmer_list, kmer_count):

def plot_ends_base_content_curve(**endBaseQual_dict):
    # head base content
    plt.clf()
    dataframe = pd.DataFrame([endBaseQual_dict['HeadBaseContent_dict']])
    A_ratio = np.array(dataframe['A'][0])/dataframe['S'][0]
    T_ratio = np.array(dataframe['T'][0])/dataframe['S'][0]
    G_ratio = np.array(dataframe['G'][0])/dataframe['S'][0]
    C_ratio = np.array(dataframe['C'][0])/dataframe['S'][0]
    numpyarray = np.array([A_ratio, T_ratio, G_ratio, C_ratio]).T
    PercentQualityFrame = pd.DataFrame(data=numpyarray, columns=['A', 'T', 'C', 'G'])
    sns.lineplot(data=PercentQualityFrame, dashes=False)
    plt.savefig('../test/read_head_base_content' + '.lineplot.png')
    # tail base content
    plt.clf()
    dataframe = pd.DataFrame([endBaseQual_dict['TailBaseContent_dict']])
    A_ratio = np.array(dataframe['A'][0])/dataframe['S'][0]
    T_ratio = np.array(dataframe['T'][0])/dataframe['S'][0]
    G_ratio = np.array(dataframe['G'][0])/dataframe['S'][0]
    C_ratio = np.array(dataframe['C'][0])/dataframe['S'][0]
    numpyarray = np.array([A_ratio, T_ratio, G_ratio, C_ratio]).T
    PercentQualityFrame = pd.DataFrame(data=numpyarray, columns=['A', 'T', 'C', 'G'])
    sns.lineplot(data=PercentQualityFrame, dashes=False).invert_xaxis()
    plt.savefig('../test/read_tail_base_content' + '.lineplot.png')

def plot_ends_base_quality_curve(**endBaseQual_dict):
    # head base content
    plt.clf()
    dataframe = pd.DataFrame([endBaseQual_dict['HeadQualContent_dict']])
    basetotal = dataframe['S'][0]
    A_quality = np.array(dataframe['A'][0])/basetotal['A']
    T_quality = np.array(dataframe['T'][0])/basetotal['T']
    G_quality = np.array(dataframe['G'][0])/basetotal['G']
    C_quality = np.array(dataframe['C'][0])/basetotal['C']
    Mean_quality = np.array( np.array(dataframe['A'][0]) + np.array(dataframe['T'][0]) + np.array(dataframe['G'][0]) + np.array(dataframe['C'][0]))/sum(basetotal.values())
    numpyarray = np.array([A_quality, T_quality, G_quality, C_quality, Mean_quality]).T
    PercentQualityFrame = pd.DataFrame(data=numpyarray, columns=['A', 'T', 'C', 'G', 'Mean'])
    sns.lineplot(data=PercentQualityFrame, dashes=False)
    plt.savefig('../test/read_head_base_quality' + '.lineplot.png')
    # tail base content
    plt.clf()
    dataframe = pd.DataFrame([endBaseQual_dict['TailQualContent_dict']])
    basetotal = dataframe['S'][0]
    A_quality = np.array(dataframe['A'][0])/basetotal['A']
    T_quality = np.array(dataframe['T'][0])/basetotal['T']
    G_quality = np.array(dataframe['G'][0])/basetotal['G']
    C_quality = np.array(dataframe['C'][0])/basetotal['C']
    Mean_quality = np.array( np.array(dataframe['A'][0]) + np.array(dataframe['T'][0]) + np.array(dataframe['G'][0]) + np.array(dataframe['C'][0]))/sum(basetotal.values())
    numpyarray = np.array([A_quality, T_quality, G_quality, C_quality, Mean_quality]).T
    PercentQualityFrame = pd.DataFrame(data=numpyarray, columns=['A', 'T', 'C', 'G', 'Mean'])
    sns.lineplot(data=PercentQualityFrame, dashes=False)
    plt.savefig('../test/read_tail_base_quality' + '.lineplot.png')

def plot_read_percent_qualtiy_curve(**allBaseQual_dict):
    plt.clf()
    percent_pos_avg_quality = np.array(allBaseQual_dict['PercentBaseQual_dict']['Q'])/ allBaseQual_dict['PercentBaseQual_dict']['S']
    numpyarray = np.array([percent_pos_avg_quality, np.arange(len(percent_pos_avg_quality))+1]).T
    PercentQualityFrame = pd.DataFrame(data=numpyarray, columns=['average_quality', 'relative_position'])
    sns.lineplot(data=PercentQualityFrame, x="relative_position", y="average_quality")
    plt.savefig('../test/read_relative_position_avg_qual' + '.lineplot.png')

fastq_json_file_path = '../test/ecoli.seq.json'

with open(fastq_json_file_path) as jsonfile:
    fq_datum_dict = json.loads(json.dumps(eval(jsonfile.read()))) # pprint: ' -> json: "
    seq_qual_dict = fq_datum_dict['seq_qual_dict']
    homopolymer_dict = fq_datum_dict['homopolymer_dict']
    endBaseQual_dict = fq_datum_dict['endBaseQual_dict']
    allBaseQual_dict = fq_datum_dict['allBaseQual_dict']
    # print(fq_datum_dict['seq_qual_dict'].keys())
    # print(fq_datum_dict['homopolymer_dict'].keys())
    # print(fq_datum_dict['endBaseQual_dict'].keys())
    # print(fq_datum_dict['allBaseQual_dict'].keys())
    plot_length_Nx_average_bar(**seq_qual_dict)
    plot_gc_content_frequency_distribution(**seq_qual_dict)
    plot_read_quality_frequency_distritution(**seq_qual_dict)
    plot_read_length_frequency_distribution(**seq_qual_dict)
    plot_read_length_cumulative_distribution(**seq_qual_dict)
    plot_length_quality_cross_scatter(**seq_qual_dict)
    plot_read_percent_qualtiy_curve(**allBaseQual_dict)
    plot_ends_base_content_curve(**endBaseQual_dict)
    plot_ends_base_quality_curve(**endBaseQual_dict)
    plot_homopolymer_frequency(**homopolymer_dict)
