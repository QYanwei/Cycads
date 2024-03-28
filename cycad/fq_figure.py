#!/usr/bin/env python3
# coding: utf-8
import re, os, sys, time
from math import floor, ceil
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import seaborn as sns
import json
import warnings

warnings.filterwarnings("ignore", "is_categorical_dtype")
warnings.filterwarnings("ignore", "use_inf_as_na")

plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['figure.dpi'] = 300
palette = dict(A="tab:red", T="tab:green", C="tab:blue", G="tab:purple", S='tab:black')
figsize = (5, 4)
hist_kw = dict(facecolor='tab:blue', edgecolor='k', linewidth=0.5)
grid_kw = dict(color='k', alpha=0.1)


def post_process_ax(ax):
    ax.spines[['right', 'top']].set_visible(False)

def plot_length_Nx_average_bar(seq_qual_dict):
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
    average_length = sum(seq_qual_dict['LEN']) / len(seq_qual_dict['LEN'])
    LengthFrame = pd.DataFrame( np.array( Nx_reads_length( seq_qual_dict['LEN']) ).astype(int).T, columns=['read_length'])
    LengthFrame.loc[: , 'nx'] = [ 'N10', 'N20', 'N30', 'N40', 'N50', 'N60', 'N70', 'N80', 'N90']
    
    fig, ax = plt.subplots(figsize=figsize)
    g = sns.barplot(LengthFrame, ax=ax, x='nx', y='read_length', color='steelblue')
    ax.axhline(average_length, color='k', ls='--', label='Mean read length')
    ax.legend(loc='upper right')
    ax.set_xlabel("")
    ax.set_ylabel("Read length (bp)")
    ax.grid(axis='y', **grid_kw)
    post_process_ax(ax)
    return fig

def plot_gc_content_frequency_distribution(seq_qual_dict):
    gc_percentage = np.array( seq_qual_dict['GC'] ).astype(float) * 100
    fig, ax = plt.subplots(figsize=figsize)
    ax.hist(gc_percentage, bins=np.arange(0, 100, 2), density=True, **hist_kw)
    ax.xaxis.set_major_formatter(mtick.PercentFormatter())
    ax.set_xlim(0, 100)
    ax.set_xlabel("Read GC content")
    ax.set_xticks(list(range(0, 101, 10)))
    ax.set_ylabel("Density")
    ax.set_yticks([])
    post_process_ax(ax)
    return fig

def plot_read_quality_frequency_distritution(seq_qual_dict):
    quality_values = np.array( seq_qual_dict['QUAL1'] )
    fig, ax = plt.subplots(figsize=figsize)
    xmin = floor(quality_values.min())
    xmax = ceil(quality_values.max())
    bins = np.arange(xmin, xmax + 1, 0.5)
    ax.hist(quality_values, bins=bins, density=True, **hist_kw)
    ax.set_xlabel("Mean per-read base quality")
    ax.set_xlim(xmin, xmax)
    ax.set_ylabel("Density")
    ax.set_yticks([])
    post_process_ax(ax)
    return fig

def plot_read_length_frequency_distribution(seq_qual_dict):
    lengths = np.array( seq_qual_dict['LEN'] )

    fig, ax = plt.subplots(figsize=figsize)
    ax.hist(lengths / 1000, bins=50, **hist_kw)
    ax.set_xlim(left=0)
    ax.set_xlabel("Read length (kb)")
    ax.set_ylabel("Density")
    ax.set_yticks([])
    post_process_ax(ax)
    return fig

def plot_read_length_cumulative_distribution(seq_qual_dict):
    lengths = np.array( seq_qual_dict['LEN'] )
    total_bases = lengths.sum()
    xs = np.arange(lengths.min(), lengths.max() + 100, 100) / 1000
    ys = [lengths[lengths >= x * 1000].sum() * 100 / total_bases for x in xs]

    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(xs, ys)
    ax.set_ylim(0, 100)
    ax.set_xlim(left=0)
    ax.set_xlabel("Min. read length (kb)")
    ax.set_ylabel("Cumulative fraction of bases")
    ax.yaxis.set_major_formatter(mtick.PercentFormatter())
    ax.grid(axis='both', **grid_kw)
    post_process_ax(ax)
    return fig

def plot_length_quality_cross_scatter(seq_qual_dict, output):
    LenQualFrame = pd.DataFrame( np.array([ seq_qual_dict['LEN'], seq_qual_dict['QUAL1'] ] ).T, columns = [ 'read_length', 'read_averageQ' ])
    plt.clf()
    sns.set_style("ticks")
    ax = sns.scatterplot(data=LenQualFrame, x="read_length", y="read_averageQ", linewidth=0, size=1, alpha=0.5, legend=None)
    bottom, top = ax.get_ylim()
    ax.set_ylim(bottom + 0.2, top - 0.2)
    # sns.scatterplot(data=LenQualFrame, x="read_length", y="read_averageQ", linewidth=0, size=2, alpha=0.3, legend=None)
    plt.savefig(output+'/report_html/read_length_quality_cross' + '.scatterplot.png')

def plot_homopolymer_frequency(homopolymer_dict, output):
    homopolymer_len_max = 10
    homopolymer_range_dict = {'A': {}, 'G': {}, 'C': {}, 'T': {} }
    for b in homopolymer_dict.keys():
        for i in homopolymer_dict[b].keys():
            if int(i) <= homopolymer_len_max:
                homopolymer_range_dict[b][int(i)] = homopolymer_dict[b][i]
            elif int(i) > homopolymer_len_max:
                if i not in homopolymer_range_dict[b].keys():
                    homopolymer_range_dict[b][homopolymer_len_max] = homopolymer_dict[b][i]
                else:
                    homopolymer_range_dict[b][homopolymer_len_max] += homopolymer_dict[b][i]
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
    plt.savefig(output+'/report_html/read_homopolymer_frequency' + '.lineplot.png')

# def plot_kmer_spectrum_heatmap(kmer_list, kmer_count):

def plot_ends_base_content_curve(endBaseQual_dict, output):
    # head base content
    plt.clf()
    dataframe = pd.DataFrame([endBaseQual_dict['HeadBaseContent_dict']])
    A_ratio = np.array(dataframe['A'][0])/dataframe['S'][0]
    T_ratio = np.array(dataframe['T'][0])/dataframe['S'][0]
    G_ratio = np.array(dataframe['G'][0])/dataframe['S'][0]
    C_ratio = np.array(dataframe['C'][0])/dataframe['S'][0]
    GC_ratio = (np.array(dataframe['G'][0]) + np.array(dataframe['C'][0]))/dataframe['S'][0]
    numpyarray = np.array([A_ratio, T_ratio, G_ratio, C_ratio, GC_ratio]).T
    PercentQualityFrame = pd.DataFrame(data=numpyarray, columns=['A', 'T', 'C', 'G', 'GC'])
    sns.lineplot(data=PercentQualityFrame, dashes=False)
    plt.savefig(output+'/report_html/read_head_base_content' + '.lineplot.png')
    # tail base content
    plt.clf()
    dataframe = pd.DataFrame([endBaseQual_dict['TailBaseContent_dict']])
    A_ratio = np.array(dataframe['A'][0])/dataframe['S'][0]
    T_ratio = np.array(dataframe['T'][0])/dataframe['S'][0]
    G_ratio = np.array(dataframe['G'][0])/dataframe['S'][0]
    C_ratio = np.array(dataframe['C'][0])/dataframe['S'][0]
    GC_ratio = (np.array(dataframe['G'][0]) + np.array(dataframe['C'][0]))/dataframe['S'][0]
    numpyarray = np.array([A_ratio, T_ratio, G_ratio, C_ratio, GC_ratio]).T
    PercentQualityFrame = pd.DataFrame(data=numpyarray, columns=['A', 'T', 'C', 'G', 'GC'])
    sns.lineplot(data=PercentQualityFrame, dashes=False).invert_xaxis()
    plt.savefig(output+'/report_html/read_tail_base_content' + '.lineplot.png')

def plot_ends_base_quality_curve(endBaseQual_dict, output):
    # head base content
    plt.clf()
    dataframe = pd.DataFrame([endBaseQual_dict['HeadQualContent_dict']])
    basetotal = dataframe['S'][0]
    A_quality = np.array(dataframe['A'][0])/np.array(basetotal['A'])
    T_quality = np.array(dataframe['T'][0])/np.array(basetotal['T'])
    G_quality = np.array(dataframe['G'][0])/np.array(basetotal['G'])
    C_quality = np.array(dataframe['C'][0])/np.array(basetotal['C'])
    base_total = np.array(basetotal['A']) + np.array(basetotal['T']) + np.array(basetotal['G']) + np.array(basetotal['C'])
    Mean_quality = np.array( np.array(dataframe['A'][0]) + np.array(dataframe['T'][0]) + np.array(dataframe['G'][0]) + np.array(dataframe['C'][0]))/np.array(base_total)
    numpyarray = np.array([A_quality, T_quality, G_quality, C_quality, Mean_quality]).T
    PercentQualityFrame = pd.DataFrame(data=numpyarray, columns=['A', 'T', 'C', 'G', 'Mean'])
    sns.lineplot(data=PercentQualityFrame, dashes=False)
    plt.savefig(output+'/report_html/read_head_base_quality' + '.lineplot.png')
    # tail base content
    plt.clf()
    dataframe = pd.DataFrame([endBaseQual_dict['TailQualContent_dict']])
    basetotal = dataframe['S'][0]
    A_quality = np.array(dataframe['A'][0])/np.array(basetotal['A'])
    T_quality = np.array(dataframe['T'][0])/np.array(basetotal['T'])
    G_quality = np.array(dataframe['G'][0])/np.array(basetotal['G'])
    C_quality = np.array(dataframe['C'][0])/np.array(basetotal['C'])
    base_total = np.array(basetotal['A']) + np.array(basetotal['T']) + np.array(basetotal['G']) + np.array(basetotal['C'])
    Mean_quality = np.array( np.array(dataframe['A'][0]) + np.array(dataframe['T'][0]) + np.array(dataframe['G'][0]) + np.array(dataframe['C'][0]))/np.array(base_total)
    numpyarray = np.array([A_quality, T_quality, G_quality, C_quality, Mean_quality]).T
    PercentQualityFrame = pd.DataFrame(data=numpyarray, columns=['A', 'T', 'C', 'G', 'Mean'])
    sns.lineplot(data=PercentQualityFrame, dashes=False)
    plt.savefig(output+'/report_html/read_tail_base_quality' + '.lineplot.png')

def plot_read_percent_qualtiy_curve(allBaseQual_dict, output):
    plt.clf()
    percent_pos_avg_quality = np.array(allBaseQual_dict['PercentBaseQual_dict']['Q'])/ allBaseQual_dict['PercentBaseQual_dict']['S']
    numpyarray = np.array([percent_pos_avg_quality, np.arange(len(percent_pos_avg_quality))+1]).T
    PercentQualityFrame = pd.DataFrame(data=numpyarray, columns=['average_quality', 'relative_position'])
    sns.lineplot(data=PercentQualityFrame, x="relative_position", y="average_quality")
    plt.savefig(output+'/report_html/read_relative_position_avg_qual' + '.lineplot.png')

def fq_figure_action(args):
    output = args["sample_name"]
    if os.path.exists(output + "/" + output +"_fq.json"):
        fastq_json_file_path = args["sample_name"]+"/" + args["sample_name"]+"_fq.json"
        with open(fastq_json_file_path) as jsonfile:
            fq_datum_dict = json.loads(json.dumps(eval(jsonfile.read()))) # pprint: ' -> json: "
            seq_qual_dict = fq_datum_dict['seq_qual_dict']
            homopolymer_dict = fq_datum_dict['homopolymer_dict']
            endBaseQual_dict = fq_datum_dict['endBaseQual_dict']
            allBaseQual_dict = fq_datum_dict['allBaseQual_dict']

        fig = plot_length_Nx_average_bar(seq_qual_dict)
        fig.savefig(output+'/report_html/read_length_biostat' + '.barplot.png')

        fig = plot_gc_content_frequency_distribution(seq_qual_dict)
        fig.savefig(output+'/report_html/read_gc_histplot' + '.barplot.png')

        fig = plot_read_quality_frequency_distritution(seq_qual_dict)
        fig.savefig(output+'/report_html/read_quality_histplot' + '.barplot.png')

        fig = plot_read_length_frequency_distribution(seq_qual_dict)
        fig.savefig(output+'/report_html/read_length_histplot_nolog' + '.barplot.png')

        fig = plot_read_length_cumulative_distribution(seq_qual_dict)
        fig.savefig(output+'/report_html/read_length_cumulative' + '.barplot.png')

        plot_length_quality_cross_scatter(seq_qual_dict, output)
        plot_read_percent_qualtiy_curve(allBaseQual_dict, output)
        plot_ends_base_content_curve(endBaseQual_dict, output)
        plot_ends_base_quality_curve(endBaseQual_dict, output)
        plot_homopolymer_frequency(homopolymer_dict, output)
    else:
        raise IOError("FASTQ stat JSON file not found!")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-name", "--sample_name", default='cycads_report', required=False, help="prefix of output file name")
    args = vars(parser.parse_args())
    fq_figure_action(args)
