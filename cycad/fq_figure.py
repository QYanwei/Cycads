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
import pickle
import warnings

warnings.filterwarnings("ignore", "is_categorical_dtype")
warnings.filterwarnings("ignore", "use_inf_as_na")

plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['figure.dpi'] = 300
palette = dict(A="tab:red", T="tab:green", C="tab:blue", G="tab:purple", S='tab:black')
figsize = (5, 4)
hist_kw = dict(facecolor='tab:blue', edgecolor='k', linewidth=0.5)
grid_kw = dict(color='k', alpha=0.1)
title_kw = dict(fontsize=10)


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
    ax.set_title("Read length Nx metrics", **title_kw)
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
    ax.set_title("Read GC content distribution", **title_kw)
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
    ax.set_title("Mean base quality distribution", **title_kw)
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
    ax.set_title("Read length distribution", **title_kw)
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
    ax.set_title("Cumulative fraction of bases vs. read length", **title_kw)
    post_process_ax(ax)
    return fig

def plot_length_quality_2d_histogram(seq_qual_dict):
    lengths = np.array(seq_qual_dict['LEN']) / 1000
    quality_values = np.array(seq_qual_dict['QUAL1'])

    fig, ax = plt.subplots(figsize=figsize)
    ax.hist2d(lengths, quality_values, bins=(50, 50), cmap='Blues')
    ax.set_xlim(lengths.min(), lengths.max())
    ax.set_ylim(quality_values.min(), quality_values.max())
    ax.set_xlabel("Read length (kb)")
    ax.set_ylabel("Mean per-read base quality")
    ax.set_title("Joint distribution of read length and mean base quality", **title_kw)
    post_process_ax(ax)
    return fig
    

def plot_homopolymer_frequency(homopolymer_dict):
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
    xs = list(range(2, homopolymer_len_max + 1))
    fig, ax = plt.subplots(figsize=figsize)

    for base in ("A", "T", "C", "G"):
        ys = [homopolymer_range_dict[base][x] for x in xs]
        color = palette[base]
        ax.plot(xs, ys, label=base)
    ax.legend(loc='upper right')
    ax.set_xlim(min(xs), max(xs))
    ax.set_xlabel("Homopolymer length")
    ax.set_yscale("log")
    ax.set_ylabel("Total number of occurances")
    ax.set_title("Homopolymer length distribution", **title_kw)
    post_process_ax(ax)

    return fig


def plot_ends_base_content_curve(endBaseQual_dict):
    for key in ("HeadBaseContent_dict", "TailBaseContent_dict"):
        fig, ax = plt.subplots(figsize=figsize)
        if key == "HeadBaseContent_dict":
            fig1 = fig
        else:
            fig2 = fig
        data = endBaseQual_dict[key]
        for base in ("A", "T", "C", "G"):
            ys = np.array(data[base]) / data['S']
            ys *= 100
            xs = np.arange(ys.shape[0])
            if key == "TailBaseContent_dict":
                xs = xs.max() - xs
            color = palette[base]
            ax.plot(xs, ys, color=color, ls='-', label=base)
        ys = (np.array(data["G"]) + np.array(data["C"])) / data['S']
        ys *= 100
        ax.plot(xs, ys, color='tab:brown', ls='-', label="G+C")

        ax.grid(axis='both', **grid_kw)
        
        if key == "HeadBaseContent_dict":
            ax.set_xlim(0, xs.max() + 1)
            ax.set_xlabel("Distance from 5' end")
            ax.legend(loc='upper right')
            ax.set_title("Base content near 5' end", **title_kw)
        else:
            ax.set_xlim(xs.max() + 1, 0)
            ax.set_xlabel("Distance from 3' end")
            ax.legend(loc='upper left')
            ax.set_title("Base content near 3' end", **title_kw)
            
        ax.yaxis.set_major_formatter(mtick.PercentFormatter())
        ax.set_ylabel("Frequency of base")
        post_process_ax(ax)
    
    return fig1, fig2

def plot_ends_base_quality_curve(endBaseQual_dict):

    for key in ("HeadQualContent_dict", "TailQualContent_dict"):
        fig, ax = plt.subplots(figsize=figsize)
        if key == "HeadQualContent_dict":
            fig1 = fig
        else:
            fig2 = fig
        data = endBaseQual_dict[key]
        for base in ("A", "T", "C", "G"):
            ys = np.array(data[base]) / np.array(data['S'][base])
            xs = np.arange(ys.shape[0])
            if key == "TailQualContent_dict":
                xs = xs.max() - xs
            color = palette[base]
            ax.plot(xs, ys, color=color, ls='-', label=base)
        ys = sum(np.array(quality_values) for key, quality_values in data.items() if key in ("A", "T", "C", "G")) / sum(np.array(base_counts) for key, base_counts in data['S'].items() if key in ("A", "T", "C", "G"))
        ax.plot(xs, ys, color='tab:brown', ls='-', label="Overall")

        ax.grid(axis='both', **grid_kw)
        
        if key == "HeadQualContent_dict":
            ax.set_xlim(0, xs.max() + 1)
            ax.set_xlabel("Distance from 5' end")
            ax.legend(loc='lower right')
            ax.set_title("Base content near 5' end", **title_kw)
        else:
            ax.set_xlim(xs.max() + 1, 0)
            ax.set_xlabel("Distance from 3' end")
            ax.legend(loc='lower left')
            ax.set_title("Base content near 3' end", **title_kw)
            
        ax.set_ylabel("Mean base quality")
        post_process_ax(ax)

    return fig1, fig2
    

def plot_read_percent_qualtiy_curve(allBaseQual_dict):
    quality_values = np.array(allBaseQual_dict['PercentBaseQual_dict']['Q'])/ allBaseQual_dict['PercentBaseQual_dict']['S']
    xs = np.arange(1, quality_values.shape[0]+1)
    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(xs, quality_values, clip_on=False, color='steelblue')
    ax.set_xlim(0, 100)
    ax.set_xlabel("Relative position in read")
    ax.set_ylabel("Mean base quality")
    ax.set_title("Mean base quality vs. relative position in read", **title_kw)
    ax.grid(axis='both', **grid_kw)
    post_process_ax(ax)
    return fig

    

def fq_figure_action(args):
    output_folder = args["sample_name"]
    pickle_path = os.path.join(output_folder, "fq.pickle")
    if not os.path.isfile(pickle_path):
        raise IOError(f"Unable to find FASTQ stats file: {pickle_path}")

    with open(pickle_path, 'rb') as picklefile:
        fq_datum_dict = pickle.load(picklefile)
        seq_qual_dict = fq_datum_dict['seq_qual_dict']
        homopolymer_dict = fq_datum_dict['homopolymer_dict']
        endBaseQual_dict = fq_datum_dict['endBaseQual_dict']
        allBaseQual_dict = fq_datum_dict['allBaseQual_dict']

    fig = plot_length_Nx_average_bar(seq_qual_dict)
    fig.savefig(output_folder+'/report_html/read_length_biostat' + '.barplot.png')

    fig = plot_gc_content_frequency_distribution(seq_qual_dict)
    fig.savefig(output_folder+'/report_html/read_gc_histplot' + '.barplot.png')

    fig = plot_read_quality_frequency_distritution(seq_qual_dict)
    fig.savefig(output_folder+'/report_html/read_quality_histplot' + '.barplot.png')

    fig = plot_read_length_frequency_distribution(seq_qual_dict)
    fig.savefig(output_folder+'/report_html/read_length_histplot_nolog' + '.barplot.png')

    fig = plot_read_length_cumulative_distribution(seq_qual_dict)
    fig.savefig(output_folder+'/report_html/read_length_cumulative' + '.barplot.png')

    fig = plot_length_quality_2d_histogram(seq_qual_dict)
    fig.savefig(output_folder+'/report_html/read_length_quality_cross' + '.scatterplot.png')

    fig = plot_read_percent_qualtiy_curve(allBaseQual_dict)
    fig.savefig(output_folder+'/report_html/read_relative_position_avg_qual' + '.lineplot.png')

    fig1, fig2 = plot_ends_base_content_curve(endBaseQual_dict)
    fig1.savefig(output_folder+'/report_html/read_head_base_content' + '.lineplot.png')
    fig2.savefig(output_folder+'/report_html/read_tail_base_content' + '.lineplot.png')

    fig1, fig2 = plot_ends_base_quality_curve(endBaseQual_dict)
    fig1.savefig(output_folder+'/report_html/read_head_base_quality' + '.lineplot.png')
    fig2.savefig(output_folder+'/report_html/read_tail_base_quality' + '.lineplot.png')

    fig = plot_homopolymer_frequency(homopolymer_dict)
    fig.savefig(output_folder+'/report_html/read_homopolymer_frequency' + '.lineplot.png')
        

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-name", "--sample_name", default='cycads_report', required=False, help="prefix of output file name")
    args = vars(parser.parse_args())
    fq_figure_action(args)
