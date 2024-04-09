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
from .plots import *

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
    
    fig, ax = plt.subplots(**figure_kw)
    g = sns.barplot(LengthFrame, ax=ax, x='nx', y='read_length', color='steelblue')
    ax.axhline(average_length, color='k', ls='--', label='Mean read length')
    ax.legend(loc='upper right')
    ax.set_xlabel("")
    ax.set_ylabel("Read length (bp)", **labelsize_kw)
    ax.grid(axis='y', **grid_kw)
    ax.tick_params(axis='both', **ticksize_kw)
    ax.set_title("Read length Nx metrics", **title_kw)
    post_process_ax(ax)
    return fig

def plot_gc_content_frequency_distribution(seq_qual_dict):
    gc_percentage = np.array( seq_qual_dict['GC'] ).astype(float) * 100
    fig, ax = plt.subplots(**figure_kw)
    ax.hist(gc_percentage, bins=np.arange(0, 100, 1), density=True, **hist_kw)
    ax.xaxis.set_major_formatter(mtick.PercentFormatter())
    ax.set_xlim(0, 100)
    ax.set_xlabel("Read GC content", **labelsize_kw)
    ax.set_xticks(list(range(0, 101, 10)))
    ax.set_ylabel("Density", **labelsize_kw)
    ax.set_yticks([])
    ax.tick_params(axis='both', **ticksize_kw)
    ax.set_title("Read GC content distribution", **title_kw)
    post_process_ax(ax)
    return fig

def plot_read_quality_frequency_distritution(seq_qual_dict):
    quality_values = np.array( seq_qual_dict['QUAL2'] )
    fig, ax = plt.subplots(**figure_kw)
    xmin = floor(quality_values.min())
    xmax = ceil(quality_values.max())
    bins = np.arange(xmin, xmax + 1, 0.5)
    ax.hist(quality_values, bins=bins, density=True, **hist_kw)
    ax.set_xlabel("Mean per-read base quality", **labelsize_kw)
    ax.set_xlim(xmin, xmax)
    ax.set_ylabel("Density", **labelsize_kw)
    ax.set_yticks([])
    ax.tick_params(axis='both', **ticksize_kw)
    ax.set_title("Mean base quality distribution", **title_kw)
    post_process_ax(ax)
    return fig

def plot_read_length_frequency_distribution(seq_qual_dict):
    lengths = np.array( seq_qual_dict['LEN'] )

    fig, ax = plt.subplots(**figure_kw)
    ax.hist(lengths/1000, bins=100, **hist_kw)
    ax.set_xlim(left=0)
    ax.set_xlabel("Read length (kb)", **labelsize_kw)
    ax.set_ylabel("Density", **labelsize_kw)
    ax.set_yticks([])
    ax.tick_params(axis='both', **ticksize_kw)
    ax.set_title("Read length distribution", **title_kw)
    post_process_ax(ax)
    return fig

def plot_read_length_cumulative_distribution(seq_qual_dict):
    lengths = np.array( seq_qual_dict['LEN'] )
    total_bases = lengths.sum()
    xs = np.arange(lengths.min(), lengths.max() + 100, 100) / 1000
    ys = [lengths[lengths >= x * 1000].sum() * 100 / total_bases for x in xs]

    fig, ax = plt.subplots(**figure_kw)
    ax.plot(xs, ys)
    ax.set_ylim(0, 100)
    ax.set_xlim(left=0)
    ax.set_xlabel("Min. read length (kb)", **labelsize_kw)
    ax.set_ylabel("Cumulative fraction of bases", **labelsize_kw)
    ax.yaxis.set_major_formatter(mtick.PercentFormatter())
    ax.grid(axis='both', **grid_kw)
    ax.tick_params(axis='both', **ticksize_kw)
    ax.set_title("Cumulative fraction of bases vs. read length", **title_kw)
    post_process_ax(ax)
    return fig

def plot_length_quality_2d_histogram(seq_qual_dict):
    lengths = np.array(seq_qual_dict['LEN']) / 1000
    quality_values = np.array(seq_qual_dict['QUAL2'])

    fig, ax = plt.subplots(**figure_kw)
    ax.hist2d(lengths, quality_values, bins=(50, 50), cmap='Blues')
    ax.set_xlim(lengths.min(), lengths.max())
    ax.set_ylim(quality_values.min(), quality_values.max())
    ax.set_xlabel("Read length (kb)", **labelsize_kw)
    ax.set_ylabel("Mean per-read base quality", **labelsize_kw)
    ax.set_title("Read length VS mean base quality", **title_kw)
    ax.tick_params(axis='both', **ticksize_kw)
    post_process_ax(ax)
    return fig
    

def plot_homopolymer_frequency(homopolymer_dict, args):
    homopolymer_len_max = args['max_homopolymer_size']
    homopolymer_len_min = args['min_homopolymer_size']
    homopolymer_range_dict = {'A': {}, 'G': {}, 'C': {}, 'T': {} }
    xs = list(range( homopolymer_len_min, homopolymer_len_max + 1))
    for b in homopolymer_dict.keys():
        for x in xs:
            if x <= homopolymer_len_max and x not in homopolymer_dict[b].keys():
                homopolymer_range_dict[b][x] = 0
            elif x <= homopolymer_len_max and x in homopolymer_dict[b].keys():
                homopolymer_range_dict[b][x] = homopolymer_dict[b][x]
            elif x > homopolymer_len_max and x in homopolymer_dict[b].keys():
                homopolymer_range_dict[b][homopolymer_len_max] += homopolymer_dict[b][x]
            else:
                pass
    fig, ax = plt.subplots(**figure_kw)
    for base in ("A", "T", "C", "G"):
        ys = [homopolymer_range_dict[base][x] for x in xs]
        color = palette[base]
        ax.plot(xs, ys, label=base)
    ax.legend(loc='upper right')
    ax.set_xlim(min(xs), max(xs))
    ax.set_xlabel("Homopolymer length", **labelsize_kw)
    ax.set_yscale("log")
    ax.set_ylabel("Total number of occurances", **labelsize_kw)
    ax.set_title("Homopolymer length distribution", **title_kw)
    ax.tick_params(axis='both', **ticksize_kw)
    post_process_ax(ax)
    return fig


def plot_ends_base_content_curve(endBaseQual_dict):
    for key in ("HeadBaseContent_dict", "TailBaseContent_dict"):
        fig, ax = plt.subplots(**figure_kw)
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
            ax.set_xlabel("Distance from 5' end", **labelsize_kw)
            ax.legend(loc='upper right')
            ax.set_title("Base content near 5' end", **title_kw)
        else:
            ax.set_xlim(xs.max() + 1, 0)
            ax.set_xlabel("Distance from 3' end", **labelsize_kw)
            ax.legend(loc='upper left')
            ax.set_title("Base content near 3' end", **title_kw)
        ax.yaxis.set_major_formatter(mtick.PercentFormatter())
        ax.set_ylabel("Frequency of base", **labelsize_kw)
        ax.tick_params(axis='both', **ticksize_kw)
        post_process_ax(ax)
    return fig1, fig2

def plot_ends_base_quality_curve(endBaseQual_dict):
    for key in ("HeadQualContent_dict", "TailQualContent_dict"):
        fig, ax = plt.subplots(**figure_kw)
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
            ax.set_xlabel("Distance from 5' end", **labelsize_kw)
            ax.legend(loc='lower right')
            ax.set_title("Base content near 5' end", **title_kw)
        else:
            ax.set_xlim(xs.max() + 1, 0)
            ax.set_xlabel("Distance from 3' end", **labelsize_kw)
            ax.legend(loc='lower left')
            ax.set_title("Base content near 3' end", **title_kw)
            
        ax.tick_params(axis='both', **ticksize_kw)
        ax.set_ylabel("Mean base quality", **labelsize_kw)
        post_process_ax(ax)
    return fig1, fig2

def plot_read_percent_qualtiy_curve(allBaseQual_dict):
    quality_values = np.array(allBaseQual_dict['PercentBaseQual_dict']['Q'])/ allBaseQual_dict['PercentBaseQual_dict']['S']
    xs = np.arange(1, quality_values.shape[0]+1)
    fig, ax = plt.subplots(**figure_kw)
    ax.plot(xs, quality_values, clip_on=False, color='steelblue')
    ax.set_xlim(0, 100)
    ax.set_xlabel("Relative position in read", **labelsize_kw)
    ax.set_ylabel("Mean base quality", **labelsize_kw)
    ax.set_title("Mean base quality in Read relative pos.", **title_kw)
    ax.grid(axis='both', **grid_kw)
    ax.tick_params(axis='both', **ticksize_kw)
    post_process_ax(ax)
    return fig

def fq_figure_action(args):
    pickle_path = args['fastq_pickle_path']
    if not os.path.isfile(pickle_path):
        raise IOError(f"Unable to find FASTQ stats file: {pickle_path}")

    with open(pickle_path, 'rb') as picklefile:
        fq_datum_dict = pickle.load(picklefile)
    seq_qual_dict = fq_datum_dict['seq_qual_dict']
    homopolymer_dict = fq_datum_dict['homopolymer_dict']
    endBaseQual_dict = fq_datum_dict['endBaseQual_dict']
    allBaseQual_dict = fq_datum_dict['allBaseQual_dict']
    output_folder = args['report_dir']
    # read_length_biostat
    fig = plot_length_Nx_average_bar(seq_qual_dict)
    fig.savefig(os.path.join(output_folder, 'read_length_biostat.barplot.png'))
    # read_gc_histplot
    fig = plot_gc_content_frequency_distribution(seq_qual_dict)
    fig.savefig(os.path.join(output_folder, 'read_gc_histplot.barplot.png'))
    # read_quality_histplot
    fig = plot_read_quality_frequency_distritution(seq_qual_dict)
    fig.savefig(os.path.join(output_folder, 'read_quality_histplot.barplot.png'))
    # read_length_histplot_nolog
    fig = plot_read_length_frequency_distribution(seq_qual_dict)
    fig.savefig(os.path.join(output_folder, 'read_length_histplot_nolog.barplot.png'))
    # read_length_cumulative
    fig = plot_read_length_cumulative_distribution(seq_qual_dict)
    fig.savefig(os.path.join(output_folder, 'read_length_cumulative.barplot.png'))
    # read_length_quality_cross
    fig = plot_length_quality_2d_histogram(seq_qual_dict)
    fig.savefig(os.path.join(output_folder, 'read_length_quality_cross.scatterplot.png'))
    # read_relative_position_avg_qual
    fig = plot_read_percent_qualtiy_curve(allBaseQual_dict)
    fig.savefig(os.path.join(output_folder, 'read_relative_position_avg_qual.lineplot.png'))
    # read_head_base_content
    fig1, fig2 = plot_ends_base_content_curve(endBaseQual_dict)
    fig1.savefig(os.path.join(output_folder, 'read_head_base_content.lineplot.png'))
    fig2.savefig(os.path.join(output_folder, 'read_tail_base_content.lineplot.png'))
    # read_head_base_quality
    fig1, fig2 = plot_ends_base_quality_curve(endBaseQual_dict)
    fig1.savefig(os.path.join(output_folder, 'read_head_base_quality.lineplot.png'))
    fig2.savefig(os.path.join(output_folder, 'read_tail_base_quality.lineplot.png'))
    # read_homopolymer_frequency
    fig = plot_homopolymer_frequency(homopolymer_dict, args)
    fig.savefig(os.path.join(output_folder, 'read_homopolymer_frequency.lineplot.png'))
