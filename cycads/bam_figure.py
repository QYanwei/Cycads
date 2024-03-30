#!/usr/bin/env python3
# coding: utf-8
import re, os, sys, time
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import pickle
import warnings
warnings.filterwarnings("ignore", "is_categorical_dtype")
warnings.filterwarnings("ignore", "use_inf_as_na")


from .plots import *

def plot_substitution_frequency(overall_aln_event_stat_dict):
    substitution = np.array(list(overall_aln_event_stat_dict['all_mis_typ_dict'].values()))
    all_mis_typ_dict = dict(zip(overall_aln_event_stat_dict['all_mis_typ_dict'], substitution * 100 / sum(substitution) ))
    substitution_type = pd.DataFrame(all_mis_typ_dict, index=[0])
    substitution_type.index = ['Substitution']
    
    fig, ax = plt.subplots(**figure_kw)
    sns.barplot(ax=ax, data=substitution_type)
    ax.grid(axis='y', **grid_kw)
    tick_labels = [col[0] + "→" + col[1] for col in substitution_type.columns]
    ticks = list(range(len(tick_labels)))
    ax.set_xticks(ticks)
    ax.set_xticklabels(tick_labels)
    ax.set_xlabel("Substitution type")
    ax.set_ylabel("Number of substitutions")
    ax.set_title("Summary of substitutions by type" , **title_kw)
    post_process_ax(ax)

    return fig



def plot_insertion_deletion_frequency(overall_aln_event_stat_dict):
    indel_len_max = 10
    indel_range_dict = {'Insertion': {}, 'Deletion':{} }
    ins_sum = sum(overall_aln_event_stat_dict['all_ins_len_dict'].values())
    del_sum = sum(overall_aln_event_stat_dict['all_del_len_dict'].values())
    for i in overall_aln_event_stat_dict['all_ins_len_dict'].keys():
        if int(i) < indel_len_max:
            indel_range_dict['Insertion'][int(i)] = overall_aln_event_stat_dict['all_ins_len_dict'][i] * 100 / ins_sum
        elif int(i) >= indel_len_max:
            if indel_len_max in indel_range_dict['Insertion'].keys():
                indel_range_dict['Insertion'][indel_len_max] += overall_aln_event_stat_dict['all_ins_len_dict'][i]
            else:
                indel_range_dict['Insertion'][indel_len_max] = overall_aln_event_stat_dict['all_ins_len_dict'][i]
        else:
            pass
    for i in overall_aln_event_stat_dict['all_del_len_dict'].keys():
        if int(i) < indel_len_max:
            indel_range_dict['Deletion'][int(i)] = overall_aln_event_stat_dict['all_del_len_dict'][i] * 100 / del_sum
        elif int(i) >= indel_len_max:
            if indel_len_max in indel_range_dict['Deletion'].keys():
                indel_range_dict['Deletion'][indel_len_max] += overall_aln_event_stat_dict['all_del_len_dict'][i]
            else:
                indel_range_dict['Deletion'][indel_len_max] = overall_aln_event_stat_dict['all_del_len_dict'][i]
        else:
            pass
    indel_range_dict['Insertion'][indel_len_max] = indel_range_dict['Insertion'][indel_len_max] * 100 / ins_sum
    indel_range_dict['Deletion'][indel_len_max] = indel_range_dict['Deletion'][indel_len_max] * 100 / del_sum

    # Insertion
    Ins_dataframe = pd.DataFrame(indel_range_dict['Insertion'], index=[0])
    Ins_dataframe.index = ['Insertion']


    # Insertion

    fig, ax = plt.subplots(**figure_kw)
    fig1 = fig
    sns.barplot(ax=ax, data=Ins_dataframe, color='steelblue')
    ax.grid(axis='y', **grid_kw)
    ax.set_xlabel("Insertion size (bp)")

    ticks = list(range(0, indel_len_max))
    tick_labels = [str(x+1) for x in ticks]
    tick_labels[-1] += "+"
    ax.set_xticks(ticks)
    ax.set_xticklabels(tick_labels)

    ax.set_ylabel("Number of insertions")
    ax.set_title("Size distribution of insertions", **title_kw)
    post_process_ax(ax)
    

    # Deletion
    Del_dataframe = pd.DataFrame(indel_range_dict['Deletion'], index=[0])
    Del_dataframe.index = ['Deletion']

    fig, ax = plt.subplots(**figure_kw)
    fig2 = fig
    sns.barplot(ax=ax, data=Del_dataframe, color='steelblue')
    ax.grid(axis='y', **grid_kw)
    ax.set_xticks(ticks)
    ax.set_xticklabels(tick_labels)
    ax.set_xlabel("Deletion size (bp)")
    ax.set_ylabel("Number of deletions")
    ax.set_title("Size distribution of insertions", **title_kw)
    post_process_ax(ax)

    return fig1, fig2

def plot_query_identity_rate_densities(query_aln_event_stat_dict):
    # identity rate
    qry_idy_rate = np.array(query_aln_event_stat_dict['qry_idy_rate'])
    qry_idy_df = pd.DataFrame([qry_idy_rate])
    qry_idy_df.index = ['Identity']

    fig, ax = plt.subplots(**figure_kw)
    sns.kdeplot(data=qry_idy_df.T, fill=True, alpha=.5, linewidth=0)
    ax.set_xlabel("Per-read identity")
    ax.set_xlim(right=1)
    ax.set_yticks([])
    ax.set_ylabel("Density")
    ax.set_title("Distribution of per-read identity", **title_kw)
    post_process_ax(ax)
    return fig

def plot_overall_homopolymer_length_event_frequency(homopolymer_aln_event_stat_dict):
    qry_hpm_event = np.array([homopolymer_aln_event_stat_dict['S'][i] for i in homopolymer_aln_event_stat_dict['S'].keys()])
    qry_hpm_event2= []
    for i in qry_hpm_event:
        lt = []
        p = 0
        if sum(i) > 0:
            for j in i /sum(i):
                p += j
                lt.append(p)
        qry_hpm_event2.append(lt)
    homopolymer_dataframe = pd.DataFrame(qry_hpm_event2)
    homopolymer_dataframe.index = range(2,len(qry_hpm_event2)+2,1)
    # homopolymer_dataframe.columns = ["≥4 bp contraction", "3 bp contraction", "2 bp contraction", "1 bp contraction",
    #                                   "others","correct", 
    #                                  "1 bp expansion", "2 bp expansion", "3 bp expansion", "≥4 bp expansion"]
    homopolymer_dataframe.columns = ["≤ -4 bp", "-3 bp", "-2 bp", "-1 bp",
                                      "others","correct", 
                                     "+1 bp", "+2 bp", "+3 bp", "≥ +4 bp"]

    cmap = plt.get_cmap("coolwarm")
    colors = [cmap(0), cmap(0.1), cmap(0.2), cmap(0.3), 
               "lightgray", "wheat",
              cmap(0.7), cmap(0.8), cmap(0.9), cmap(0.999)]
    
    fig, ax = plt.subplots(**figure_kw)
    xs = homopolymer_dataframe.index
    columns = homopolymer_dataframe.columns
    for j, col in enumerate(columns):
        if j == 0:
            y0 = [0 for _ in xs]
        
        y1 = homopolymer_dataframe[col]
        color = colors[j]
        ax.fill_between(xs, y0, y1, label=col, color=color)
        y0 = y1
    ax.legend(loc='upper left')
    ax.set_xlim(xs.min(), xs.max())
    ax.set_ylim(0, 1)
    ax.set_xlabel("Homopolymer length (bp)")
    ax.set_ylabel("Fraction of homopolymers")
    ax.set_title("Summary of homopolymer errors" , **title_kw)

    post_process_ax(ax)
    return fig


def plot_overall_alignment_frequency(overall_aln_event_sum_dict):
    all_event = overall_aln_event_sum_dict['identity'] + overall_aln_event_sum_dict['substitution'] + overall_aln_event_sum_dict['contraction'] + overall_aln_event_sum_dict['expansion']
    all_dif_event = overall_aln_event_sum_dict['substitution'] + overall_aln_event_sum_dict['contraction'] + overall_aln_event_sum_dict['expansion']
    all_dif = all_dif_event / all_event *100
    all_mis = overall_aln_event_sum_dict['substitution'] / all_event *100
    all_del = overall_aln_event_sum_dict['contraction'] / all_event *100
    all_ins = overall_aln_event_sum_dict['expansion'] / all_event *100


    all_aln_rates = [all_dif, all_mis, all_del, all_ins]

    labels = ['Overall', 'Mismatch', 'Deletion', 'Insertion']
    x = np.arange(len(labels))  # the label locations
    width = 0.35  # the width of the bars
    fig, ax = plt.subplots(**figure_kw)
    ax.bar(x, all_aln_rates, width, label='Overall error rate', color='steelblue')
    ax.grid(axis='y', **grid_kw)
    ax.set_ylabel('Error rate (%)')
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.legend(loc='upper right')
    ax.set_title("Error rates by type" , **title_kw)
    post_process_ax(ax)

    return fig

def bam_figure_action(args):
    output_foler = os.path.join(args["output_dir"], args["sample_name"])
    bam_pickle_file_path = args['bam_pickle_path']
    with open(bam_pickle_file_path, 'rb') as picklefile:
        bam_datum_dict = pickle.load(picklefile)

    fig = plot_overall_homopolymer_length_event_frequency(bam_datum_dict['homopolymer_aln_event_stat_dict'])
    fig.savefig(f'{output_foler}/report_html/query_homopolymer_length_event.lineplot.png')
    
    fig1, fig2 = plot_insertion_deletion_frequency(bam_datum_dict['overall_aln_event_stat_dict'])
    fig1.savefig(f'{output_foler}/report_html/query_insertion_frequency.barplot.png')
    fig2.savefig(f'{output_foler}/report_html/query_deletion_frequency.barplot.png')

    fig = plot_overall_alignment_frequency(bam_datum_dict['overall_aln_event_sum_dict'])
    fig.savefig(f'{output_foler}/report_html/query_all_error_item.barplot.png')
    
    fig = plot_substitution_frequency(bam_datum_dict['overall_aln_event_stat_dict'])
    fig.savefig(f'{output_foler}/report_html/query_all_substitution_errors.barplot.png')

    fig = plot_query_identity_rate_densities(bam_datum_dict['query_aln_event_stat_dict'])
    fig.savefig(f'{output_foler}/report_html/query_events_curve_idy.displot.png')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output_dir", default='./', required=False, help="Output direcotry")
    parser.add_argument("-n", "--sample_name", default='cycads_report', required=False, help="prefix of output file name")

    args = vars(parser.parse_args())
    bam_figure_action(args)

