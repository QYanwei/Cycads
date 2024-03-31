import os,re,sys,time
import datetime
import argparse
import pickle
import jinja2
from statistics import mean
# import html template
def import_jinja_template(template_file): # TODO: enable custom template for advanced users
    # get template_fq_report.html
    if template_file:
        try:
            with open(template_file) as fp:
                template = jinja2.Template(fp.read())
                return template
        except (FileNotFoundError, IOError, jinja2.exceptions.TemplateNotFound, jinja2.exceptions.TemplateSyntaxError):
            print("File not found, non-readable or invalid\n")
def N50(list_read_length):
    # Calculate the total length of the sequences
    total_length = sum(list_read_length)
    list_read_length.sort(reverse=True)
    base_sum = 0
    for length in list_read_length:
        base_sum += length
        if base_sum/total_length == 0.5:
           return length
# fq report data
def generate_fq_report_strings(args, output_folder):
    table_name = ["FileName", "TotalRead", "TotalBases", "GC%", "AvgLen", "MaxLen", "MinLen", "N50", "AvgReadQ", "MaxReadQ", "MinReadQ"]
    #table_list = ["TB2000B609-202403200954240_read.fq.gz", "873663", "138045423", "54.542", "158.008", "31187", "1", "27.5", "7", "0"]
    fq_sum = args['fastq_summary_path']
    with open(fq_sum) as txtfile:
        txtfile.readline()
        fq_info = txtfile.readline().strip().split()
    table_list = fq_info[:7]
    fq_pickle = args['fastq_pickle_path']
    with open(fq_pickle, 'rb') as picklefile:
        fq_datum_dict = pickle.load(picklefile)
        seq_qual_dict = fq_datum_dict['seq_qual_dict']
        table_list.append(N50(seq_qual_dict["LEN"]))
        table_list.append(mean(seq_qual_dict["QUAL2"]))
        table_list.append(max(seq_qual_dict["QUAL2"]))
        table_list.append(min(seq_qual_dict["QUAL2"]))
    fq_table_string = "<thead><tr>"
    for i in table_name:
        fq_table_string += "<th>{}</th>".format(i)
    fq_table_string += "</tr><tr>"
    for i in table_list:
        fq_table_string += "<th>{}</th>".format(i)
    fq_table_string += "</tr><thead>"
    # plots string create
    plots_list = [["read_length_quality_cross.scatterplot.png",
                   "read_length_histplot_nolog.barplot.png"],
                  ["read_relative_position_avg_qual.lineplot.png",
                   "read_length_cumulative.barplot.png"],
                  ["read_gc_histplot.barplot.png",
                   "read_quality_histplot.barplot.png"],
                  ["read_length_biostat.barplot.png",
                   "read_homopolymer_frequency.lineplot.png"],
                  ["read_head_base_quality.lineplot.png",
                   "read_tail_base_quality.lineplot.png"],
                  ["read_head_base_content.lineplot.png",
                   "read_tail_base_content.lineplot.png"]]
    plots_src = [["Read length vs quality scatter",
                  "Read relative position average Q score"],
                 ["Read length distribution",
                  "Read length cumulative curve"],
                 ["Read GC content histgram",
                  "Read quality score"],
                 ["Read Nx and average length bar",
                  "Read homopolymer frequency"],
                 ["Read head base quality curves",
                  "Read tail base quality curves"],
                 ["Read head base content curves",
                  "Read tail base content curves"]]
    fq_plots_string = ""
    for i in range(6):
        fq_plots_string += "<tr>"
        for j in range(2):
            fq_plots_string += "<td><figure>"
            fq_plots_string += "<img src=\"./{}\" width=\"auto\" height=\"480\" />".format(plots_list[i][j])
            #fq_plots_string += "<figcaption class=\"figure-caption text-center\">{}</figcaption>".format(plots_src[i][j])
            fq_plots_string += "</td>"
        fq_plots_string += "</tr>"
    return fq_table_string, fq_plots_string
    
# bam report data
def generate_bam_report_string(args, output_folder):
    table_list = [args["sample_name"]]
    table_name = ["SampleName", "TotalReads", "MappedReads", "Indentity(%)","TotalErr(%)", "MismatchErr(%)", "InsertionErr(%)", "DeletionErr(%)", "HomopolymerErr(%)"]
    bam_pickle = args['bam_pickle_path']
    with open(bam_pickle, 'rb') as picklefile:
        bam_datum_dict = pickle.load(picklefile)
        overall_event_dict = bam_datum_dict['overall_aln_event_sum_dict']
        all_evt = overall_event_dict['substitution'] + overall_event_dict['contraction'] + overall_event_dict['expansion'] + overall_event_dict['identity']
        all_idy = round(overall_event_dict['identity'] / all_evt *100, 2)
        all_dif = round(100 - all_idy, 2)
        all_mis = round(overall_event_dict['substitution'] / all_evt *100, 2)
        all_del = round(overall_event_dict['contraction'] / all_evt *100, 2)
        all_ins = round(overall_event_dict['expansion'] / all_evt *100, 2)
        hpm_evt = overall_event_dict['hpm_substitution'] + overall_event_dict['hpm_contraction'] + overall_event_dict['hpm_expansion']
        hpm_dif = round(hpm_evt/ all_evt * 100, 2)
#        non_hpm_dif = round(all_dif - hpm_dif, 2)
        table_list.append(overall_event_dict['total_reads'])
        table_list.append(overall_event_dict['mapped_reads'])
        table_list.extend([ all_idy, all_dif, all_mis, all_ins, all_del, hpm_dif])
    bam_table_string = "<thead><tr>"
    for i in table_name:
        bam_table_string += "<th>{}</th>".format(i)
    bam_table_string += "</tr><tr>"
    for i in table_list:
        bam_table_string += "<th>{}</th>".format(i)
    bam_table_string += "</tr><thead>"

    plots_list = [["query_all_error_item.barplot.png",
                   "query_all_substitution_errors.barplot.png"],
                  ["query_deletion_frequency.barplot.png",
                   "query_insertion_frequency.barplot.png"],
                  ["query_events_curve_idy.displot.png",
                  "query_homopolymer_length_event.lineplot.png"]]
    plots_src = [["Sequencing error rate",
                  "Substitution type error frequency"],
                 ["Deletion frequency",
                  "Insertion frequency"],
                 ["Read identity rate distribution",
                 "Homopolymer error type with its length"]]
    bam_plots_string = ""
    for i in range(3):
        bam_plots_string += "<tr>"
        for j in range(2):
            bam_plots_string += "<td><figure>"
            bam_plots_string += "<img src=\"./{}\" width=\"auto\" height=\"480\" />".format(plots_list[i][j])
            bam_plots_string += "<figcaption class=\"figure-caption text-center\">{}</figcaption>".format(plots_src[i][j])
            bam_plots_string += "</td>"
        bam_plots_string += "</tr>"
    return bam_table_string, bam_plots_string

def generate_report_html(args, output_folder, flag):
    output_path = os.path.join(output_folder, "summary.html")

    # set the basic report infortion
    report_title = "CycloneSEQ quality reporter"
    report_subtitle = "Created on {} with {} {}".format(datetime.datetime.now().strftime("%d/%m/%y"), "Cycads", "0.3.0")
    pwd_config_file = os.path.realpath(__file__)
    if flag == 0:
        fq_table_string, fq_plots_string = generate_fq_report_strings(args, output_folder)
        template_file = "config/template_fq_report.j2"
        template = import_jinja_template(template_file)
        rendering = template.render(
           fq_table=fq_table_string,
           fq_plots=fq_plots_string,
           report_title=report_title,
           report_subtitle=report_subtitle)
        # Write to HTML file
        with open(output_path, "w") as fp:
            fp.write(rendering)
        
    elif flag == 1:
        bam_table_string, bam_plots_string = generate_bam_report_string(args, output_folder)
        template_file = "config/template_bam_report.j2"
        template = import_jinja_template(template_file)
        rendering = template.render(
            bam_table=bam_table_string,
            bam_plots=bam_plots_string,
            report_title=report_title,
            report_subtitle=report_subtitle)
        # Write to HTML file
        with open(output_path, "w") as fp:
            fp.write(rendering)
        
    elif flag == 2:
        fq_table_string, fq_plots_string = generate_fq_report_strings(args, output_folder)
        bam_table_string, bam_plots_string = generate_bam_report_string(args, output_folder)
        template_file = "config/template_all_report.j2"
        template = import_jinja_template(template_file)
        rendering = template.render(
            fq_table=fq_table_string,
            fq_plots=fq_plots_string,
            bam_table=bam_table_string,
            bam_plots=bam_plots_string,
            report_title=report_title,
            report_subtitle=report_subtitle)
        # Write to HTML file
        with open(output_path, "w") as fp:
            fp.write(rendering)

    abs_output_path = os.path.abspath(output_path)
    print(f"HTML report written to {abs_output_path}")

def generate_html(args):
    fq_pickle = args['fastq_pickle_path']
    bam_pickle = args['bam_pickle_path']
    output_folder = args['report_dir']
    if os.path.isfile(fq_pickle) and os.path.isfile(bam_pickle):
        generate_report_html(args, output_folder, 2)
    elif os.path.isfile(fq_pickle) and not os.path.isfile(bam_pickle):
        generate_report_html(args, output_folder, 0)
    elif not os.path.isfile(fq_pickle) and os.path.isfile(bam_pickle):
        generate_report_html(args, output_folder, 1)
    else:
        raise IOError(f"Unable to find {fq_pickle} or {bam_pickle}")


# if __name__ == '__main__':
#      parser = argparse.ArgumentParser()
#      parser.add_argument("-o", "--output_dir", default='./', required=False, help="Output direcotry")
#      parser.add_argument("-n", "--sample_name", default='cycads_report', required=False, help="prefix of output file name")
#      args = vars(parser.parse_args())
#      generate_html(args)
