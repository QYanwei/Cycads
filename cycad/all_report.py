import os,re,sys,time
import datetime
import argparse
import jinja2

def import_jinja_template(template_file):
    # get template_fq_report.html
    if template_file:
        print("Try to load provided jinja template file\n")
        try:
            with open(template_file) as fp:
                template = jinja2.Template(fp.read())
                return template
        except (FileNotFoundError, IOError, jinja2.exceptions.TemplateNotFound, jinja2.exceptions.TemplateSyntaxError):
            print("File not found, non-readable or invalid\n")


def generate_fq_report_html(args):
    #fq
    table_list = ["TB2000B609-202403200954240_read.fq.gz", "873663", "138045423", "54.542", "158.008", "31187", "1",
                  "27.5", "7"]
    fq_table_string = ""
    for i in table_list:
        fq_table_string += "<th>{}</th>".format(i)
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
            fq_plots_string += "<img src=\"./{}\" width=\"600\" height=\"580\" />".format(plots_list[i][j])
            fq_plots_string += "<figcaption class=\"figure-caption text-center\">{}</figcaption>".format(plots_src[i][j])
            fq_plots_string += "</td>"
        fq_plots_string += "</tr>"
    # bam
    table_list = ["TB2000B609-202403200954240_read.bam", "2399", "3939982", "95.02", "3.98", "1.53", "1.32", "1.13",
                  "2.3", "1.68"]
    bam_table_string = ""
    for i in table_list:
        bam_table_string += "<th>{}</th>".format(i)
    plots_list = [["query_all_error_item.barplot.png",
                   "query_insertion_frequency.barplot.png"],
                  ["query_deletion_frequency.barplot.png",
                   "query_all_substitution_errors.barplot.png"],
                  ["query_events_curve_idy.dispplot.png",
                   "query_events_curve_dif.dispplot.png"],
                  ["query_homopolymer_length_event.lineplot.png", ""]]
    plots_src = [["Sequencing error rate",
                  "Insertion frequency"],
                 ["Deletion frequency",
                  "Substitution type error frequency"],
                 ["Read identity rate distribution",
                  "Read error rate distribution"],
                 ["Homopolymer error type with its length", ""]]
    bam_plots_string = ""
    for i in range(4):
        bam_plots_string += "<tr>"
        for j in range(2):
            bam_plots_string += "<td><figure>"
            bam_plots_string += "<img src=\"./{}\" width=\"600\" height=\"580\" />".format(plots_list[i][j])
            bam_plots_string += "<figcaption class=\"figure-caption text-center\">{}</figcaption>".format(plots_src[i][j])
            bam_plots_string += "</td>"
        bam_plots_string += "</tr>"

    report_title = "CycloneSEQ quality reporter"
    report_subtitle = "Created on {} with {} {}".format(datetime.datetime.now().strftime("%d/%m/%y"), "Cycads",
                                                          "0.3.0")
    pwd_config_file = os.path.realpath(__file__)
    template_file = '/'.join(pwd_config_file.split('/')[:-1]) + "../config/template_all_report.j2"
    template = import_jinja_template(template_file)
    rendering = template.render(
        fq_table=fq_table_string,
        fq_plots=fq_plots_string,
        bam_table=bam_table_string,
        bam_plots=bam_plots_string,
        report_title=report_title,
        report_subtitle=report_subtitle)
    
    # Write to HTML file
    outfile = "../" + args["sample_name"] + "/" + args["sample_name"] + "_qc_report.html"
    with open(outfile, "w") as fp:
        fp.write(rendering)


args = {"sample_name": "test"}
generate_fq_report_html(args)
# if __name__ == '__main__':
#     parser = argparse.ArgumentParser()
#     parser.add_argument("-name", "--sample_name", default='cycads_report', required=False, help="prefix of output file name")
#     args = vars(parser.parse_args())
#     generate_fq_report_html(args)