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
    
    table_list = ["TB2000B609-202403200954240_read.fq.gz","873663","138045423","54.542","158.008","31187","1","27.5","7"]
    table_string = ""
    for i in table_list:
        table_string += "<th>{}</th>".format(i)

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

    plots_string = ""
    for i in range(6):
        plots_string += "<tr>"
        for j in range(2):
            plots_string += "<td><figure>"
            plots_string += "<img src=\"./{}\" width=\"600\" height=\"580\" />".format(plots_list[i][j])
            plots_string += "<figcaption class=\"figure-caption text-center\">{}</figcaption>".format(plots_src[i][j])
            plots_string += "</td>"
        plots_string += "</tr>"
    report_title="CycloneSEQ quality reporter"
    report_subtitle="Created on {} with {} {}".format( datetime.datetime.now().strftime("%d/%m/%y"), "Cycads", "0.3.0")
    pwd_config_file = os.path.realpath(__file__)
    template_file = '/'.join(pwd_config_file.split('/')[:-1]) + "../config/template_fq_report.j2"
    template = import_jinja_template(template_file)
    rendering = template.render(
        table=table_string,
        plots=plots_string,
        report_title=report_title,
        report_subtitle=report_subtitle)
    
    # Write to HTML file
    outfile = "../"+args["sample_name"] + "/" + args["sample_name"] + "_quality_report.html"
    with open(outfile, "w") as fp:
        fp.write(rendering)

args = {"sample_name":"test"}
generate_fq_report_html(args)
# if __name__ == '__main__':
#     parser = argparse.ArgumentParser()
#     parser.add_argument("-name", "--sample_name", default='cycads_report', required=False, help="prefix of output file name")
#     args = vars(parser.parse_args())
#     generate_fq_report_html(args)