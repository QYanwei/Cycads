import os, re, sys, time
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
    table_list = ["TB2000B609-202403200954240_read.bam","2399","3939982","95.02","3.98","1.53","1.32","1.13","2.3","1.68"]
    table_string = ""
    for i in table_list:
        table_string += "<th>{}</th>".format(i)
    
    plots_list = [["query_all_error_item.barplot.png",
                   "query_insertion_frequency.barplot.png"],
                  ["query_deletion_frequency.barplot.png",
                   "query_all_substitution_errors.barplot.png"],
                  ["query_events_curve_idy.dispplot.png",
                   "query_events_curve_dif.dispplot.png"],
                  ["query_homopolymer_length_event.lineplot.png",""]]
    plots_src = [["Sequencing error rate",
                  "Insertion frequency"],
                 ["Deletion frequency",
                  "Substitution type error frequency"],
                 ["Read identity rate distribution",
                  "Read error rate distribution"],
                 ["Homopolymer error type with its length", ""]]
    
    plots_string = ""
    for i in range(4):
        plots_string += "<tr>"
        for j in range(2):
            plots_string += "<td><figure>"
            plots_string += "<img src=\"./{}\" width=\"600\" height=\"580\" />".format(plots_list[i][j])
            plots_string += "<figcaption class=\"figure-caption text-center\">{}</figcaption>".format(plots_src[i][j])
            plots_string += "</td>"
        plots_string += "</tr>"
    report_title = "CycloneSEQ quality reporter"
    report_subtitle = "Created on {} with {} {}".format(datetime.datetime.now().strftime("%d/%m/%y"), "Cycads",
                                                          "0.3.0")
    pwd_config_file = os.path.realpath(__file__)
    template_file = '/'.join(pwd_config_file.split('/')[:-1]) + "../config/template_bam_report.j2"
    template = import_jinja_template(template_file)
    rendering = template.render(
        table=table_string,
        plots=plots_string,
        report_title=report_title,
        report_subtitle=report_subtitle)
    
    # Write to HTML file
    outfile = "../" + args["sample_name"] + "/" + args["sample_name"] + "_error_report.html"
    with open(outfile, "w") as fp:
        fp.write(rendering)


args = {"sample_name": "test"}
generate_fq_report_html(args)
# if __name__ == '__main__':
#     parser = argparse.ArgumentParser()
#     parser.add_argument("-name", "--sample_name", default='cycads_report', required=False, help="prefix of output file name")
#     args = vars(parser.parse_args())
#     generate_fq_report_html(args)