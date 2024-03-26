#!/usr/bin/env python3
# -*- coding: utf-8 -*-


# Basic library imports
import json
import datetime
import os

#import plotly.express as px
import pandas as pd

# Third party imports
#import plotly.offline as py
import jinja2

def plotting(config_dict):
    plots = list()
    titles = list()
    for method_name, method_args in config_dict.items():
        try:
            # Store plot title for HTML title and remove from data passed to plotly
            plot_title = method_args["plot_title"]
            method_args["plot_title"] = ""
            
            # Get method and generate plot
            method = getattr(plotter, method_name)
            fig = method(**method_args)
            plot = py.plot(
                fig,
                output_type='div',
                include_plotlyjs=False,
                image_width='',
                image_height='',
                show_link=False,
                auto_open=False)
            
            plots.append(plot)
            titles.append(plot_title)
        except AttributeError as E:
            print("\t\t{} is not a valid plotting method".format(method_name))
    return plots, titles
# Set a subtitle for the HTML report

# Define source files list
src_files = ""
fq_json_file = "../config/fq_config.json"
bam_json_file = "../config/bam_config.json"
for json_file_list, name in ((bam_json_file, "Alignments"), (fq_json_file, "Sequences")):
    if json_file_list:
        src_files += "<h4>Source {} files</h4><ul>".format(name)
        for f in json_file_list:
            f = os.path.abspath(f)
            src_files += "<li>{}</li>".format(f)
        src_files += "</ul>"

package_name=""
package_version=""
report_subtitle="Generated on {} with {} {}".format( datetime.datetime.now().strftime("%d/%m/%y"), package_name, package_version)

report_title="Cycads report"

# Load HTML template for Jinja
template_file="../config/template.html"
with open(template_file) as fp:
    template = jinja2.Template(fp.read())

# Render plots, Rendering plots in d3js
config_dict = {}
plots, titles=plotting(config_dict)
rendering = template.render(
    plots=plots,
    titles=titles,
    plotlyjs=py.get_plotlyjs(),
    report_title=report_title,
    report_subtitle=report_subtitle,
    src_files=src_files
)

# Write to HTML file
outfile='../test/cycads_report.html'
print("Writing to HTML file")

with open(outfile, "w") as fp:
    fp.write(rendering)
