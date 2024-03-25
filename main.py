
#!/usr/bin/env python3
# coding: utf-8

# Copyright (C) 2021, Yanwei Qi.
# qiyanwei1@genomics.cn, or qiyanweii@icloud.com.

# Cycads is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# Cycads is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
@name: Cycads long reads quality analyser
@description: All-in-one parsing fasta format file
@author: Yanwei Qi
@email: qiyanwei1@genomics.cn
@last modified by: qiyanwei

change log:
    2024/03/25  Coding completed.
    2024/02/18  Sketch finished.
    2024/01/18  Project created.
"""

import sys
import argparse

# from cycad import fq_datum
# from cycad import fq_figure
# from cycad import fq_report
# from cycad import fq_filter
# from cycad import bam_datum
# from cycad import bam_figure
# from cycad import bam_report

# basic stent function
## check file
def file_checkin(infile, func):
    if os.path.exists(infile):
        print(func + ': check...' + infile + ' file existed.')
        return 1
    else:
        print(func + ': check...' + infile + ' file not found.')
        return 0
##/

def version():
    version_file = open('%sVERSION.txt' % './')
    return version_file.readline().strip()

def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.

def print_helpdoc():
    help_message = """
              ........:::=== Cycads V%s ===:::........
    ============================================================================
                Quality control & Data filtering & Error analysis
                             for Long-read sequencing
    ============================================================================
    """ % version()
    print(help_message)





# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    ############################################## initialize subparsers ###############################################

    # initialize the options parser
    parser = argparse.ArgumentParser(description='Cycads:Long reads sequencing quality analyzer')
    parser.add_argument("-i",    "--fastq",     required=True,  help="sequences data.fq/fq.gz")
    parser.add_argument("-bam",  "--alignment", required=False, help="alignment data.bam")
    parser.add_argument("-ref",  "--reference", required=False, help="alignment data.fasta")
    parser.add_argument("-list", "--fqlist",    required=False, help="fastq path list.txt")
    parser.add_argument("-name", "--sample_name", default='cycad', required=False, help="prefix of output file name")
    
    # for quality control.

    parser.add_argument("-P", "--platform", required=False, help="cyclone")
    parser.add_argument("-K", "--kmersize", type=int, required=False, help="5")
    parser.add_argument("-BQH", "--base_quality_head", type=int, default=200, required=False, help="200")
    parser.add_argument("-BQT", "--base_quality_tail", type=int, default=200, required=False, help="200")
    parser.add_argument("-BCH", "--base_content_head", type=int, default=200, required=False, help="200")
    parser.add_argument("-BCT", "--base_content_tail", type=int, default=200, required=False, help="200")
    parser.add_argument("-hpmin", "--homopolymer_min_length", type=int, required=False, help="2")
    parser.add_argument("-hpmax", "--homopolymer_max_length", type=int, required=False, help="16")
    # for data filtering.
    parser.add_argument("-filter", "--filtering", required=False, help="get clean data")
    parser.add_argument("-Qmin", "--minium_quality",  type=float, default='10',   required=False, help="filter low quality reads")
    parser.add_argument("-Lmin", "--minium_length",   type=int,   default='1000', required=False, help="filter short reads")
    parser.add_argument("-Lmax", "--maxium_length",   type=int,   default='1000', required=False, help="filter large reads")
    parser.add_argument("-Hcut", "--cut_head_length", type=int,   default='200',  required=False, help="trim head sequences")
    parser.add_argument("-Tcut", "--cut_tail_length", type=int,   default='200', required=False,  help="trim head sequences")
    parser.add_argument("-Htrim", "--trim_head_homopolymer", type=int, default='200', required=False, help="trim head sequences homopolymer")
    parser.add_argument("-Ttrim", "--trim_tail_homopolymer", type=int, default='200', required=False, help="trim head sequences homopolymer")
    parser.add_argument("-Dx", "--downsampling_dx", type=str, default='10', required=False, help="depth")
    parser.add_argument("-Dg", "--downsampling_gs", type=str, default='4m', required=False, help="genome size")
    
    # for error analysis.
    
    # parser.add_argument("-hpmin", "--homopolymer_min_length", type=int, required=False, help="2")
    # parser.add_argument("-hpmax", "--homopolymer_max_length", type=int, required=False, help="16")
    ############################## parse provided arguments and run corresponding function #############################
    
    # get and check options
    #    args = None
    args = parser.parse_args()
    if (len(sys.argv) == 1) or (sys.argv[1] == '-h') or (sys.argv[1] == '-help') or (sys.argv[1] == '--help'):
        print_helpdoc()
        sys.exit(0)
    else:
        args = vars(parser.parse_args())
    print(args)
    
    

