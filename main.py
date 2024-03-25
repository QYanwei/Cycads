
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

from cycad import fq_datum
from cycad import fq_figure
from cycad import fq_report
from cycad import fq_filter
from cycad import bam_datum
from cycad import bam_figure
from cycad import bam_report

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
    version_file = open('%s/VERSION' % './VERSION')
    return version_file.readline().strip()

def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print_hi('PyCharm')
    version()


