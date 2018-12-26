#!/usr/local/bin/python
"""
***********************************************
- PROGRAM: draw_volcano_plots.py
- CONTACT: Gaurav Jain(gauravj49@gmail.com)
***********************************************
"""
print (__doc__)

# Built in modules
import argparse
import os.path
import sys

# 3rd party modules
import textwrap
import re
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib as mp
#mp.use('Agg') # to use matplotlib without X11
import matplotlib.pyplot as plt
import subprocess
import binascii as bi
import scipy.stats as stats
from collections import *
from numpy import nanmean

# for looping files in a dir
import glob

# user defined modules
from gjainLIB import *      # import all the functions from the Gaurav`s python library

### for color scale
from  matplotlib import colors
from itertools import cycle, islice # barplot colors
import seaborn as sns

# For clustering
from scipy import stats
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as dist
import numpy.ma as ma
import time

################ USER CONFIGURATION ###################
matplotlib.rcParams['pdf.fonttype'] = 'truetype'
matplotlib.rcParams['pdf.fonttype'] = 42
pd.set_option('display.width', 1000)
plt.style.use("classic")
#plt.rcParams['axes.edgecolor'] = "#777777"
#######################################################

def main():
    # Get input options
    args                 = check_options()
    input_data_file1     = args.input_data_file1
    input_data_file2     = args.input_data_file2
    output_file          = args.output_file
    indexName            = args.indexName
    filter_smallrna      = args.filter_smallrna     # if set, filter out low expressed smallrnas

    # Get the input data in a dataframe
    inputDF1 = pd.read_csv(input_data_file1, comment='#', delimiter="\t", index_col=indexName, engine='python')
    inputDF2 = pd.read_csv(input_data_file2, comment='#', delimiter="\t", index_col=indexName, engine='python')

    # Get index names in a list
    inputIndexNames1 = inputDF1.index.tolist()
    inputIndexNames2 = inputDF2.index.tolist()

    # Get the base colors for all the points
    inputDF1['colors'] = ["black" for i in xrange(inputDF1.shape[0])]
    inputDF2['colors'] = ["red" for i in xrange(inputDF2.shape[0])]

     # Get the lables
    label1 = "{0}".format(get_file_info(input_data_file1)[1].replace("_DE_RESULTS_common_gene_names", "").replace("_", " "))
    label2 = "{0}".format(get_file_info(input_data_file2)[1].replace("_DE_RESULTS_common_gene_names", "").replace("_", " "))

    # Get data for x-axis
    x1 = inputDF1['log2FoldChange']
    x2 = inputDF2['log2FoldChange']

    # Get data for y-axis
    y1 = -np.log10(inputDF1['pvalue'])
    y2 = -np.log10(inputDF2['pvalue'])

    ymin = min(min(y1), min(y2))
    ymax = max(max(y1), max(y2))

    ylim  = [60, ymax+10]
    ylim2 = [-2, 10]
    ylimratio  = (ylim[1]-ylim[0])/(ylim2[1]-ylim2[0]+ylim[1]-ylim[0])
    ylim2ratio = (ylim2[1]-ylim2[0])/(ylim2[1]-ylim2[0]+ylim[1]-ylim[0])

    # Help: http://matplotlib.org/examples/pylab_examples/broken_axis.html
    # Let's 'break' or 'cut-out' the y-axis into two portions and use the:
    # 	- top    (axtop) for the outliers, 
    # 	- bottom (axbot) for the details of the majority of our data
    import matplotlib.gridspec as gridspec
    gs  = gridspec.GridSpec(2, 1, height_ratios=[ylim2ratio, ylimratio])
    fig = plt.figure()
    ax  = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])

    # Plot the same volcano plot on both axes
    # Plot bottom plots first so that legend is plotted at the top
    ax.scatter( x1,y1, c=inputDF1['colors'], edgecolor='k', linewidth=0.2, label=label1)
    ax.scatter( x2,y2, c=inputDF2['colors'], edgecolor='k', linewidth=0.2, label=label2)
    ax2.scatter(x1,y1, c=inputDF1['colors'], edgecolor='k', linewidth=0.2, label=label1)
    ax2.scatter(x2,y2, c=inputDF2['colors'], edgecolor='k', linewidth=0.2, label=label2)
    ax.set_ylim(ylim)
    ax2.set_ylim(ylim2)
    plt.subplots_adjust(hspace=0.03)
    
    ax.spines['bottom'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax.xaxis.tick_top()
    ax.tick_params(labeltop='off')
    ax2.xaxis.tick_bottom()
    
    # Add the title
    plt.suptitle("{0}".format(get_file_info(output_file)[1].replace("_", " ")))

    ax2.set_xlabel('Log2FoldChange')
    ax2.set_ylabel('-log10(pvalue)')
    ax2.yaxis.set_label_coords(0.05, 0.5, transform=fig.transFigure)
    
    kwargs = dict(color='k', clip_on=False)
    xlim = ax.get_xlim()
    xlim = [-1.2, 1.0]
    dx = .02*(xlim[1]-xlim[0])
    dy = .01*(ylim[1]-ylim[0])/ylim2ratio
    ax.plot((xlim[0]-dx,xlim[0]+dx), (ylim[0]-dy,ylim[0]+dy), **kwargs)
    ax.plot((xlim[1]-dx,xlim[1]+dx), (ylim[0]-dy,ylim[0]+dy), **kwargs)
    dy = .01*(ylim2[1]-ylim2[0])/ylimratio
    ax2.plot((xlim[0]-dx,xlim[0]+dx), (ylim2[1]-dy,ylim2[1]+dy), **kwargs)
    ax2.plot((xlim[1]-dx,xlim[1]+dx), (ylim2[1]-dy,ylim2[1]+dy), **kwargs)
    ax.set_xlim(xlim)
    ax2.set_xlim(xlim)
    ax.legend(loc='upper center', ncol=2, fontsize = 'x-small')

    # Add the title and legend
    plt.suptitle("{0}".format(get_file_info(output_file)[1].replace("_", " ")))

    # Get output pdf dir
    pdfdir = "{0}/pdf".format(get_file_info(output_file)[0])
    create_dir(pdfdir)

    # Save the plots
    print "- Saving heatmap in output file:\n\t- {0}_heatmap.png\n\t- {1}/{2}_heatmap.pdf".format(get_file_info(output_file)[3], pdfdir, get_file_info(output_file)[1])
    plt.savefig("{0}_combined_volcano_plot.png"    .format(get_file_info(output_file)[3]), dpi=300, bbox_inches='tight')
    plt.savefig("{0}/{1}_combined_volcano_plot.pdf".format(pdfdir, get_file_info(output_file)[1]) , bbox_inches='tight')
    plt.close()
    
################ USER DEFINED FUNCTIONS ###################
def check_options():
    ''' Checks the options to the program '''

    # Create parser object
    parser = argparse.ArgumentParser(add_help=True,
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=textwrap.dedent('''\
        ----------------- SAMPLE USAGE ------------------
        - python draw_combined_volcano_plots.py -af=results_DEseq2/A_over_C_DE_RESULTS_common_gene_names.txt -bf=results_DEseq2/B_over_D_DE_RESULTS_common_gene_names.txt -of=diagnosis_plots/A_over_C_combined_with_B_over_D_volcano_plots.txt -ix=GeneName
        -------------------------------------------------
        CONTACT: 
        	Gaurav Jain
        	gauravj49@gmail.com
        -------------------------------------------------
        '''))

    # Add arguments 
    parser.add_argument("-af", metavar="--inpfl1", help="*Input data file1. Ex: deseq2 results file1", dest="input_data_file1", type=str, required=True)
    parser.add_argument("-bf", metavar="--inpfl2", help="*Input data file2. Ex: deseq2 results file2", dest="input_data_file2", type=str, required=True)
    parser.add_argument("-of", metavar="--opfile", help="*Output file name" , dest="output_file", type=str, required=True)
    parser.add_argument("-ix", metavar="--idxnme", help=" Index column name\n - (Default = GeneName)", dest="indexName", type=str, default="GeneName")
    parser.add_argument('-nf', "--filter_smallrna", help=" if set, filter any feature (based on low expression) from the  inputDF\n - (Default = False)", action='store_true', default=False)

    # Print the help message only if no arguments are supplied
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    # Save the STDOUT output in a log file
    if parser.parse_args().output_file:
        logdir = "{0}/logs".format(get_file_info(parser.parse_args().output_file)[0])
        create_dir(logdir)
        logfile = "{0}/{1}.log".format(logdir, get_file_info(parser.parse_args().output_file)[1])
    else:
        logdir  = "{0}/logs".format(os.getcwd())
        create_dir(logdir)
        logfile = "{0}/{1}.log".format(logdir,get_file_info(sys.argv[0])[1])

    logf = open(logfile, 'w')
    sys.stdout = Log(logf, sys.stdout)

    # Parse command line with parse_args and store it in an object
    args = parser.parse_args()
    print_initial_arguments(parser)
    return args

# main function
if __name__=="__main__":
      main()
