#!/usr/local/bin/python
"""
***********************************************
- PROGRAM: merge_illumina_samples.py
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
import subprocess
import math
from collections import *

# User defined modules
# Import all the functions from the Gaurav`s python library
from gjainLIB import *      

# For looping files in a dir
import glob
from os import listdir
from os.path import isfile, join

################ USER CONFIGURATION ###################
#######################################################

def main():
    # Get input options
    args             = check_options()
    input_dir        = args.input_dir.rstrip('\/')
    output_dir       = args.output_dir
    del_unmerged_dir = args.del_unmerged_dir

    # Create output directory if not present
    create_dir(output_dir)

    # Get unmerged fastQ file names
    unmerged_files_list = [ f for f in listdir(input_dir) if isfile(join(input_dir,f))]
    #unmerged_files_list = list(glob.iglob(input_dir))
    
    # Get the dictionary of merged file names
    # Exp1_hs_smallrna_sr_Gaurav_D_6_AGTCAA_L008_R1_001.fastq.gz
    # Exp1_hs_smallrna_sr_Gaurav_D_6_AGTCAA_L008_R1_002.fastq.gz
    # Exp1_hs_smallrna_sr_Gaurav_D_6_AGTCAA_L008_R1_003.fastq.gz
    merged_names = defaultdict(list)
    for f in unmerged_files_list:
        merged_names["_".join(f.split('_')[:-4])].append(f)
        
    # Cat all the samples with the same sample name
    # Exp1_hs_smallrna_sr_Gaurav_D_6.fastq.gz
    i=1
    for k,v in merged_names.iteritems():
        print "\n{0}) {1}".format(i,k)
        for subf in v:
            print "\t-{0}".format(subf)
        if len(v) > 1:
            cmd = "cat {0}/{1}* > {2}/{1}.fastq.gz".format(input_dir, k, output_dir)
        else:
            cmd = "mv {3}/{0} {2}/{1}.fastq.gz".format(v[0], k, output_dir, input_dir)
        
        # Run the command
        os.system(cmd)
        i += 1

    # if del_unmerged_dir is set, then delete the directory
    if del_unmerged_dir:
        print "\n- Deleting \"{0}\" directory\n".format(input_dir)
        os.system("rm -rf {0}".format(input_dir))

################ USER DEFINED FUNCTIONS ###################
def read_fasta(fp):
    ''' Reads the fasta file and returns the name and seq '''
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))

def check_options():
    ''' Checks the options to the program '''

    # Create parser object
    parser = argparse.ArgumentParser(add_help=True,
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=textwrap.dedent('''\
        ----------------- SAMPLE USAGE ------------------
        - python scripts/merge_illumina_samples.py -id=input/raw_fastq/unmerged -od=input/raw_fastq/merged
        - python scripts/merge_illumina_samples.py -id=input/raw_fastq/unmerged -od=input/raw_fastq/merged -du
        -------------------------------------------------
        CONTACT: 
        	Gaurav Jain
        	gauravj49@gmail.com
        -------------------------------------------------
        '''))

    # Add arguments 
    parser.add_argument("-id", metavar="--indir" , help="*Input directory containing the unmerged fastQ files", dest="input_dir", type=str, required=True)
    parser.add_argument("-od", metavar="--otdir" , help="*Output directory to save the merged fastQ files", dest="output_dir", type=str, required=True)
    parser.add_argument("-du",         "--dudir" , help="if set, delete unmerged fastq files dir \n(Default: False)", action='store_true', dest="del_unmerged_dir", default=False)
    
    # Print the help message only if no arguments are supplied
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    # Save the STDOUT output in a log file
    if parser.parse_args().output_dir:
        logdir = "{0}/logs".format(parser.parse_args().output_dir)
        create_dir(logdir)
        logfile = "{0}/merged_samples{1}.log".format(logdir, get_file_info(parser.parse_args().output_dir)[1])
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
