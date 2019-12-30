#!/usr/local/bin/python
"""
***********************************************
- PROGRAM: demultiplex_fastq_py2.py
- CONTACT: Gaurav Jain(gauravj49@gmail.com)

- DECRIPTION: The program mainly demultiplexes a fastq file if it was not demultiplexed before by using the indexes provided in a file. 
***********************************************
"""
print (__doc__)

# Built in modules
import argparse
import os.path
import sys
import io

# 3rd party modules
import re
import csv
import math
import itertools
import operator
import gzip
import textwrap
import glob
from collections import *

################ USER CONFIGURATION ###################
#######################################################

def main():
    # Get input options
    args               = check_options()
    input_sample_sheet = args.input_sample_sheet
    input_fastq_file   = args.input_fastq_file
    output_dir         = args.output_dir
    index_mismatches   = args.index_mismatches
    
    ##########################################

    # Get the output directory
    if output_dir == None:
        output_dir = "{0}/demultiplexed".format(get_file_info(input_fastq_file)[0])
    create_dir(output_dir)

    # Get the dictionary of indexes with sample names
    samidx_dict = get_sample_index(input_sample_sheet)

    # Get the dictionary of open filehandles for output files
    flhidx_dict = get_open_fh_index(samidx_dict, output_dir)
    undfh       = open("{0}/undetermined.fastq".format(output_dir), 'w')


    # Demultiplex the fastq file using the indexes and save them in separate files
    for seqID, seq, qualityID, quality in read_fastq(input_fastq_file):
        index = seqID.split(":")[-1]
        if index_mismatches == 0:
            if index in samidx_dict.iterkeys():
                flhidx_dict[index].write("{0}\n".format("\n".join([seqID, seq, qualityID, quality])))
        else:
            idxhit = 0  # Total number of hits for each index
            hitidx = "" # Index with only one hit
            for ssindex in samidx_dict.iterkeys():
                if index_mismatches >= hamdist(index, ssindex):
                    hitidx  = ssindex
                    idxhit += 1
            if idxhit == 1:
                # This is used to avoid assiging a read with same hamming distance for indexes
                flhidx_dict[hitidx].write("{0}\n".format(
                    "\n".join([seqID, seq, qualityID, quality])))
            else:
                undfh.write("{0}\n".format("\n".join([seqID, seq, qualityID, quality])))

    # Close all the open file handles
    close_fh_index(flhidx_dict)

    # Gzip all the demultiplexed files
    for f in glob.glob("{0}/*.fastq".format(output_dir)):
        os.system("gzip {0}".format(f))

    # Print output files
    print "Final output files are:"
    for f in glob.glob("{0}/*.gz".format(output_dir)):
        print "\t- {0}".format(f)

################ USER DEFINED FUNCTIONS ###################
def get_sample_index(input_sample_sheet):
    ''' Read a csv file with sample name as value and indexes as key'''
    # Inititalize the dictionary
    d = defaultdict()
    # Get the information in a dictionary
    with open(input_sample_sheet,'rU') as fg:
        ssreader = csv.reader(fg)
        for r in ssreader:
            d[r[1]] = r[0]
    return d

def read_fastq(input_fastq_file):
    ''' 
        Reads the fastq file and returns the name and seq 
        Source: https://stackoverflow.com/questions/39150965/fastest-way-to-read-a-fastq-with-scikit-bio
    '''
    with open_file_handle(input_fastq_file, 'r') as fg:
        # # Input fastq file
        # 1              2   3          4 5    6     7    8 9 10 11
        # @HWI-BRUNOP20X:765:A419SRABXX:1:1251:11801:6250 1:N:0:CCGCGGTT
        # GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
        # +
        # !''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65

        # 1 : Instrument name
        # 2 : Run number
        # 3 : Flowcell id
        # 4 : Lane number
        # 5 : Tile number
        # 6 : x_position
        # 7 : y_position
        # 8 : Read
        # 9 : is_filtered_flag
        # 10: Control number
        # 11: Index
        for seqID, seq, qualityID, quality in itertools.izip_longest(*[fg] * 4):
            if any(line is None for line in (seqID, seq, qualityID, quality)):
                raise Exception("Invalid fastq file with number of lines does not match to the multiple of four.")
            if not seqID.startswith('@'):
                raise Exception("Invalid seqID:{0}".format(seqID))
            if qualityID != '+\n':
                raise Exception("Invalid qualityID:{0}".format(qualityID))
            if quality == '\n':
                raise Exception("Fastq record is missing quality scores.")

            yield (seqID.strip('\n'), seq.strip('\n'), qualityID.strip('\n'), quality.strip('\n'))

def open_file_handle(input_file, mode='rU'):
    ''' Opens a file depending on whether it's gzipped or not '''
    if get_file_info(input_file)[2] == 'gz':
        return gzip.open(input_file, mode)
    else:
        return open(input_file, mode)

def hamdist(str1, str2):
    '''
        Calculate the Hamming distance (or number of differences) between two strings of the same length
        Source: http://code.activestate.com/recipes/499304-hamming-distance/
    '''
    assert len(str1) == len(str2)
    ne = operator.ne
    return sum(imap(ne, str1, str2))

def get_open_fh_index(inputdict, output_dir):
    ''' Return a dictionary of open file handles of indxes'''
    # Inititalize the dictionary
    d = defaultdict()
    # Get the information in a dictionary
    for index, sample_name in inputdict.iteritems():
        output_file = "{0}/{1}.fastq".format(output_dir,sample_name)
        d[index] = open(output_file, 'w')
    return d

def close_fh_index(inputdict):
    ''' Closes open file handles of indxes'''
    for fh in inputdict.itervalues():
        fh.close()

def get_file_info(file_name_with_path):
    ''' Get path, basename, ext, path+basename and basename+ext of a file '''

    # get the path of the file
    input_file_path = os.path.dirname(os.path.abspath(file_name_with_path))

    # get the filename without path but with extension
    input_file_name = os.path.basename(file_name_with_path)

    # get basename of the file (no path or extension)
    input_file_basename = os.path.splitext(input_file_name)[0]

    # get extension of the file
    input_file_extension = os.path.splitext(input_file_name)[1]
    # remove "."
    input_file_extension = re.sub(r"^\.","",input_file_extension)

    # get the path+filename
    path_filename = input_file_path + "/" + input_file_basename

    return (input_file_path, input_file_basename, input_file_extension, path_filename, input_file_name)

def create_dir(dirname):
    '''Creates a directory if not exists. Returns relevant status messages on error'''
    try:
        os.makedirs(dirname)
    except OSError:
        if os.path.exists(dirname):
            # We are nearly safe
            pass
        else:
            # There was an error on creation, so make sure we know about it
            raise

def useful_lines(fh):
    '''Skips all the comments (right now #) from the file and returns the generator object of lines without comments.
    It remembers the state where it left. So it will after the last executed line '''
    for line in (l.strip() for l in fh if not l.startswith('#')):
        yield line

def print_initial_arguments(parser):
    ''' Prints all the initial arguments entered '''

    print "\n------Input Arguments ------"
    opts = vars(parser.parse_args())
    maxl = max(len(text) for text in opts.keys())
    for k,v in opts.items():
        print "%-*s = %s" % (maxl, k, v)
    print "-"*29, "\n"

# class for Logging
class Log(object):
    def __init__(self, *files):
        self.files = files

    def write(self, obj):
        for f in self.files:
            f.write(obj)

    def flush(self):
        pass

def check_options():
    ''' Checks the options to the program '''

    # Create parser object
    parser = argparse.ArgumentParser(add_help=True,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent('''\
        ----------------- SAMPLE USAGE ------------------
        - python demultiplex_fastq_py2.py -if=demultiplex.fastq.gz -is=samples.csv
        - python demultiplex_fastq_py2.py -if=input/fastq/experimentX_sample1.fastq.gz -is=input/sampleSheets/experimentX.csv -od=output/demultiplexed/experimentX
        -------------------------------------------------
        CONTACT: 
        	Gaurav Jain
        	gaurvj49@gmail.com
        -------------------------------------------------
        '''))

    # Add arguments 
    parser.add_argument("-is", metavar="--insamf", help="*Input sample sheet file" , dest="input_sample_sheet" , type=str, required=True)
    parser.add_argument("-if", metavar="--infile", help="*Input fastq file"        , dest="input_fastq_file"   , type=str, required=True)
    parser.add_argument("-od", metavar="--optdir", help=" Output directory"        , dest="output_dir"         , type=str)
    parser.add_argument("-im", metavar="--idxmis", help=" Number of mismatches in the index (Default: 0)"      , dest="index_mismatches" , type=int, default = 0)
    
    # Print the help message only if no arguments are supplied
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    # Save the STDOUT output in a log file
    if parser.parse_args().output_dir != None:
        logdir = "{0}/logs".format(parser.parse_args().output_dir)
        create_dir(logdir)
        logfile = "{0}/{1}.log".format(logdir, get_file_info(parser.parse_args().input_fastq_file)[1])
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
