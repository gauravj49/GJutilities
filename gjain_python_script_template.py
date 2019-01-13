#!/usr/local/bin/python
"""
***********************************************
- PROGRAM: gjain_python_script_template.py
- CONTACT: Gaurav Jain (gaurav.jain@dzne.edu)
- LICENSE: Copyright (C) <year>  <name of author>
           This program is free software: you can redistribute it and/or modify
           it under the terms of the GNU General Public License as published by
           the Free Software Foundation, either version 3 of the License, or
           (at your option) any later version.

           This program is distributed in the hope that it will be useful,
           but WITHOUT ANY WARRANTY; without even the implied warranty of
           MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
           GNU General Public License for more details.

           You should have received a copy of the GNU General Public License
           along with this program.  If not, see <http://www.gnu.org/licenses/>.
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

# For looping files in a dir
import glob
from os import listdir
from os.path import isfile, join

# user defined modules
#from gjainLIB import *      # import all the functions from the Gaurav`s python library

################ USER CONFIGURATION ###################
#######################################################

def main():
    # Get input options
    args               = check_options()
    my_argument1       = args.my_argument1
    my_argument2       = args.my_argument2
    my_argument3       = args.my_argument3
    some_flag          = args.smflag
    my_argument4       = args.my_argument4
    
    # Your main function code

    
################ USER DEFINED FUNCTIONS ###################

def check_options():
    ''' Checks the options to the program '''

    # Create parser object
    parser = argparse.ArgumentParser(add_help=True,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent('''\
        ----------------- SAMPLE USAGE ------------------
        - python gjain_python_script_template.py -a1="my_first_argument" 
        - python gjain_python_script_template.py -a1="my_first_argument" -a2="my_second_argument" -a3=49
        - python gjain_python_script_template.py -a1="my_first_argument" -a2="my_second_argument" -a3=49 -fg
        - python gjain_python_script_template.py -a1="my_first_argument" -a2="my_second_argument" -a3=49 -fg -a4='r'
        -------------------------------------------------
        CONTACT: 
        	Gaurav Jain
        	gaurav.jain@dzne.de
        -------------------------------------------------
        '''))

    # Add arguments 
    parser.add_argument("-a1", metavar="--myarg1" , help="*My argument one"  , dest="my_argument1"  , type=str, required=True)
    parser.add_argument("-a2", metavar="--myarg2" , help=" My argument two"  , dest="my_argument2"  , type=str)
    parser.add_argument("-a3", metavar="--myarg3" , help=" My argument three", dest="my_argument3"  , type=int)
    parser.add_argument("-fg",         "--smflag" , help=" This is some flag. if set, perform some additional tasks", action='store_true', default=False)
    parser.add_argument('-a4',         "--selval" , help=" The argument that can have only certain options\n - Please enter:\n\t -sc='n' for NONE\n\t -sc='r' for ROW\n\t -sc='c' for COLUMN", dest="my_argument4", default='n', action='store', choices=['n', 'r', 'c'])
    
    # Print the help message only if no arguments are supplied
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    # Save the STDOUT output in a log file
    if parser.parse_args().my_argument1:
        logdir = "{0}/logs".format(get_file_info(parser.parse_args().my_argument1)[0])
        create_dir(logdir)
        logfile = "{0}/{1}.log".format(logdir, get_file_info(parser.parse_args().my_argument1)[1])
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
