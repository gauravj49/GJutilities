
GJ-utilities
====================================================

## Overview:
GJutilities is a collection of scripts useful in the analysis and management of:
* Initial sequencing data from Illumina high throughput NGS machines
* Wrappers for running and parsing commonly used command line tools like bedtools, samtools etc. 
* Various data visualization techniques

## Table of contents
<!--ts-->
1. [Merge Illumina Samples](#merge_illumina_samplespy)



<!--te-->

## Description
For various scripts:

1. #### merge_illumina_samples.py 
	```
	[gjain@gjvbx GJutilities ] $ python merge_illumina_samples.py 

	***********************************************
	- PROGRAM: merge_illumina_samples.py
	- CONTACT: Gaurav Jain(gauravj49@gmail.com)
	***********************************************

	usage: merge_illumina_samples.py [-h] -id --indir -od --otdir [-du]

	optional arguments:
	-h, --help    show this help message and exit
	-id --indir   *Input directory containing the unmerged fastQ files
	-od --otdir   *Output directory to save the merged fastQ files
	-du, --dudir  if set, delete unmerged fastq files dir 
					(Default: False)

	----------------- SAMPLE USAGE ------------------
	- python scripts/merge_illumina_samples.py -id=input/raw_fastq/unmerged -od=input/raw_fastq/merged
	- python scripts/merge_illumina_samples.py -id=input/raw_fastq/unmerged -od=input/raw_fastq/merged -du
	-------------------------------------------------
	CONTACT: 
		Gaurav Jain
		gauravj49@gmail.com
	-------------------------------------------------
	```
	* [top ▴](#table-of-contents)


