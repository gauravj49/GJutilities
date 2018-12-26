
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
1. [Combined Volcano Plots](#draw_volcano_plotspy)


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


1. #### draw_volcano_plots.py
	```
	[gjain@gjvbx GJutilities ] $ python draw_combined_volcano_plots.py 

	***********************************************
	- PROGRAM: draw_volcano_plots.py
	- CONTACT: Gaurav Jain(gauravj49@gmail.com)
	***********************************************

	usage: draw_combined_volcano_plots.py [-h] -af --inpfl1 -bf --inpfl2 -of
										--opfile [-ix --idxnme] [-nf]

	optional arguments:
	-h, --help              Show this help message and exit
	-af --inpfl1          * Input data file1. Ex: deseq2 results file1
	-bf --inpfl2          * Input data file2. Ex: deseq2 results file2
	-of --opfile          * Output file name
	-ix --idxnme            Index column name
						    - (Default = GeneName)
	-nf, --filter_smallrna
							if set, filter any feature (based on low expression) from the  inputDF
							- (Default = False)

	----------------- SAMPLE USAGE ------------------
	- python draw_combined_volcano_plots.py -af=results_DEseq2/A_over_C_DE_RESULTS_common_gene_names.txt -bf=results_DEseq2/B_over_D_DE_RESULTS_common_gene_names.txt -of=diagnosis_plots/A_over_C_combined_with_B_over_D_volcano_plots.txt -ix=GeneName
	-------------------------------------------------
	CONTACT: 
		Gaurav Jain
		gauravj49@gmail.com
	-------------------------------------------------
	```

	* Example output plot
	  ![a_over_c_combined_with_b_over_d_volcano_plots_combined_volcano_plot](https://user-images.githubusercontent.com/10153240/50427639-e8d2bc00-08ad-11e9-983b-6ab146a019f9.png)  
	  
	* [top ▴](#table-of-contents)



