
GJainLIB
====================================================

## Overview:
gjainLIB: A custom python module containing:
* Optimal python functions commonly used in my python scripts
* Mostly collected from stackoverflow and modified to improve performance

## Table of contents
<!--ts-->
1. [get_file_info](#get-file-info)
1. [create_dir](#create_dir)
1. [get_fileslist_in_dir](#get_fileslist_in_dir)
1. [print_initial_arguments](#print_initial_arguments)
1. [useful_lines](#useful_lines)
1. [get_unique_list](#get_unique_list)
1. [get_line_counts](#get_line_counts)
1. [get_line_counts_pygen](#get_line_counts_pygen)
1. [get_temp_file](#get_temp_file)
1. [Log](#Log)


<!--te-->

## Description
The module contains functions and classes. They are:

1. #### ``get_file_info():``  
	* Get path, basename, ext, path+basename and basename+ext of a file
	* Input parameters:
		* ``file_name_with_path: str``  
			Filename with or without absolute/relative path
	* Returns: 
		* input_file_path
		* input_file_basename
		* input_file_extension
		* path_filename
		* input_file_name
	* [top ▴](#Table-of-contents)

1. #### ``create_dir():`` 
	* Creates a directory if not exists
	* Input parameters:
		* ``dirname: str``  
			Folder/directory name with or without absolute/relative path
	* Returns: 
		* Relevant status messages on error
	* [top ▴](#Table-of-contents)

1. #### ``get_fileslist_in_dir():``  
	* Get the list of files inside a directory 
	* If extension is passed then return files for that extension
	* Input parameters:
		* ``dirname: str``  
			Folder/directory name with or without absolute/relative path
		* ``exts: str, list, optional``
			List of extensions to filter out for filetype
	* Returns: 
		* specific files if extensions are passed
		* all files if extensions are not passed
	* [top ▴](#Table-of-contents)

1. #### ``print_initial_arguments():``
	* Print the script name and all the initial arguments entered in the script
	* Input parameters:
		* ``parser: object``  
			Instance of ``ArgumentParser`` object
	* Returns:
		* None
	* [top ▴](#Table-of-contents)

1. #### ``useful_lines():`` 
	* Skips all the comments (right now #) from the file and returns the generator object of lines without comments
	* It remembers the state where it left. So it will start after the last executed line
	* Input parameters:
		* ``fh: file``  
			File handle of the file
	* Returns: 
		* generator for the last read line in the file
	* [top ▴](#Table-of-contents)

1. #### ``get_unique_list():`` 
	* Get unique element in a list with the order of elements preserved
	* Input parameters:
		* ``seq: list``  
			File handle of the file
		* ``idfun: function, optional``  
			A function that can transform the list values
	* Returns: 
		* generator for the unique and ordered preserved list
	* [top ▴](#Table-of-contents)

1. #### ``get_line_counts():``  
	* Get the number of lines in a file efficiently
	* Implemented with generators expressions in-line using itertools
	* Input parameters:
		* ``filename: str``  
			Filename with or without absolute/relative path
	* Returns: 
		* Total number of lines in the file
	* [top ▴](#Table-of-contents)

1. #### ``get_line_counts_pygen():``
	* Get the number of lines in a file efficiently
	* Implemented with generators expressions without using itertools
	* Works for files as big as 100+ Gb
	* Input parameters:
		* ``filename: str``  
			Filename with or without absolute/relative path
	* Returns: 
		* Total number of lines in the file
	* [top ▴](#Table-of-contents)

1. #### ``get_temp_file():``
	* Create a temporary file
	* The file will be deleted at program exit
	* Input parameters:
		* ``suffix: str, optional``  
			Filename with or without absolute/relative path
		* ``prefix: str, optional``
			Temporary file extension
		* ``dir: str, optional``
			Folder/directory name with or without absolute/relative path where the temporary file will be created
	* Returns: 
		* Filename with absolute path
	* [top ▴](#Table-of-contents)

1. #### ``Log():``  
	* A class for saving the STDOUT output in a log file 
	* Input parameters:
		* ``logf: file``  
			File handle of the log file
	* Returns: 
		* Write STDOUT in the log file
		* Return None
	* [top ▴](#Table-of-contents)





