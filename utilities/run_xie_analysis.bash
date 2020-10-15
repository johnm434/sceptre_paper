#!/bin/bash

# This script runs the Xie data analysis. It is assumed that this script is being executed from within the utilities directory.

# Set the machine.
machine=local

# Obtain the filepaths to the code and "offsite" directories
code_dir=$(bash get_file_paths.bash $machine code)
offsite_dir=$(bash get_file_paths.bash $machine data_results)

echo Check the availability of the required packages
Rscript $code_dir"/analysis_drivers_xie/"check_packages_0.R $code_dir

echo Initialize the offsite directory structure.
# Rscript $code_dir"/analysis_drivers_xie/"check_directory_structure_1.R $code_dir $offsite_dir

echo Download the data. Note that one of the downloads must be done manually.
# Rscript $code_dir"/analysis_drivers_xie/"download_data_2.R $code_dir $offsite_dir

echo Pre-process the data.
# Rscript $code_dir"/analysis_drivers_xie/"pre_process_data_3.R $code_dir $offsite_dir

echo Construct model covariate matrix and perform quality control.
# Rscript $code_dir"/analysis_drivers_xie/"build_model_covariate_matrix_4.R $code_dir $offsite_dir
