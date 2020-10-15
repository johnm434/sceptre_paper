#!/bin/bash

# This script runs the Xie data analysis. It is assumed that this script is being executed from within the utilities directory.
# Note: the R packages purr, stringr, rhdf5, openxlsx, and ravel must be installed for this analysis to work.

# Set the machine.
machine=local

echo Build and install the sceptre package.
# bash build_and_install_package.bash sceptre $machine

# Obtain the filepaths to the code and "offsite" directories
code_dir=$(bash get_file_paths.bash $machine code)
offsite_dir=$(bash get_file_paths.bash $machine data_results)

echo Initialize the offsite directory structure.
# Rscript $code_dir"/analysis_drivers_xie/"check_directory_structure_1.R $code_dir $offsite_dir

echo Download the data.
# Rscript $code_dir"/analysis_drivers_xie/"download_data_2.R $code_dir $offsite_dir

echo pre-process the data
# Rscript $code_dir"/analysis_drivers_xie/"pre_process_data_3.R $code_dir $offsite_dir

echo Construct model covariate matrix and perform quality control.
# Rscript $code_dir"/analysis_drivers_xie/"build_model_covariate_matrix_4.R $code_dir $offsite_dir
