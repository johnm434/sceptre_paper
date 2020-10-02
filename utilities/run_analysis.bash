#!/bin/bash

# This bash script runs the entire analysis. The code is commented to increase ease of adoption and use.
# Note: These additional packages must be installed for the entire analysis pipeline to work: stringr, monocle. These packages are not required to independently use the sceptre package.

# Set the machine.
machine=local

# Build and install the sceptre and katsevich2020 packages
# bash build_and_install_package.bash sceptre $machine
# bash build_and_install_package.bash katsevich2020 $machine

# Obtain the filepaths to the code and "offsite" directories
code_dir=$(bash get_file_paths.bash $machine code)
offsite_dir=$(bash get_file_paths.bash $machine data_results)

# Initialize the offsite directory structure
Rscript $code_dir"/analysis_drivers/"check_directory_structure.R $offsite_dir

# Download the data
Rscript $code_dir"/analysis_drivers/"download_raw_data.R $offsite_dir

# Pre-process the data
Rscript $code_dir"/analysis_drivers/"pre_process_data.R $offsite_dir
