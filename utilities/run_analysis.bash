#!/bin/bash

# This bash script runs the entire analysis. The code is commented to increase ease of adoption and use. 

# Set the machine.
machine=local

# Build and install the sceptre and katsevich2020 packages
# bash build_and_install_package.bash sceptre $machine
# bash build_and_install_package.bash katsevich2020 $machine

# Obtain the filepaths to the code and "offsite" directories
code_dir=$(bash get_file_paths.bash $machine code)
offsite_dir=$(bash get_file_paths.bash $machine data_results)

# Initialize the offsite directory structure
Rscript check_directory_structure.R $offsite_dir

# Download the data
