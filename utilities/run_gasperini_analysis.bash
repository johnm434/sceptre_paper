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
# Rscript $code_dir"/analysis_drivers/"check_directory_structure_1.R $code_dir $offsite_dir

# Download the data
# Rscript $code_dir"/analysis_drivers/"download_raw_data_2.R $code_dir $offsite_dir

echo Pre-process the data.
# Rscript $code_dir"/analysis_drivers/"pre_process_data_3.R $code_dir $offsite_dir

echo Perform quality control.
# Rscript $code_dir"/analysis_drivers/"perform_quality_control_4.R

echo Create the precomputation and results dictionaries.
pod_sizes=$(Rscript $code_dir"/analysis_drivers/"create_file_dictionaries_5.R $code_dir $offsite_dir)
n_gene_pods="$(echo $pod_sizes | cut -d' ' -f1)"
n_gRNA_pods="$(echo $pod_sizes | cut -d' ' -f2)"
n_pair_pods="$(echo $pod_sizes | cut -d' ' -f3)"
echo number of gene pods: $n_gene_pods
echo number of gRNA pods: $n_gRNA_pods
echo number of gRNA-gene pair pods: $n_pair_pods

echo Run gene precomputation.

echo Run gRNA precomputation.
