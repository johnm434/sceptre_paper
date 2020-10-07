#!/bin/bash

# This bash script runs the Gasperini data analysis. The code is commented to increase ease of adoption and use. It is assumed that this bash file is being executed from within the utilities directory.
# Note: stringr and monocle must be installed for this analysis to work.

# Set the machine.
machine=local

echo Build and install the sceptre package.
# bash build_and_install_package.bash sceptre $machine

# Obtain the filepaths to the code and "offsite" directories
code_dir=$(bash get_file_paths.bash $machine code)
offsite_dir=$(bash get_file_paths.bash $machine data_results)

echo Initialize the offsite directory structure.
# Rscript $code_dir"/analysis_drivers_gasp/"check_directory_structure_1.R $code_dir $offsite_dir

echo Download the data.
# Rscript $code_dir"/analysis_drivers_gasp/"download_raw_data_2.R $code_dir $offsite_dir

echo Pre-process the data.
# Rscript $code_dir"/analysis_drivers_gasp/"pre_process_data_3.R $code_dir $offsite_dir

echo Perform quality control.
# Rscript $code_dir"/analysis_drivers_gasp/"perform_quality_control_4.R $code_dir $offsite_dir

echo Create the precomputation and results dictionaries.
pod_sizes=$(Rscript $code_dir"/analysis_drivers_gasp/"create_file_dictionaries_5.R $code_dir $offsite_dir)
n_gene_pods="$(echo $pod_sizes | cut -d' ' -f1)"
n_gRNA_pods="$(echo $pod_sizes | cut -d' ' -f2)"
n_pair_pods="$(echo $pod_sizes | cut -d' ' -f3)"
echo number of gene pods: $n_gene_pods
echo number of gRNA pods: $n_gRNA_pods
echo number of gRNA-gene pair pods: $n_pair_pods

echo Run gene precomputation across all gene pods.
seq 1 $n_gene_pods | xargs -I{} -n 1 -P 2 Rscript $code_dir"/analysis_drivers_gasp/"run_gene_precomputation_6.R $code_dir $offsite_dir {} &

echo Run gRNA precomputation across all gRNA pods.
seq 1 $n_gRNA_pods | xargs -I{} -n 1 -P 2 Rscript $code_dir"/analysis_drivers_gasp/"run_gRNA_precomputation_7.R $code_dir $offsite_dir {} &
wait

echo Run gRNA-gene pair analysis across all pair pods.
seq 1 $n_pair_pods | xargs -I{} -n 1 -P 3 Rscript $code_dir"/analysis_drivers_gasp/"run_pair_analysis_at_scale_8.R $code_dir $offsite_dir {} &
wait

echo Collect and save results.
Rscript $code_dir"/analysis_drivers/"aggregate_results_9.R $code_dir $offsite_dir
