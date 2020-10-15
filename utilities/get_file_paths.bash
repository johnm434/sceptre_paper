#!/bin/bash

# This script takes as an argument (i) the name of a machine ("local", "hydra", "ubergenno", "uberduo", "bridges") and (ii) the desired directory ("code" or "data_results"). It outputs the full file path of that directory on that machine.

# Define machine_name, code_dir, and data_results_dir arrays.
machine_name=(local hydra ubergenno uberduo bridges)
code_dir=(/Users/timbarry/Box/SCEPTRE/sceptre_paper /home/tbarry2/SCEPTRE/sceptre_paper ubergenno_code_dir uberduo_code_dir bridges_code_dir)
data_results_dir=(/Volumes/tims_new_drive/research/sceptre_files /raid6/Tim/sceptre_offsite_dir ubergenno_data_dir uberduo_data_dir bridges_data_dir)

# find the index of the selected machine
for ((index=0; index<${#machine_name[@]}; index++))
  do
		if [ ${machine_name[$index]} == $1 ]
		then
		  machine_idx=$index
		fi
  done

# echo the correct directory
if [ $2 == "code" ]
then
  echo ${code_dir[$machine_idx]}
fi

if [ $2 == "data_results" ]
then
  echo ${data_results_dir[$machine_idx]}
fi
