#!/bin/bash

# this script assumes execution is occurring in the sceptre_paper directory (i.e., the command is bash bash_scripts/build_and_install_package.bash).

# Obtain the version of the package through an awk call on DESCRIPTION
version=$(awk '/Version/ {print $2}' "sceptre/DESCRIPTION")

# build the package
R CMD build sceptre

# install the package
R CMD INSTALL "sceptre_$version.tar.gz"

# remove the tar.gz file
rm "sceptre_$version.tar.gz"
