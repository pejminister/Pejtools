#!/bin/sh

# This is a simple version of linux commmand "split" except that it keeps the header line
# The first parameter is the inputfile name, and the second one is the number of lines in each slice
# Example: Pej_Split_Header Myfile.txt 500
# NOTE: This file assumes the first line is a header line!
#=================================
# Pejman, March 2018
#=================================

tail -n +2 $1 | split -l $2 -a 6 - $1"_slice_"
for file in $1"_slice_"*
do
    head -n 1 $1 > tmp_file
    cat $file >> tmp_file
    mv tmp_file $file
done

mkdir $1"_sliced"
mv $1"_slice_"* $1"_sliced/"
