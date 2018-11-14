# This Extracts the first 5 columns (chr and position, ID, and alleles) and a given field from the info column of a VCF file.
# Example usage:
# ./Pej_VCF_getinfo.sh Myfile.vcf SEVERE_GENE
# this would get the filed "SEVERE_GENE" out of the info column
#=================================
# Pejman, April 2015
#=================================

#!/bin/sh
awk -F'\t' -v Key="$2" -v OFS="\t" '{
                           n = match($8,Key"=[^;]+") ;
                           if (n != 0) {
                              info=substr($8,RSTART+length(Key)+1,RLENGTH-length(Key)-1);
                              print  $1,$2,$3,$4,$5,info
                           }
                       }' $1
