#!/bin/sh

# this is a variation of Pej_Chrushit.sh, for cases that the file is far from sorted by the desired column, in which case Pej_Chrushit.sh gets very slow, this version sorts the file first and then crushes it. 
# This crushes a file by a given column:
# Chrushit Myfile 1 # crush myfile based on 1st column
# the closer to sorted the column is the faster this becomes.
# NOTE: This file assumes the first line is a header line!
#=================================
# Pejman, March 2015
#=================================

DIR=$(dirname $1)
Fbname=$(basename $1)
OutFldr=$DIR'/'$Fbname'_Crushed_by_Column'$2'/'


if [ -d  $OutFldr ]; then
    echo Nothing done!
    echo The directory $OutFldr already exists, delete it and rerun the script.
    exit -1
fi

mkdir $OutFldr
sort-header -k$2,$2 $1 | tawk -v Col=$2 -v OF=$OutFldr '{
if (NR==1)
{
    Header=$0
    next
}
CurrF=(OF $Col ".cr");
if(CurrF!=OldF)
{
    close(OldF);
    OldF=CurrF;
    if ( system( "[ -f "CurrF" ] " ) !=0)
    {
       #This is a new file
      print Header >CurrF   
    }  
}
print >>CurrF;

}' 
echo ALL Done!


