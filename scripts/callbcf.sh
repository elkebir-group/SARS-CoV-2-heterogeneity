#!/bin/bash

if [ $# == 0 ]
then
  echo "Usage: $0 <output file name>"
  exit 1
fi

filelist=" "

for file in *.bam
do
    filelist="${filelist} ${file}"
done

echo "bcftools mpileup -f ../reference.fasta -q 20 -Q 20 -a DP,AD,ADF,ADR,SP --skip-indels $filelist > $1"  
