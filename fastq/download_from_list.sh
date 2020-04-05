#!/bin/bash

if [ $# != 1 ] 
then
  echo "Usage: $0 <SRA ACC list file>"
  exit 1
fi

acc_file=$1

num_entries=$( wc -l $acc_file | rev | cut -d " " -f 2 | rev )

echo "number of enteries in the accession file is " $num_entries

counter=0

if [ ! -d ../metadata ]
then
	mkdir ../metadata
fi 

for i in $(seq 1 $num_entries)
do 
    line=$( sed -n -e "$i"p $acc_file )


    if [ ! -f $line.fastq ]
    then
    	echo "downloading $line.fastq"
        fastq-dump $line
        echo "$line.fastq downloaded"
        counter=$(( $counter + 1 ))
    fi

	if [ ! -f ../metadata/$line.csv ]
	then
		echo "donwloading metadata for $line" 
		wget -O ../metadata/$line.csv "http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=${line}"
	fi

done

echo "downloaded $counter new files"

