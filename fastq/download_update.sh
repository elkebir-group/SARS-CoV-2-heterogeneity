#!/bin/bash

# download from SRA
IDS=$(curl 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=txid2697049\[Organism%3Anoexp\]&retmax=1000' | grep "<Id>" | grep -o "[[:digit:]]\+" | tr "\n" ",")
curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&amp;id=$IDS" | grep -o "SRR[[:digit:]]\+" | sort | uniq > sra.txt

acc_file="sra.txt"

num_entries=$( wc -l $acc_file | rev | cut -d " " -f 2 | rev )

echo "number of enteries in the accession file is " $num_entries

counter=0

for i in $(seq 1 $num_entries)
do 
    line=$( sed -n -e "$i"p $acc_file )

    echo downloading $line.fastq

#    if [ ! -f $line.fastq ]
#    then
#        echo "$line not found"
#    fi

    if [ -f $line.fastq ]
    then
        echo "$line.fastq file already exists"
    else
        echo "$line not found"
        fastq-dump $line
        echo "$line.fastq downloaded"
        counter=$(( $counter + 1 ))
    fi

done

echo "downloaded $counter new files"
