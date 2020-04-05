#!/bin/bash

if [ ! -d illumina_sequences ]
then
	mkdir illumina_sequences
fi

echo "identifying the illumina sequences"

for f in metadata/*.csv
do
	if grep -Fqi "Illumina" $f
	then
		#echo "found in $f"
		seqName=$( basename $f .csv )
		if [ ! -f illumina_sequences/$seqName.fastq ]
		then
			ln -s ../fastq/${seqName}.fastq illumina_sequences/.
		fi
	fi
done

# align

echo "alignging sam files"

if [ ! -d sam ]
then
    mkdir sam
fi

if [ -f aligning.txt ]
then
	rm -rf aligning.txt
fi

counter=0
for f in illumina_sequences/*.fastq
do
    seqName=$( basename $f .fastq)
    #echo "looking at $seqName"
    if [ ! -f sorted_bam/$seqName.sorted.bam ]
    then
        echo "aligning $seqName"
        echo "bwa mem -t 4 reference.fasta $f > sam/$(basename $f .fastq).sam" >> aligning.txt
		counter=$(( $counter + 1 ))
    fi
done

if [ $counter -gt 0 ]
then
	parallel -j12 --bar :::: aligning.txt
fi

# sort

echo "sorting sam files"

if [ -f sorting.txt ]
then
    rm -rf sorting.txt
fi

if [ ! -d sorted_sam ]
then
    mkdir sorted_sam
fi

counter=0
for f in illumina_sequences/*.fastq
do
    seqName=$( basename $f .fastq )
    if [ ! -f sorted_bam/${seqName}.sorted.bam ]
    then
        #echo "sorting  $( basename $f .sam)"
        echo "samtools sort sam/${seqName}.sam -o sorted_sam/${seqName}.sorted.sam" >>  sorting.txt
		counter=$(( $counter + 1 ))
    fi
done

if [ $counter -gt 0 ]
then
	parallel -j12 --bar :::: sorting.txt
fi

# convert sam2bam

echo "converting sam files to bam files"

if [ -f converting.txt ]
then
    rm -rf converting.txt
fi

if [ ! -d sorted_bam ]
then
	mkdir sorted_bam
fi

counter=0
for f in illumina_sequences/*.fastq
do
    seqName=$( basename $f .fastq )
    if [ ! -f sorted_bam/${seqName}.sorted.bam ]
    then
        #echo "converting $( basename $f .sorted.sam)"
        echo "samtools view -S -b sorted_sam/${seqName}.sorted.sam > sorted_bam/${seqName}.sorted.bam" >> converting.txt
		counter=$(( $counter + 1 ))
    fi
done

if [ $counter -gt 0 ]
then
	parallel -j12 --bar :::: converting.txt
fi

# de-duplication

echo "marking duplicate reads"

if [ -f marking.txt ]
then
	rm -rf marking.txt
fi

if [ ! -d sorted_dedup_bam ]
then
	mkdir sorted_dedup_bam
fi

counter=0
for f in illumina_sequences/*.fastq
do
	seqName=$( basename $f .fastq )
	if [ ! -f sorted_dedup_bam/${seqName}.sorted.dedup.bam ]
	then
		echo "gatk --java-options \"-Djava.io.tmpdir=/scratch/tmp -Xmx12G\" MarkDuplicates -I sorted_bam/${seqName}.sorted.bam -O sorted_dedup_bam/${seqName}.sorted.dedup.bam -M sorted_dedup_bam/${seqName}.dedup.metrics.txt" >> marking.txt
		counter=$(( $counter + 1 ))
	fi
done

if [ $counter -gt 0 ]
then
	parallel -j10 --bar :::: marking.txt
fi

# index bam files

echo "indexing bam files"

if [ -f indexing.txt ]
then
	rm -rf indexing.txt
fi

counter=0
for f in illumina_sequences/*.fastq
do
	seqName=$( basename $f .fastq )	
	if [ ! -f sorted_dedup_bam/${seqName}.sorted.dedup.bam.bai ]
	then
		echo "samtools index sorted_dedup_bam/$seqName.sorted.dedup.bam" >> indexing.txt
		counter=$(( $counter + 1 ))
	fi
done

if [ $counter -gt 0 ]
then
	parallel -j10 --bar :::: indexing.txt
fi
