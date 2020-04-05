#!/bin/bash

if [ $# == 0 ] 
then
  echo "Usage: $0 <output folder> (<mean depth threshold>, default=50) (<coverage breadth threshold>, default=20)"
  exit 1
fi

seqFolder="sorted_dedup_bam"
outFolder=$1

depth_threshold=50
cov_threshold=20

if [ $# -ge 2 ]
then
	depth_threshold=$2
fi

if [ $# -ge 3 ]
then
	cov_threshold=$3
fi

if [ -f depth.txt ]
then
	rm -rf depth.txt
fi

if [ ! -d $outFolder ]
then
	mkdir $outFolder
fi

if [ ! -d $seqFolder/depth_output ]
then
	mkdir $seqFolder/depth_output
fi

counter=0
for f in $seqFolder/*.bam
do
	seqName=$( basename ${f} .sorted.dedup.bam )
	if [ ! -f $seqFolder/depth_output/$seqName.depth ]
	then
		echo "samtools depth -a $f > $seqFolder/depth_output/$seqName.depth" >> depth.txt
		counter=$(( $counter  + 1 ))
	fi
done

if [ $counter -gt 0 ]
then
	parallel -j12 --bar :::: depth.txt
fi

statFile="$seqFolder/all_stats_file.csv"

if [ ! -f $statFile ]
then
	echo "accNum,mean_depth,coverage" > $statFile
fi

for f in $seqFolder/*.bam
do

	seqName=$( basename ${f} .sorted.dedup.bam )
	
	depthFile="${seqFolder}/depth_output/${seqName}.depth"

	read mean_depth <<< $( cat $depthFile | awk '{c++; s+=$3}END{print s/c}' )
	read coverage <<< $( cat $depthFile | awk '{c++; if($3>0) total+=1}END{print (total/c)*100}' )
				
	if ! grep -Fq $seqName $statFile
	then

		echo $seqName,$mean_depth,$coverage >> $statFile

	fi

	if (( $( echo "$mean_depth > $depth_threshold" | bc -l ) )) && (( $( echo "$coverage > $cov_threshold" | bc -l ) )) && [[ ! -f $outFolder/$seqName.sorted.dedup.bam ]]
	then
		ln -s ../$f $outFolder/.
		ln -s ../$f.bai $outFolder/.
	fi

done

for f in $outFolder/*.bam
do
	echo $( basename $f .sorted.dedup.bam )
done > filtered_list.txt	
