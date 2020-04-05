# SARS-CoV-2-heterogeneity

This repository has python and bash scripts to perform SNV calling on bulk-sequencing SARS-CoV-2 data on the SRA database.
The steps in the pipeline are as follows

	0. Download the fastq from the SRA database and identify the ILLUMINA sequences
	1. Map reads to the reference genome using bwa
	2. Sort the reads using samtools and convert to bam files
	3. Mark duplicates reads using gatk
	4. Filter sequences based on some mean depth and breath coverage threshold 
	5. Call SNVs using bcftools
	6. Download metadata using pysradb (https://github.com/saketkc/pysradb)		
	7. Filter and analyze the called SNVs

# Downloading Script

This script must be run in the fastq folder

Usage: ./download_from_list.sh \<SRA ACC list file\>

	$ cd fastq
	$ ./download_from_list.sh SRR_Acc_List_April2.txt

The metadata will be stored in ../metadata folder.

# Preprocessing Script

Thie script should be run in the SARS-CoV-2-heterogeneity folder.

Usage: ./preprocess_illumina.sh

	$ cd ..
	$ ./preprocess_illumina.sh

It will generate the following folders

	1. illumina_seqeunces : It will contain soft link to the illumina sequences
	2. sam : It will contain the sam files
	3. sorted_sam : The sorted sam files
	4. sorted_bam : The sorted bam files
	5. sorted_dedup_bam : The sorted bam files with marked duplicates and its indexing

# Filtering Script

This script must also be from SARS-CoV-2-heterogeneity folder.

Usage: ./filter_illumina \<output folder\> (\<mean depth threshold\>, default=50) (\<coverage breadth threshold\>, default=20) 

	$ ./filter_illumina filtered_illumina_sequences

It will generate the folder filtered_illumina_sequences which will contain soft link to the sequences that pass the threshold

# SNV calling

Usage: ./callbcf.sh \<output file name\>

	$ cd filtered_illumina_sequences
	$ ln -s ../scripts/callbcf.sh .
	$ ./callbcf.sh snv_file.vcf

# Downloading Metadata

	$ cd ..
	$ for f in filtered_illumina_sequences/*.bam; do echo $( basename $f .sorted.dedup.bam ); done > filtered_list_April2.txt
	$ pysradb metadata $(<filtered_list_April2.txt) --saveto filtered_list_metadata_April2.txt --detailed

# Postprocessing

This requires a fasta file with the consensus sequences (consensus.fasta).

	usage: process_vcf.py [-h] [--vcfFile VCFFILE] [--outputDir OUTPUTDIR]
			      [--metadata METADATA] [--consensus CONSENSUS]
			      [--qualitythreshold QUALITYTHRESHOLD]

	optional arguments:
	  -h, --help            show this help message and exit
	  --vcfFile VCFFILE, -i VCFFILE
				input alternate allele count file in space separated
				format
	  --outputDir OUTPUTDIR, -o OUTPUTDIR
				output directory for the sample summary files
	  --metadata METADATA, -m METADATA
				metadata file in csv format
	  --consensus CONSENSUS, -c CONSENSUS
				consensus seqeuncesin fasta format
	  --qualitythreshold QUALITYTHRESHOLD, -q QUALITYTHRESHOLD
				threshold for the phred quality score


Example:

	$ cd ../script
	$ python process.py -i ../filtered_illumina_sequences/snv_file.vcf -o ../results -m ../filtered_list_metadata_April2.txt -c ../consensus.fasta -q 20 
