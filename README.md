# SARS-CoV-2-heterogeneity

This repository has python and bash scripts to perform SNV calling on bulk-sequencing SARS-CoV-2 data on the SRA database.
The steps in the pipeline are as follows
	0. Download the fastq from the SRA database
	1. Map reads to the reference genome using bwa
	2. Sort the reads using samtools and convert to bam files
	3. Mark duplicates reads using gatk
	4. Call SNVs using bcftools
