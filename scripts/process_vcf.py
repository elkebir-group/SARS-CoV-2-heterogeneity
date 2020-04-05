#!/usr/bin/python
from __future__ import print_function

import os, sys, glob
import math
import numpy as np
import pandas as pd
#import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import allel

parser = argparse.ArgumentParser()
#parser.add_argument('--refFile','-r', help='input reference allele count file in space separated format',default='../results/ref_file.out')
#parser.add_argument('--altFile','-a', help='input alternate allele count file in space separated format',default='../results/alt_file.out')
parser.add_argument('--vcfFile','-i', help='input alternate allele count file in space separated format',default='../results/covid.vcf')
parser.add_argument('--outputDir','-o', help='output directory for the sample summary files',default='../test/results')
parser.add_argument('--metadata','-m', help='metadata file in csv format',default='../results/metadata.csv')
parser.add_argument('--consensus','-c',help='consensus seqeuncesin fasta format',default='../restuls/consensus.fasta')

args = parser.parse_args()

print('reading the vcf file')

callset = allel.read_vcf('../test/covid_annotated_dedup_April2.vcf', fields='*')

sampleNames=[element.split('.')[0] for element in callset['samples']]

data_ref = {}

data_ref['pos'] = list(callset['variants/POS'])
data_ref['ref'] = list(callset['variants/REF'])
data_ref['alt'] = list([alt[0] for alt in callset['variants/ALT']])
for i in range(len(sampleNames)):
    sample = sampleNames[i]
    data_ref[sample] = list([ad[i][0] for ad in callset['calldata/AD']])
        
df_ref = pd.DataFrame(data_ref).set_index(['pos','ref','alt'])

data_alt = {}

data_alt['pos'] = list(callset['variants/POS'])
data_alt['ref'] = list(callset['variants/REF'])
data_alt['alt'] = list([alt[0] for alt in callset['variants/ALT']])
for i in range(len(sampleNames)):
    sample = sampleNames[i]
    data_alt[sample] = list([ad[i][1] for ad in callset['calldata/AD']])
        
df_alt = pd.DataFrame(data_alt).set_index(['pos','ref','alt'])


print('removing sequences with the same sample accession numbers')

mData = args.metadata
df_meta = pd.read_csv(mData, sep='\t', index_col=['run_accession'])

sampleAccessions = set()
idx = 0
for column in df_ref:
    name = df_meta.loc[column]['sample_accession']
    if name in sampleAccessions:
        print('removing ',str(column),' with sample name ',str(name))
        idx += 1
        del df_ref[column]
        del df_alt[column]
    else:
        sampleAccessions.add(name)

print('removed ',str(idx),' samples')


df_vaf=(df_alt/(df_alt+df_ref))
df_dep = df_alt + df_ref

samples = list(df_vaf.columns)
positions = list(df_vaf.index)

print('number of samples is ', str(len(samples)))
print(samples)

outDir = args.outputDir + '/'

if not os.path.exists(outDir):
    os.mkdir(outDir)


subclonal_pos = set()
shared_count = {}

if not os.path.exists(outDir + "/sample_summary.tsv"):
    print('writing summary files')
    with open(outDir + "/sample_summary.tsv", "w") as fsum:
        fsum.write("\t".join(["sample", "mutations", "type", "depth"]) + "\n")
        for sample in samples:
            clonal = 0
            subclonal = 0
            depth = []
            with open(outDir + "/%s.tsv" % sample, "w") as f:
                f.write("\t".join(["pos", "ref_allele", "alt_allele", "ref", "alt", "vaf", "tot", "type"]) + "\n")
                for pos in positions:
                    vaf = df_vaf.loc[pos][sample]
                    ref = df_ref.loc[pos][sample]
                    alt = df_alt.loc[pos][sample]
                    dep = df_dep.loc[pos][sample]
                    depth.append(dep)
                    t = "none"
                    if 0.05 <= vaf <= 0.9 and alt >= 5:
                        t = "subclonal"
                        subclonal += 1
                        subclonal_pos.add(pos)
                        shared_count
                    elif vaf > 0.9:
                        t = "clonal"
                        clonal += 1
                    if t != "none":
                        if pos not in shared_count:
                            shared_count[pos] = 0
                        shared_count[pos] += 1
                        
                    f.write("\t".join(map(str, [pos[0], pos[1], pos[2], ref, alt, vaf, dep, t])) + "\n")
                    
            fsum.write("\t".join(map(str, [sample, clonal, "clonal", np.mean(depth)])) + "\n")
            fsum.write("\t".join(map(str, [sample, subclonal, "subclonal", np.mean(depth)])) + "\n")
else:
    print('finding subclonal mutations in samples')
    for sample in samples:
        for pos in positions:
            vaf = df_vaf.loc[pos][sample]
            alt = df_alt.loc[pos][sample]
            t = "none"
            if 0.05 <= vaf <= 0.9 and alt >= 5:
                subclonal_pos.add(pos)
                t = "subclonal"
            elif vaf > 0.9:
                t = "clonal"
                
            if t != "none":
                if pos not in shared_count:
                    shared_count[pos] = 0
                shared_count[pos] += 1

df_sum = pd.read_csv(outDir + "/sample_summary.tsv", sep="\t")

if not os.path.exists(outDir + "/subclonal.pdf"):
    print('plotting number of subclonal mutations per sample')    
    sns.barplot(palette=[sns.color_palette()[1]]+[sns.color_palette()[0]]*15, x="mutations", y="sample", data=df_sum[df_sum["type"] == "subclonal"].groupby("mutations").count().reset_index())
    plt.gca().set_xlabel("number of subclonal mutations")
    plt.gca().set_ylabel("number of samples")
    plt.savefig(outDir + "/subclonal.pdf")

if not os.path.exists(outDir + "/subclonal_mutations.pdf"):
    print('making subclonal mutation distribution plot')
    df_plot = pd.DataFrame([(pos, shared_count[pos]) for pos in subclonal_pos], columns=['pos', 'samples']).groupby('samples').count().reset_index()
    sns.barplot(palette=[sns.color_palette()[2]] + [sns.color_palette()[3]]*11,
                data= df_plot, 
                x="samples", y="pos")
    plt.gca().set_xlabel("number of samples")
    plt.gca().set_ylabel("number of subclonal mutations")
    plt.savefig(outDir + "/subclonal_mutations.pdf")

print('picking the samples that have high depth and subclonal mutations')
variant_positions_with_clonal_reduced = set([])
heterogeneous_samples_with_clonal_reduced = []
idx = 1
for sample in samples:
    #print  df_meta.loc[sample]['geo_loc_name'] 
    df = pd.read_csv(outDir + "/%s.tsv" % sample, sep="\t")
    df_res = df[(df["vaf"] >= 0.1) & (df["tot"] >= 50)]
    if len(df_res) > 0:
        heterogeneous_samples_with_clonal_reduced.append(sample)
        #print(idx, sample, list(df_res["pos"]))
        variant_positions_with_clonal_reduced |= set(df_res["pos"])
        idx += 1

df_het_vaf_with_clonal_reduced = df_vaf[heterogeneous_samples_with_clonal_reduced].reset_index()
df_het_vaf_with_clonal_reduced = df_het_vaf_with_clonal_reduced[df["pos"].isin(variant_positions_with_clonal_reduced)].set_index(["pos", "ref", "alt"])

print('making VAF=0 if the number of reads supporting alternate allele is less than 5')
for index in df_het_vaf_with_clonal_reduced.index:
    for column in df_het_vaf_with_clonal_reduced.columns:
        if df_alt[column][index] < 5 & df_dep[column][index] != 0:
            #print(column, index)
            df_het_vaf_with_clonal_reduced.loc[index,column] = 0

df_het_vaf_with_clonal_reduced.columns = [str(col_name) + '|' + str(df_meta.loc[col_name]['geo_loc_name']) for col_name in df_het_vaf_with_clonal_reduced.columns]
            
if not os.path.exists(outDir + "/VAF_heatmap_with_clonal_reduced.pdf"):
    print('making VAF heat map')
    cg = sns.clustermap(df_het_vaf_with_clonal_reduced, yticklabels=1, xticklabels=1, annot=True, fmt=".1f", row_cluster=False, col_cluster=False)
    cg.ax_row_dendrogram.set_visible(False)
    cg.ax_col_dendrogram.set_visible(False)
    cg.cax.set_visible(False)
    # cg.ax_row_dendrogram.set_xlim([0,0])
    plt.gcf().set_size_inches(50,50)
    plt.tight_layout()
    plt.savefig(outDir + "/VAF_heatmap_with_clonal_reduced.pdf")

print('reading consensus sequences')
sequences = {}

with open(args.consensus,'r') as input:
    for line in input:
        if line.startswith('>'):
            data = line.strip().split('\n')
            fastaName = line[1:-1]
            flag = 1
        if not line.startswith('>'):
            data = line.strip().split('\n')
            if flag == 1:
                sequences[fastaName] = "".join(data)              
                flag = 0
            else:
                sequences[fastaName] = sequences[fastaName] + "".join(data)                

reference = sequences["Wuhan-Hu-1/2019"]

# generate the codon table
print('generating codon table')
bases = "TCAG"
codons = [a + b + c for a in bases for b in bases for c in bases]
amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
codon_table = dict(zip(codons, amino_acids))
    

print('final results are being generated')
df = df_het_vaf_with_clonal_reduced

NS = []
AA = []
geneList = []
nSRAsubclonal = []
nSRAclonal = []
nSRA = []
nConsensus = []

for ind in df.index:
    pos = ind[0]
    ref = ind[1]
    alt = ind[2]
    
    if 266 < pos < 13468:
        gene = "ORF1a"
        start = 266
    elif 13468 < pos < 21555:
        gene="ORF1b"             
        start = 13468
    elif 21563 < pos < 25384:
        gene="S"
        start = 21563
    elif 25393 < pos < 26220:
        gene="ORF3a"
        start = 25393
    elif 26245 < pos < 26472:
        gene="E"
        start = 26245
    elif 26523 < pos < 27191:
        gene="M"     
        start = 26523
    elif 27202 < pos < 27387:
        gene="ORF6"
        start = 27202
    elif 27394 < pos < 27759:
        gene="ORF7a"
        start = 27394
    elif 27756 < pos < 27887:
        gene="ORF7b"
        start = 27756
    elif 27894 < pos < 28259:
        gene="ORF8"
        start = 27894
    elif 28274 < pos < 29533:
        gene="N"
        start = 28274
    elif 28284 < pos < 28577:
        gene="ORF9b"
        start = 28284
    elif 29558 < pos < 29674:
        gene="ORF10"
        start = 29558
    elif 28734 < pos < 28955:
        gene="ORF14"
        start = 28734
    else:
        gene="noncoding"

    geneList.append(gene)
        
    if gene != "noncoding":
        offset = (pos - start)%3
        aa_pos = math.floor((pos - start)/3) + 1
        a = pos - offset
        b = pos + (1 - offset)
        c = pos + (2 - offset)

        refCodon = reference[a-1] + reference[b-1] + reference[c-1]

        if offset == 0:
            modCodon = alt + reference[b-1] + reference[c-1]
        elif offset == 1:
            modCodon = reference[a-1] + alt + reference[c-1]
        else:
            modCodon = reference[a-1] + reference[b-1] + alt

        if codon_table[refCodon] != codon_table[modCodon]:
            synResult = 'N'
        else:
            synResult = 'S'

        AA.append(codon_table[refCodon] + str(aa_pos) + codon_table[modCodon])

        NS.append(synResult)
    else:
        AA.append('N/A')
        NS.append('N/A')
        
    numConsensus = 0
    for sample in sequences:
        allele = sequences[sample][pos - 1]
        if allele == alt:
            numConsensus += 1
    
    nConsensus.append(numConsensus)

    numSubclonal = 0
    numClonal = 0
    for column in df.columns:
        vaf = df[column][ind]
        if vaf >= 0.1 and vaf <= 0.9:
            numSubclonal += 1
        elif vaf >= 0.9:
            numClonal += 1
    
    nSRAsubclonal.append(numSubclonal)
    nSRAclonal.append(numClonal)
    nSRA.append(numSubclonal + numClonal)

nCols = len(df.columns)
    
# adding columns to the vaf dataframe
#N/S, AA, nSRAsubclonal, nSRAclonal, nSRA, nConsensus    
df["gene"] = geneList
df["N/S"] = NS
df["AA"] = AA
df["nSRAsubclonal"] = nSRAsubclonal
df["nSRAclonal"] = nSRAclonal
df["nSRA"] = nSRA
df["nConsensus"] = nConsensus

cols = df.columns.to_list()
newcols = cols[nCols:] + cols[:nCols]
df = df[newcols]

print('writing final vaf results')
df.to_csv(outDir + "/final_vaf_results.csv", na_rep='NULL')

# changing dataframe to get ref values
for column in list(df.columns)[7:]:
    seqName=column.split('|')[0]
    for index in df.index:
        #df[column][index] = df_ref[seqName][index]
        df.loc[index,column] = int(df_ref[seqName][index])

print('writing final ref results')
df.to_csv(outDir + "/final_ref_results.csv", na_rep='NULL')

# changing dataframe to get alt values
for column in list(df.columns)[7:]:
    seqName=column.split('|')[0]
    for index in df.index:
        #df[column][index] = df_alt[seqName][index]
        df.loc[index,column] = int(df_alt[seqName][index])

print('writing final alt results')
df.to_csv(outDir + "/final_alt_results.csv", na_rep='NULL')
