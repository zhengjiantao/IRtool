#!/usr/bin/env python
# encoding: utf-8

import os
import sys
import csv
import pandas as pd
import numpy as np
#from tqdm import *
from functools import reduce
from collections import defaultdict,Counter

def ProduceIntronGTF(genome_gtf, gene_isoform_intron_map):
    """
    append intron transcript to gtf (filter by gene map), current finish single intron transcript
    todo: intron combination transcript
    output: intron.gtf, genome_extend.gtf, intron_to_transcript.txt
    """
    gtf = pd.read_csv(genome_gtf, sep="\t", header=None)
    gtf[[0,1,2,8]] = gtf[[0,1,2,8]].astype(str)
    valid_gene = gtf[8].apply(lambda x: True if x.split("gene_id \"")[1].split("\";")[0] in gene_isoform_intron_map else False)
    valid_gtf = gtf[valid_gene]
    valid_gtf = valid_gtf[valid_gtf[2].apply(lambda x: True if x in ["gene","transcript","exon"] else False)]

    gene_idx = list(valid_gtf[valid_gtf[2] == "gene"].index)
    gene_idx += [valid_gtf.index[-1]+1]
    all_gtf = pd.DataFrame() # valid_gtf + intron_gtf
    intron_gtf = pd.DataFrame() # only intron transcript gtf
    intron_to_transcript = pd.DataFrame() # intron-transcript-gene

    #for gpre_idx, gnext_idx in tqdm(zip(gene_idx[:-1], gene_idx[1:])):
    for gpre_idx, gnext_idx in zip(gene_idx[:-1], gene_idx[1:]):
        gene_fragment = valid_gtf.loc[gpre_idx:gnext_idx-1]
        ensg = gene_fragment.loc[gpre_idx][8].split("gene_id \"")[1].split("\";")[0]
        retained_introns = gene_isoform_intron_map[ensg]
        all_gtf = all_gtf.append(gene_fragment.loc[gpre_idx]) # add gene row
        intron_gtf = intron_gtf.append(gene_fragment.loc[gpre_idx])

        transcript_idx = list(gene_fragment[gene_fragment[2] == "transcript"].index)
        transcript_idx += [gene_fragment.index[-1]+1]
        for tpre_idx, tnext_idx in zip(transcript_idx[:-1], transcript_idx[1:]):
            transcript_fragment = gene_fragment.loc[tpre_idx:tnext_idx-1].sort_values(by=3) # sort by start coordinate
            enst = transcript_fragment.loc[tpre_idx][8].split("transcript_id \"")[1].split("\";")[0]
            all_gtf = all_gtf.append(transcript_fragment) # add old transcript
            count = 0 # new transcript id

            #### find exon which overlap coordinate with intron, create new transcript and append
            coordinates = transcript_fragment.loc[:,3:4]
            transcript_coordinate = coordinates.iloc[0].values.tolist()
            for intron in retained_introns:
                new_transcript_fragment = pd.DataFrame()
                new_transcript_fragment = new_transcript_fragment.append(transcript_fragment.loc[tpre_idx]) # add transcript row
                create_flag = False
                # check intron whether in transcript
                if transcript_coordinate[0]<=intron[1] and intron[2]<=transcript_coordinate[1]-1: # 1 base
                    # find overlap exon: 5' overlap or 3' overlap (3' exon 5')
                    for exon in coordinates.iloc[1:].itertuples():
                        df_index = exon[0]
                        if create_flag:
                            new_transcript_fragment = new_transcript_fragment.append(transcript_fragment.loc[df_index])
                            continue
                        if intron[1]==exon[2] and exon[1]!=exon[2]: # intron start==exon end, 5' overlap
                            create_flag = True
                            extend_exon = transcript_fragment.loc[df_index]
                            extend_exon.iloc[4] = intron[2]
                            new_transcript_fragment = new_transcript_fragment.append(extend_exon)
                        elif intron[2]==exon[1]-1 and exon[1]!=exon[2]: # 3' overlap
                            create_flag = True
                            extend_exon = transcript_fragment.loc[df_index]
                            extend_exon.iloc[3] = intron[1]+1
                            new_transcript_fragment = new_transcript_fragment.append(extend_exon)
                        else:
                            new_transcript_fragment = new_transcript_fragment.append(transcript_fragment.loc[df_index])

                if create_flag:
                    count += 1
                    # need change trancript_id: ENST -> ENST_count, todo: change trancript_name
                    if "ENSMUST" in new_transcript_fragment.iloc[0,8]:
                        t = new_transcript_fragment.iloc[0,8].index('ENSMUST')+11+7
                    else:
                        t = new_transcript_fragment.iloc[0,8].index('ENST')+11+4
                    new_transcript_fragment[8] = new_transcript_fragment[8].apply(lambda x: "{}_ir{}{}".format(x[:t],count,x[t:]))
                    new_enst = "{}_ir{}".format(enst, count)

                    all_gtf = all_gtf.append(new_transcript_fragment) # add new intron transcript
                    intron_gtf = intron_gtf.append(new_transcript_fragment)
                    intron_to_transcript = intron_to_transcript.append(pd.Series([reduce(lambda a,b: str(a)+"-"+str(b), intron), new_enst, ensg]), ignore_index=True)

    all_gtf[[3,4]] = all_gtf[[3,4]].astype(int)
    intron_gtf[[3,4]] = intron_gtf[[3,4]].astype(int)
    return all_gtf, intron_gtf, intron_to_transcript

def ReadIreadResult(iread_res):
    '''read iread result, and return a map: gene-intron list(chr,start,end)'''
    print(sample_id)
    iread_data = pd.read_csv(iread_res,sep="\t")
    retained_intron = iread_data[iread_data["retention_at_given_cutoff"]=="yes"]
    print("{}: the number of retained intron is {}".format(sample_id, retained_intron.shape[0]))

    gene_isoform_intron_map = defaultdict(list)
    for items in retained_intron.itertuples():
        intron_info = items[1]
        chrom,start,end,ensg,_ = intron_info.split("-")
        gene_isoform_intron_map[ensg].append([str(chrom),int(start),int(end)])

    max_introns = 0
    max_gene = ""
    for gene, introns in gene_isoform_intron_map.items():
        if len(introns) > max_introns:
            max_introns = len(introns)
            max_gene = gene
    print("{}: the number of genes is {}, the gene with the mose introns is {} with {}".format(sample_id, len(gene_isoform_intron_map), max_gene, max_introns))
    return gene_isoform_intron_map

print(__name__)
if __name__ == "__main__":
    genome_gtf = sys.argv[1]
    sample_id = sys.argv[2]
    method_res = sample_id+"Aligned.sortedByCoord.out.ir.txt"
    out_dir = sys.argv[3]
    print("processing {}".format(sample_id))
    gene_isoform_intron_map = ReadIreadResult(method_res)
    all_gtf, intron_gtf, intron_to_transcript = ProduceIntronGTF(genome_gtf, gene_isoform_intron_map)

    all_gtf_file = os.path.join(out_dir,sys.argv[2]+"genome_extend.gtf")
    intron_gtf_file = os.path.join(out_dir,sys.argv[2]+"intron.gtf")
    intron_to_transcript_file = os.path.join(out_dir,sys.argv[2]+"intron_to_transcript.txt")

    all_gtf.to_csv(all_gtf_file,sep="\t",header=None,index=None,quoting=csv.QUOTE_NONE)
    intron_gtf.to_csv(intron_gtf_file,sep="\t",header=None,index=None,quoting=csv.QUOTE_NONE)
    intron_to_transcript.to_csv(intron_to_transcript_file,sep="\t",header=None,index=None)
