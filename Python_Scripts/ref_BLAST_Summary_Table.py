#!/bin/python
#blast_output /Users/brycematlock/Documents/UPDx_28s_Workflow/Results/Experiment_10-11-23-13_55_23/sequence_analysis_10112023-13_55_23/Negative_Reads/ref_blast/UPDx12neg1_S12_L001..parasite_BLAST2_results.modified
#cluster_summary /Users/brycematlock/Documents/UPDx_28s_Workflow/Results/Experiment_10-11-23-13_55_23/FINAL_run_summary_10112023-13_55_23/Negatives/UPDx12neg1_S12_L001..summary/Cluster_Summary2_Table.csv
#Complete_Summary_Table /Users/brycematlock/Documents/negcomplete.csv
#BLAST_Summary_Table /Users/brycematlock/Documents/negblast.csv

import re
import argparse
import pandas as pd


parser = argparse.ArgumentParser()
parser.add_argument('--blast_output', type=str, required=True, help='The full reference blast output that needs to be analyzed')
parser.add_argument('--cluster_summary', type=argparse.FileType('r'), required=True, help='The full cluster summary table that needs to be merged')
parser.add_argument('--Complete_Summary_Table', type=argparse.FileType('a'), required=True, help='The full blast summary table that needs to be created')
parser.add_argument('--BLAST_summary_table', type=argparse.FileType('a'), help='The new parasite BLAST output table that needs to be created')
args = parser.parse_args()

SeqID_Library = open(args.blast_output)
SeqID_Contents = SeqID_Library.read()
pattern = re.compile(r'>([A-Z0-9:-]+)\t[0-9.]+\t[A-Za-z0-9-_.():{}[\]<>*,]+\t[A-Z-]+.*[0-9a-z-.]+\t[0-9]+')
matches = pattern.findall(SeqID_Contents)
seqid_list = []
for match in matches:
    seqid_list.append(match)
    SeqID_DF = pd.DataFrame(seqid_list, columns=['SeqID'])


Pident_Library = open(args.blast_output)
Pident_Contents = Pident_Library.read()
pattern = re.compile(r'>[A-Z0-9:-]+\t([0-9.]+)\t[A-Za-z0-9-_.():{}[\]<>*,]+\t[A-Z-]+.*[0-9a-z-.]+\t[0-9]+')
matches = pattern.findall(Pident_Contents)
pident_list = []
for match in matches:
    pident_list.append(match)
    Pident_DF = pd.DataFrame(pident_list, columns=['PIDENT'])
    
    

Temporary_Table1 = pd.concat([SeqID_DF.reset_index(drop=True), Pident_DF.reset_index(drop=True)], axis=1)

SseqID_Library = open(args.blast_output)
SseqID_Contents = SseqID_Library.read()
pattern = re.compile(r'>[A-Z0-9:-]+\t[0-9.]+\t([A-Za-z0-9-_.():{}[\]<>*,]+)\t[A-Z-]+.*[0-9a-z-.]+\t[0-9]+')
matches = pattern.findall(SseqID_Contents)
sseqid_list = []
for match in matches:
    sseqid_list.append(match)
    SseqID_DF = pd.DataFrame(sseqid_list, columns=['SseqID'])

Temporary_Table2 = pd.concat([Temporary_Table1.reset_index(drop=True), SseqID_DF.reset_index(drop=True)], axis=1)


Sequence_Library = open(args.blast_output)
Sequence_Contents = Sequence_Library.read()
pattern = re.compile(r'>[A-Z0-9:-]+\t[0-9.]+\t[A-Za-z0-9-_.():{}[\]<>*,]+\t([A-Z-]+).*[0-9a-z-.]+\t[0-9]+')
matches = pattern.findall(Sequence_Contents)
sequence_list = []
for match in matches:
    sequence_list.append(match)
    Sequence_DF = pd.DataFrame(sequence_list, columns=['BLAST_Sequence'])

Temporary_Table3 = pd.concat([Temporary_Table2.reset_index(drop=True), Sequence_DF.reset_index(drop=True)], axis=1)

Evalue_Library = open(args.blast_output)
Evalue_Contents = Evalue_Library.read()
pattern = re.compile(r'>[A-Z0-9:-]+\t[0-9.]+\t[A-Za-z0-9-_.():{}[\]<>*,]+\t[A-Z-]+\t([0-9a-z-.]+)\t[0-9]+')
matches = pattern.findall(Evalue_Contents)
evalue_list = []
for match in matches:
    evalue_list.append(match)
    Evalue_DF = pd.DataFrame(evalue_list, columns=['E-Value'])

Temporary_Table4 = pd.concat([Temporary_Table3.reset_index(drop=True), Evalue_DF.reset_index(drop=True)], axis=1)

Qcovg_Library = open(args.blast_output)
Qcovg_Contents = Qcovg_Library.read()
pattern = re.compile(r'>[A-Z0-9:-]+\t[0-9.]+\t[A-Za-z0-9-_.():{}[\]<>*,]+\t[A-Z-]+\t[0-9a-z-.]+\t([0-9]+)')
matches = pattern.findall(Qcovg_Contents)
qcovg_list = []
for match in matches:
    qcovg_list.append(match)
    Qcovg_DF = pd.DataFrame(qcovg_list, columns=['Query_Coverage'])


BLAST_Summary_Table = pd.concat([Temporary_Table4.reset_index(drop=True), Qcovg_DF.reset_index(drop=True)], axis=1)

#filtered_BLAST_Summary_Table1= BLAST_Summary_Table.loc[BLAST_Summary_Table['SseqID'].str.contains("Fungi") == False]


BLAST_Summary_Table['BLAST_Sequence'] = BLAST_Summary_Table['BLAST_Sequence'].str.wrap(15)

BLAST_Summary_Table['Query_Coverage'] = BLAST_Summary_Table['Query_Coverage'].astype('int')



BLAST_Summary_Table.to_csv(args.BLAST_summary_table)

Cluster_Summary_Table = pd.read_csv(args.cluster_summary)

Complete_Summary_Table = pd.merge(Cluster_Summary_Table, BLAST_Summary_Table, on='SeqID')
#Complete_Summary_Table = Complete_Summary_Table1.sort_values(['PIDENT'], ascending=False)

#filtered_Complete_Summary_Table = Complete_Summary_Table.loc[Complete_Summary_Table['SseqID'].str.contains("Fungi") == False]

Complete_Summary_Table.to_csv(args.Complete_Summary_Table)


