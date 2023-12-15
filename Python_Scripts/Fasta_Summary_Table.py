#!/bin/python
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--specimen_name', type=str, required=True, help='The name of the specimen for which the fasta tables are being created')
parser.add_argument('--raw_fasta_txt_1', type=argparse.FileType('r'), required=True, help='The raw quality check text file for read 1 that needs to be converted to a table')
parser.add_argument('--raw_fasta_table_1', type=argparse.FileType('a'), required=True, help='The raw quality check table for read 1 that needs to be created')
parser.add_argument('--raw_fasta_txt_2', type=argparse.FileType('r'), required=True, help='The raw quality check text file for read 2 that needs to be converted to a table')
parser.add_argument('--raw_fasta_table_2', type=argparse.FileType('a'), required=True, help='The raw quality check table for read 2 that needs to be created')
parser.add_argument('--clean_fasta_txt_1', type=argparse.FileType('r'), required=True, help='The clean quality check text file for read 1 that needs to be converted to a table')
parser.add_argument('--clean_fasta_table_1', type=argparse.FileType('a'), required=True, help='The clean quality check table for read 1 that needs to be created')
parser.add_argument('--clean_fasta_txt_2', type=argparse.FileType('r'), required=True, help='The clean quality check text file for read 2 that needs to be converted to a table')
parser.add_argument('--clean_fasta_table_2', type=argparse.FileType('a'), required=True, help='The clean quality check table for read 2 that needs to be created')
args = parser.parse_args()

colnames = ['FQC Stat', 'Result']
title_row_raw_read1 = pd.DataFrame({'FQC Stat' : args.specimen_name, 'Result' : 'Raw-Read1'}, index=[0])
title_row_raw_read2 = pd.DataFrame({'FQC Stat' : args.specimen_name, 'Result' : 'Raw-Read2'}, index=[0])
title_row_clean_read1 = pd.DataFrame({'FQC Stat' : args.specimen_name, 'Result' : 'Clean-Read1'}, index=[0])
title_row_clean_read2 = pd.DataFrame({'FQC Stat' : args.specimen_name, 'Result' : 'Clean-Read2'}, index=[0])

raw_fasta_table_1_temp = pd.read_csv(args.raw_fasta_txt_1, sep=' ', names=colnames)
raw_fasta_table_1 = pd.concat([title_row_raw_read1, raw_fasta_table_1_temp])
raw_fasta_table_1.to_csv(args.raw_fasta_table_1, index=False)


raw_fasta_table_2_temp = pd.read_csv(args.raw_fasta_txt_2, sep=' ', names=colnames)
raw_fasta_table_2 = pd.concat([title_row_raw_read2, raw_fasta_table_1_temp])
raw_fasta_table_2.to_csv(args.raw_fasta_table_2, index=False)


clean_fasta_table_1_temp = pd.read_csv(args.clean_fasta_txt_1, sep=' ', names=colnames)
clean_fasta_table_1 = pd.concat([title_row_clean_read1, raw_fasta_table_1_temp])
clean_fasta_table_1.to_csv(args.clean_fasta_table_1, index=False)


clean_fasta_table_2_temp = pd.read_csv(args.clean_fasta_txt_2, sep=' ', names=colnames)
clean_fasta_table_2 = pd.concat([title_row_clean_read2, raw_fasta_table_1_temp])
clean_fasta_table_2.to_csv(args.clean_fasta_table_2, index=False)
