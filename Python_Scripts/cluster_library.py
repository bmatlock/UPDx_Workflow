#!/bin/python
import re
import argparse
import pandas as pd




parser = argparse.ArgumentParser()
parser.add_argument('--clstr_file', type=argparse.FileType('r'), required=True, help='The clstr file to be analyzed')
parser.add_argument('--fasta_file', type=argparse.FileType('r'), required=True, help='The fasta file to be analyzed')
parser.add_argument('--summary_table', type=argparse.FileType('a'), required=True, help='Name of the count summary file to add the representative sequences to')
args = parser.parse_args()
Library = args.clstr_file.readlines()
Sequences = args.fasta_file.readlines()
LibraryIndex = 0
while LibraryIndex < len(Library):
    FC = re.search(r'^(>Cluster.*)', Library[LibraryIndex])
    if FC != None:
        firstcol = FC.group().replace('>','')
        args.summary_table.write(firstcol+',')
        LibraryIndex = LibraryIndex + 1
        while LibraryIndex < len(Library):
            SC = re.search(r'^(.*\*)', Library[LibraryIndex])
            if SC != None:
                seccol = SC.group()
                seccol = re.sub(r'.*>','',seccol)
                seccol = re.sub(r'\.\.\. \*','',seccol)
                args.summary_table.write(seccol+'\n')
                break
            LibraryIndex = LibraryIndex + 1
    LibraryIndex = LibraryIndex + 1
CountIndex = 0
args.clstr_file.close()









