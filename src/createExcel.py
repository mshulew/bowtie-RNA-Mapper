#!/usr/bin/env python3
# coding: utf-8

"""
combined tsv files into a single Excel spreadsheet
"""

import sys
import pandas as pd

if __name__ == "__main__":

# parse arguments
    counts_filename = sys.argv[1]
    rpkm_filename = sys.argv[2]
    tpm_filename = sys.argv[3]
    
# create dataframes
    counts_df = pd.read_csv(counts_filename, sep='\t', index_col=[0])
    rpkm_df = pd.read_csv(rpkm_filename, sep='\t', index_col=[0])
    tpm_df = pd.read_csv(tpm_filename, sep='\t', index_col=[0])
    
# save to Excel
    writer = pd.ExcelWriter('consolidated_countsRPKMTPM.xlsx', engine='xlsxwriter')
    counts_df.to_excel(writer, sheet_name='counts')
    rpkm_df.to_excel(writer, sheet_name='rpkm')
    tpm_df.to_excel(writer, sheet_name='tpm')
    writer.save()
