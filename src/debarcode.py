#!/usr/bin/env python3
# coding: utf-8

"""
transfer barcode from R2 to R1
usage: paste <(zcat R1.fastq.gz ) <(zcat R2.fastq.gz ) python3.6 /opt/biorad/src/debarcode.py
usage: ./fastqmerge.sh fastq_repo/postverification/V46_S4_L001_R{1,2}_001.fastq.gz | python3.6 bowtie_aligner/src/debarcode.py \
| tee >(grep 'Reads' > debarcode_stats.txt) | grep -ve 'Reads' > debarcode.fastq
"""

import sys
from collections import Counter

if __name__ == "__main__":
    start, end = 0, 8

# get input data
    input_file_data = sys.stdin.read().strip().split('\n')
     
# process file data
    bad_counter = 0
    good_counter = 0
    for a in range(0,len(input_file_data),4):
        barcode = input_file_data[a+1].split('\t')[1][start:end]
        missing = 8 - len(barcode)
        # If the length - number of N's is less than 7, throw it out
        # basically allows an edit distance of 1
        if missing + Counter(barcode).get("N", 0) <= 1:
            good_counter += 1
            print('@{}_{}'.format(barcode,input_file_data[a].split('\t')[0][1:]))
            print(input_file_data[a+1].split('\t')[0])
            print(input_file_data[a+2].split('\t')[0])
            print(input_file_data[a+3].split('\t')[0])
        else:
            bad_counter += 1
            
# generate output log
    print('Total Reads\t{}'.format(str(good_counter + bad_counter)))
    print('Good Reads\t{}'.format(str(good_counter)))
    print('Bad Reads\t{}'.format(str(bad_counter)))
