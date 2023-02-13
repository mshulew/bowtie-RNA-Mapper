#!/usr/bin/env python3
# coding: utf-8

"""
split merged fastq into separate files
usage: usage: cat mergedfastq | python3.6 splitmerge.py [r1 output file] [r2 output file]
"""

import sys
from collections import Counter

if __name__ == "__main__":
    r1_filename = sys.argv[1]
    r2_filename = sys.argv[2]
    
    r1_data = []
    r2_data = []

# get input data
    input_file_data = sys.stdin.read().strip().split('\n')
    
        
# process file data
    if len(input_file_data) > 1:
        for entry in input_file_data:
            r1_data.append(entry.split('\t')[0])
            r2_data.append(entry.split('\t')[1])

# write to file
    output_file = open(r1_filename, 'w')
    for entry in r1_data:
        output_file.write('{}\n'.format(entry))
    output_file.close()
    
    output_file = open(r2_filename, 'w')
    for entry in r2_data:
        output_file.write('{}\n'.format(entry))
    output_file.close()
    
            
            
              
     
