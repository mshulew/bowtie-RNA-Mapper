#!/usr/bin/env python3
# coding: utf-8

"""
split singleline merged fastq into separate files
"""

import sys

if __name__ == "__main__":
    r1_filename = sys.argv[1]
    if len(sys.argv) > 2:
        r2_filename = sys.argv[2]
    
    r1_data = []
    r2_data = []

# get input data
    input_file_data = sys.stdin.read().strip().split('\n')
    
# process file data
    if len(input_file_data) > 1:
        for entry in input_file_data:
            r1_data.append(entry.split('\t')[0])
            r1_data.append(entry.split('\t')[1])
            r1_data.append(entry.split('\t')[2])
            r1_data.append(entry.split('\t')[3])
            if len(sys.argv) > 2:
                r2_data.append(entry.split('\t')[4])
                r2_data.append(entry.split('\t')[5])
                r2_data.append(entry.split('\t')[6])
                r2_data.append(entry.split('\t')[7])

# write to file
    output_file = open(r1_filename, 'w')
    for entry in r1_data:
        output_file.write('{}\n'.format(entry))
    output_file.close()
    
    if len(sys.argv) > 2:
        output_file = open(r2_filename, 'w')
        for entry in r2_data:
            output_file.write('{}\n'.format(entry))
        output_file.close()
    
            
            
              
     
