#!/usr/bin/env python3
# coding: utf-8

"""
debarcode for SEQuoia Express
handles barcode at R1 or R2 -- update to use only R2
SE or PE reads
"""

import sys
from collections import Counter

if __name__ == "__main__":

# get arguments
    umi_location = sys.argv[1]

# get input data
    input_file_data = sys.stdin.read().strip().split('\n')

# move barcode from R1 to read name
# remove 1st base + barcode from R1
    good_counter = 0
    bad_counter = 0

    for a in range(0,len(input_file_data),4):
        if umi_location == 'r1':
            barcode = input_file_data[a+1].split('\t')[0][1:9]
        elif umi_location == 'r2':
            barcode = input_file_data[a+1].split('\t')[1][0:8]

        missing = 8 - len(barcode)

# if the length - number of N's is less than 7, throw out read
# basically allows an edit distance of 1
        if missing + Counter(barcode).get("N", 0) <= 1:
            good_counter += 1

# SE read
            if len(input_file_data[0].split('\t')) == 1:
                if umi_location == 'r1':
                    print('@{}_{}'.format(barcode,input_file_data[a].split('\t')[0][1:]))
                    print(input_file_data[a+1].split('\t')[0])
                    print(input_file_data[a+2].split('\t')[0])
                    print(input_file_data[a+3].split('\t')[0])

# PE reads
            if len(input_file_data[0].split('\t')) == 2:
                print('@{}_{}'.format(barcode,input_file_data[a].split('\t')[0][1:]))
                print('@{}_{}'.format(barcode,input_file_data[a].split('\t')[1][1:]))
                print(input_file_data[a+1].split('\t')[0])
                print(input_file_data[a+1].split('\t')[1][8:])
                print(input_file_data[a+2].split('\t')[0])
                print(input_file_data[a+2].split('\t')[1])
                print(input_file_data[a+3].split('\t')[0])
                print(input_file_data[a+3].split('\t')[1][8:])

        else:
            bad_counter += 1

# generate output log
    print('Total Reads\t{}'.format(str(good_counter + bad_counter)))
    print('Good Reads\t{}'.format(str(good_counter)))
    print('Bad Reads\t{}'.format(str(bad_counter)))
