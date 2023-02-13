#!/usr/bin/env python3
# coding: utf-8

"""
generates stats on RNAs detected and reads mapping to an RNA
"""

import sys

if __name__ == "__main__":

    gene_count_filename = sys.argv[1]
    threshold_read_num = int(sys.argv[2])
    output_filename = sys.argv[3]
    
# import gene count data
    gene_count_data = []
    with open(gene_count_filename, 'r') as gene_count_file:
        for line in gene_count_file:
            if line[0] != '#' and 'Geneid' not in line:
                split_line = line.strip('\n').split('\t')
                gene_name = split_line[0]
                chromosome = split_line[1]
                start_pos = split_line[2]
                stop_pos = split_line[3]
                gene_strand = split_line[4]
                gene_length = int(split_line[5])
                read_counts = int(split_line[6])
                rpkm = float(split_line[7])
                tpm = float(split_line[8])
                gene_count_data.append([gene_name,chromosome,start_pos,stop_pos,gene_strand,gene_length,read_counts,rpkm,tpm])
                
# count: total number of RNAs with >0 reads
#        total number of reads mapped to an RNA
#        total number of RNAs with >= threshold reads
#        total number of reads mapped to an RNA with >= threshold reads
    total_rnas = 0
    total_reads_mapped = 0
    total_rnas_meeting_threshold = 0
    total_reads_mapped_to_rna_meeting_threshold = 0
    
    for entry in gene_count_data:
        if entry[6] > 0:
            total_rnas += 1
            total_reads_mapped += entry[6]
            if entry[6] >= 5:
                total_rnas_meeting_threshold += 1
                total_reads_mapped_to_rna_meeting_threshold += entry[6]
                
    # write to file
    output_file = open(output_filename, 'w')
    output_file.write('Total number of reads mapping to an RNA\t{}\n'.format(total_reads_mapped))
    output_file.write('Total number of RNAs detected\t{}\n'.format(total_rnas))
    output_file.write('Total number of reads mapping to an RNA with >={} reads\t{}\n'.format(str(threshold_read_num),total_reads_mapped_to_rna_meeting_threshold))
    output_file.write('Total number of RNAs with >={} reads detected\t{}\n'.format(str(threshold_read_num),total_rnas_meeting_threshold))
    output_file.close()   
