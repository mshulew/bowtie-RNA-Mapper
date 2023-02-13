#!/usr/bin/env python3
# coding: utf-8

"""
makes gtf file from fasta file containing miRNa sequences
trims 3' As in order to be compatible with poly A trimming of reads
version 1.0 December 27, 2019
"""

import sys

if __name__ == "__main__":

# parse arguments
    if len(sys.argv) == 3:
        input_filename = sys.argv[1]
        output_filename = sys.argv[2]
    else:
        print('fasta to gtf missing input and/or output filenames')
        exit()
        
# import miRNA fasta sequences and create a gtf file
    mirna_gtf_data = []
    with open(input_filename, 'r') as mirna_fasta_file:
        for line in mirna_fasta_file:
            if line[0] == '>':
                mirna_name = line.strip('\n')[1:]
            else:
                mirna_sequence = line.strip('\n')
                
# trim 3 prime As for gtf file
                terminal_a_length = 0
                for a in range(5,0,-1):
                    if mirna_sequence[-a:] == 'A'*a:
                        terminal_a_length = a
                        break
                        
                sequence_for_gtf = mirna_sequence[:len(mirna_sequence) - terminal_a_length]
                
                mirna_gtf_annotation = 'gene_id "{}"; transcript_id "{}";'.format(mirna_name,mirna_name) 
                mirna_gtf_data.append([mirna_name,'miRNA','exon','1',str(len(sequence_for_gtf)),'0','+','.',mirna_gtf_annotation])
                
# write gtf file
    outputfile = open(output_filename, 'w')
    for entry in mirna_gtf_data:
        outputfile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(entry[0],
                                                                       entry[1],
                                                                       entry[2],
                                                                       entry[3],
                                                                       entry[4],
                                                                       entry[5],
                                                                       entry[6],
                                                                       entry[7],
                                                                       entry[8]))
    outputfile.close()
        
    
    
