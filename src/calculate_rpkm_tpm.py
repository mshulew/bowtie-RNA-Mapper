#!/usr/bin/env python3
# coding: utf-8

"""
calculates RPKM and RPM from featureCounts output file
"""

import sys

if __name__ == "__main__":

    version_num = 0.1

# parse arguments
    input_filename = sys.argv[1]
    output_filename = sys.argv[2]

# import featureCounts output file
    input_data = []
    with open(input_filename, 'r') as input_file:
        for line in input_file:
            if line[0] != '#' and 'Geneid' not in line:
                split_line = line.strip('\n').split('\t')
                gene_name = split_line[0]
                chromosome = split_line[1]
                start_pos = split_line[2]
                stop_pos = split_line[3]
                gene_strand = split_line[4]
                gene_length = int(split_line[5])
                read_counts = int(split_line[6])
                input_data.append([gene_name,chromosome,start_pos,stop_pos,gene_strand,gene_length,read_counts])
                
# calculate rpkm
# step 1: count all reads and divide by 1,000,000 (= per million scaling factor)
    total_reads = 0
    for gene in input_data:
        total_reads += gene[6]
    per_million_scaling_factor = total_reads/1000000
    
# step 2: divide each read count by scaling factor to calculate RPM (reads per million)
# step 3: divid RPM by length of gene in kilobase (= RPKM)
    for a in range(0, len(input_data)):
        if per_million_scaling_factor != 0:
            input_data[a].append(input_data[a][6]/per_million_scaling_factor)
            input_data[a].append(input_data[a][7]/(input_data[a][5]/1000))
        else:
            input_data[a].append(0)
            input_data[a].append(0)
     
# Calculate TPM
# step 1: divide read counts by length of each gene in kilobase (= RPK)
        input_data[a].append(input_data[a][6]/(input_data[a][5]/1000))
        
# step 2: add all RPK values and divide by 1,000,000 (= per million scaling factor)
    total_rpk = 0
    for gene in input_data:
        total_rpk += gene[9]
    per_million_scaling_factor = total_rpk/1000000
    
# step 3: divide rpk values by per million scaling factor
    for b in range(0, len(input_data)):
        if per_million_scaling_factor != 0:
            input_data[b].append(input_data[b][9]/per_million_scaling_factor)
        else: 
            input_data[b].append(0)
            
# write to file
    outputfile = open(output_filename, 'w')
    outputfile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('Geneid','Chr','Start','End','Strand','Length','counts','rpkm','tpm'))
    for gene in input_data:
        outputfile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(gene[0],
                                                                       gene[1],
                                                                       gene[2],
                                                                       gene[3],
                                                                       gene[4],
                                                                       gene[5],
                                                                       gene[6],
                                                                       gene[8],
                                                                       gene[10]))
    outputfile.close()    

        
        
      
        
                
        
