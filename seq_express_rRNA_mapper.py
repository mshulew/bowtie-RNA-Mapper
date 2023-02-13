#!/usr/bin/env python3
# coding: utf-8

"""
rRNA Mapper for SEQuoia-Express
"""
import sys
import os
import subprocess
import shutil

if __name__ == "__main__":

    version_num = 1.0

    print('*'*100)
    print('SEQuoia-Express rRNA Mapper Launcher')
    print('2020 Mark Shulewiz Bio-Rad Laboratories, Inc.')
    print('version ' + str(version_num))
    print('*'*100)

# parse arguments
    input_filenames = []
    output_dir = os.getcwd()
    trimming = True
    ref_dir = '~/efs/bowtie-aligner'
    skipumi = False
    umiType = 'R2'
    reversestrand = False
    fracoverlap = '0'
    genome = 'hg38'
    singleend = False

    if len(sys.argv) == 2:
        if sys.argv[1].lower() == '--h' or sys.argv[1].lower() == '--help':
            print('Maps FASTQ reads to human rRNA sequences')
            print('Usage: ' + os.path.basename(__file__) + ' --I [input filename] --O [output directory name]')
            print('Options: --umiType --genome --reverseStrand --noTrim --skipUmi')  
            print('--umiType for UMI location (R1, R2 or none; default is R1)')
            print('--skipUmi to skip UMI deduplication')
            print('--SE for single end reads; default is paired end reads')
            print('--genome to select species (hg38, mm10, rnor6; default is hg38)')
            print('--noTrim to turn off poly A trimming and removal of first 5 prime base')
            print('--reverseStrand for kits that generate reverse complement reads')
            print('--fracOverlap = minimum fraction of overlap between reads and genes; default is 0')
            print('--removeReads = generate fastq file(s) with assigned reads removed')
            print('*'*100)
            exit()       
    elif len(sys.argv) > 2:
        input_trigger = False
        for a in range(0, len(sys.argv)):
            if input_trigger == True:
                if sys.argv[a][0:2] == '--':
                    input_trigger = False
                else:
                    input_filenames.append(sys.argv[a])
            if input_trigger == False:
                if a < len(sys.argv) - 1:
                    if sys.argv[a].lower() == '--i':
                        input_trigger = True
                    elif sys.argv[a].lower() == '--o':
                        output_dir = sys.argv[a + 1]
                    elif sys.argv[a] == '--refs':
                        ref_dir = sys.argv[a + 1]
                    elif sys.argv[a] == '--fracOverlap':
                        fracoverlap = sys.argv[a + 1]
                    elif sys.argv[a]== '--genome':
                        genome = sys.argv[a + 1]
                    elif sys.argv[a] == '--umiType':
                        umiType = sys.argv[a + 1]
                if sys.argv[a] == '--reverseStrand':
                    reversestrand = True
                if sys.argv[a] == '--noTrim':
                    trimming = False
                if sys.argv[a] == '--skipUmi':
                    skipumi = True
                if sys.argv[a] == '--SE':
                    singleend = True

    else:
        print('missing parameters/options')
        exit()

# copy input files into a temp directory
    if os.path.isdir(output_dir):
        shutil.rmtree(output_dir)
    tmp_dir = output_dir + '/tmp'
    os.makedirs(tmp_dir,exist_ok=True)

    for filename in input_filenames:
        if os.path.exists(filename):
            base_filename = filename.split('/')[(len(filename.split('/'))-1)]
            shutil.copy(filename,tmp_dir + '/' + base_filename)
        else:
            print('file not found: {}'.format(filename))
            exit()

# create options for command line
    if umiType == 'A':
        umiType = 'none'
    elif umiType == 'B':
        umiType = 'R1'
    elif umiType == 'C':
        umiType = 'R2'

    if umiType == 'none':
        skipumi = True
    else:
        umiType = umiType.upper()
    cl_options = ' --fracOverlap ' + fracoverlap +  ' --genome ' + genome + ' --umiType ' + umiType
    if skipumi:
        cl_options = cl_options + ' --skipUmi'
    if reversestrand:
        cl_options = cl_options + ' --reverseStrand'
    if not trimming:
        cl_options = cl_options + ' --noTrim' 
    if singleend:
        cl_options = cl_options + ' --SE' 

    cmd = "nextflow run {} -profile docker --inDir {} --outDir {} {}".format(ref_dir + '/bowtie-express.nf',tmp_dir,output_dir,cl_options)
    process = subprocess.Popen(cmd, shell=True, executable='/bin/bash')
    process.wait()

# delete temporary directory
    shutil.rmtree(tmp_dir)
