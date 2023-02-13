#!/usr/bin/env python3
# coding: utf-8

"""
allows tRNAmapper to accept multiple samples & consolidates output into a single report
launches tRNAmapper through nextflow 
"""
import sys
import os
import subprocess
import shutil

if __name__ == "__main__":

    version_num = 1.0

    print('*'*100)
    print('tRNA Mapper Launcher')
    print('2020 Mark Shulewiz Bio-Rad Laboratories, Inc.')
    print('version ' + str(version_num))
    print('*'*100)

# parse arguments
    input_filenames = []
    output_dir = os.getcwd()
    trimming = True
    ref_dir = '~/efs/bowtie-aligner'
    skipumi = False
    reversestrand = False
    removereads = False
    fracoverlap = '0'
    
    if len(sys.argv) == 2:
        if sys.argv[1] == '--h' or sys.argv[1] == '--help':
            print('*'*100)
            print(os.path.basename(__file__) + ' version ' + str(version_num))
            print('Mark Shulewitz 2019')
            print('Bio-Rad Laboratories, Inc.')
            print('Usage: ' + os.path.basename(__file__) + ' --I [input files] --O [output directory]')
            print('--O: output directory (optional; default is current directory)')
            print('--skipUmi (do not use UMI deduplication)')
            print('--fracOverlap = minimum fraction of overlap between reads and genes; default is 0')
            print('--removeReads (generate fastq file(s) with assigned reads removed)')
            print('--reverseStrand (map to reverse strand)')
            print('--noTrim (do not trim)')
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
                    elif sys.argv[a].lower() == '--refs':
                        ref_dir = sys.argv[a + 1]
                    elif sys.argv[a].lower() == '--fracoverlap':
                        fracoverlap = sys.argv[a + 1]
                if sys.argv[a].lower() == '--reversestrand':
                    reversestrand = True
                if sys.argv[a].lower() == '--notrim':
                    trimming = False
                if sys.argv[a].lower() == '--skipumi':
                    skipumi = True
                if sys.argv[a].lower() == '--removereads':
                    removereads = False
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
    cl_options = ' --fracOverlap ' + fracoverlap
    if skipumi:
        cl_options = cl_options + ' --skipUmi'
    if reversestrand:
        cl_options = cl_options + ' --reversestrand'
    if not trimming:
        cl_options = cl_options + ' --notrim'
    if removereads:
        cl_options = cl_options + ' --removeReads'
        
    cmd = "nextflow run {} -profile docker --inDir '{}' --outDir {} --pipeline tRNA {}".format(ref_dir + '/main.nf',tmp_dir,output_dir,cl_options)
    process = subprocess.Popen(cmd, shell=True, executable='/bin/bash')
    process.wait()
    
# delete temporary directory
    shutil.rmtree(tmp_dir)
