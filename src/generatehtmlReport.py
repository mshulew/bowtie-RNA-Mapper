#!/usr/bin/env python3
# coding: utf-8

"""
Move output files into separate directories for each sample
Perform Rmarkdown for each sample to generate report
version 1.0 September 4, 2020
"""

import sys
import os
import subprocess
import shutil

if __name__ == "__main__":

# parse arguments
    input_dir = sys.argv[1]
    rrna_reference_key = sys.argv[2]
    pipeline = sys.argv[3]

# make temp dir
    os.makedirs('./tmp')
# iterate through input dir and transfer files for each sample into a separate directory
    for foldername in os.listdir(input_dir):
        for filename in os.listdir(input_dir + '/' + foldername):
            sample_id = filename.split('_')[0]
# make sample directory if not already made
            if not os.path.isdir('tmp/' + sample_id):
                os.makedirs('tmp/' + sample_id)
# make folder directory if not already make
            if not os.path.isdir('tmp/' + sample_id + '/' + foldername):
                os.makedirs('tmp/' + sample_id + '/' + foldername)
            cmd = 'cp {} {}'.format(input_dir + '/' + foldername + '/' + filename, 'tmp/' + sample_id + '/' + foldername + '/') 
            process = subprocess.Popen(cmd, shell=True, executable='/bin/bash')
            process.wait()
# iterate through sample folders and perform Rmarkdown - I use a bash script to run rmarkdown
    for foldername in os.listdir('tmp'):
        cmd = "bash /opt/biorad/src/rmarkdown.sh {} {} {}".format(foldername, rrna_reference_key, pipeline)
        process = subprocess.Popen(cmd, shell=True, executable='/bin/bash')
        process.wait()
# remove temp dir
    shutil.rmtree('./tmp')
        

            
                               
                
    
    
