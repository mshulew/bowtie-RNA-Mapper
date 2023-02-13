#!/usr/bin/env python3
# coding: utf-8

# transfer UMI barcode from query name to tag

import pysam
import sys


if __name__ == "__main__":
    samfile=pysam.AlignmentFile(sys.argv[1], "rb")
    outputfile=sys.argv[2]

    outputsam = pysam.AlignmentFile(outputfile, "wb", template=samfile)
    for read in samfile:
        spl = sorted(read.query_name.split("_"), key=len)
        bc = spl[0]
        rd = spl[1]
        read.set_tag("XU", bc)
        read.query_name = rd
        outputsam.write(read)
    outputsam.close()
        

    
	
