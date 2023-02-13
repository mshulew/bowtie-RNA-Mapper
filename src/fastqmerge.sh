#!/bin/bash
# merges R1 & R2 files

if [[ "$1" =~ \.gz$ ]] && [[ "$2" =~ \.gz$ ]]; then 
        paste <(zcat "$1") <(zcat "$2")
else
        paste "$1" "$2"
fi
       
