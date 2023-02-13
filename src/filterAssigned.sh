#!/bin/bash
# filters Assigned reads from BAM file

samtools view $1 | grep "Assigned" | cut -f 1 | sed -e 's/^.*_//' > "$2"_assigned.txt
if ! [ -s "$2"_assigned.txt ]; then
  printf "NO_ASSIGNED_READS_FOUND" > "$2"_assigned.txt
fi
