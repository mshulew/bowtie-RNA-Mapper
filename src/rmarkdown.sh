#!/bin/bash
# runs R markdown

inDir=$1
reference_key=$2
pipeline=$3

Rscript -e 'rmarkdown::render(input = "htmlReport.Rmd", params = list(inDir = "./tmp/'$inDir'", pipeline = "'$pipeline'", rrna_reference_key = "'$reference_key'", sample_index = "'$inDir'"), clean = TRUE, output_format = "html_document", output_file = "'$inDir'_htmlReport.html")'
