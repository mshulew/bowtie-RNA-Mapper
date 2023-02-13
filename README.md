# bowtie-aligner
- Maps reads from fastq files against a fasta file of RNAs  
- Can be used to map to miRNAs, rRNAs or any other RNA with a known sequence  
- Accepts compressed (.fastq.gz) and uncompressed (.fastq) files  
- When using UMI deduplication, must supply R1 & R2 file  
- When skipping UMI deduplication, will ignore R2 file  

## rRNA mapper
map reads to rRNAs  

Usage:  
required command: rRNAmapper --I {list of fastq file(s) to analyze} --O {output directory}  
optional: --reverseStrand --noTrim --allowUmi  
          reverseStrand = map to reverse strand  
          noTrim = skip trimming step  
          allowUmi = allow umi deduplication (requires R1 & R2 files)  
          
## miRNA mapper
map reads to miRNAs    

Usage:  
required command: miRNAmapper --I {list of fastq file(s) to analyze} --O {output directory}  
optional: --spikeType (miltenyi, ercc, miltenyiercc) --reverseStrand --noTrim --skiUmi  
          spikeType = default is natural, can change to miltenyi, ercc or both  
          reverseStrand = map to reverse strand  
          noTrim = skip trimming step  
          skipUmi = skip UMi deduplication (default is for UMI deduplication)  

## to build docker container and push to Docker Hub
git clone https://github.com/mshulewitz/bowtie_aligner.git  
cd bowtie_aligner  
docker build -t mshulew/bowtie-aligner .  
docker push mshulew/bowtie-aligner:latest

## to pull docker container
docker pull mshulew/bowtie-aligner:latest

## to execute pipelines with Nextflow
nextflow run bowtie-aligner/main.nf -profile docker --max_cpus 16
--reads {path to input files}
--bowtie_fa {path to fasta file for bowtie index}
--pipeline (miRNA or rRNA)
--outDir {output path} 

Example command line to test without UMI deduplication:  
nextflow run bowtie_aligner/main.nf -profile docker --reads 'fastq_repo/short/Kit4RNA3.fastq' --bowtie_fa bowtie-aligner/ref/natural.fa --outDir test --pipeline miRNA --skipUmi  

Example command line to test with UMI deduplication:
nextflow run bowtie_aligner/main.nf -profile docker --max_cpus 16 --reads 'fastq_repo/postverification/V46_S4_L001_R{1,2}_001.fastq' --bowtie_fa mirna-search/natural.fa --outDir test --pipeline miRNA  

## For full help with Nextflow:  
nextflow run bowtie-aligner/main.nf --help

