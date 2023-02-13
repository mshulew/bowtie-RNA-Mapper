#!/usr/bin/env nextflow

"""
bowtie aligner pipeline
2020 Mark Shulewitz Bio-Rad Laboratories, Inc.

Processes SEQuoia Express fastq files
SE or PE

"""

def helpMessage() {
  log.info Header()
  log.info """
  Usage:

  The typical command for running the pipeline is as follows:

  nextflow run bowtie-aligner/main.nf -profile docker \
  --inDir{path to folder with fastq files}
  --outDir {path to output directory}

  Required:
    --inDir                 Path to directory containing input fastq file(s)
    --outDir                Path ot output directory

  Optional arguments:
    --genome                reference genome for rRNA pipeline (hg38, mm10, rnor6; default is hg38)
    --skipUmi               Skip UMI deduplication
    --SE                    Single end reads (default is for paired end reads)
    --umiType               UMI location (R1, R2 or none; default is R1)
    --reverseStrand         Reads analyzed are reverse complement of coding strand
    --noTrim                Do not trim reads
    --count_threshold       Minimum number of reads that must map to an RNA; default is 5
    --fracOverlap           Minimum fractions of overlap between reads and genes; default is 0
    --max_cpus              Maximum number of CPUs available; default = 16
    --max_memory            Maximum memory available; default = 60
    """.stripIndent()
}

// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}

// Housekeeping
allowedGenomes = ["hg38","mm10","rnor6","hg38-45"]

if (params.genome != "NONE" && !allowedGenomes.contains(params.genome)) {
    log.error "$params.genome is not in $allowedGenomes"
    exit 1
}

//Set fasta file based on pipeline and spikeType

if (params.genome == "hg38"){
  bowtie_fa = file("${baseDir}/ref/hs_rRNA_references.fa", checkIfExists: true)
  rrna_reference_key = file("${baseDir}/ref/hs_rRNA_reference_key.tsv", checkIfExists: true)
}
else if (params.genome == "hg38-45"){
  bowtie_fa = file("${baseDir}/ref/hs_45SrRNA_references.fa", checkIfExists: true)
  rrna_reference_key = file("${baseDir}/ref/hs_45SrRNA_reference_key.tsv", checkIfExists: true)
}
else if (params.genome == "mm10"){
  bowtie_fa = file("${baseDir}/ref/mm_rRNA_references.fa", checkIfExists: true)
  rrna_reference_key = file("${baseDir}/ref/mm_rRNA_reference_key.tsv", checkIfExists: true)
}
else if (params.genome == "rnor6"){
  bowtie_fa = file("${baseDir}/ref/rnor_rRNA_references.fa", checkIfExists: true)
  rrna_reference_key = file("${baseDir}/ref/rnor_rRNA_reference_key.tsv", checkIfExists: true)
}

def summary = [:]
summary['Run Name'] = workflow.runName
if(workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Config Profile'] = workflow.profile
summary['Input directory'] = params.inDir
summary['Genome'] = params.genome
summary['fasta file'] = bowtie_fa
summary['Output dir'] = params.outDir
summary['skipUmi'] = params.skipUmi
summary['SE reads?'] = params.SE
summary['UMI location'] = params.umiType
if(params.reverseStrand) summary['Strand'] = 'Reverse'
else summary['Strand'] = 'Forward'
if(params.noTrim) summary['Trimming'] = 'None'
summary['Fraction Overlap'] = params.fracOverlap
summary['Cores'] = params.max_cpus
summary['Memory'] = params.max_memory
log.info Header()
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "----------------------------------------------------"

// Add input files to channels

if (!params.SE) {
    Channel
        .fromFilePairs("$params.inDir/*{R1,R2}*", size: 2, checkIfExists: true)
        .ifEmpty { exit 1, "R1 and/or R2 file(s) can not be found" }
        .into { pipelineinfo_ch; fastqc_ch; input_reads_ch }
} else {
    Channel
        .fromFilePairs("$params.inDir/*R1*", size: 1, checkIfExists: true)
        .ifEmpty { exit 1, "R1 file can not be found" }
        .into { pipelineinfo_ch; fastqc_ch; input_reads_ch }
}

process pipelineInfo {
    publishDir "${params.outDir}/$sample_id", mode: 'copy'
  
    input:
    set val(sample_id), file(reads) from pipelineinfo_ch
  
    output:
    file '*_pipeline_info.tsv' into pipelineinfo_info_ch
  
    script:
    sample_id = sample_id.split( '_' )[ 0 ]
    reads_string = reads.toString()
    """
    printf "Metric\tValue\n" > ${sample_id}_pipeline_info.tsv
    printf "Pipeline\t${params.pipeline}\n" >> ${sample_id}_pipeline_info.tsv
    printf "Name\t${workflow.manifest.name}\n" >> ${sample_id}_pipeline_info.tsv
    printf "Description\t${workflow.manifest.description}\n" >> ${sample_id}_pipeline_info.tsv
    printf "Author\t${workflow.manifest.author}\n" >> ${sample_id}_pipeline_info.tsv
    printf "Version\t${workflow.manifest.version}\n" >> ${sample_id}_pipeline_info.tsv
    printf "Input reads\t${params.inDir}\n" >> ${sample_id}_pipeline_info.tsv
    printf "Reads processed\t${reads_string}\n" >> ${sample_id}_pipeline_info.tsv
    printf "Output directory\t${params.outDir}\n" >> ${sample_id}_pipeline_info.tsv
    printf "Reference fasta file\t${bowtie_fa}\n" >> ${sample_id}_pipeline_info.tsv
    printf "UMI location\t${params.umiType}\n" >> ${sample_id}_pipeline_info.tsv
    printf "Single end reads\t${params.SE}\n" >> ${sample_id}_pipeline_info.tsv
    printf "Skip UMI\t${params.skipUmi}\n" >> ${sample_id}_pipeline_info.tsv
    printf "Fraction overlap\t${params.fracOverlap}\n" >> ${sample_id}_pipeline_info.tsv
    printf "Reverse Strand\t${params.reverseStrand}\n" >> ${sample_id}_pipeline_info.tsv
    printf "No trimming\t${params.noTrim}\n" >> ${sample_id}_pipeline_info.tsv
    """   
}

process fastQc {
    tag "fastqc on $sample_id"
    publishDir "${params.outDir}/$sample_id/fastqc", mode: 'copy'

    input:
    set val(sample_id), file(reads) from fastqc_ch

    output:
    file '*_fastqc.{zip,html}' into fastqc_report_ch
  
    script:
    sample_id = sample_id.split( '_' )[ 0 ]
    """
    fastqc $reads
    """
}

if (!params.skipUmi) {
    process debarcode {
        tag "debarcode on $sample_id"
        label 'mid_memory'
        publishDir "${params.outDir}/$sample_id/debarcode", mode: 'copy', pattern: "*.txt"
    
        input:
        set sample_id, file(reads) from input_reads_ch
        
        output:
        set val(sample_id), file('*_debarcoded.fastq') into debarcoded_ch
        file '*_debarcode_stats.txt' into debarcode_stats_ch
        
        script:
        sample_id = sample_id.split( '_' )[ 0 ]
        if(params.umiType=='R1') {
            if(!params.SE) {
                (r1, r2) = reads
                """
                paste <(zcat $r1 | sed -e 's/\r\$//g') <(zcat $r2 | sed -e 's/\r\$//g') | \
                python3.6 /opt/biorad/src/express_debarcode.py r1 | \
                tee >(grep 'Reads' > ${sample_id}_debarcode_stats.txt) | grep -ve 'Reads' | \
                paste - - | python3.6 /opt/biorad/src/splitmerge.py ${sample_id}_r1_debarcoded.fastq ${sample_id}_r2_debarcoded.fastq
                """
            } else {
                """
                zcat $reads | sed -e 's/\r\$//g' | python3.6 /opt/biorad/src/express_debarcode.py r1 | \
                tee >(grep 'Reads' > debarcode_stats.txt) | grep -ve 'Reads' > ${sample_id}_r1_debarcoded.fastq
                """
            }
        } else if(params.umiType == 'R2') {
            (r1, r2) = reads
            """
            paste <(zcat $r1 | sed -e 's/\r\$//g') <(zcat $r2 | sed -e 's/\r\$//g') | \
            python3.6 /opt/biorad/src/express_debarcode.py r2 | \
            tee >(grep 'Reads' > ${sample_id}_debarcode_stats.txt) | grep -ve 'Reads' | \
            paste - - | python3.6 /opt/biorad/src/splitmerge.py ${sample_id}_r1_debarcoded.fastq ${sample_id}_r2_debarcoded.fastq
            """
        }
    }
} else {
    debarcoded_ch = input_reads_ch
    debarcode_stats_ch = Channel.empty()
}

process cutadapt {
    tag "cutadapt on $sample_id"
    label 'mid_cpu'
    publishDir "${params.outDir}/$sample_id/cutadapt", mode: 'copy', pattern: "*.{fastq,log}"
  
    input:
    set val(sample_id), file(reads) from debarcoded_ch
  
    output:
    set val(sample_id), file('*_trimmed.fastq') into trimmed_file_ch
    file '*_cutadapt.log' into trimmed_log_ch
  
    script:
    cutting = "-u 9"
  
    if(params.noTrim) {
        cutting = ""
    } else {
        if(!params.SE){
            if(params.umiType == 'R2'){
                cutting = "-u 1 -U 8"
            } else if(params.umiType == 'none'){
                cutting = "-u 1"
            }
            
        } else {
            if(params.umiType == 'none'){
                cutting = "-u 1"
            }
        }
    }
    
    if (params.skipUmi){
        sample_id = sample_id.split( '_' )[ 0 ]
    }
    
    if(!params.SE){
        (r1,r2) = reads 
        """
        cutadapt ${cutting} -m 15 -j $task.cpus -q 0,0 -o ${sample_id}_R1_trimmed.fastq -p ${sample_id}_R2_trimmed.fastq $r1 $r2 1> ${sample_id}_cutadapt.log
        gzip --keep *_trimmed.fastq
        """
    } else {
        """
        cutadapt ${cutting} -m 15 -j $task.cpus -q 0,0 -o ${sample_id}_R1_trimmed.fastq ${reads} 1> ${sample_id}_cutadapt.log
        gzip --keep *_trimmed.fastq
        """
    }
}

process bowtieIndex {
  publishDir "${params.outDir}/bowtieIndex", mode: 'copy'

  input:
  file bowtie_fa from bowtie_fa
  
  output:
  file 'bowtie_index.*' into bowtie_index_ch
  file 'bowtie.gtf' into bowtie_gtf_ch
  
  script:
  """
  bowtie-build ${bowtie_fa} bowtie_index
  python3.6 /opt/biorad/src/gtf_from_fasta.py ${bowtie_fa} bowtie.gtf
  """
}

process bowtie {
  tag "bowtie on $sample_id"
  label 'mid_cpu'
  publishDir "${params.outDir}/$sample_id/bowtie", mode: 'copy'
  
  input:
  set val(sample_id), file(trimmed) from trimmed_file_ch
  file bowtie_index from bowtie_index_ch
 
  output:
  set val(sample_id), file('*_bowtie.bam') into bowtie_ch
  file '*_bowtie.log' into bowtie_log_ch
  file '*_bowtie.bam'
  
  script:
  
  if(!params.SE){
  (r1,r2) = trimmed
  """
  bowtie bowtie_index -1 ${r1} -2 ${r2} --threads $task.cpus --sam bowtie.sam --maxins 1000 2> ${sample_id}_bowtie.log
  samtools view -Sb bowtie.sam > temp.bam
  samtools sort temp.bam -o ${sample_id}_bowtie.bam
  sambamba index -t $task.cpus ${sample_id}_bowtie.bam
  """
  } else {
  """
  bowtie bowtie_index ${trimmed} --threads $task.cpus --sam bowtie.sam 2> ${sample_id}_bowtie.log
  samtools view -Sb bowtie.sam > temp.bam
  samtools sort temp.bam -o ${sample_id}_bowtie.bam
  sambamba index -t $task.cpus ${sample_id}_bowtie.bam
  """
  }
}

if (!params.skipUmi) {
  process umiTagging {
    tag "umi_tagging on $sample_id"
    label 'mid_cpu'
    publishDir "${params.outDir}/$sample_id/umiTagging", mode: 'copy'
        
    input:
    set val(sample_id), file(bams) from bowtie_ch
        
    output:
    set val(sample_id), file('bowtie.tagged.bam*') into umitagged_ch
        
    script:
    (bam, bai) = bams
    """
    python3.6 /opt/biorad/src/tagFullBamFile.py ${bam} bowtie.tagged.bam
    sambamba index -t $task.cpus bowtie.tagged.bam
    """
  }
  process deduplication {
    label 'mid_cpu'
    label 'high_memory'
    tag "deduplication on $sample_id"
    publishDir "${params.outDir}/$sample_id/deduplication", mode: 'copy'
        
    input:
    set val(sample_id), file(bams) from umitagged_ch
        
    output:
    set val(sample_id), file('bowtie.deduplicated.out.bam*') into deduplicated_ch
    file '*_dedup.log' into dedup_log_ch
        
    script:
    (bam, bai) = bams
    """
    mkdir -p ./deduplicated
    umi_tools dedup -I $bam --paired --output-stats=./deduplicated \
    --method unique --log ./${sample_id}_dedup.log \
    --extract-umi-method=tag --umi-tag=XU > ./bowtie.deduplicated.out.bam
    sambamba index -t $task.cpus ./bowtie.deduplicated.out.bam
    printf "unique_input_reads: " >> ./${sample_id}_dedup.log; samtools view $bam | cut -f1 | sort -u | wc -l >> ./${sample_id}_dedup.log
    printf "unique_output_reads: " >> ./${sample_id}_dedup.log; samtools view ./bowtie.deduplicated.out.bam | cut -f1 | sort -u | wc -l >> ./${sample_id}_dedup.log
    """
  }
} else {
    deduplicated_ch = bowtie_ch
    dedup_log_ch = Channel.empty()
}

process featureCounts {
  label 'mid_cpu'
  tag "featureCounts on $sample_id"
  publishDir "${params.outDir}/$sample_id/featureCounts", mode: 'copy'
  
  input:
  set val(sample_id), file(bams) from deduplicated_ch
  file bowtiegtf from bowtie_gtf_ch
  
  output:
  set val (sample_id), file('*_gene_counts') into counts_ch
  file '*.featureCounts.bam'
  file '*_gene_counts.summary' into counts_summary_ch
  
  script:
  (bam, bai) = bams
  strand = params.reverseStrand ? "-s 2" : "-s 1"
  """
  featureCounts -T $task.cpus --primary -M -t exon -g gene_id $strand -p -Q 1 \
    --fracOverlap $params.fracOverlap --fracOverlapFeature $params.fracOverlap \
    -a $bowtiegtf \
    -o ./${sample_id}_gene_counts \
    -R BAM ${bam}
  """
}
  
process calcRPKMTPM {
  tag "calcRPKMTPM on $sample_id"
  publishDir "${params.outDir}/$sample_id/calcRPKMTPM", mode: 'copy'
    
  input:
  set val(sample_id), file(counts) from counts_ch
    
  output:
  file '*_RNAstats.tsv' into rnastats_ch
  file '*_gene_counts_rpkmtpm.tsv' into rpkm_tpm_ch, countconsolidator_ch
    
  script:
  """
  python3.6 /opt/biorad/src/calculate_rpkm_tpm.py ${counts} ./${sample_id}_gene_counts_rpkmtpm.tsv
  python3.6 /opt/biorad/src/RNAstats.py ${sample_id}_gene_counts_rpkmtpm.tsv $params.count_threshold ${sample_id}_RNAstats.tsv
  """
}

process generateReport {
  publishDir"${params.outDir}/reports", mode: 'copy'
    
  input:
  file('out/pipelineinfo/*') from pipelineinfo_info_ch.collect()
  file(fastqc: 'out/fastqc/*') from fastqc_report_ch.collect()
  file('out/debarcode/*') from debarcode_stats_ch.collect().ifEmpty([])
  file('out/cutadapt/*') from trimmed_log_ch.collect()
  file('out/bowtie/*') from bowtie_log_ch.collect()
  file('out/deduplication/*') from dedup_log_ch.collect().ifEmpty([])
  file('out/counts/*') from counts_summary_ch.collect()
  file('out/counts/*') from rpkm_tpm_ch.collect()
  file('out/counts/*') from rnastats_ch.collect()
  file rrna_reference_key from rrna_reference_key
    
  output:
  file '*.html'
  file '*.tsv' into consolidatedreports_ch
    
  script:
  """
  cp /opt/biorad/src/htmlReport.Rmd .
  python3.6 /opt/biorad/src/generatehtmlReport.py out $rrna_reference_key $params.pipeline
  """ 
}

process consolidateReports {
    publishDir"${params.outDir}", mode: 'copy'
    
    input:
    file('reports/*') from consolidatedreports_ch.collect()
    
    output:
    file '*.html'
    file '*.tsv'
    
    script:
    """
    cp /opt/biorad/src/consolidatedReport.Rmd .
    Rscript -e 'rmarkdown::render(input = "consolidatedReport.Rmd", params = list(inputPath = "./reports"), clean = TRUE, output_format = "html_document", output_file = "consolidatedReport.html")'
    """
}

process consolidateCounts {
  publishDir "${params.outDir}", mode: 'copy', pattern: '*.xlsx'
  publishDir "${params.outDir}/tsv_files", mode: 'copy', pattern: '*.tsv'
  
  input:
  file('counts/*') from countconsolidator_ch.collect()
  
  output:
  file '*.xlsx'
  file '*.tsv'
  
  script:
  """
  Rscript /opt/biorad/src/countConsolidator.R './counts'
  python3.6 /opt/biorad/src/createExcel.py consolidated_counts.tsv consolidated_rpkm.tsv consolidated_tpm.tsv
  """
}

def Header() {
    return """
    Fasta RNA analysis pipeline
    """.stripIndent()
    }
