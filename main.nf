#!/usr/bin/env nextflow

"""
bowtie aligner pipeline v. 2.1
2021 Mark Shulewitz Bio-Rad Laboratories, Inc.

Process multiple sets of a file in a single pipeline

"""

def helpMessage() {
  log.info Header()
  log.info """
  Usage:

  The typical command for running the pipeline is as follows:

  nextflow run bowtie-aligner/main.nf -profile docker \
  --inDir {path to folder with fastq files}
  --outDir {path to output directory}
  --pipeline {miRNA/rRNA/tRNA}

  Required:
    --inDir                 Path to directory containing input fastq file(s)
    --outDir                Path ot output directory
    --pipeline              Type of RNA to analyze (miRNA or rRNA; default is rRNA)

  Optional arguments:
    --genome                for rRNA pipeline (hg38, mm10, rnor6; default is hg38)
    --skipUmi               Skip UMI deduplication
    --spikeType             miltenyi, ercc or ercc/miltenyi (default in NONE)
    --libraryType           SEQuoia or Takara; default = SEQuoia
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
acceptablePipelines = ["miRNA","rRNA","tRNA"]
allowedSpikes = ["NONE","ercc", "miltenyi", "ercc/miltenyi"]
allowedGenomes = ["hg38","mm10","rnor6"]
allowedLibraries  = ["SEQuoia","Takara"]

if (params.pipeline == 'NONE') {
  log.error "Analytical pipeline not selected from $acceptablePipelines"
  exit 1
} else {
  if ( !acceptablePipelines.contains(params.pipeline) ) {
  log.error "$params.pipeline not an accepted pipeline ($acceptablePipelines)"
  exit 1
  }
}

if (params.spikeType != "NONE" && !allowedSpikes.contains(params.spikeType)) {
    log.error "$params.spikeType is not in $allowedSpikes"
    exit 1
}

if (params.genome != "NONE" && !allowedGenomes.contains(params.genome)) {
    log.error "$params.genome is not in $allowedGenomes"
    exit 1
}

//Set fasta file based on pipeline and spikeType

if (params.genome == "hg38"){
  rrna_reference_key = file("${baseDir}/ref/hs_rRNA_reference_key.tsv", checkIfExists: true)
}
else if (params.genome == "mm10"){
  rrna_reference_key = file("${baseDir}/ref/mm_rRNA_reference_key.tsv", checkIfExists: true)
}
else if (params.genome == "rnor6"){
  rrna_reference_key = file("${baseDir}/ref/rnor_rRNA_reference_key.tsv", checkIfExists: true)
}

if (params.pipeline == "rRNA"){
  if (params.genome == "hg38"){
    bowtie_fa = file("${baseDir}/ref/hs_rRNA_references.fa", checkIfExists: true)
  }
  else if (params.genome == "mm10"){
    bowtie_fa = file("${baseDir}/ref/mm_rRNA_references.fa", checkIfExists: true)
  }
  else if (params.genome == "rnor6"){
    bowtie_fa = file("${baseDir}/ref/rnor_rRNA_references.fa", checkIfExists: true)
  }
}  

else if (params.pipeline == "miRNA") {
  if (params.spikeType == "NONE"){
    bowtie_fa = file("${baseDir}/ref/natural.fa", checkIfExists: true)
  }
  else if (params.spikeType == "miltenyi"){
    bowtie_fa = file("${baseDir}/ref/miltenyi_and_natural.fa", checkIfExists: true)
  }
  else if (params.spikeType == "ercc"){
    bowtie_fa = file("${baseDir}/ref/ercc_and_natural.fa", checkIfExists: true)
  }
  else {
    bowtie_fa = file("${baseDir}/ref/ercc_miltenyi_and_natural.fa", checkIfExists: true)
  }   
}

else if (params.pipeline == "tRNA"){
  bowtie_fa = file("${baseDir}/ref/human_tRNAs.fa", checkIfExists: true)
}


// if using pre-generated bowtie indices use:
// bowtie_index = file("${params.bowtie_index}.fa")
// bowtie_indices = Channel.fromPath( "${params.bowtie_index}*.ebwt" ).toList()
// if( !bowtie_index.exists() ) exit 1, "Reference genome Bowtie not found: ${params.bowtie_index}"

// bowtie_gtf = file("${params.bowtie_index}.gtf")

def summary = [:]
summary['Run Name'] = workflow.runName
if(workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Config Profile'] = workflow.profile
summary['Input directory'] = params.inDir
summary['Pipeline'] = params.pipeline
if (params.pipeline == "miRNA") summary['Spike'] = params.spikeType
else summary['Spike'] = 'N/A'
summary['Genome'] = params.genome
summary['fasta file'] = bowtie_fa
summary['Output dir'] = params.outDir
summary['skipUmi'] = params.skipUmi
summary['Cores'] = params.max_cpus
summary['Memory'] = params.max_memory
if(params.reverseStrand) summary['Strand'] = 'Reverse'
else summary['Strand'] = 'Forward'
if(params.noTrim) summary['Trimming'] = 'None'
summary['Fraction Overlap'] = params.fracOverlap
summary['Library Type']     = params.libraryType
log.info Header()
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "----------------------------------------------------"

// Add input files to channels

if (!params.skipUmi) {
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
  printf "Genome\t${params.genome}\n" >> ${sample_id}_pipeline_info.tsv
  printf "Spike type\t${params.spikeType}\n" >> ${sample_id}_pipeline_info.tsv
  printf "Library type\t${params.libraryType}\n" >> ${sample_id}_pipeline_info.tsv
  printf "Reference fasta file\t${bowtie_fa}\n" >> ${sample_id}_pipeline_info.tsv
  printf "Skip UMI\t${params.skipUmi}\n" >> ${sample_id}_pipeline_info.tsv
  printf "Fraction overlap\t${params.fracOverlap}\n" >> ${sample_id}_pipeline_info.tsv
  printf "Reverse Strand\t${params.reverseStrand}\n" >> ${sample_id}_pipeline_info.tsv
  printf "No trimming\t${params.noTrim}\n" >> ${sample_id}_pipeline_info.tsv
  """   
}

process fastQc {
  tag "fastqc on $sample_id"
  label 'mid_memory'
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
        label 'high_memory'
        publishDir "${params.outDir}/$sample_id/debarcode", mode: 'copy'
        
        input:
        set val(sample_id), file(reads) from input_reads_ch
        
        output:
        set val(sample_id), file('*_R1.fastq.gz') into debarcoded_ch
        file '*_debarcode_stats.txt' into debarcode_stats_ch
        
        script:
        sample_id = sample_id.split( '_' )[ 0 ]
        """
        bash /opt/biorad/src/fastqmerge.sh $reads | \
        python3.6 /opt/biorad/src/debarcode.py \
        | tee >(grep 'Reads' > ${sample_id}_debarcode_stats.txt) | grep -ve 'Reads' | gzip > ${sample_id}_debarcoded_R1.fastq.gz
        """
    }
} else {
    debarcoded_ch = input_reads_ch
    debarcode_stats_ch = Channel.empty()
}

process cutadapt {
  tag "cutadapt on $sample_id"
  label 'mid_cpu'
  publishDir "${params.outDir}/$sample_id/cutadapt", mode: 'copy', pattern: "*.{fastq.gz,log}"
  
  input:
  set val(sample_id), file(reads) from debarcoded_ch
  
  output:
  set val(sample_id), file('*_trimmed.fastq') into trimmed_file_ch
  file '*_cutadapt.log' into trimmed_log_ch
  
  script:
  if (params.skipUmi){
      sample_id = sample_id.split( '_' )[ 0 ]
  }
  if (!params.noTrim) {
  if (params.libraryType != 'Takara'){
  """
  cutadapt -u 1 -a 'AAAAAAAAAA' -m 15 -j $task.cpus -q 0,0 -o ${sample_id}_trimmed.fastq ${reads} 1> ${sample_id}_cutadapt.log
  gzip --keep *_trimmed.fastq
  """
  } else { // trimming for Takara kit
  """
  cutadapt -u 3 -a 'AAAAAAAAAA' -m 15 -j $task.cpus -q 0,0 -o ${sample_id}_trimmed.fastq ${reads} 1> ${sample_id}_cutadapt.log
  gzip --keep *_trimmed.fastq
  """
  }
  } else {
  """
  cutadapt -m 15 -j $task.cpus -q 0,0 -o ${sample_id}_trimmed.fastq ${reads} 1> ${sample_id}_cutadapt.log
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
  """
  bowtie bowtie_index ${trimmed} --threads $task.cpus --sam bowtie.sam 2> ${sample_id}_bowtie.log
  samtools view -Sb bowtie.sam > temp.bam
  samtools sort temp.bam -o ${sample_id}_bowtie.bam
  sambamba index -t $task.cpus ${sample_id}_bowtie.bam
  """
}

if (!params.skipUmi) {
  process umiTagging {
    tag "umi tagging on $sample_id"
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
    tag "deduplication on $sample_id"
    label 'mid_cpu'
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
    umi_tools dedup -I $bam \
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
  tag "featureCounts on $sample_id"
  label 'mid_cpu'
  publishDir "${params.outDir}/$sample_id/featureCounts", mode: 'copy'
  
  input:
  set val(sample_id), file(bams) from deduplicated_ch
  file bowtiegtf from bowtie_gtf_ch
  
  output:
  set val (sample_id), file('*_gene_counts') into counts_ch
  file '*.featureCounts.bam*'
  file '*_gene_counts.summary' into counts_summary_ch
  
  script:
  (bam, bai) = bams
  strand = params.reverseStrand ? "-s 2" : "-s 1"
  """
  featureCounts -T $task.cpus --primary -M -t exon -g gene_id $strand -Q 1 \
    --fracOverlap $params.fracOverlap --fracOverlapFeature $params.fracOverlap \
    -a $bowtiegtf \
    -o ./${sample_id}_gene_counts \
    -R BAM ${bam}
  sambamba index -t $task.cpus ${bam}.featureCounts.bam
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
