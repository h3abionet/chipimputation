#!/usr/bin/env nextflow

/*
========================================================================================
=                                 h3achipimputation                                    =
========================================================================================
 h3achipimputation imputation pipeline.
----------------------------------------------------------------------------------------
 @Authors

----------------------------------------------------------------------------------------
 @Homepage / @Documentation
  https://github.com/h3abionet/chipimputation
----------------------------------------------------------------------------------------

----------------------------------------------------------------------------------------

================================================================================
=                           C O N F I G U R A T I O N                          =
================================================================================
*/

// Show help message
if (params.help){
    helpMessage()
    exit 0
}

// Configurable variables
params.name = false
params.email = false
params.plaintext_email = false

output_docs = file("$baseDir/docs/output.md")

// Has the run name been specified by the user?
// this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
    custom_runName = workflow.runName
}

// check if study genotype files exist
target_datasets = []
if(params.target_datasets) {
    params.target_datasets.each { target ->
        if (!file(target.value).exists()) exit 1, "Target VCF file ${target.value} not found. Please check your config file."
        target_datasets << [target.key, file(target.value)]
    }
}

// Validate eagle map file for phasing step and create channel if file exists
eagle_genetic_map = params.eagle_genetic_map
if(params.eagle_genetic_map) {
    if (!file(eagle_genetic_map).exists()) {
        System.err.println "MAP file ${params.eagle_genetic_map} not found. Please check your config file."
        exit 1
    }
}

// Validate reference genome
if(params.reference_genome) {
    if (!file(params.reference_genome).exists()) {
        System.err.println "Reference genome (reference_genome) file ${params.reference_genome} not found. Please check your config file."
        exit 1
    }
}

// Create channel for the study data from VCF files
Channel
        .from(target_datasets)
        .set{ target_datasets }

// TODO Be able to run just on a specified chunk


// Header log info
log.info """
=======================================================

h3achipimputation v${params.version}"

======================================================="""
def summary = [:]
summary['Pipeline Name']    = 'h3achipimputation'
summary['Pipeline Version'] = params.version
summary['Run Name']         = custom_runName ?: workflow.runName
summary['Target datasets']  = params.target_datasets.values().join(', ')
summary['Reference panels']  = params.ref_panels.keySet().join(', ')
summary['Max Memory']       = params.max_memory
summary['Max CPUs']         = params.max_cpus
summary['Max Time']         = params.max_time
summary['Output dir']       = params.outDir
summary['Result dir']       = params.resultDir
summary['Working dir']      = workflow.workDir
summary['Current path']     = "$PWD"
summary['Container Engine'] = workflow.containerEngine
summary['Git info']         = "${workflow.repository} - ${workflow.revision} [${workflow.commitId}]"
summary['Command line']     = workflow.commandLine
if(workflow.containerEngine) {
    summary['Container'] = workflow.container
    summary['Current home'] = "$HOME"
    summary['Current user'] = "$USER"
    summary['Current path'] = "$PWD"
    summary['Working dir'] = workflow.workDir
    summary['Output dir'] = params.outDir
    summary['Script dir'] = workflow.projectDir
    summary['Config Profile'] = workflow.profile
}

if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="


def create_workflow_summary(summary) {

    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'h3achipimputation-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'h3achipimputation Workflow Summary'
    section_href: 'https://github.com/h3abionet/chipimputation'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
    """.stripIndent()

    return yaml_file
}


/*
 * STEP 1: Parse software version numbers
 */
//process get_software_versions {
//    tag "get_software_versions"
//    output:
//        file("software_versions_mqc.yaml") into software_versions_yaml
//    script:
//        """
//        echo $params.version > v_pipeline.txt
//        echo $workflow.nextflow.version > v_nextflow.txt
//        minimac4 --version > v_minimac4.txt
//        eagle --version > v_eagle.txt
//        bcftools --version > v_bcftools.txt
//        ${params.plink} --version > v_${params.plink}.txt
//        scrape_software_versions.py > software_versions_mqc.yaml
//        """
//}


/*
 * STEP 2 - Check user's provided chromosomes vs those in map file
 */
target_datasets.into{ target_datasets; target_datasets_check }
process check_chromosome {
    tag "check_chromosome_${target_name}"
    input:
        set val(target_name), file(target_vcfFile) from target_datasets_check
    output:
        file(chromFile) into check_chromosome
        set val(target_name), file(mapFile) into mapFile_cha
    script:
        base = file(target_vcfFile.baseName).baseName
        chromFile = "${base}_chromosomes.txt"
        mapFile = "${base}.map"
        """
        zcat ${target_vcfFile} | grep -v "^#" | awk -F' ' '{print \$1}' | sort -n | uniq >  ${chromFile}
        zcat ${target_vcfFile} | grep -v "^#" | awk -F' ' '{print \$1"\t"\$2"\t"\$3"\t"\$4"\t"\$5}' | sort -n | uniq > ${mapFile}
        """
}

// Check if specified chromosomes exist in VCF file
check_chromosome.into{ check_chromosome; check_chromosome1 }
chromosomes_ = []
check_chromosome1.toSortedList().val.each{ check_file ->
    chromosomes_ = chromosomes_ + file(check_file).readLines()
}
chromosomes_ = chromosomes_.unique().collect { it as int }.sort()
// Chromosomes from the dataset VCF file
if (params.chromosomes == '' || params.chromosomes == 'ALL'){
    chromosomes = chromosomes_
}
else{
    not_chrs = []
    in_chrs = []
    params.chromosomes.split(',').each { chrm ->
        if (!(chrm.toInteger() in chromosomes_)){
            not_chrs << chrm
        }
        else{
            in_chrs << chrm
        }
    }
    if (in_chrs.isEmpty()){
        if (chromosomes_.isEmpty()) {
            System.err.println "|-- No Chromosome(s) found not in target(s) dataset(s)! The pipeline will exit."
            exit 1
        }
        else{
            chromosomes = chromosomes_
        }
    }

    if (!(not_chrs.isEmpty())){
        System.err.println "|-- Chromosome(s) ${not_chrs.join(', ')} not in target datasets and will be ignored."
        if (in_chrs.isEmpty()){
            chromosomes = chromosomes_
        }
        else {
            chromosomes = in_chrs
        }
    }
    else {
        chromosomes = in_chrs
    }
}
println "|-- Chromosomes used: ${chromosomes.join(', ')}"


// check if ref files exist
params.ref_panels.each { ref ->
    chromosomes.each { chrm ->
        m3vcf = sprintf(params.ref_panels[ref.key].m3vcfFile, chrm)
        vcf = sprintf(params.ref_panels[ref.key].vcfFile, chrm)
        if(!file(m3vcf).exists()) exit 1, "File ${m3vcf} not found. Please check your config file."
        if(!file(vcf).exists()) exit 1, "File ${vcf} not found. Please check your config file."
    }
}

/*
 * STEP 3: QC
*/
target_datasets.into{ target_datasets; target_datasets_qc }
process check_mismatch {
    tag "check_mismatch_${target_name}"
//    publishDir "${params.outDir}", overwrite: true, mode:'symlink'
    input:
        set val(target_name), file(target_vcfFile), file(reference_genome) from target_datasets_qc.combine([file(params.reference_genome)])
    output:
        set val(target_name), file(target_vcfFile), file("${base}_checkRef_warn.log"), file("${base}_checkRef_summary.log") into check_mismatch
    script:
        base = file(target_vcfFile.baseName).baseName
        """
        samtools faidx ${reference_genome} 
        nblines=\$(zcat ${target_vcfFile} | wc -l)
        if (( \$nblines > 1 ))
        then
            bcftools norm --check-ref w \
                -f ${reference_genome} \
                ${target_vcfFile} \
                -Oz -o /dev/null
            cp .command.err ${base}_checkRef_warn.log
            bcftools +fixref \
                ${target_vcfFile} \
                -- \
                -f ${reference_genome} \
                2>&1 | tee "${base}_checkRef_summary.log"
            rm -f ${base}_clean_mind.*
        fi

        """
}

check_mismatch.into{ check_mismatch; check_mismatch_1 }
check_mismatch_noMis = Channel.create()
check_mismatch_1.toSortedList().val.each{ target_name, target_vcfFile, warn, sumary ->
    mismatch = 0
    file(warn).readLines().each{ it ->
        if(it.contains("REF_MISMATCH")){
            mismatch += 1
        }
    }
    if ( mismatch != 0 ) {
        System.err.println "|-- ${mismatch} ref mismatch sites found in '${target_name}' dataset! The pipeline will exit."
        exit 1
    }
    else{
        check_mismatch_noMis << [ target_name, target_vcfFile, warn, sumary ]
    }
}
check_mismatch_noMis.close()

/*
 * STEP 4 - Identify chromosomes and start/stop positions per chromosome and generate chunks
*/
mapFile_cha.into{ mapFile_cha; mapFile_cha_chunks }
process generate_chunks {
    tag "generate_chunks_${target_name}"
//    publishDir "${params.outDir}", overwrite: true, mode:'copy'
    echo true
    input:
        set val(target_name), file(mapFile), chromosomes from mapFile_cha_chunks.combine([chromosomes.join(',')])
    output:
        set val(target_name), file(chunkFile) into generate_chunks
    script:
        if(params.chunk){chunk = params.chunk}else{chunk=''} // To impute only a chunk like 1000000-1100000
        chunkFile = "chunks.txt"
        chunk_size = params.chunk_size
        template "generate_chunks.py"
}


/*
 * STEP 5: QC
*/
check_mismatch_noMis.into{ check_mismatch_noMis; check_mismatch_noMis_1 }
process target_qc {
    tag "target_qc_${target_name}"
    input:
        set val(target_name), file(target_vcfFile), file(mismatch_warn), file(mismatch_summary) from check_mismatch_noMis_1
    output:
        set val(target_name), file("${base}_clean.vcf.gz") into target_qc
    script:
        base = file(target_vcfFile.baseName).baseName
        """
        bcftools view \
            -i 'ALT="."' ${target_vcfFile} | \
        bcftools query \
            -f '%CHROM  %POS  %REF  %ALT\\n' \
            > ${base}_noALT.snp
        bcftools view \
            -e 'ALT="."' ${target_vcfFile} \
            -Oz -o ${base}_noALT.vcf.gz
        ${params.plink} --vcf ${base}_noALT.vcf.gz \
            --set-missing-var-ids --rm-dup force-first \
            --recode vcf \
            --out ${base}_clean_mind
        ${params.plink} --vcf ${base}_clean_mind.vcf \
            --exclude ${base}_clean_mind.dupvar \
            --recode vcf \
            --out ${base}_clean
        bgzip -f ${base}_clean.vcf
        """
}
// TODO include number of each filtered SNPs from QC in the report


"""
Split VCF per chromosomes
"""
target_qc.into{ target_qc; target_qc_1 }
generate_chunks.into{ generate_chunks; generate_chunks_1 }
all_chunks = generate_chunks_1.toSortedList().val
def transform_chunk = { target_name, target_vcfFile ->
    chunks_datas = []
    all_chunks.each{ target_name_, chunk_file ->
        chunks = file(chunk_file).text.split()
        chunks.each{ chunk_data ->
            data = chunk_data.split(',')
            chrm = data[0]
            chunk_start = data[1]
            chunk_end = data[2]
            if (target_name == target_name_) {
                chunks_datas << [chrm, chunk_start, chunk_end, target_name, file(target_vcfFile)]
            }
        }
    }
    return chunks_datas
}
target_qc_chunk = target_qc_1
        .flatMap{ it -> transform_chunk(it) }

/*
 * STEP 6:
*/
process split_target_to_chunk {
    tag "split_${target_name}_${chrm}:${chunk_start}-${chunk_end}"
    input:
        set chrm, chunk_start, chunk_end, target_name, file(target_vcfFile) from target_qc_chunk
    output:
        set chrm, chunk_start, chunk_end, target_name, file(target_vcfFile_chunk) into split_vcf_to_chrm
    script:
        base = file(target_vcfFile.baseName).baseName
        target_vcfFile_chunk = "${base}.chr${chrm}_${chunk_start}-${chunk_end}.vcf.gz"
        """
        bcftools index --tbi -f ${target_vcfFile}
        bcftools view \
            --regions ${chrm}:${chunk_start}-${chunk_end} \
            -m2 -M2 -v snps \
            ${target_vcfFile} \
            -Oz -o ${target_vcfFile_chunk}
        """
}

split_vcf_to_chrm.into{ split_vcf_to_chrm; split_vcf_to_chrm_1 }
def transform_qc_chunk = { chrm, chunk_start, chunk_end, target_name, target_vcfFile ->
    chunks_datas = []
    params.ref_panels.each { ref ->
        ref_m3vcf = sprintf(params.ref_panels[ref.key].m3vcfFile, chrm)
        ref_vcf = sprintf(params.ref_panels[ref.key].vcfFile, chrm)
        chunks_datas << [chrm, chunk_start, chunk_end, target_name, file(target_vcfFile), ref.key, file(ref_vcf), file(ref_m3vcf), file(params.eagle_genetic_map)]
    }
    return chunks_datas
}

target_qc_chunk_ref = split_vcf_to_chrm_1
        .flatMap{ it -> transform_qc_chunk(it) }


/*
 * STEP 7: Phase each chunk using eagle
*/
split_vcf_to_chrm.into{ split_vcf_to_chrm; split_vcf_to_chrm_1 }
process phase_target_chunk {
    tag "phase_${target_name}_${chrm}:${chunk_start}-${chunk_end}_${ref_name}"
    input:
        set chrm, chunk_start, chunk_end, target_name, file(target_vcfFile_chunk), ref_name, file(ref_vcf), file(ref_m3vcf), file(eagle_genetic_map) from target_qc_chunk_ref
    output:
        set chrm, chunk_start, chunk_end, target_name, file("${file_out}.vcf.gz"), ref_name, file(ref_vcf), file(ref_m3vcf) into phase_target
    script:
        file_out = "${file(target_vcfFile_chunk.baseName).baseName}_${ref_name}-phased"
        """
        nblines=\$(zcat ${target_vcfFile_chunk} | grep -v '^#' | wc -l)
        if (( \$nblines > 0 ))
        then
            bcftools index --tbi -f ${ref_vcf}
            bcftools index --tbi -f ${target_vcfFile_chunk}
            eagle \
                --vcfTarget=${target_vcfFile_chunk} \
                --geneticMapFile=${eagle_genetic_map} \
                --vcfRef=${ref_vcf} \
                --vcfOutFormat=z \
                --noImpMissing \
                --chrom=${chrm} \
                --bpStart=${chunk_start} \
                --bpEnd=${chunk_end} \
                --bpFlanking=${params.buffer_size} \
                --outPrefix=${file_out} 2>&1 | tee ${file_out}.log
            if [ ! -f "${file_out}.vcf.gz" ]; then
                touch ${file_out}.vcf && bgzip -f ${file_out}.vcf
            fi
        else
            touch ${file_out}.vcf && bgzip -f ${file_out}.vcf
        fi
        """
}


/*
 * STEP 8:
*/
process impute_target {
    tag "imp_${target_name}_${chrm}:${chunk_start}-${chunk_end}_${ref_name}"
    publishDir "${params.resultDir}/impute/${ref_name}/${target_name}/${chrm}", overwrite: true, mode:'symlink'
    publishDir "${params.resultDir}/impute/${target_name}/${ref_name}/${chrm}", overwrite: true, mode:'symlink'
    input:
        set chrm, chunk_start, chunk_end, target_name, file(target_phased_vcfFile), ref_name, file(ref_vcf), file(ref_m3vcf) from phase_target
    output:
        set chrm, chunk_start, chunk_end, target_name, ref_name, file("${base}_imputed.dose.vcf.gz"), file("${base}_imputed.info") into impute_target
    shell:
        base = "${file(target_phased_vcfFile.baseName).baseName}"
        """
        nblines=\$(zcat ${target_phased_vcfFile} | grep -v '^#' | wc -l)
        if (( \$nblines > 0 ))
        then
            minimac4 \
                --refHaps ${ref_m3vcf} \
                --haps ${target_phased_vcfFile} \
                --format GT \
                --allTypedSites \
                --chr ${chrm} --start ${chunk_start} --end ${chunk_end} --window ${params.buffer_size} \
                --prefix ${base}_imputed
        else
             touch ${base}_imputed.dose.vcf && bgzip -f ${base}_imputed.dose.vcf
             touch ${base}_imputed.info
        fi
        """
}

def helpMessage() {
    log.info"""
    =========================================
    h3achipimputation v${params.version}
    =========================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run h3abionet/chipimputation --reads '*_R{1,2}.fastq.gz' -profile standard,docker

    Mandatory arguments (Must be specified in the configuration file, and must be surrounded with quotes):
      --target_datasets             Path to input study data (Can be one ou multiple for multiple runs)
      --genome                      Human reference genome for checking REF mismatch
      --ref_panels                  Reference panels to impute to (Can be one ou multiple for multiple runs)
      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: standard, conda, docker, singularity, test

    Other options:
      --outDir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      --name                        Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
      --project_name                Project name. If not specified, target file name will be used as project name
    """.stripIndent()
}
