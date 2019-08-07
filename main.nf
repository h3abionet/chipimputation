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
        if (!file(target.value).exists() && !file(target.value).isFile()) exit 1, "|-- ERROR: Target VCF file ${target.value} not found. Please check your config file."
        target_datasets << [target.key, file(target.value)]
    }
}
else{
    params.target_datasets.each { target ->
        System.err.println "|-- ERROR: Target VCF file ${target.value} not found. Please check your config file."
        exit 1
    }
}

// Validate eagle map file for phasing step and create channel if file exists
if(params.eagle_genetic_map) {
    if (!file(params.eagle_genetic_map).exists() && !file(params.eagle_genetic_map).isFile()) {
        System.err.println "|-- ERROR: MAP file ${params.eagle_genetic_map} not found. Please check your config file."
        exit 1
    }
}
else{
    System.err.println "|-- ERROR: MAP file ${params.eagle_genetic_map} not found. Please check your config file."
    exit 1
}

// Validate reference genome
if(params.reference_genome) {
    if ((!file(params.reference_genome).exists() && !file(params.reference_genome).isFile())) {
        System.err.println "|-- ERROR: Reference genome file ${params.reference_genome} not found. Please check your config file."
        exit 1
    }
}
else{
    System.err.println "|-- ERROR: Reference genome file ${params.reference_genome} not found. Please check your config file."
    exit 1
}

// Create channel for the study data from VCF files
Channel
        .from(target_datasets)
        .set{ target_datasets }


// Header log info
log.info """
=======================================================
h3achipimputation v${params.version}"
======================================================="""
def summary = [:]
summary['Pipeline Name']    = 'h3achipimputation'
summary['Pipeline version'] = params.version
summary['Run Name']         = custom_runName ?: workflow.runName
summary['Target datasets']  = params.target_datasets.values().join(', ')
summary['Reference panels']  = params.ref_panels.keySet().join(', ')
summary['Max Memory']       = params.max_memory
summary['Max CPUs']         = params.max_cpus
summary['Max Time']         = params.max_time
summary['Output dir']       = params.outDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['Current path']     = "$PWD"
summary['Git info']         = "${workflow.repository} - ${workflow.revision} [${workflow.commitId}]"
summary['Command line']     = workflow.commandLine
if(workflow.containerEngine) {
    summary['Container Engine'] = workflow.containerEngine
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
process check_chromosome {
    tag "check_chromosome_${target_name}"
    input:
    set target_name, file(target_vcfFile) from target_datasets
    output:
    set target_name, file(chromFile) into check_chromosome
    set target_name, file(target_vcfFile), file(mapFile) into mapFile_cha,mapFile_cha_1
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
chromosomes_ = [:]
chromosomes_['ALL'] = []
valid_chrms = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]
not_chrs = []
in_chrs = []
notValid_chrs = []
check_chromosome.toSortedList().val.each{ target_name, check_file ->
    chromosomes_[target_name] = file(check_file).readLines().unique().collect { it as int }.sort()
    chromosomes_[target_name].each { chrm ->
        if(!(chrm in chromosomes_['ALL'])) {
            if (chrm.toInteger() in valid_chrms){
                chromosomes_['ALL'] << chrm.toInteger()
            }
            else{
                notValid_chrs << chrm.toInteger()
            }
        }
    }
}
if (params.chromosomes == '' || params.chromosomes == 'ALL'){
    chromosomes = chromosomes_['ALL']
}
else{
    params.chromosomes.split(',').each { chrm ->
        chrm = chrm.toInteger()
        if (!(chrm in chromosomes_['ALL'])){
            not_chrs << chrm
        }
        else{
            in_chrs << chrm
        }
    }
    if (in_chrs.isEmpty()){
        System.err.println "|-- ERROR- No Chromosome(s) found not in target(s) dataset(s)! The pipeline will exit."
        exit 1
    }

    if (!(not_chrs.isEmpty())){
        System.err.println "|-- WARN- Chromosome(s) ${not_chrs.join(', ')} not in target datasets and will be ignored."
    }
    chromosomes = in_chrs

}
// Ignore invalid chromosome in VCF
if (!(notValid_chrs.isEmpty())){
    System.err.println "|-- ERROR- Chromosome(s) ${notValid_chrs.join(', ')} not valid chromosomes. Check your VCF file and remove invalid chromosomes! The pipeline will exit."
    exit 1
}
ignore_chrms = [:]
toImpute_chrms = [:]
mapFile_cha_1.toSortedList().val.each { target_name, target_vcfFile, mapFile ->
    chromosomes_[target_name].each{ chrm ->
        chrm = chrm.toInteger()
        if(!(chrm in chromosomes)){
            if(!(target_name in ignore_chrms)){
                ignore_chrms[target_name] = []
            }
            ignore_chrms[target_name] << chrm
        }
        else{
            if(!(target_name in toImpute_chrms)){
                toImpute_chrms[target_name] = []
            }
            toImpute_chrms[target_name] << chrm
        }
    }
}

targets_toImpute_list = []
mapFile_cha.toSortedList().val.each { target_name, target_vcfFile, mapFile ->
    if(target_name in toImpute_chrms){
        targets_toImpute_list << [ target_name, target_vcfFile, mapFile, file(params.reference_genome) ]
    }
    else{
        System.err.println "|-- WARN- Dataset ${target_name} does not contain the specified chromosome(s) ${chromosomes.join(', ')} and will be ignored."
    }
}
targets_toImpute = Channel.from(targets_toImpute_list)


//def mapFile_cha_ = { target_name, target_vcfFile, mapFile ->
//    targets_toImpute_list = []
//    if(target_name in toImpute_chrms){
//        targets_toImpute_list << [ target_name, target_vcfFile, mapFile, file(params.reference_genome) ]
//    }
//    else{
//        System.err.println "|-- WARN- Dataset ${target_name} does not contain the specified chromosome(s) ${chromosomes.join(', ')} and will be ignored."
//    }
//    return targets_toImpute_list
//}
//targets_toImpute = mapFile_cha.collect{ it -> mapFile_cha_(it) }.view()

println "|-- Chromosomes used: ${chromosomes.join(', ')}"
if(params.chunk){
    println "|-- Chunks to impute: ${(params.chunk.split(',')).join(', ')}"
}

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
process check_mismatch {
    tag "check_mismatch_${target_name}"
    label "medium"
    publishDir "${params.outDir}/reports/${target_name}", overwrite: true, mode:'copy', pattern: "*checkRef_*.log*"
    input:
    set target_name, file(target_vcfFile), file(mapFile), file(reference_genome) from targets_toImpute
    output:
    set target_name, file(target_vcfFile), file(mapFile), file("${base}_checkRef_warn.log"), file("${base}_checkRef_summary.log") into check_mismatch
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

check_mismatch_noMis = Channel.create()
check_mismatch.toSortedList().val.each{ target_name, target_vcfFile, mapFile, warn, sumary ->
    mismatch = 0
    // TODO use summary instead, print mismatch, non-biallelic, non-ACGT
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
        check_mismatch_noMis << [ target_name, target_vcfFile, mapFile, warn, sumary, toImpute_chrms[target_name]]
    }
}
check_mismatch_noMis.close()
check_mismatch_noMis.into{ check_mismatch_noMis; check_mismatch_noMis_1 }

/*
 * STEP 4 - Identify chromosomes and start/stop positions per chromosome and generate chunks
*/
process generate_chunks {
    tag "generate_chunks_${target_name}_${chrms[0]}_${chrms[-1]}"
    publishDir "${params.outDir}/reports/${target_name}", overwrite: true, mode:'copy'
    label "small"
    input:
    set target_name, file(target_vcfFile), file(mapFile), file(mismatch_warn), file(mismatch_summary), chrms from check_mismatch_noMis
    output:
    set target_name, file(chunkFile) into generate_chunks
    script:
    if(params.chunk){chunk = params.chunk} else{chunk=''}
    chromosomes = chrms.join(',')
    chunkFile = "chunks.txt"
    chunk_size = params.chunk_size
    template "generate_chunks.py"
}


/*
 * STEP 5: QC
*/
process target_qc {
    tag "target_qc_${target_name}_${chrms[0]}_${chrms[-1]}"
    label "medium"
    publishDir "${params.outDir}/qc/${target_name}", overwrite: true, mode:'copy', pattern: "*clean.vcf.gz*"
    input:
    set target_name, file(target_vcfFile), file(mapFile), file(mismatch_warn), file(mismatch_summary), chrms from check_mismatch_noMis_1
    output:
    set target_name, file("${base}_clean.vcf.gz") into target_qc
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
        bcftools norm \
            --rm-dup both \
            ${base}_noALT.vcf.gz \
            -Oz -o ${base}_clean.vcf.gz
        """
}


"""
Split VCF per chromosomes
"""
generate_chunks.into{ generate_chunks; generate_chunks_1 }
all_chunks = generate_chunks_1.toSortedList().val
all_chunks.each{ target_name_, chunk_file ->
    chunks = file(chunk_file).text.split()
    if(chunks.size() == 0){
        System.err.println "|-- ERROR- No valid chunks (${(params.chunk.split(',')).join(', ')}) in not specified chromosomes (${chromosomes.join(', ')}). Check your VCF file and correct your chunks for specified chromosomes! The pipeline will exit."
        exit 1
    }
}

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
target_qc_chunk = target_qc
        .flatMap{ it -> transform_chunk(it) }


/*
 * STEP 6:
*/
process split_target_to_chunk {
    tag "split_${target_name}_${chrm}:${chunk_start}-${chunk_end}"
    label "medium"
    input:
    set chrm, chunk_start, chunk_end, target_name, file(target_vcfFile) from target_qc_chunk
    output:
    set chrm, chunk_start, chunk_end, target_name, file(target_vcfFile_chunk) into split_vcf_to_chrm
    script:
    base = file(target_vcfFile.baseName).baseName
    target_vcfFile_chunk = "${base}.chr${chrm}_${chunk_start}-${chunk_end}.vcf.gz"
    start = chunk_start - params.buffer_size
    if(chunk_start.toInteger() - params.buffer_size.toInteger() <= 0){ end = 1 }
    end = chunk_end.toInteger() + params.buffer_size.toInteger()
    """
        bcftools index --tbi -f ${target_vcfFile}
        bcftools view \
            --regions ${chrm}:${start}-${end} \
            -m2 -M2 -v snps \
            ${target_vcfFile} \
            -Oz -o ${target_vcfFile_chunk}
        """
}

def transform_qc_chunk = { chrm, chunk_start, chunk_end, target_name, target_vcfFile ->
    chunks_datas = []
    params.ref_panels.each { ref ->
        ref_m3vcf = sprintf(params.ref_panels[ref.key].m3vcfFile, chrm)
        ref_vcf = sprintf(params.ref_panels[ref.key].vcfFile, chrm)
        chunks_datas << [chrm, chunk_start, chunk_end, target_name, file(target_vcfFile), ref.key, file(ref_vcf), file(ref_m3vcf), file(params.eagle_genetic_map)]
    }
    return chunks_datas
}

target_qc_chunk_ref = split_vcf_to_chrm
        .flatMap{ it -> transform_qc_chunk(it) }


/*
 * STEP 7: Phase each chunk using eagle
*/
process phase_target_chunk {
    tag "phase_${target_name}_${chrm}:${chunk_start}-${chunk_end}_${ref_name}"
    label "bigmem"
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
    label "bigmem"
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
                --format GT,DS \
                --allTypedSites \
                --minRatio ${params.minRatio} \
                --chr ${chrm} --start ${chunk_start} --end ${chunk_end} --window ${params.buffer_size} \
                --prefix ${base}_imputed
        else
             touch ${base}_imputed.dose.vcf && bgzip -f ${base}_imputed.dose.vcf
             touch ${base}_imputed.info
        fi
        """
}


'''
Combine output
'''

// Create a dataflow instance of all impute results
imputeCombine = [:]
infoCombine = [:]
infoCombine_all = [:]
impute_target_list = impute_target.toSortedList().val
impute_target_list.each{ chrm, chunk_start, chunk_end, target_name, ref_name, impute, info ->
    ref_vcf = file(sprintf(params.ref_panels[ref_name].vcfFile, chrm))
    id = target_name +"__"+ ref_name +"__"+ chrm
    if(!(id in imputeCombine)){
        imputeCombine[id] = [target_name, ref_name, ref_vcf, chrm, []]
    }
    imputeCombine[id][4] << impute
    if(!(id in infoCombine)){
        infoCombine[id] = [target_name, ref_name, ref_vcf, chrm, []]
    }
    infoCombine[id][4] << info
    id1 = target_name +"__"+ ref_name
    if(!(id1 in infoCombine_all)){
        infoCombine_all[id1] = [target_name, ref_name, ref_vcf, []]
    }
    infoCombine_all[id1][3] << info
}


"""
Combine impute chunks to chromosomes
"""
process combineImpute {
    //maxForks 1 // TODO: this is only because bcftools sort is using a common TMPFOLDER
    tag "impComb_${target_name}_${ref_name}_${chrm}"
    publishDir "${params.outDir}/imputed/${ref_name}", overwrite: true, mode:'symlink', pattern: '*imputed.gz'
    label "bigmem"
    input:
    set target_name, ref_name, file(ref_vcf), chrm, file(imputed_files) from imputeCombine.values()
    output:
    set target_name, ref_name, file(ref_vcf), chrm, file(comb_impute) into combineImpute
    script:
    comb_impute = "${target_name}_${ref_name}_chr${chrm}.imputed.gz"
    """
        bcftools concat \
            ${imputed_files} \
            -Oz -o ${target_name}.tmp.vcf.gz
        ## Recalculate AC, AN, AF
        bcftools +fill-tags ${target_name}.tmp.vcf.gz -Oz -o ${target_name}.tmp1.vcf.gz -- -t AC,AN,AF,MAF
        bcftools sort ${target_name}.tmp1.vcf.gz -T . -Oz -o ${comb_impute}
        rm ${target_name}.tmp*.vcf.gz
        """
}


"""
Combine impute info chunks to chromosomes
"""
process combineInfo {
    tag "infoComb_${target_name}_${ref_name}_${chrm}"
    publishDir "${params.outDir}/imputed/${ref_name}", overwrite: true, mode:'copy', pattern: '*imputed_info'
    label "medium"
    input:
    set target_name, ref_name, file(ref_vcf), chrm, file(info_files) from infoCombine.values()
    output:
    set target_name, ref_name, file(ref_vcf), chrm, file(comb_info) into combineInfo_frq
    script:
    comb_info = "${target_name}_${ref_name}_chr${chrm}.imputed_info"
    """
        head -n1 ${info_files[0]} > ${comb_info}
        tail -q -n +2 ${info_files.join(' ')} >> ${comb_info}
        """
}


"""
Combine all impute info chunks by dataset
"""
process combineInfo_all {
    tag "infoComb_${target_name}_${ref_name}_${chrms}"
    publishDir "${params.outDir}/imputed/${ref_name}", overwrite: true, mode:'copy', pattern: '*imputed_info'
    label "medium"
    input:
    set target_name, ref_name, file(ref_vcf), file(info_files) from infoCombine_all.values()
    output:
    set target_name, ref_name, file(ref_vcf), file(comb_info) into combineInfo_all,combineInfo_all_frq
    script:
    chrms = chromosomes_[target_name][0]+"-"+chromosomes_[target_name][-1]
    comb_info = "${target_name}_${ref_name}_chrs${chrms}.imputed_info"
    """
        head -n1 ${info_files[0]} > ${comb_info}
        tail -q -n +2 ${info_files.join(' ')} >> ${comb_info}
        """
}


"""
Generating report
"""
combineInfo_all_list = combineInfo_all.toSortedList().val
target_infos = [:] // Grouping by target
ref_infos = [:] // Grouping by ref
ref_panels = params.ref_panels.keySet().join('_')
target_names = params.target_datasets.keySet().join('_')
combineInfo_all_list.each{ target_name, ref_name, ref_vcf, comb_info ->
    if(!(target_name in target_infos)){
        target_infos[target_name] = [ target_name, ref_panels, []]
    }
    target_infos[target_name][2] << ref_name+"=="+comb_info
    if(!(ref_name in ref_infos)){
        ref_infos[ref_name] = [ ref_name, target_name, []]
    }
    ref_infos[ref_name][2] << target_name+"=="+comb_info
}


"""
Filtering all reference panels by maf for a dataset
"""
//TODO generate filtered info by reference panels.
process filter_info_target {
    tag "filter_${target_name}_${ref_panels}_${chrms}"
    publishDir "${params.outDir}/reports/${ref_panels}", overwrite: true, mode:'copy', pattern: "${comb_info}*"
    label "medium"
    input:
    set target_name, ref_panels, ref_infos from target_infos.values()
    output:
    set target_name, ref_panels, file("${well_out}.tsv") into target_info_Well
    set target_name, ref_panels, file("${acc_out}.tsv") into target_info_Acc,target_info_Acc_1
    script:
    chrms = chromosomes_[target_name][0]+"-"+chromosomes_[target_name][-1]
    comb_info = "${target_name}_${ref_panels}_${chrms}.imputed_info"
    well_out = "${comb_info}_well_imputed"
    acc_out = "${comb_info}_accuracy"
    infos = ref_infos.join(',')
    impute_info_cutoff = params.impute_info_cutoff
    template "filter_info_minimac.py"
}


"""
Report 1: Well imputed all reference panels by maf for a dataset
"""
//TODO do this by chromosomes for each dataset
target_info_Well.into{ target_info_Well; target_info_Well_1}
process report_well_imputed_target {
    tag "report_wellImputed_${target_name}_${ref_panels}_${chrms}"
    publishDir "${params.outDir}/reports/${ref_panels}", overwrite: true, mode:'copy'
    label "medium"
    input:
    set target_name, ref_panels, file(inWell_imputed) from target_info_Well_1
    output:
    set target_name, ref_panels, file("${outWell_imputed}.tsv"), file("${outWell_imputed}_summary.tsv") into report_well_imputed_target
    script:
    chrms = chromosomes_[target_name][0]+"-"+chromosomes_[target_name][-1]
    outWell_imputed = "${target_name}_${ref_panels}_${chrms}.imputed_info_performance_by_maf_report"
    group = "REF_PANEL"
    template "report_well_imputed.py"
}


"""
Plot performance all reference panels by maf for a dataset
"""
process plot_performance_target{
    tag "plot_performance_dataset_${target_name}_${ref_panels}_${chrms}"
    publishDir "${params.outDir}/plots/${ref_panels}", overwrite: true, mode:'copy'
    input:
    set target_name, ref_panels, file(well_imputed_report), file(well_imputed_report_summary) from report_well_imputed_target
    output:
    set target_name, ref_panels, file(plot_by_maf) into plot_performance_target
    script:
    plot_by_maf = "${well_imputed_report.baseName}.tiff"
    chrms = chromosomes_[target_name][0]+"-"+chromosomes_[target_name][-1]
    report = well_imputed_report
    group = "REF_PANEL"
    xlab = "MAF bins"
    ylab = "Number of well imputed SNPs"
    template "plot_results_by_maf.R"
}

"""
Repor 2: Accuracy all reference panels by maf for a dataset
"""
process report_accuracy_target {
    tag "report_acc_${target_name}_${ref_panels}_${chrms}"
    publishDir "${params.outDir}/reports/${ref_panels}/", overwrite: true, mode:'copy'
    label "medium"
    input:
    set target_name, ref_panels, file(inSNP_acc) from target_info_Acc_1
    output:
    set target_name, ref_panels, file(outSNP_acc) into report_SNP_acc_target
    script:
    chrms = chromosomes_[target_name][0]+"-"+chromosomes_[target_name][-1]
    outSNP_acc = "${target_name}_${ref_panels}_${chrms}.imputed_info_report_accuracy.tsv"
    group = "REF_PANEL"
    template "report_accuracy_by_maf.py"
}


"""
Plot accuracy all reference panels by maf for a dataset
"""
process plot_accuracy_target{
    tag "plot_accuracy_dataset_${target_name}_${ref_panels}_${chrms}"
    publishDir "${params.outDir}/plots/${ref_panels}", overwrite: true, mode:'copy'
    input:
    set target_name, ref_panels, file(accuracy_report) from report_SNP_acc_target
    output:
    set target_name, ref_panels, file(plot_by_maf) into plot_accuracy_target
    script:
    plot_by_maf = "${accuracy_report.baseName}_accuracy_by_maf.tiff"
    chrms = chromosomes_[target_name][0]+"-"+chromosomes_[target_name][-1]
    report = accuracy_report
    group = "REF_PANEL"
    xlab = "MAF bins"
    ylab = "Concordance rate"
    template "plot_results_by_maf.R"
}



"""
Filtering all targets by maf for a reference panel
"""
process filter_info_ref {
    tag "filter_${ref_name}_${target_names}_${chrms}"
    label "bigmem"
    publishDir "${params.outDir}/reports/${ref_name}", overwrite: true, mode:'copy', pattern: "${comb_info}*"
    input:
    set ref_name, target_names, target_infos from ref_infos.values()
    output:
    set ref_name, target_names, file("${well_out}.tsv") into ref_info_Well
    set ref_name, target_names, file("${acc_out}.tsv") into ref_info_Acc
    script:
    chrms = chromosomes[0]+"-"+chromosomes[-1]
    comb_info = "${ref_name}_${target_names}_${chrms}.imputed_info"
    well_out = "${comb_info}_well_imputed"
    acc_out = "${comb_info}_accuracy"
    infos = file(target_infos.join(','))
    impute_info_cutoff = params.impute_info_cutoff
    template "filter_info_minimac.py"
}


"""
Report: Well imputed all targets by maf for a reference panel
"""
process report_well_imputed_ref {
    tag "report_wellImputed_${ref_name}_${target_names}_${chrms}"
    publishDir "${params.outDir}/reports/${ref_name}", overwrite: true, mode:'copy'
    label "medium"
    input:
    set ref_name, target_names, file(inWell_imputed) from ref_info_Well
    output:
    set ref_name, target_names, file("${outWell_imputed}.tsv"), file("${outWell_imputed}_summary.tsv") into report_well_imputed_ref
    script:
    chrms = chromosomes[0]+"-"+chromosomes[-1]
    outWell_imputed = "${ref_name}_${target_names}_${chrms}.imputed_info_report_well_imputed"
    group = "DATASET"
    template "report_well_imputed.py"
}


"""
Plot performance all targets by maf for a reference panel
"""
process plot_performance_ref{
    tag "plot_performance_dataset_${ref_name}_${target_names}_${chrms}"
    publishDir "${params.outDir}/plots/${ref_name}", overwrite: true, mode:'copy'
    input:
    set ref_name, target_names, file(well_imputed_report), file(well_imputed_report_summary) from report_well_imputed_ref
    output:
    set ref_name, target_names, file(plot_by_maf) into plot_performance_ref
    script:
    plot_by_maf = "${well_imputed_report.baseName}_performance_by_maf.tiff"
    chrms = chromosomes[0]+"-"+chromosomes[-1]
    report = well_imputed_report
    group = "DATASET"
    xlab = "MAF bins"
    ylab = "Number of well imputed SNPs"
    template "plot_results_by_maf.R"
}


"""
Repor 2: Accuracy all targets by maf for a reference panel
"""
process report_accuracy_ref {
    tag "report_acc_${ref_name}_${target_names}_${chrms}"
    publishDir "${params.outDir}/reports/${ref_name}/", overwrite: true, mode:'copy'
    label "medium"
    input:
    set ref_name, target_names, file(inSNP_acc) from ref_info_Acc
    output:
    set ref_name, target_names, file(outSNP_acc) into report_SNP_acc_ref
    script:
    chrms = chromosomes[0]+"-"+chromosomes[-1]
    outSNP_acc = "${ref_name}_${target_names}_${chrms}.imputed_info_report_accuracy.tsv"
    group = "DATASET"
    template "report_accuracy_by_maf.py"
}

"""
Plot accuracy all reference panels by maf for a dataset
"""
process plot_accuracy_ref{
    tag "plot_accuracy_dataset_${ref_name}_${target_names}_${chrms}"
    publishDir "${params.outDir}/plots/${ref_name}", overwrite: true, mode:'copy'
    input:
    set ref_name, target_names, file(accuracy_report) from report_SNP_acc_ref
    output:
    set ref_name, target_names, file(plot_by_maf) into plot_accuracy_ref
    script:
    plot_by_maf = "${accuracy_report.baseName}_by_maf.tiff"
    chrms = chromosomes[0]+"-"+chromosomes[-1]
    report = accuracy_report
    group = "REF_PANEL"
    xlab = "MAF bins"
    ylab = "Concordance rate"
    template "plot_results_by_maf.R"
}



"""
Step: generate allele frequency
"""
process generate_frequency {
    tag "frq_${target_name}_${ref_name}_${chrm}"
    publishDir "${params.outDir}/frqs/${ref_name}", overwrite: true, mode:'copy', pattern: '*frq'
    label "medium"
    input:
    set target_name, ref_name, file(ref_vcf), chrm, file(impute_vcf) from combineImpute
    output:
    set target_name, ref_name, file(ref_vcf), chrm, file(dataset_frq), file(ref_frq) into frq_dataset,frq_dataset_info
    script:
    ref_frq = "${file(ref_vcf.baseName).baseName}.frq"
    dataset_frq = "${file(impute_vcf.baseName).baseName}.frq"
    """
        # For datastet
        echo -e 'CHR\tPOS\tSNP\tREF\tALT\tAF' > ${dataset_frq}
        bcftools query -f '%CHROM\t%POS\t%CHROM\\_%POS\\_%REF\\_%ALT\t%REF\t%ALT\t%INFO/AF\\n' ${impute_vcf} >> ${dataset_frq}
        # For the reference panel
        echo -e 'CHR\tPOS\tSNP\tREF\tALT\tAF' > ${ref_frq}
        bcftools +fill-tags ${ref_vcf} -Oz -o ${ref_name}_AF.vcf.gz -- -t AF
        bcftools query -f '%CHROM\t%POS\t%CHROM\\_%POS\\_%REF\\_%ALT\t%REF\t%ALT\t%INFO/AF\\n' ${ref_name}_AF.vcf.gz >> ${ref_frq}
        """
}


"""
Plot number of imputed SNPs over the mean r2 for all reference panels
"""
combineInfo_frq_ = combineInfo_frq.combine(frq_dataset_info, by:[0,1,3]).map{it -> [it[0], it[1], it[2], it[4], it[6], it[7]]}
combineInfo_frq_.into{ combineInfo_frq; combineInfo_frq_comp }
process plot_r2_SNPpos {
    tag "plot_r2_SNPpos_${target_name}_${ref_name}_${chrm}"
    publishDir "${params.outDir}/plots/${ref_name}/r2_SNPpos", overwrite: true, mode:'copy'
    label "medium"
    input:
    set target_name, ref_name, chrm, file(target_info), file(target_frq), file(ref_frq) from combineInfo_frq
    output:
    set target_name, ref_name, file(output) into plot_r2_SNPpos
    script:
    info = target_info
    target = target_frq
    output = "${target_name}_${ref_name}_${chrm}_r2_SNPpos.png"
    template "r2_pos_plot.R"
}


"""
Plot frequency of imputed SNPs against SNP frequencies in reference panels
"""
process plot_freq_comparison {
    tag "plot_freq_comparison_${target_name}_${ref_name}_${chrm}"
    publishDir "${params.outDir}/plots/${ref_name}/freq_comparison", overwrite: true, mode:'copy'
    label "medium"
    input:
    set target_name, ref_name, chrm, file(target_info), file(target_frq), file(ref_frq) from combineInfo_frq_comp
    output:
    set target_name, ref_name, file(outputcolor) into plot_freq_comparison
    script:
    info = target_info
    target = target_frq
    frq = ref_frq
    //output = "${target_name}_${ref_name}_${chrm}_freq_comparison.png"
    outputcolor = "${target_name}_${ref_name}_${chrm}_freq_comparison_color.png"
    template "AF_comparison.R"
}


"""
Plot number of imputed SNPs over the mean r2 for all reference panels
"""
process plot_r2_SNPcount {
    tag "plot_r2_SNPcount_${target_name}_${ref_panels}_${chrms}"
    publishDir "${params.outDir}/plots/${ref_panels}", overwrite: true, mode:'copy'
    label "medium"
    input:
    set target_name, ref_panels, infos from target_infos.values()
    output:
    set target_name, ref_panels, file(plot_out) into plot_r2_SNPcount
    script:
    chrms = chromosomes_[target_name][0]+"-"+chromosomes_[target_name][-1]
    plot_out = "${target_name}_${ref_panels}_${chrms}_r2_SNPcount.png"
    infos = infos.join(',')
    impute_info_cutoff = params.impute_info_cutoff
    template "r2_Frequency_plot.R"
}


"""
Plot histograms of number of imputed SNPs over the mean r2 for all reference panels
"""
process plot_hist_r2_SNPcount {
    tag "plot_hist_r2_SNPcount_${target_name}_${ref_panels}_${chrms}"
    publishDir "${params.outDir}/plots/${ref_panels}/", overwrite: true, mode:'copy'
    label "medium"
    input:
    set target_name, ref_panels, infos from target_infos.values()
    output:
    set target_name, ref_panels, file(plot_out) into plot_hist_r2_SNPcount
    script:
    chrms = chromosomes_[target_name][0]+"-"+chromosomes_[target_name][-1]
    plot_out = "${target_name}_${ref_panels}_${chrms}_r2_SNPcount_hist.png"
    infos = infos.join(',')
    impute_info_cutoff = params.impute_info_cutoff
    template "r2_Frequency_plot_histogram.R"
}


"""
Plot MAF of imputed SNPs over r2 for all references
"""
process plot_MAF_r2 {
    tag "plot_MAF_r2_${target_name}_${ref_panels}_${chrms}"
    publishDir "${params.outDir}/plots/${ref_panels}", overwrite: true, mode:'copy'
    label "medium"
    input:
    set target_name, ref_panels, infos from target_infos.values()
    output:
    set target_name, ref_panels, file(plot_out) into plot_MAF_r2
    script:
    chrms = chromosomes_[target_name][0]+"-"+chromosomes_[target_name][-1]
    plot_out = "${target_name}_${ref_panels}_${chrms}_MAF_r2.png"
    infos = infos.join(',')
    impute_info_cutoff = params.impute_info_cutoff
    template "Frequency_r2_MAF_plot.R"
}

/*
 * Completion e-mail notification
 */
workflow.onComplete {
    // Set up the e-mail variables
    def subject = "[h3abionet/chipimputation] Successful: $workflow.runName"
    if(!workflow.success){
        subject = "[h3abionet/chipimputation] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if(workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if(workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    if(workflow.container) email_fields['summary']['Docker image'] = workflow.container
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // TODO nf-core: If not using MultiQC, strip out this code (including params.maxMultiqcEmailFileSize)
    // On success try attach the multiqc report
    def mqc_report = null
    try {
        if (workflow.success) {
            mqc_report = multiqc_report.getVal()
            if (mqc_report.getClass() == ArrayList){
                log.warn "[h3abionet/chipimputation] Found multiple reports from process 'multiqc', will use only one"
                mqc_report = mqc_report[0]
            }
        }
    } catch (all) {
        log.warn "[h3abionet/chipimputation] Could not attach MultiQC report to summary email"
    }

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir", mqcFile: mqc_report/*, mqcMaxSize: params.maxMultiqcEmailFileSize.toBytes() */]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // // Send the HTML e-mail
    // if (params.email) {
    //     try {
    //       if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
    //       // Try to send HTML e-mail using sendmail
    //       [ 'sendmail', '-t' ].execute() << sendmail_html
    //       log.info "[h3abionet/chipimputation] Sent summary e-mail to $params.email (sendmail)"
    //     } catch (all) {
    //       // Catch failures and try with plaintext
    //       [ 'mail', '-s', subject, params.email ].execute() << email_txt
    //       log.info "[h3abionet/chipimputation] Sent summary e-mail to $params.email (mail)"
    //     }
    // }

    // Write summary e-mail HTML to a file
    //def output_d = new File( "${params.outDir}/pipeline_info/" )
    def output_d = new File( "${params.outDir}/" )
    if( !output_d.exists() ) {
        output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";

    if (workflow.stats.ignoredCountFmt > 0 && workflow.success) {
        log.info "${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}"
        log.info "${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCountFmt} ${c_reset}"
        log.info "${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCountFmt} ${c_reset}"
    }

    if(workflow.success){
        log.info "${c_purple}[h3abionet/chipimputation]${c_green} Pipeline completed successfully${c_reset}"
    } else {
        checkHostname()
        log.info "${c_purple}[h3abionet/chipimputation]${c_red} Pipeline completed with errors${c_reset}"
    }

    // Copy the test config file to the current directory if test profile
    if ('test' in workflow.profile.split(',')) {
        println workflow.projectDir
        println workflow.scriptFile
        println workflow.profile
        log.info "${c_purple}[h3abionet/chipimputation]${c_green} Pipeline completed successfully${c_reset}"
    }

}

def checkHostname(){
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if(params.hostnames){
        def hostname = "hostname".execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if(hostname.contains(hname) && !workflow.profile.contains(prof)){
                    log.error "====================================================\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "============================================================"
                }
            }
        }
    }
}

def helpMessage() {
    log.info"""
    =========================================
    h3achipimputation v${params.version}
    =========================================
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run h3abionet/chipimputation -profile standard,docker
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

// println
