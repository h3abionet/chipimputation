#!/usr/bin/env nextflow

/*
 * Authors:
 *      Mamana Mbiyavanga
 *
 *  On behalf of the H3ABionet Consortium
 *  2017
 *
 *
 * Description  : Nextflow pipeline for ...
 *
*/

//---- General definitions --------------------------------------------------//

version = '1.0'

println "=================================================="


def helpMessage() {
    log.info"""
        ${version}
    """.stripIndent()
}

// Show help emssage
//if (params.help){
//    helpMessage()
//    exit 0
//}

// Check that Nextflow version is up to date enough
// try / throw / catch works for NF versions < 0.25 when this was implemented
nf_required_version = '0.25.0'
try {
    if( ! nextflow.version.matches(">= $nf_required_version") ){
        throw GroovyException('Nextflow version too old')
    }
} catch (all) {
    log.error "====================================================\n" +
            "  Nextflow version $nf_required_version required! You are running v$workflow.nextflow.version.\n" +
            "  Pipeline execution will continue, but things may break.\n" +
            "  Please run `nextflow self-update` to update Nextflow.\n" +
            "============================================================"
}

println "|-- Project : $workflow.projectDir"
println "|-- Git info: $workflow.repository - $workflow.revision [$workflow.commitId]"
println "|-- Command line: $workflow.commandLine"
println "|-- Datasets: ${file(params.bedFile).getName()}, ${file(params.bimFile).getName()}, ${file(params.famFile).getName()}"

// Check if chromosomes
if (params.chromosomes == '' || params.chromosomes == 'ALL'){
    chromosomes = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22"
}
chromosomes = params.chromosomes.split(',')
chromosomes_ = []
bimFile = file(params.bimFile)
// TODO do this in a process to control for memory
bimFile.eachLine{ line ->
    chrm = line.split()[0]
    if ( !chromosomes_.contains(chrm) ){
        chromosomes_ << chrm
    }
}

not_chrs = []
chromosomes.each { chromosome ->
    if (!(chromosome in chromosomes_)){
        not_chrs << chromosome
    }
}
if (!(not_chrs.isEmpty())){
    println "|-- Chromosome(s) ${not_chrs.join(', ')} not in dataset ${params.bedFile} and will be ignored."
    chromosomes = chromosomes_
}
println "|-- Chromosomes used: ${chromosomes.join(', ')}"


//// Help functions
datas = [:]
def chunk_split(map_file, chunk=params.chunk_size, chromosomes) {
    '''
    Return: chunck files in the output folder
    '''
    chromosomes = (chromosomes.split(',')).collect{ it.toInteger() }
    data = map_file.readLines().collect{ [it.split()[0], it.split()[3]] }
    data.each{ dat ->
        chromosome = dat[0].toInteger()
        myPos = dat[1].toInteger()
        if( chromosome != '' && chromosome in chromosomes){
            if ( !(chromosome in datas.keySet()) ){
                datas[chromosome] = []
            }
            datas[chromosome] << dat[1].toInteger()
        }
    }
    data = [:]
    datas.keySet().each { chromosome ->
        if ( !(chromosome in data.keySet()) ){
            data[chromosome] = []
        }
        chunk_ = chunk.toInteger()
        max_ = datas[chromosome].max()
        min_ = datas[chromosome].min() - (datas[chromosome].min() % 10) + 1
        myPos = (min_..max_).step(chunk_)
        myPos.each { pos ->
            start_ = pos
            end_ = start_ + chunk_ - 1
            data[chromosome] << [start_, end_]
        }
    }
    return data
}


// Create channel for the study data from ped and map files and check if files exist
bed_data = Channel
        .fromPath(params.bedFile)
        .ifEmpty { exit 1, "BED file not found: ${params.bedFile}" }
fam_data = Channel.fromPath(params.famFile)
bim_data = Channel.fromPath(params.bimFile)
bim_data.into{bim_data; bim_data_all}
ped_map_data = bed_data
        .merge(fam_data){ o,e -> [o,e] }
        .merge(bim_data){ o,e -> o + [e] }


// TODO check if files (study data and reference) exist if( !fasta.exists() ) exit 1, "Fasta file not found: ${params.fasta}"

// TODO Be able to run everything on a specified chunk

"""
Create a channel with bed/fam/bim for each chromosome
"""
ped_map_data.into{ ped_map_data; ped_map_data_chan }
ped_map_data_all = ped_map_data_chan.flatMap { bed_File, fam_File, bim_File ->
    def results = []
    chromosomes.each { chromosome ->
        results.push( [ chromosome, bed_File, fam_File, bim_File] )
    }
    return results
}


"""
Identify chromosomes and start/stop positions per chromosome and generate chunks
"""
ped_map_data_all.into{ ped_map_data_all; ped_map_data__identifyChromosomes; ped_map_data__identifyChromosomes_1}
bim_data_all.into{bim_data_all; bim_data_chunks}

process generate_chunks {
    tag "generate_chunks"
    input:
        val bim_data from bim_data_chunks
        val chromosome from chromosomes.join(',')
    output:
        val chunks into identifyChromosomes, identifyChromosomes_all, identifyChromosomes_prehase
    script:
        // This will always be run locally
        chunks = chunk_split(bim_data, params.chunk_size, chromosomes = chromosome)
        """
        """
}


"""
Plink bed/fam/bim per chromosomes
"""
ped_map_data_all.into{ped_map_data_all; ped_map_data_1}
process plink_to_chrm {
    tag "plink2chrm_${chromosome}"
    memory { 4.GB * task.attempt }
    clusterOptions  "-l nodes=1:ppn=${task.cpus}:series600"
    publishDir "${params.output_dir}/qc/${chromosome}", overwrite: true, mode:'symlink'

    input:
        set val(chromosome), file(data_bed), file(data_fam), file(bim_data) from ped_map_data_1
    output:
        set val(chromosome), file("${data_bed.baseName}.chr${chromosome}_clean.bed"), file("${data_bed.baseName}.chr${chromosome}_clean.bim"), file("${data_bed.baseName}.chr${chromosome}_clean.fam") into plink_to_chrm_all
    script:
        /*
        2- Exclude samples with missing more than 5% of genotype calls
        3- Two versions of you dataset, one with minor allele frequencies (MAFs) greater than 5% and one with MAFs less than 5%
        */
        """
        plink2 --bfile ${data_bed.baseName} \
            --chr ${chromosome} --allow-no-sex --make-bed \
            --out ${data_bed.baseName}.chr${chromosome}
        plink2 --bfile ${data_bed.baseName}.chr${chromosome} \
            --mind ${params.cut_mind} --allow-no-sex --recode \
            --out ${data_bed.baseName}.chr${chromosome}_clean_mind  
        plink2 --file ${data_bed.baseName}.chr${chromosome}_clean_mind \
            --allow-no-sex --list-duplicate-vars ids-only suppress-first \
            --out ${data_bed.baseName}.chr${chromosome}_clean_mind
        plink2 --file ${data_bed.baseName}.chr${chromosome}_clean_mind \
            --allow-no-sex \
            --exclude ${data_bed.baseName}.chr${chromosome}_clean_mind.dupvar \
            --make-bed --out ${data_bed.baseName}.chr${chromosome}_clean
        rm -f ${data_bed.baseName}.chr${chromosome}.* \
            ${data_bed.baseName}.chr${chromosome}_clean_mind.* \
            ${data_bed.baseName}.chr${chromosome}_clean_mind.dupvar
    """
}


"""
check the study genotypes
"""
plink_to_chrm_all.into{plink_to_chrm_all; plink_to_chrm__checkGenotypes}
process checkGenotypes {
    tag "checkGenotypes_${chromosome}"
    memory { 4.GB * task.attempt }
    validExitStatus 0,1,2
    errorStrategy 'ignore'
    publishDir "${params.output_dir}/qc/checkGenotypes_shapeit/${chromosome}", overwrite: true, mode:'symlink'

    input:
        set val(chromosome), file(bedFile), file(bimFile), file(famFile) from plink_to_chrm__checkGenotypes
    output:
        set val(chromosome), file(bedFile), file(bimFile), file(famFile), file("${bedFile.baseName}.alignments.log"), file("${bedFile.baseName}.alignments.snp.strand.exclude") into checkGenotypes_all
    script:
        ref_hapFile = file( sprintf(params.ref_hapFile, chromosome) )
        ref_legendFile = file( sprintf(params.ref_legendFile, chromosome) )
        ref_sampleFile = file( params.ref_sampleFile )
        """
        shapeit -check \
            --input-bed ${bedFile} ${bimFile} ${famFile} \
            --input-ref ${ref_hapFile} ${ref_legendFile} ${ref_sampleFile} \
            --output-log ${bedFile.baseName}.alignments.log
        """
}


"""
Step 3: Pre-phase each chromosome using shapeit
"""
checkGenotypes_all.into{checkGenotypes_all; checkGenotypes_prephase}

process prephase {
    tag "prephase_${chromosome}"
    memory { 4.GB * task.attempt }
    cpus { 2 * task.attempt }
    time { 4.h * task.attempt }
    clusterOptions  "-l nodes=1:ppn=${task.cpus}:series600"
    publishDir "${params.output_dir}/prephased/${chromosome}", overwrite: true, mode:'symlink'

    input:
        set val(chromosome), file(bedFile), file(bimFile), file(famFile), file(shapeit_check_log), file(shapeit_check_snp_strand_exclude) from checkGenotypes_prephase
    output:
        set val(chromosome), file("${bedFile.baseName}.phased.haps"), file("${bedFile.baseName}.phased.sample") into prephase_all
    script:
        ref_hapFile = file( sprintf(params.ref_hapFile, chromosome) )
        ref_legendFile = file( sprintf(params.ref_legendFile, chromosome) )
        ref_sampleFile = file( params.ref_sampleFile )
        """
        shapeit -phase \
            --input-bed ${bedFile} ${bimFile} ${famFile} \
            --input-ref ${ref_hapFile} ${ref_legendFile} ${ref_sampleFile} \
            --exclude-snp ${shapeit_check_snp_strand_exclude} \
            -T ${task.cpus} \
            -O ${bedFile.baseName}.phased
        """
//        ref_mapFile     = file( sprintf(params.ref_mapFile, chromosome) )
//        eagle \
//            --bfile=${data_ped.baseName} \
//            --geneticMapFile=${ref_mapFile} \
//            --chrom=${chromosome} \
//            --outPrefix=${data_ped.baseName}.phased \
//            --numThreads=${task.cpus} \
//            2>&1 | tee ${data_ped.baseName}.phased.log


}
prephase_all.into{ prephase_all; prephase_sub; prephase__chunk; prephase__impute }
prephase_sub.subscribe {
    println "|--- Finished ${it[0]}, ${it[2]}"
}

'''
Impute using the prephased genotypes with impute2
'''

chunk_data = identifyChromosomes_all.toSortedList().val[0]
chunk_prephased = prephase__chunk.flatMap {
    chromosome      = it[0]
    haps            = it[1]
    sample          = it[2]
    ref_hapFile     = file( sprintf(params.ref_hapFile, chromosome) )
    ref_legendFile  = file( sprintf(params.ref_legendFile, chromosome) )
    ref_mapFile     = file( sprintf(params.ref_mapFile, chromosome) )
    ref_sampleFile  = file( params.ref_sampleFile )
    res = []
    chunk_data[chromosome.toInteger()].each { chunks ->
        res << tuple(chromosome, chunks[0], chunks[1], haps, sample, ref_hapFile, ref_legendFile, ref_mapFile, ref_sampleFile)
    }
    res
}
process impute {
    tag "imp_chr${chromosome}_${chunkStart}-${chunkEnd}"
    memory { 4.GB * task.attempt }
    cpus { 2 * task.attempt }
    clusterOptions  "-l nodes=1:ppn=${task.cpus}:series600" // Specific to UCT HPC
    publishDir "${params.impute_result}/impute/${chromosome}", overwrite: true, mode:'symlink'

    input:
        set val(chromosome), val(chunkStart), val(chunkEnd), file(haps), file(sample), file(ref_hapFile), file(ref_legendFile), file(ref_mapFile), file(ref_sampleFile) from chunk_prephased
    output:
        set val(chromosome), file("${outfile}.imputed.gz"), file("${outfile}.imputed_haps.gz"), file("${outfile}.imputed_info"), file("${outfile}.imputed_summary") into impute_all
    shell:
        outfile = "${haps.baseName}_${chunkStart}-${chunkEnd}"
        """
        impute2 \
            -use_prephased_g \
            -known_haps_g ${haps} \
            -h ${ref_hapFile}  \
            -l ${ref_legendFile}  \
            -m ${ref_mapFile}  \
            -int ${chunkStart} ${chunkEnd} \
            -Ne 15000 \
            -buffer 250 \
            -phase \
            -o_gz \
            -o ${outfile}.imputed || true

        ## Sometimes there are no (type2) SNP's in a region
        if [ ! -f "${outfile}.imputed.gz" ]; then
            if grep 'ERROR: There are no type 2 SNPs after applying the command-line settings for this run' ${outfile}.imputed_summary || \
                grep 'Your current command-line settings imply that there will not be any SNPs in the output file, so IMPUTE2 will not perform any analysis or print output files.' ${outfile}.imputed_summary || \
                grep 'There are no SNPs in the imputation interval' ${outfile}.imputed_summary; then
                touch ${outfile}.imputed
                bgzip -f ${outfile}.imputed
                touch ${outfile}.imputed_haps
                bgzip -f ${outfile}.imputed_haps
                touch ${outfile}.imputed_info
            fi
        fi
        """
}
impute_all.into{ impute_all; impute_sub }
impute_sub.subscribe {
    println "|--- Finished ${it[0]}, ${it[-1]}"
}


'''
Combine output
'''
impute_all.into{impute_all; impute__imputeCombine_cha}

// Create a dataflow instance of all impute results
imputeCombine_impute = []
imputeCombine_info = []
impute__imputeCombine_cha_list = impute__imputeCombine_cha.toSortedList().val
impute__imputeCombine_cha_list.each{ chromosome, impute, haps, info, summary ->
    imputeCombine_impute << [chromosome, impute]
    imputeCombine_info << [chromosome, info]
}
// Create channels
imputeCombine_impute_cha = Channel
        .from(imputeCombine_impute)
        .groupTuple()
imputeCombine_info_cha = Channel
        .from(imputeCombine_info)
        .groupTuple()

process imputeCombine {
    tag "impComb_chr${chromosome}"
    memory { 2.GB * task.attempt }
    clusterOptions  "-l nodes=1:ppn=${task.cpus}:series600" // Specific to UCT HPC
    publishDir "${params.impute_result}/combined", overwrite: true, mode:'copy'
    input:
        set val(chromosome), file(imputed_files) from imputeCombine_impute_cha
        set val(chromo), file(info_files) from imputeCombine_info_cha
    output:
        set val(chromosome), file(comb_impute), file(comb_info) into imputeCombine_all
    script:
        comb_impute = "${file(params.bedFile).getBaseName()}_chr${chromosome}.imputed.gz"
        comb_info = "${file(params.bedFile).getBaseName()}_chr${chromosome}.imputed_info"
        """
        zcat ${imputed_files.join(' ')} | bgzip -c > ${comb_impute}
        head -n 1 ${info_files[0]} > ${comb_info}
        tail -q -n +2 ${info_files.join(' ')}>> ${comb_info}
        """
}
imputeCombine_all.into { imputeCombine_all; imputeCombine_sub; imputeCombine_toPlink }
imputeCombine_sub.subscribe{
    println "|--- Finished ${it.join(', ')}"
}


prephase_all.into { prephase_all; prephase_toPlink}
process imputeToPlink {
    tag "toPlink_chr${chromosome}"
    memory { 2.GB * task.attempt }
    publishDir "${params.impute_result}/plink", overwrite: true, mode:'copy'
    
    input:
        set val(chromosome), file(chromosome_imputed_gz), file(chromosome_imputed_info) from imputeCombine_toPlink
        set val(chrom), file(prephased_haps), file(prephased_sample) from prephase_toPlink
    output:
        set val(chromosome), file("${chromosome_imputed_gz.baseName}.bed"), file("${chromosome_imputed_gz.baseName}.bim"), file("${chromosome_imputed_gz.baseName}.fam") into imputeToPlink
    script:
        """
        gunzip -c ${chromosome_imputed_gz} > ${chromosome_imputed_gz.baseName}
        plink2 \
            --gen ${chromosome_imputed_gz.baseName} \
            --sample ${prephased_sample} \
            --oxford-single-chr ${chromosome} \
            --oxford-pheno-name plink_pheno \
            --hard-call-threshold 0.1 \
            --make-bed --out ${chromosome_imputed_gz.baseName} || true
        rm -f ${chromosome_imputed_gz.baseName}
        """
}

workflow.onComplete {
    def subject = 'My pipeline execution'
    def recipient = 'mypandos@gmail.com'

    ['mail', '-s', subject, recipient].execute() << """

    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    Error report: ${workflow.errorReport ?: '-'}
    """
}

