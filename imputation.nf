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

println "=================================================="

chromosomes = params.chromosomes.split(',')
chunk_size = params.chunk_size[0..-7]
ref_dir = file(params.ref_dir)

println "Project : $workflow.projectDir"
println "Git info: $workflow.repository - $workflow.revision [$workflow.commitId]"
println "Cmd line: $workflow.commandLine"
println "Chromosomes used: ${chromosomes.join(', ')}"

//// Help functions
datas = []
def chunk_split(map_file, chromosome, chunk=params.chunk_size) {
    '''
    Return: chunck files in the output folder
    '''
    data = map_file.readLines().collect{[it.split()[0], it.split()[3]]}
    data.each{ dat ->
        if( dat[0].toInteger() == chromosome.toInteger() && dat[0] != ''){
            datas << dat[1].toInteger()
        }
    }
    chunk_ = chunk.toInteger()
    max_ = datas.max()
    min_ = datas.min() - (datas.min() % 10) + 1
    myPos = (min_..max_).step(chunk_)
    chr_chunks = []
    for(pos in myPos){
        start_ = pos
        end_ = start_ + chunk_ - 1
        chr_chunks << "${start_} ${end_}\n"
    }
    return chr_chunks.join()
}

// Create channel for the study data from ped and map files
bed_data = Channel.fromPath(params.bedFile)
fam_data = Channel.fromPath(params.famFile)
bim_data = Channel.fromPath(params.bimFile)
ped_map_data = bed_data
                .merge(fam_data){ o,e -> [o,e] }
                .merge(bim_data){ o,e -> o + [e] }

// TODO check if files exist
// TODO check if chromosomes specified in config are in data
// TODO Be able to run everything on a specified chunk

"""
Create a channel with bed/fam/bim for each chromosome
"""
ped_map_data.into{ped_map_data;ped_map_data_chan}
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
ped_map_data_all.into{ped_map_data_all; ped_map_data__identifyChromosomes; ped_map_data__identifyChromosomes_1}
process generate_chunks {
    tag "generate_chunks_${chromosome}"
    input:
        set val(chromosome), file(data_bed), file(data_fam), val(bim_data) from ped_map_data__identifyChromosomes
    output:
        set val(chromosome), val(chunks) into identifyChromosomes, identifyChromosomes_all, identifyChromosomes_prehase
    exec:
        // This will always be run locally
        chunks = chunk_split(bim_data, chromosome, params.chunk_size)
}


"""
Plink bed/fam/bim per chromosomes
"""
ped_map_data_all.into{ped_map_data_all; ped_map_data__plink_to_chrm}
process plink_to_chrm {
    tag "plink2chrm_${chromosome}"
    publishDir "${params.output_dir}/qc/${chromosome}", overwrite: true, mode:'symlink'

    input:
        set val(chromosome), file(data_bed), file(data_fam), file(bim_data) from ped_map_data__plink_to_chrm
    output:
        set val(chromosome), file("${data_bed.baseName}.chr${chromosome}_clean.ped"), file("${data_bed.baseName}.chr${chromosome}_clean.map") into plink_to_chrm_all
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
            --recode --out ${data_bed.baseName}.chr${chromosome}_clean
        rm -f ${data_bed.baseName}.chr${chromosome}.* \
            ${data_bed.baseName}.chr${chromosome}_clean_mind.* \
            ${data_bed.baseName}.chr${chromosome}_clean_mind.dupvar
    """
}

"""
check the study genotypes
"""
plink_to_chrm_all.into{plink_to_chrm__checkGenotypes}
process checkGenotypes {
    tag "checkGenotypes_${chromosome}"
    validExitStatus 0,1,2
    errorStrategy 'ignore'
    publishDir "${params.output_dir}/checkGenotypes_shapeit/${chromosome}", overwrite: true, mode:'symlink'

    input:
        set val(chromosome), file(pedFile), file(mapFile) from plink_to_chrm__checkGenotypes
    output:
        set val(chromosome), file(pedFile), file(mapFile), file("${pedFile.baseName}.alignments.log"), file("${pedFile.baseName}.alignments.snp.strand.exclude") into checkGenotypes_all
    script:
        ref_hapFile = file( sprintf(params.ref_hapFile, chromosome) )
        ref_legendFile = file( sprintf(params.ref_legendFile, chromosome) )
        ref_sampleFile = file( params.ref_sampleFile )

        """
        shapeit -check \
            --input-ped ${pedFile} ${mapFile} \
            --input-ref ${ref_hapFile} ${ref_legendFile} ${ref_sampleFile} \
            --output-log ${pedFile.baseName}.alignments.log
        """
}


"""
Step 3: Pre-phase each chromosome using shapeit
"""
checkGenotypes_all.into{checkGenotypes_all; checkGenotypes_prephase}

process prephase {
    tag "prephase_${chromosome}"
    memory { 2.GB * task.attempt }
    publishDir "${params.output_dir}/prephased/${chromosome}", overwrite: true, mode:'symlink'

    input:
        set val(chromosome), file(data_ped), file(data_map), file(shapeit_check_log), file(shapeit_check_snp_strand_exclude) from checkGenotypes_prephase
        set val(chrom), val(chunks) from identifyChromosomes_prehase
    output:
        set val(chromosome), file("${data_ped.baseName}.prephased.haps"), file("${data_ped.baseName}.prephased.sample") into prephase_all, prephase__chunk, prephase__impute
    script:
        ref_hapFile = file( sprintf(params.ref_hapFile, chromosome) )
        ref_legendFile = file( sprintf(params.ref_legendFile, chromosome) )
        ref_sampleFile = file(params.ref_sampleFile)
        """
        shapeit -phase \
            --input-ped ${data_ped} ${data_map} \
            --input-ref ${ref_hapFile} ${ref_legendFile} ${ref_sampleFile} \
            --exclude-snp ${shapeit_check_snp_strand_exclude} \
            -T ${task.attempt} \
            -O ${data_ped.baseName}.prephased
        """
}


'''
Impute using the prephased genotypes with impute2
'''
chunk_data = identifyChromosomes_all.toSortedList().val[0][1..-1]
chunk_prephased = prephase__chunk.flatMap {
    chromosome      = it[0]
    haps            = it[1]
    sample          = it[2]
    ref_hapFile     = file( sprintf(params.ref_hapFile, chromosome) )
    ref_legendFile  = file( sprintf(params.ref_legendFile, chromosome) )
    ref_mapFile     = file( sprintf(params.ref_mapFile, chromosome) )
    ref_sampleFile  = file( params.ref_sampleFile )
    res = []
    for (chunks in chunk_data) {
        chunks = chunks.split()
        res << tuple(chromosome, chunks[0], chunks[1], haps, sample, ref_hapFile, ref_legendFile, ref_mapFile, ref_sampleFile)
    }
    res
}

process impute {
    tag "imp_chr${chromosome}_${chunkStart}-${chunkEnd}"
    memory { 2.GB * task.attempt }
    publishDir "${params.impute_result}/impute/${chromosome}", overwrite: true, mode:'symlink'

    input:
        set val(chromosome), val(chunkStart), val(chunkEnd), file(haps), file(sample), file(ref_hapFile), file(ref_legendFile), file(ref_mapFile), file(ref_sampleFile) from chunk_prephased
    output:
        set val(chromosome), file("${outfile}.imputed.gz"), file("${outfile}.imputed_haps.gz"), file("${outfile}.imputed_info"), file("${outfile}.imputed_summary") into impute_all, impute__sub
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


'''
Combine output
'''
impute_all.into{impute_all; impute__imputeCombine_cha}

// Create a dataflow instance of all impute results
imputeCombine_impute = []
imputeCombine_info = []
impute__imputeCombine_cha_list = impute__imputeCombine_cha.toSortedList().val
impute__imputeCombine_cha_list.each{chromosome, impute, haps, info, summary ->
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
    publishDir "${params.impute_result}/combined", overwrite: true, mode:'symlink'
    input:
        set val(chromosome), file(imputed_files) from imputeCombine_impute_cha
        set val(chromo), file(info_files) from imputeCombine_info_cha
    output:
        set val(chromosome), file(comb_impute), file(comb_info) into imputeCombine_all
    script:
        comb_impute = "${file(params.bedFile).getBaseName()}_chr${chromosome}.imputed.gz"
        comb_info = "${file(params.bedFile).getBaseName()}_chr${chromosome}.imputed_info"
        """
        zcat $imputed_files | bgzip -c > ${comb_impute}
        cat $info_files > ${comb_info}
        """
}


imputeCombine_all.into { imputeCombine_all; imputeCombine_toPlink }
prephase_all.into { prephase_all; prephase_toPlink}
println prephase_all.toSortedList().val.each { chromosome, prephased_haps, prephased_sample ->

}
process imputeToPlink {
    tag "toPlink_chr${chromosome}"
    memory { 2.GB * task.attempt }
    publishDir "${params.impute_result}/plink", overwrite: true, mode:'symlink'
    
    input:
        set val(chromosome), file(chromosome_imputed_gz) from imputeCombine_toPlink
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
            --make-bed --out ${chromosome_imputed_gz.baseName}
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

