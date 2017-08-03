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
// TODO
// Be able to run everything on a specified chunk

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
        set val(chromosome), val(chunks) into identifyChromosomes_all, identifyChromosomes_prehase
    script:
        chunks = chunk_split(bim_data, chromosome, params.chunk_size)
        """
        sleep 1
        """
}


"""
Plink bed/fam/bim per chromosomes
"""
ped_map_data_all.into{ped_map_data_all; ped_map_data__plink_to_chrm}
process plink_to_chrm {
    tag "plink2chrm_${chromosome}"
    publishDir "${params.output_dir}/data/qc/${chromosome}", overwrite: true, mode:'symlink'

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

///
// check the study genotypes
///
plink_to_chrm_all.into{plink_to_chrm__checkGenotypes}
process checkGenotypes {
    tag "checkGenotypes_${chromosome}"
    validExitStatus 0,1,2
    errorStrategy 'ignore'
    publishDir "${params.output_dir}/data/checkGenotypes_shapeit/${chromosome}", overwrite: true, mode:'symlink'

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
    memory { 2.GB * task.cpus }
    publishDir "${params.output_dir}/data/prephased/${chromosome}", overwrite: true, mode:'symlink'

    input:
        set val(chromosome), file(data_ped), file(data_map), file(shapeit_check_log), file(shapeit_check_snp_strand_exclude) from checkGenotypes_prephase
        set val(chrom), val(chunks) from identifyChromosomes_prehase
    output:
        set val(chromosome), file("${data_ped.baseName}.prephased.haps"), file("${data_ped.baseName}.prephased.sample") into prephase_all, prephase__chunk, prephase__impute
    script:
        ref_hapFile = file( sprintf(params.ref_hapFile, chromosome) )
        ref_legendFile = file( sprintf(params.ref_legendFile, chromosome) )
        ref_sampleFile = file( params.ref_sampleFile )
        """
        shapeit -phase \
            --input-ped ${data_ped} ${data_map} \
            --input-ref ${ref_hapFile} ${ref_legendFile} ${ref_sampleFile} \
            --exclude-snp ${shapeit_check_snp_strand_exclude} \
            -T ${task.cpus} \
            -O ${data_ped.baseName}.prephased
        """
}
prephase_all.subscribe{
    println "|-- Finished for ${it[1][-1]}, ${it[2][-1]}"
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
//    chunk_data = identifyChromosomes_all.toSortedList().val[0][1..-1]
    for (chunks in chunk_data) {
        chunks = chunks.split()
        res << tuple(chromosome, chunks[0], chunks[1], haps, sample, ref_hapFile, ref_legendFile, ref_mapFile, ref_sampleFile)
    }
    res
}

process impute {
    tag "imp_chr${chromosome}_${chunkStart}-${chunkEnd}"
    memory { 2.GB * task.cpus }
    publishDir "${params.impute_result}/${chromosome}", overwrite: true, mode:'symlink'

    input:
        set val(chromosome), val(chunkStart), val(chunkEnd), file(haps), file(sample), file(ref_hapFile), file(ref_legendFile), file(ref_mapFile), file(ref_sampleFile) from chunk_prephased
    output:
        set val(chromosome), file("chr${chromosome}_${chunkStart}-${chunkEnd}.imputed.gz"), file("chr${chromosome}_${chunkStart}-${chunkEnd}.imputed_haps.gz"), file("chr${chromosome}_${chunkStart}-${chunkEnd}.imputed_info"), file("chr${chromosome}_${chunkStart}-${chunkEnd}.imputed_summary") into impute_all, impute__sub
    shell:
        outfile = "chr${chromosome}_${chunkStart}-${chunkEnd}"
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
impute__sub.subscribe{
    println "|-- Finished for ${it[-1][-1]}"
}


'''
Combine output
'''
impute_all.into{impute_all; impute__imputeCombine_cha}

// Create a dataflow instance of all impute results
imputeCombine_map = []
impute__imputeCombine_cha_list = impute__imputeCombine_cha.toSortedList().val
impute__imputeCombine_cha_list.each{chromosome, impute, haps, info, summary ->
    imputeCombine_map << [chromosome, impute]
}
// Create a channel
imputeCombine_cha = Channel
        .from(imputeCombine_map)
        .groupTuple()

imputeCombine_cha.into { imputeCombine_cha; imputeCombine_cha1 }
process imputeCombine {
    tag "impComb_chr${chromosome}"
    cpus { 2 * task.attempt }
    memory { 2.GB * task.cpus }
    publishDir "${params.impute_result}/combined", overwrite: true, mode:'symlink'

    input:
        set val(chromosome), file(imputed_files) from imputeCombine_cha
    output:
        set val(chromosome), file("chr${chromosome}.imputed.gz") into imputeCombine_all
    script:
        """
        zcat $imputed_files | bgzip -c > chr${chromosome}.imputed.gz
        """
}
imputeCombine_all.subscribe{
    println "|-- Finished for ${it[1][-1]}"
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

