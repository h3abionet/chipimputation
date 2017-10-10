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
//
//
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
    chromosomes = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22".split(',')
}
else{
    chromosomes = params.chromosomes.split(',')
}

// Help functions

// Create channel for the study data from ped and map files and check if files exist
bed_data = Channel
        .fromPath(params.bedFile)
        .ifEmpty { exit 1, "BED file not found: ${params.bedFile}" }
fam_data = Channel
        .fromPath(params.famFile)
        .ifEmpty { exit 1, "FAM file not found: ${params.famFile}" }
bim_data = Channel
        .fromPath(params.bimFile)
        .ifEmpty { exit 1, "BIM file not found: ${params.bimFile}" }

// TODO check if files (study data and reference) exist if( !fasta.exists() ) exit 1, "Fasta file not found: ${params.fasta}"
// TODO Be able to run everything on a specified chunk


"""
Check user's provided chromosomes vs those in map file
"""
bim_data.into{bim_data; bim_data_check}
process check_chromosome {
    tag "check_chromosome_${bim_data.baseName}"
    memory { 8.GB * task.attempt }
    cpus { 4 * task.attempt }
    time { 4.h * task.attempt }
    clusterOptions  "-l nodes=1:ppn=${task.cpus}:series600"
    publishDir "${params.output_dir}", overwrite: true, mode:'symlink'
    
    input:
        file bim_data from bim_data_check
    output:
        file chromFile into check_chromosome
    script:
        chromFile = "${bim_data.baseName}_chromosomes.txt"
        """
        awk '{print \$1}' ${bim_data}  | sort | uniq >  ${chromFile}
        """
}
chromosomes_ = file(check_chromosome.toSortedList().val[0]).readLines() // Chromosomes from the dataset map file
not_chrs = []
in_chrs = []
chromosomes.each { chromosome ->
    if (!(chromosome in chromosomes_)){
        not_chrs << chromosome
    }
    else{
        in_chrs << chromosome
    }
}
old_chrs = chromosomes
if (in_chrs.isEmpty()){
    println "|-- Chromosome(s) ${old_chrs.join(', ')} not in dataset ${params.bedFile} and the pipeline will exit."
    exit 1
}
if (!(not_chrs.isEmpty())){
    println "|-- Chromosome(s) ${not_chrs.join(', ')} not in dataset ${params.bedFile} and will be ignored."
    chromosomes = in_chrs
}
println "|-- Chromosomes used: ${chromosomes.join(', ')}"


"""
Identify chromosomes and start/stop positions per chromosome and generate chunks
"""
bim_data.into{ bim_data; bim_data_chunks }

process generate_chunks {
    tag "generate_chunks_${bim_data.baseName}"
    memory { 8.GB * task.attempt }
    cpus { 4 * task.attempt }
    time { 4.h * task.attempt }
    clusterOptions  "-l nodes=1:ppn=${task.cpus}:series600"
    publishDir "${params.output_dir}", overwrite: true, mode:'copy'
    
    input:
        file bim_data from bim_data_chunks
    output:
        file chunkFile into generate_chunks
    script:
        chunkFile = "${bim_data.baseName}_chunks.txt"
        """
        python ${params.scripts}/generate_chunks.py ${bim_data} ${chunkFile} ${params.chunk_size}
        """
}
generate_chunks.into{generate_chunks; generate_chunks_sub; generate_chunks_all}
generate_chunks_sub.subscribe {
    println "|--- Finished ${it}"
}


"""
Plink bed/fam/bim per chromosomes
"""
bed_data.into{ bed_data; bed_data_plink }
bim_data.into{ bed_data; bim_data_plink }
fam_data.into{ bed_data; fam_data_plink }
generate_chunks.into{ generate_chunks; generate_chunks_plink}
process plink_to_chrm {
    tag "plink2chrm_${chromosome}_${data_bed.baseName}"
    memory { 4.GB * task.attempt }
    clusterOptions  "-l nodes=1:ppn=${task.cpus}:series600"
    publishDir "${params.output_dir}/qc/${chromosome}", overwrite: true, mode:'symlink'
    
    input:
        each chromosome from chromosomes
        file data_bed from bed_data_plink
        file data_bim from bim_data_plink
        file data_fam from fam_data_plink
        file chunkFile from generate_chunks_plink
    output:
        set val(chromosome), file("${data_bed.baseName}.chr${chromosome}_clean.bed"), file("${data_bed.baseName}.chr${chromosome}_clean.bim"), file("${data_bed.baseName}.chr${chromosome}_clean.fam") into plink_to_chrm
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
            --set-missing-var-ids @:# \
            --exclude ${data_bed.baseName}.chr${chromosome}_clean_mind.dupvar \
            --make-bed --out ${data_bed.baseName}.chr${chromosome}_clean
        rm -f ${data_bed.baseName}.chr${chromosome}.* \
            ${data_bed.baseName}.chr${chromosome}_clean_mind.* \
            ${data_bed.baseName}.chr${chromosome}_clean_mind.dupvar
    """
}
plink_to_chrm.into{ plink_to_chrm; plink_to_chrm_sub }
plink_to_chrm_sub.subscribe {
    println "|--- Finished ${it[0]}, ${it[2]}"
}

//"""
//check the study genotypes
//"""
//plink_to_chrm.into{plink_to_chrm; plink_to_chrm_checkGenotypes}
//process checkGenotypes {
//    tag "checkGenotypes_${chromosome}_${bedFile.baseName}"
//    memory { 4.GB * task.attempt }
//    validExitStatus 0,1,2
//    errorStrategy 'ignore'
//    publishDir "${params.output_dir}/qc/checkGenotypes_shapeit/${chromosome}", overwrite: true, mode:'symlink'
//    input:
//        set val(chromosome), file(bedFile), file(bimFile), file(famFile) from plink_to_chrm_checkGenotypes
//    output:
//        set val(chromosome), file(bedFile), file(bimFile), file(famFile), file("${bedFile.baseName}.alignments.log"), file("${bedFile.baseName}.alignments.snp.strand.exclude") into checkGenotypes
//    script:
//        ref_hapFile_1 = file( sprintf(params.ref_hapFile_1, chromosome) )
//        ref_legendFile_1 = file( sprintf(params.ref_legendFile_1, chromosome) )
//        ref_sampleFile_1 = file( params.ref_sampleFile_1 )
//        """
//        shapeit -check \
//            --input-bed ${bedFile} ${bimFile} ${famFile} \
//            --input-ref ${ref_hapFile_1} ${ref_legendFile_1} ${ref_sampleFile_1} \
//            --output-log ${bedFile.baseName}.alignments.log
//        """
//}


"""
Step 3: Pre-phase each chromosome using shapeit
"""
plink_to_chrm.into{plink_to_chrm; plink_to_chrm_1}
process prephase {
    tag "prephase_${chromosome}_${bedFile.baseName}"
    memory { 8.GB * task.attempt }
    cpus { 4 * task.attempt }
    time { 4.h * task.attempt }
    clusterOptions  "-l nodes=1:ppn=${task.cpus}:series600"
    publishDir "${params.output_dir}/prephased/${chromosome}", overwrite: true, mode:'symlink'
    input:
        set val(chromosome), file(bedFile), file(bimFile), file(famFile) from plink_to_chrm_1
    output:
        set val(chromosome), file("${vcf_out}.haps.gz"), file("${vcf_out}.sample") into prephase_all
    script:
        ref_hapFile_1 = file( sprintf(params.ref_hapFile_1, chromosome) )
        ref_legendFile_1 = file( sprintf(params.ref_legendFile_1, chromosome) )
        ref_sampleFile_1 = file( params.ref_sampleFile_1 )
        vcf_out = "${bedFile.baseName}.phased"
        """
        eagle \
            --bfile=${bedFile.baseName} \
            --geneticMapFile=${params.eagle_genetic_map} \
            --chrom=${chromosome} \
            --genoErrProb 0.003 --pbwtOnly \
            --allowRefAltSwap \
            --outPrefix=${vcf_out} 2>&1 | tee ${vcf_out}.log
        """

}
prephase_all.into{ prephase_all; prephase_sub; prephase__chunk; prephase__impute }
prephase_sub.subscribe {
    println "|--- Finished ${it[0]}, ${it[1]}"
}
if ( params.ref_1_name != null){
    if ( !('ref_2_name' in params.keySet()) || params.ref_2_name.length() == 0){
        '''
        Impute using the prephased genotypes with impute2 with 1 reference panel
        '''
        chunk_data = file(generate_chunks_all.toSortedList().val[0]).readLines()
        chunk_prephased = prephase__chunk.flatMap { chromosome, haps, sample ->
            ref_hapFile_1     = file( sprintf(params.ref_hapFile_1, chromosome) )
            ref_legendFile_1  = file( sprintf(params.ref_legendFile_1, chromosome) )
            ref_mapFile_1     = file( sprintf(params.ref_mapFile_1, chromosome) )
            ref_sampleFile_1  = file( params.ref_sampleFile_1 )
            tmp = []
            res = []
            chunk_data.each { chunk_dat ->
                chunks = chunk_dat.split()
                if (chromosome.toInteger() == chunks[0].toInteger()){
                    if(!(chunk_dat in tmp)){
                        res << [chromosome, chunks[1], chunks[2], haps, sample, ref_hapFile_1, ref_legendFile_1, ref_mapFile_1, ref_sampleFile_1]
                        tmp << chunk_dat
                    }
                }
            }
            res
        }
        
        process impute {
            tag "imp_chr${chromosome}_${chunkStart}-${chunkEnd}_${haps.baseName}"
            memory { 4.GB * task.attempt }
            cpus { 2 * task.attempt }
            clusterOptions  "-l nodes=1:ppn=${task.cpus}:series600" // Specific to UCT HPC
            publishDir "${params.impute_result}/impute/${chromosome}", overwrite: true, mode:'symlink'
        
            input:
                set val(chromosome), val(chunkStart), val(chunkEnd), file(haps), file(sample), file(ref_hapFile_1), file(ref_legendFile_1), file(ref_mapFile_1), file(ref_sampleFile_1) from chunk_prephased
            output:
                set val(chromosome), file("${outfile}.imputed.gz"), file("${outfile}.imputed_haps.gz"), file("${outfile}.imputed_info"), file("${outfile}.imputed_summary") into impute_all
            shell:
                outfile = "${haps.baseName}_${chunkStart}-${chunkEnd}"
                """
                impute2 \
                    -use_prephased_g \
                    -known_haps_g ${haps} \
                    -h ${ref_hapFile_1}  \
                    -l ${ref_legendFile_1}  \
                    -m ${ref_mapFile_1}  \
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
    }
    else{
        tmp = []
        refs_data = []
        chunk_data = file(generate_chunks_all.toSortedList().val[0]).readLines()
        chunk_prephased = prephase__chunk.toSortedList().val
        chunk_prephased.each { chromosome, study_haps, study_sample ->
            chunk_data.each { chunk_dat ->
                chunks = chunk_dat.split()
                chrm = chunks[0]
                chunkStart = chunks[1]
                chunkEnd = chunks[2]
                ref_hapFile_1 = file(sprintf(params.ref_hapFile_1, chrm))
                ref_legendFile_1 = file(sprintf(params.ref_legendFile_1, chrm))
                ref_mapFile_1 = file(sprintf(params.ref_mapFile_1, chrm))
                ref_sampleFile_1 = file(params.ref_sampleFile_1)
                ref_hapFile_2 = file(sprintf(params.ref_hapFile_2, chrm))
                ref_legendFile_2 = file(sprintf(params.ref_legendFile_2, chrm))
                ref_mapFile_2 = file(sprintf(params.ref_mapFile_2, chrm))
                ref_sampleFile_2 = file(params.ref_sampleFile_2)
                if ( chrm in chromosomes ) {
                    if ( chrm.toInteger() == chromosome.toInteger() ){
                        if (!(chunk_dat in tmp)) {
                            refs_data << [chrm, chunkStart, chunkEnd, study_haps, study_sample, ref_hapFile_1, ref_legendFile_1, ref_mapFile_1, ref_sampleFile_1, ref_hapFile_2, ref_legendFile_2, ref_mapFile_2, ref_sampleFile_2]
                            tmp << chunk_dat
                        }
                    }
                }
            }
        }
        process cross_impute_2refs {
            tag "cross_impute_2refs_chr${chromosome}_${chunkStart}-${chunkEnd}_${params.ref_1_name}_${params.ref_2_name}"
            memory { 6.GB * task.attempt }
            cpus { 3 * task.attempt }
            clusterOptions "-l nodes=1:ppn=${task.cpus}:series600" // Specific to UCT HPC
            publishDir "${params.impute_result}/cross_impute_${params.ref_1_name}_${params.ref_2_name}/${chromosome}", overwrite: true, mode: 'symlink'
            input:
                set val(chromosome), val(chunkStart), val(chunkEnd), file(study_haps), file(study_sample), file(ref_hapFile_1), file(ref_legendFile_1), file(ref_mapFile_1), file(ref_sampleFile_1), file(ref_hapFile_2), file(ref_legendFile_2), file(ref_mapFile_2), file(ref_sampleFile_2) from refs_data
            output:
                set val(chromosome), val(chunkStart), val(chunkEnd), file(study_haps), file(study_sample), file(haps_file), file(legend_file), file(ref_mapFile_1) into cross_impute_2refs
            shell:
                outfile = "${params.ref_1_name}_${params.ref_2_name}_chr${chromosome}_${chunkStart}-${chunkEnd}.impute2"
                haps_file = "${outfile}.hap.gz"
                legend_file = "${outfile}.legend.gz"
                buffer = (params.chunk_size.toInteger()/4000).toInteger() // in kb
                """
                impute2 \
                    -merge_ref_panels \
                    -m ${ref_mapFile_1} \
                    -h ${ref_hapFile_1} ${ref_hapFile_2} \
                    -l ${ref_legendFile_1} ${ref_legendFile_2} \
                    -Ne ${params.NE} \
                    -burnin ${params.impute_burnin} \
                    -iter ${params.impute_iter} \
                    -buffer ${buffer} \
                    -int ${chunkStart} ${chunkEnd} \
                    -include_buffer_in_output \
                    -merge_ref_panels_output_ref ${outfile} \
                    -fill_holes \
                    -no_remove \
                    -o ${outfile} \
                    -o_gz \
                    || true
                #rm ${ref_hapFile_1.baseName} ${ref_hapFile_2.baseName} ${params.ref_1_name}_${ref_legendFile_1.baseName} ${params.ref_2_name}_${ref_legendFile_2.baseName}
                """
        }
    }
}

//
//'''
//Combine output
//'''
//impute_all.into{impute_all; impute__imputeCombine_cha}
//
//// Create a dataflow instance of all impute results
//imputeCombine_impute = []
//imputeCombine_info = []
//impute__imputeCombine_cha_list = impute__imputeCombine_cha.toSortedList().val
//impute__imputeCombine_cha_list.each{ chromosome, impute, haps, info, summary ->
//    imputeCombine_impute << [chromosome, impute]
//    imputeCombine_info << [chromosome, info]
//}
//// Create channels
//imputeCombine_impute_cha = Channel
//        .from(imputeCombine_impute)
//        .groupTuple()
//imputeCombine_info_cha = Channel
//        .from(imputeCombine_info)
//        .groupTuple()
//
//process imputeCombine {
//    tag "impComb_chr${chromosome}"
//    memory { 2.GB * task.attempt }
//    clusterOptions  "-l nodes=1:ppn=${task.cpus}:series600" // Specific to UCT HPC
//    publishDir "${params.impute_result}/combined", overwrite: true, mode:'copy'
//    input:
//        set val(chromosome), file(imputed_files) from imputeCombine_impute_cha
//        set val(chromo), file(info_files) from imputeCombine_info_cha
//    output:
//        set val(chromosome), file(comb_impute), file(comb_info) into imputeCombine_all
//    script:
//        comb_impute = "${file(params.bedFile).getBaseName()}_chr${chromosome}.imputed.gz"
//        comb_info = "${file(params.bedFile).getBaseName()}_chr${chromosome}.imputed_info"
//        """
//        zcat ${imputed_files.join(' ')} | bgzip -c > ${comb_impute}
//        head -n 1 ${info_files[0]} > ${comb_info}
//        tail -q -n +2 ${info_files.join(' ')}>> ${comb_info}
//        """
//}
//imputeCombine_all.into { imputeCombine_all; imputeCombine_sub; imputeCombine_toPlink }
//imputeCombine_sub.subscribe{
//    println "|--- Finished ${it.join(', ')}"
//}
//
//
//prephase_all.into { prephase_all; prephase_toPlink}
//process imputeToPlink {
//    tag "toPlink_chr${chromosome}"
//    memory { 2.GB * task.attempt }
//    publishDir "${params.impute_result}/plink", overwrite: true, mode:'copy'
//
//    input:
//        set val(chromosome), file(chromosome_imputed_gz), file(chromosome_imputed_info) from imputeCombine_toPlink
//        set val(chrom), file(prephased_haps), file(prephased_sample) from prephase_toPlink
//    output:
//        set val(chromosome), file("${chromosome_imputed_gz.baseName}.bed"), file("${chromosome_imputed_gz.baseName}.bim"), file("${chromosome_imputed_gz.baseName}.fam") into imputeToPlink
//    script:
//        """
//        gunzip -c ${chromosome_imputed_gz} > ${chromosome_imputed_gz.baseName}
//        plink2 \
//            --gen ${chromosome_imputed_gz.baseName} \
//            --sample ${prephased_sample} \
//            --oxford-single-chr ${chromosome} \
//            --hard-call-threshold 0.1 \
//            --make-bed --out ${chromosome_imputed_gz.baseName} || true
//        rm -f ${chromosome_imputed_gz.baseName}
//        """
//        //--oxford-pheno-name plink_pheno
//}
//
//workflow.onComplete {
//    def subject = 'My pipeline execution'
//    def recipient = 'mamana.mbiyavanga@uct.ac.za'
//
//    ['mail', '-s', subject, recipient].execute() << """
//
//    Pipeline execution summary
//    ---------------------------
//    Completed at: ${workflow.complete}
//    Duration    : ${workflow.duration}
//    Success     : ${workflow.success}
//    workDir     : ${workflow.workDir}
//    exit status : ${workflow.exitStatus}
//    Error report: ${workflow.errorReport ?: '-'}
//    """
//}
