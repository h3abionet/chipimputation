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
nf_required_version = '0.31.1'
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
//println "|-- Git info: $workflow.repository - $workflow.revision [$workflow.commitId]"
println "|-- Command line: $workflow.commandLine"
println "|-- Datasets: ${params.target_datasets.values().join(', ')}"
//println "|-- Datasets: ${file(params.vcfFile).getName()}"

// Help functions

// check if study genotype files exist
target_datasets = []
params.target_datasets.each { target ->
    if (!file(target.value).exists()) exit 1, "VCF file ${target.value} not found. Please check your config file."
    target_datasets << [target.key, file(target.value)]
}
if (!file(params.eagle_genetic_map).exists()) exit 1, "MAP file ${params.eagle_genetic_map} not found. Please check your config file."

//// Create channel for the study data from ped and map files
target_datasets = Channel
        .from(target_datasets)

// TODO Be able to run just on a specified chunk


"""
Check user's provided chromosomes vs those in map file
"""
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

check_chromosome.into{ check_chromosome; check_chromosome1 }
// Check if specified chromosomes exist in bim file or VCF
chromosomes_ = []
check_chromosome1.toSortedList().val.each{ check_file ->
    chromosomes_ = chromosomes_ + file(check_file).readLines()
}
chromosomes_ = chromosomes_.unique().collect { it as int }.sort()
// Chromosomes from the dataset map file
// Check if chromosomes
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
            println "|-- No Chromosome(s) found not in target(s) dataset(s)! The pipeline will exit."
            exit 1
        }
        else{
            chromosomes = chromosomes_
        }
    }

    if (!(not_chrs.isEmpty())){
        println "|-- Chromosome(s) ${not_chrs.join(', ')} not in target datasets and will be ignored."
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


"""
Identify chromosomes and start/stop positions per chromosome and generate chunks
"""
mapFile_cha.into{ mapFile_cha; mapFile_cha_chunks }
process generate_chunks {
    tag "generate_chunks_${target_name}"
//    publishDir "${params.output_dir}", overwrite: true, mode:'copy'
    echo true
    input:
        set val(target_name), file(mapFile) from mapFile_cha_chunks
    output:
        set val(target_name), file(chunkFile) into generate_chunks
    script:
        chunkFile = "${mapFile.baseName}_chunks.txt"
        chunk_size = params.chunk_size
        chromosomes = chromosomes.join(',')
        template "generate_chunks.py"
}


"""
QC
"""
target_datasets.into{ target_datasets; target_datasets_qc }
process check_mismach {
    tag "check_mismach_${target_name}"
//    publishDir "${params.output_dir}", overwrite: true, mode:'symlink'
    input:
        set val(target_name), file(target_vcfFile) from target_datasets_qc
    output:
        set val(target_name), file(target_vcfFile), file("${base}_checkRef_warn.log"), file("${base}_checkRef_summary.log") into check_mismach
    script:
        base = file(target_vcfFile.baseName).baseName
        """
        nblines=\$(zcat ${target_vcfFile} | wc -l)
        if (( \$nblines > 1 ))
        then
            bcftools norm --check-ref w \
                -f ${params.reference_genome} \
                ${target_vcfFile} \
                -Ou -o /dev/null
            cp .command.err ${base}_checkRef_warn.log
            bcftools +fixref \
                ${target_vcfFile} \
                -- \
                -f ${params.reference_genome} \
                2>&1 | tee "${base}_checkRef_summary.log"
            rm -f ${base}_clean_mind.*
        fi

        """
}

check_mismach.into{ check_mismach; check_mismach_1 }
check_mismach_noMis = Channel.create()
check_mismach_1.toSortedList().val.each{ target_name, target_vcfFile, warn, summary ->
    mismach = file(warn).readLines().size()-1
    if ( mismach != 0 ) {
        System.err.println "|-- ${mismach} ref mismach sites found in '${target_name}' dataset! The pipeline will exit."
        exit 1
    }
    else{
        check_mismach_noMis << [ target_name, target_vcfFile, warn, summary ]
    }
}
check_mismach_noMis.close()

"""
QC
"""
check_mismach_noMis.into{ check_mismach_noMis; check_mismach_noMis_1 }
process target_qc {
    tag "target_qc_${target_name}"
//    publishDir "${params.output_dir}", overwrite: true, mode:'symlink'
    input:
        set val(target_name), file(target_vcfFile), file(mismach_warn), file(mismach_summary) from check_mismach_noMis_1
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
        plink2 --vcf ${base}_noALT.vcf.gz \
            --keep-allele-order \
            --list-duplicate-vars ids-only suppress-first \
            --allow-no-sex \
            --recode vcf-iid \
            --out ${base}_clean_mind
        plink2 --vcf ${base}_clean_mind.vcf \
            --keep-allele-order \
            --exclude ${base}_clean_mind.dupvar \
            --recode vcf-iid \
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
                params.ref_panels.each { ref ->
                    ref_m3vcf = sprintf(params.ref_panels[ref.key].m3vcfFile, chrm)
                    ref_vcf = sprintf(params.ref_panels[ref.key].vcfFile, chrm)
                    chunks_datas << [chrm, chunk_start, chunk_end, target_name, file(target_vcfFile), ref.key, file(ref_vcf), file(ref_m3vcf)]
                }
            }
        }
    }
    return chunks_datas
}
target_qc_chunk = target_qc_1
        .flatMap{ it -> transform_chunk(it) }

process split_target_to_chunk {
    tag "split_${target_name}_${chrm}:${chunk_start}-${chunk_end}"
    publishDir "${params.output_dir}/qc/${chrm}/chunks", overwrite: true, mode:'symlink'
    input:
        set chrm, chunk_start, chunk_end, target_name, file(target_vcfFile), ref_name, file(ref_vcf), file(ref_m3vcf) from target_qc_chunk
    output:
        set chrm, chunk_start, chunk_end, target_name, file(target_vcfFile_chunk), ref_name, file(ref_vcf), file(ref_m3vcf) into split_vcf_to_chrm
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


"""
Phase each chromosome using eagle
"""
split_vcf_to_chrm.into{ split_vcf_to_chrm; split_vcf_to_chrm_1 }
process phase_target_chunk {
    tag "phase_${target_name}_${chrm}:${chunk_start}-${chunk_end}_${ref_name}"
    input:
        set chrm, chunk_start, chunk_end, target_name, file(target_vcfFile_chunk), ref_name, file(ref_vcf), file(ref_m3vcf) from split_vcf_to_chrm_1
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
                --geneticMapFile=${params.eagle_genetic_map} \
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
// TODO Record failed phasing chunks in report

//phase_data.into{ phase_data; phase_data_1; phase_data_2 }
//def combine_chunk_data_with_ref = { chromosome, chunk_start, chunk_end, study_haps, study_sample ->
//    res = []
//    ref_hapFile     = file( sprintf(params.ref_1.hapFile, chromosome) )
//    ref_legendFile  = file( sprintf(params.ref_1.legendFile, chromosome) )
//    ref_mapFile     = file( sprintf(params.ref_1.mapFile, chromosome) )
//    ref_sampleFile  = file( params.ref_1.sampleFile )
//    res << [chromosome, chunk_start, chunk_end, study_haps, study_sample, ref_hapFile, ref_legendFile, ref_mapFile, ref_sampleFile]
//    return res
//}
//def combine_chunk_data_with_2refs = { chrm, chunk_start, chunk_end, study_haps, study_sample ->
//    res = []
//    ref_1_hapFile = file(sprintf(params.ref_1.hapFile, chrm))
//    ref_1_legendFile = file(sprintf(params.ref_1.legendFile, chrm))
//    ref_1_mapFile = file(sprintf(params.ref_1.mapFile, chrm))
//    ref_1_sampleFile = file(params.ref_1.sampleFile)
//    ref_2_hapFile = file(sprintf(params.ref_2.hapFile, chrm))
//    ref_2_legendFile = file(sprintf(params.ref_2.legendFile, chrm))
//    ref_2_mapFile = file(sprintf(params.ref_2.mapFile, chrm))
//    ref_2_sampleFile = file(params.ref_2.sampleFile)
//    res << [chrm, chunk_start, chunk_end, study_haps, study_sample, ref_1_hapFile, ref_1_legendFile, ref_2_mapFile, ref_1_sampleFile, ref_2_hapFile, ref_2_legendFile, ref_2_sampleFile]
//    return res
//}
//if ( params.ref_1.name != null){
//    if ( !('ref_2' in params.keySet()) || params.ref_2.name.length() == 0){
//        '''
//        Impute using the prephased genotypes with impute2 with 1 reference panel
//        '''
//        chunk_prephased = phase_data_1.
//                flatMap { it -> combine_chunk_data_with_ref(it) }
//    }
//    else{
//        cross_refs_data = phase_data_1.
//                flatMap { it -> combine_chunk_data_with_2refs(it) }
//        process cross_impute_2refs {
//            tag "cross_impute_2refs_chr${chromosome}_${chunkStart}-${chunkEnd}_${params.ref_1.name}_${params.ref_2.name}"
//            memory { 15.GB * task.attempt }
//            time { 10.h * task.attempt }
//            publishDir "${params.impute_result}/cross_impute_${params.ref_1.name}_${params.ref_2.name}/${chromosome}", overwrite: true, mode: 'symlink'
//            input:
//                set val(chromosome), val(chunkStart), val(chunkEnd), file(study_haps), file(study_sample), file(ref_1_hapFile), file(ref_1_legendFile), file(ref_2_mapFile), file(ref_1_sampleFile), file(ref_2_hapFile), file(ref_2_legendFile), file(ref_2_sampleFile) from cross_refs_data
//            output:
//                set val(chromosome), val(chunkStart), val(chunkEnd), file(study_haps), file(study_sample), file(haps_file), file(legend_file), file(ref_2_mapFile), file(sample_file) into cross_impute_2refs
//            shell:
//                outfile = "${params.ref_1.name}_${params.ref_2.name}_chr${chromosome}_${chunkStart}-${chunkEnd}.impute2"
//                haps_file = "${outfile}.hap.gz"
//                legend_file = "${outfile}.legend.gz"
//                sample_file = "${outfile}.sample"
//                """
//                gunzip -c ${ref_1_hapFile} > ${ref_1_hapFile.baseName}
//                gunzip -c ${ref_2_hapFile} > ${ref_2_hapFile.baseName}
//                gunzip -c ${ref_1_legendFile} > ${ref_1_legendFile.baseName}
//                gunzip -c ${ref_2_legendFile} > ${ref_2_legendFile.baseName}
//                impute2 \
//                    -merge_ref_panels \
//                    -m ${ref_2_mapFile} \
//                    -h ${ref_1_hapFile.baseName} ${ref_2_hapFile.baseName} \
//                    -l ${ref_1_legendFile.baseName} ${ref_2_legendFile.baseName} \
//                    -Ne ${params.NE} \
//                    -burnin ${params.impute_burnin} \
//                    -iter ${params.impute_iter} \
//                    -buffer ${params.buffer_size} \
//                    -int ${chunkStart} ${chunkEnd} \
//                    -include_buffer_in_output \
//                    -merge_ref_panels_output_ref ${outfile} \
//                    -fill_holes \
//                    -no_remove \
//                    -o ${outfile} \
//                    -o_gz \
//                    || true
//                head -n 1 ${ref_1_sampleFile} > ${sample_file}
//                tail -q -n +2 ${ref_1_sampleFile} ${ref_2_sampleFile} >> ${sample_file}
//                rm -f ${ref_1_hapFile.baseName} ${ref_2_hapFile.baseName} ${ref_1_legendFile.baseName} ${ref_2_legendFile.baseName}
//                ## Sometimes there are no (type2) SNP's in a region
//                if [ ! -f "${outfile}.hap.gz" ]; then
//                    if grep 'ERROR: There are no type 2 SNPs after applying the command-line settings for this run' ${outfile}_summary || \
//                        grep 'Your current command-line settings imply that there will not be any SNPs in the output file, so IMPUTE2 will not perform any analysis or print output files.' ${outfile}_summary || \
//                        grep 'There are no SNPs in the imputation interval' ${outfile}_summary; then
//                        touch ${outfile}.hap
//                        bgzip -f ${outfile}.hap.gz
//                        touch ${outfile}.legend
//                        bgzip -f ${outfile}.legend.gz
//                        touch ${outfile}.sample
//                    fi
//                fi
//                """
//        }
//        cross_impute_2refs.into{ cross_impute_2refs; chunk_prephased}
//    }

process impute_target {
    tag "imp_${target_name}_${chrm}:${chunk_start}-${chunk_end}_${ref_name}"
    publishDir "${params.impute_result}/impute/${chromosome}", overwrite: true, mode:'symlink'
    input:
        set chrm, chunk_start, chunk_end, target_name, file(target_phased_vcfFile), ref_name, file(ref_vcf), file(ref_m3vcf) from phase_target
    output:
        set chrm, chunk_start, chunk_end, target_name, ref_name, file("${base}_imputed.dose.vcf.gz"), file("${base}_imputed.info") into impute_target
    shell:
        base = "${file(target_phased_vcfFile.baseName).baseName}"
        """
        minimac4 \
            --refHaps ${ref_m3vcf} \
            --haps ${target_phased_vcfFile} \
            --format GT \
            --allTypedSites \
            --chr ${chrm} --start ${chunk_start} --end ${chunk_end} --window ${params.flanking_region} \
            --prefix ${base}_imputed
        """
}

//impute2 \
//            -use_prephased_g \
//            -known_haps_g ${study_haps} \
//            -h ${ref_hapFile}  \
//            -l ${ref_legendFile}  \
//            -m ${ref_mapFile}  \
//            -int ${chunkStart} ${chunkEnd} \
//            -Ne 15000 \
//            -buffer ${params.buffer_size} \
//            -phase \
//            -o_gz \
//            -o ${outfile}.imputed || true
//
//## Sometimes there are no (type2) SNP's in a region
//if [ ! -f "${outfile}.imputed.gz" ]; then
//if grep 'ERROR: There are no type 2 SNPs after applying the command-line settings for this run' ${outfile}.imputed_summary || \
//                grep 'Your current command-line settings imply that there will not be any SNPs in the output file, so IMPUTE2 will not perform any analysis or print output files.' ${outfile}.imputed_summary || \
//                grep 'There are no SNPs in the imputation interval' ${outfile}.imputed_summary; then
//touch ${outfile}.imputed
//bgzip -f ${outfile}.imputed
//touch ${outfile}.imputed_haps
//bgzip -f ${outfile}.imputed_haps
//touch ${outfile}.imputed_info
//fi
//fi
//}
//
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
//    publishDir "${params.impute_result}/combined", overwrite: true, mode:'copy'
//    input:
//        set chromosome, file(imputed_files) from imputeCombine_impute_cha
//        set chromo, file(info_files) from imputeCombine_info_cha
//    output:
//        set chromosome, file(comb_impute) into imputeCombine
//        set chromosome, file(comb_info) into infoCombine
//    script:
//        comb_impute = "${file(params.bedFile).getBaseName()}_chr${chromosome}.imputed.gz"
//        comb_info = "${file(params.bedFile).getBaseName()}_chr${chromosome}.imputed_info"
//        """
//        zcat ${imputed_files.join(' ')} | bgzip -c > ${comb_impute}
//        head -n 1 ${info_files[0]} > ${comb_info}
//        tail -q -n +2 ${info_files.join(' ')}>> ${comb_info}
//        """
//}
//
//
//imputeCombine.into { imputeCombine; imputeCombine_1 }
//process imputeToPlink {
//    tag "toPlink_chr${chromosome}"
//    memory { 2.GB * task.attempt }
//    publishDir "${params.impute_result}/plink", overwrite: true, mode:'copy'
//    input:
//        set chromosome, file(chromosome_imputed_gz) from imputeCombine_1
//        set chrom, chunk_start, chunk_end, file(prephased_haps), file(prephased_sample) from phase_data_2
//    output:
//        set chromosome, file("${chromosome_imputed_gz.baseName}.bed"), file("${chromosome_imputed_gz.baseName}.bim"), file("${chromosome_imputed_gz.baseName}.fam") into imputeToPlink
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
//}
//
//
//"""
//Generating report
//"""
//infoCombine.into { infoCombine; infoCombine_1 }
//infoCombine_list = infoCombine_1.toSortedList().val
//chrm_infos = []
//infoCombine_list.each{ chrm, comb_info ->
//    chrm_infos << chrm+"=="+comb_info
//}
//process filter_info {
//    tag "filter_${projectName}_${chrms}"
//    memory { 2.GB * task.attempt }
//    publishDir "${params.output_dir}/INFOS/${projectName}", overwrite: true, mode:'copy'
//    input:
//        val(projectName) from params.project_name
//    output:
//        set val(projectName), file(well_out) into info_Well
//        set val(projectName), file(acc_out) into info_Acc
//    script:
//        chrms = chromosomes[0]+"-"+chromosomes[-1]
//        comb_info = "${projectName}_${chrms}.imputed_info"
//        well_out = "${comb_info}_well_imputed"
//        acc_out = "${comb_info}_accuracy"
//        infos = chrm_infos.join(',')
//        """
//        python2.7 ${params.homedir}/scripts/report.py \
//            --infoFiles ${infos} \
//            --outWell_imputed ${well_out} \
//            --outSNP_acc ${acc_out} \
//            --infoCutoff ${params.impute_info_cutoff}
//        """
//}
//
//
//"""
//Report 1: Well imputed
//"""
//info_Well.into{ info_Well; info_Well_1}
//process report_well_imputed {
//    tag "report_wellImputed_${projectName}_${chrms}"
//    memory { 2.GB * task.attempt }
//    publishDir "${params.output_dir}/REPORTS/${projectName}", overwrite: true, mode:'copy'
//    input:
//        set val(projectName), file(well_in) from info_Well_1
//    output:
//        set val(projectName), file(well_report_out) into report_well_imputed
//    script:
//        chrms = chromosomes[0]+"-"+chromosomes[-1]
//        comb_info = "${projectName}_${chrms}.imputed_info"
//        well_report_out = "${comb_info}_report_well_imputed.tsv"
//        """
//        python2.7 -u ${params.homedir}/scripts/report.py \
//            --inWell_imputed ${well_in} \
//            --outWell_imputed ${well_report_out} \
//            --infoCutoff ${params.impute_info_cutoff}
//        """
//}
//
//
//"""
//Repor 2: Accuracy
//"""
//info_Acc.into{ info_Acc; info_Acc_2}
//process report_SNP_acc {
//    tag "report_SNP_acc_${projectName}_${chrms}"
//    memory { 2.GB * task.attempt }
//    publishDir "${params.output_dir}/REPORTS/${projectName}", overwrite: true, mode:'copy'
//    input:
//        set val(projectName), file(acc_in) from info_Acc_2
//    output:
//        set val(projectName), file(SNP_acc_out) into report_SNP_acc
//    script:
//        chrms = chromosomes[0]+"-"+chromosomes[-1]
//        comb_info = "${projectName}_${chrms}.imputed_info"
//        SNP_acc_out = "${comb_info}_report_SNP_acc.tsv"
//        """
//        python2.7 -u ${params.homedir}/scripts/report.py \
//            --inSNP_acc ${acc_in} \
//            --outSNP_acc ${SNP_acc_out} \
//            --infoCutoff ${params.impute_info_cutoff}
//        """
//}
