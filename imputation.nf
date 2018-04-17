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
//println "|-- Git info: $workflow.repository - $workflow.revision [$workflow.commitId]"
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


// check if study genotype files exist
if(!file(params.bedFile).exists()) exit 1, "BED file ${params.bedFile} not found. Please check your config file."
if(!file(params.famFile).exists()) exit 1, "FAM file ${params.famFile} not found. Please check your config file."
if(!file(params.bimFile).exists()) exit 1, "BIM file ${params.bimFile} not found. Please check your config file."
if(!file(params.eagle_genetic_map).exists()) exit 1, "MAP file ${params.eagle_genetic_map} not found. Please check your config file."

// check if ref files exist
for (chrm in chromosomes){
    // For Ref 1 // Must exist
    if(!file(sprintf(params.ref_1.hapFile, chrm)).exists()) exit 1, "File ${sprintf(params.ref_1.hapFile, chrm)} not found. Please check your config file."
    if(!file(sprintf(params.ref_1.legendFile, chrm)).exists()) exit 1, "File ${sprintf(params.ref_1.legendFile, chrm)} not found. Please check your config file."
    if(!file(sprintf(params.ref_1.mapFile, chrm)).exists()) exit 1, "File ${sprintf(params.ref_1.mapFile, chrm)} not found. Please check your config file."
    if(!file(sprintf(params.ref_1.sampleFile, chrm)).exists()) exit 1, "File ${sprintf(params.ref_1.sampleFile, chrm)} not found. Please check your config file."
    // If Ref 2 exist
    if ('ref_2' in params.keySet()){
        if(!file(sprintf(params.ref_2.hapFile, chrm)).exists()) exit 1, "File ${sprintf(params.ref_2.hapFile, chrm)} not found. Please check your config file."
        if(!file(sprintf(params.ref_2.legendFile, chrm)).exists()) exit 1, "File ${sprintf(params.ref_2.legendFile, chrm)} not found. Please check your config file."
        if(!file(sprintf(params.ref_2.mapFile, chrm)).exists()) exit 1, "File ${sprintf(params.ref_2.mapFile, chrm)} not found. Please check your config file."
        if(!file(sprintf(params.ref_2.sampleFile, chrm)).exists()) exit 1, "File ${sprintf(params.ref_2.sampleFile, chrm)} not found. Please check your config file."
    }
}

// Create channel for the study data from ped and map files
bed_data = Channel
        .fromPath(params.bedFile)
fam_data = Channel
        .fromPath(params.famFile)
bim_data = Channel
        .fromPath(params.bimFile)
// TODO Be able to run just on a specified chunk


"""
Check user's provided chromosomes vs those in map file
"""
bim_data.into{ bim_data; bim_data_check }
process check_chromosome {
    tag "check_chromosome_${bim_data.baseName}"
    memory { 8.GB * task.attempt }
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
    println "|-- Chromosome(s) ${old_chrs.join(', ')} not in dataset ${params.bimFile} and the pipeline will exit."
    exit 1
}
if (!(not_chrs.isEmpty())){
    println "|-- Chromosome(s) ${not_chrs.join(', ')} not in dataset ${params.bimFile} and will be ignored."
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


"""
Basic QC
"""
bed_data.into{ bed_data; bed_data_plink }
bim_data.into{ bim_data; bim_data_plink }
fam_data.into{ fam_data; fam_data_plink }
generate_chunks.into{ generate_chunks; generate_chunks_plink}
process qc_data {
    tag "qc_data_${data_bed.baseName}"
    memory { 4.GB * task.attempt }
    input:
        file data_bed from bed_data_plink
        file data_bim from bim_data_plink
        file data_fam from fam_data_plink
    output:
        file("${base}_qc.vcf.gz") into qc_data
    script:
        base = "${data_bed.baseName}"
        """
        # 1. Filter
        #   - Missingness per individual > ${params.cut_mind}
        #   - Missingness per marker > ${params.cut_geno}
        #   - Hardy-Weinberg equilibrium < ${params.cut_hwe}
        ${params.plink} --bfile ${base} \
            --mind ${params.cut_mind} \
            --geno ${params.cut_geno} \
            --hwe ${params.cut_hwe} \
            --allow-no-sex --recode \
            --out ${base}_clean_mind
        # 2. Remove duplicate sample (remove first)
        ${params.plink} --file ${base}_clean_mind \
            --allow-no-sex --list-duplicate-vars ids-only suppress-first \
            --out ${base}_clean_mind
        ${params.plink} --file ${base}_clean_mind \
            --allow-no-sex \
            --set-missing-var-ids @:# \
            --keep-allele-order \
            --double-id \
            --exclude ${base}_clean_mind.dupvar \
            --recode vcf --out ${base}_qc
        bgzip -f ${base}_qc.vcf
        rm -f ${base}.* \
            ${base}_clean_mind.* \
            ${base}_clean_mind.dupvar
    """
}

"""
Generate sample file from fam
"""
qc_data.into{ qc_data; qc_data_2}
process sample_from_fam {
    tag "sample_from_fam_${base}"
    memory { 2.GB * task.attempt }
    input:
        file(vcf) from qc_data_2
    output:
        file("${base}.sample") into sample_from_fam
    script:
        base = "${file(vcf.baseName).baseName}"
        """
        ${params.plink} --vcf ${vcf} \
            --allow-no-sex \
            --double-id \
            --make-bed --out ${base}
        echo "ID_1\tID_2\tmissing\n0\t0\t0" > ${base}.sample
        awk '{print \$2"\t"\$2"\t"0}' ${base}.fam  >> ${base}.sample
        rm -f ${base}.fam ${base}.bim ${base}.bed ${base}.log
        """
}


"""
Split vcf data in chromosomes
"""
qc_data.into{ qc_data; qc_data_1 }
generate_chunks.into{ generate_chunks; generate_chunks_plink}
process qc_data_to_chrm {
    tag "vcf2chrm_${chromosome}_${base}"
    memory { 4.GB * task.attempt }
    input:
        each chromosome from chromosomes
        file(data_vcf) from qc_data_1
    output:
        set val(chromosome), file("${base}.chr${chromosome}.vcf.gz") into qc_data_to_chrm
    script:
        base = "${file(data_vcf.baseName).baseName}"
        """
        # 1. Extract single chromosomes
        ${params.plink} --vcf ${data_vcf} \
            --chr ${chromosome} --allow-no-sex \
            --keep-allele-order \
            --double-id \
            --recode vcf --out ${base}.chr${chromosome}
        bgzip -f ${base}.chr${chromosome}.vcf
        """
}


qc_data_to_chrm.into{ qc_data_to_chrm; qc_data_to_chrm_1}
generate_chunks.into{ generate_chunks; generate_chunks_1; generate_chunks_all }
all_chunks = generate_chunks_1.toSortedList().val

def transform_chunk = { chrm, vcfFile ->
    chunks_datas = []
    all_chunks.each{ chunk_file ->
        chunks = file(chunk_file).text.split()
        chunks.each{ chunk_data ->
            data = chunk_data.split(',')
            chrm_ = data[0]
            chunk_start = data[1]
            chunk_end = data[2]
            if (chrm ==chrm_){
                chunks_datas << [chrm, chunk_start, chunk_end, file(vcfFile)]
            }
        }
    }
    return chunks_datas
}
chunk_vcfFile = qc_data_to_chrm_1
        .flatMap{ it -> transform_chunk(it) }

chunk_vcfFile.into{ chunk_vcfFile; chunk_vcfFile_1}
process chunk_vcf_data {
    tag "chunk_vcf_${chromosome}_${chunk_start}_${chunk_end}"
    memory { 4.GB * task.attempt }
    input:
        set chromosome, chunk_start, chunk_end, file(vcfFile) from chunk_vcfFile_1
    output:
        set chromosome, chunk_start, chunk_end, file(vcfFile_out) into chunk_vcf_data
    script:
        vcfFile_out = "${file(vcfFile.baseName).baseName}_${chunk_start}_${chunk_end}.vcf.gz"
        buffer_size = params.chunk_size.toInteger()/2
        """
        vcftools \
            --gzvcf ${vcfFile} \
            --chr ${chromosome} \
            --from-bp ${chunk_start.toInteger() - buffer_size.toInteger()} --to-bp ${chunk_end.toInteger() + buffer_size.toInteger()} \
            --recode --stdout \
            | bgzip -c > ${vcfFile_out}
        """
}

"""
Pre-phase each chromosome using Eagle
"""
chunk_vcf_data.into{ chunk_vcf_data; chunk_vcf_data_1 }
process phase_data {
    tag "prephase_${chromosome}_${chunk_start}_${chunk_end}"
    memory { 8.GB * task.attempt }
    input:
        set chromosome, chunk_start, chunk_end, file(vcfFile) from chunk_vcf_data
    output:
        set chromosome, chunk_start, chunk_end, file("${file_out}.haps.gz"), file("${file_out}.sample"), file(vcfFile) into phase_data
    script:
        file_out = "${file(vcfFile.baseName).baseName}.phased"
        """
        nblines=\$(zcat ${vcfFile} | grep -v '^#' | wc -l)
        if (( \$nblines > 0 ))
        then
            ${params.plink} \
                --vcf ${vcfFile} \
                --set-missing-var-ids @:# \
                --double-id --recode --make-bed \
                --out ${file_out} || true
            eagle \
                --bfile=${file_out} \
                --geneticMapFile=${params.eagle_genetic_map} \
                --chrom=${chromosome} \
                --genoErrProb 0.003 --pbwtOnly \
                --allowRefAltSwap \
                --maxMissingPerIndiv 1 \
                --outPrefix=${file_out} 2>&1 | tee ${file_out}.log
            if [ ! -f "${file_out}.haps.gz" ]; then
                touch ${file_out}.haps && bgzip -f ${file_out}.haps
                touch ${file_out}.sample
            fi
            rm -f ${file_out}.bim ${file_out}.bed ${file_out}.fam ${file_out}.ped ${file_out}.map ${file_out}.log
        else
            touch ${file_out}.haps && bgzip -f ${file_out}.haps
            touch ${file_out}.sample
        fi
        """
}


"""
Impute with one (ref1) or two reference panels (ref1 and ref2) according to the configuration files
"""
phase_data.into{ phase_data; phase_data_1; phase_data_2 }
def combine_chunk_data_with_ref = { chromosome, chunk_start, chunk_end, study_haps, study_sample, study_vcf ->
    res = []
    ref_hapFile     = file( sprintf(params.ref_1.hapFile, chromosome) )
    ref_legendFile  = file( sprintf(params.ref_1.legendFile, chromosome) )
    ref_mapFile     = file( sprintf(params.ref_1.mapFile, chromosome) )
    ref_sampleFile  = file( params.ref_1.sampleFile )
    res << [chromosome, chunk_start, chunk_end, study_haps, study_sample, study_vcf, ref_hapFile, ref_legendFile, ref_mapFile, ref_sampleFile]
    return res
}
def combine_chunk_data_with_2refs = { chrm, chunk_start, chunk_end, study_haps, study_sample ->
    res = []
    ref_1_hapFile = file(sprintf(params.ref_1.hapFile, chrm))
    ref_1_legendFile = file(sprintf(params.ref_1.legendFile, chrm))
    ref_1_mapFile = file(sprintf(params.ref_1.mapFile, chrm))
    ref_1_sampleFile = file(params.ref_1.sampleFile)
    ref_2_hapFile = file(sprintf(params.ref_2.hapFile, chrm))
    ref_2_legendFile = file(sprintf(params.ref_2.legendFile, chrm))
    ref_2_mapFile = file(sprintf(params.ref_2.mapFile, chrm))
    ref_2_sampleFile = file(params.ref_2.sampleFile)
    res << [chrm, chunk_start, chunk_end, study_haps, study_sample, study_vcf, ref_1_hapFile, ref_1_legendFile, ref_2_mapFile, ref_1_sampleFile, ref_2_hapFile, ref_2_legendFile, ref_2_sampleFile]
    return res
}
if ( params.ref_1.name != null){
    if ( !('ref_2' in params.keySet()) || params.ref_2.name.length() == 0){
        '''
        Impute using the prephased genotypes with impute2 with 1 reference panel
        '''
        chunk_prephased = phase_data_1.
                flatMap { it -> combine_chunk_data_with_ref(it) }
    }
    else{
        cross_refs_data = phase_data_1.
                flatMap { it -> combine_chunk_data_with_2refs(it) }
        process cross_impute_2refs {
            tag "cross_impute_2refs_chr${chromosome}_${chunkStart}-${chunkEnd}_${params.ref_1.name}_${params.ref_2.name}"
            memory { 15.GB * task.attempt }
            time { 10.h * task.attempt }
            input:
                set val(chromosome), val(chunkStart), val(chunkEnd), file(study_haps), file(study_sample), file(study_vcf), file(ref_1_hapFile), file(ref_1_legendFile), file(ref_2_mapFile), file(ref_1_sampleFile), file(ref_2_hapFile), file(ref_2_legendFile), file(ref_2_sampleFile) from cross_refs_data
            output:
                set val(chromosome), val(chunkStart), val(chunkEnd), file(study_haps), file(study_sample), file(haps_file), file(legend_file), file(ref_2_mapFile), file(sample_file) into cross_impute_2refs
            shell:
                outfile = "${params.ref_1.name}_${params.ref_2.name}_chr${chromosome}_${chunkStart}-${chunkEnd}.impute2"
                haps_file = "${outfile}.hap.gz"
                legend_file = "${outfile}.legend.gz"
                sample_file = "${outfile}.sample"
                buffer = (params.chunk_size.toInteger()/2000).toInteger() // in kb
                """
                gunzip -c ${ref_1_hapFile} > ${ref_1_hapFile.baseName}
                gunzip -c ${ref_2_hapFile} > ${ref_2_hapFile.baseName}
                gunzip -c ${ref_1_legendFile} > ${ref_1_legendFile.baseName}
                gunzip -c ${ref_2_legendFile} > ${ref_2_legendFile.baseName}
                impute2 \
                    -merge_ref_panels \
                    -m ${ref_2_mapFile} \
                    -h ${ref_1_hapFile.baseName} ${ref_2_hapFile.baseName} \
                    -l ${ref_1_legendFile.baseName} ${ref_2_legendFile.baseName} \
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
                head -n 1 ${ref_1_sampleFile} > ${sample_file}
                tail -q -n +2 ${ref_1_sampleFile} ${ref_2_sampleFile} >> ${sample_file}
                rm -f ${ref_1_hapFile.baseName} ${ref_2_hapFile.baseName} ${ref_1_legendFile.baseName} ${ref_2_legendFile.baseName}
                ## Sometimes there are no (type2) SNP's in a region
                if [ ! -f "${haps_file}" ]; then
                    nblines=\$(zcat ${study_vcf} | grep -v '^#' | wc -l)
                    if (( \$nblines >= 1 )); then
                        bcftools convert -gensample \
                            ${study_vcf} \
                            -Oz -o ${file(study_vcf.baseName).baseName}.gen.gz
                        cp ${file(study_vcf.baseName).baseName}.gen.gz
                    else
                        if grep 'ERROR: There are no type 2 SNPs after applying the command-line settings for this run' ${outfile}_summary || \
                            grep 'Your current command-line settings imply that there will not be any SNPs in the output file, so IMPUTE2 will not perform any analysis or print output files.' ${outfile}_summary || \
                            grep 'There are no SNPs in the imputation interval' ${outfile}_summary; then
                            touch ${outfile}.hap
                            bgzip -f ${outfile}.hap
                            touch ${outfile}.legend
                            bgzip -f ${outfile}.legend
                            touch ${outfile}.sample
                        fi
                    fi
                fi
                """
        }
        cross_impute_2refs.into{ cross_impute_2refs; chunk_prephased}
    }
    process impute {
        tag "imp_${chromosome}_${chunkStart}-${chunkEnd}"
        memory { 8.GB * task.attempt }
        input:
            set val(chromosome), val(chunkStart), val(chunkEnd), file(study_haps), file(study_sample), file(study_vcf), file(ref_hapFile), file(ref_legendFile), file(ref_mapFile), file(ref_sampleFile) from chunk_prephased
        output:
            set val(chromosome), val(chunkStart), val(chunkEnd), file("${outfile}.imputed.gz"), file("${outfile}.imputed_info"), file("${outfile}.imputed_summary"), file(study_vcf), file(study_sample) into impute_all
        shell:
            outfile = "${file(study_haps.baseName).baseName}"
            """
            impute2 \
                -use_prephased_g \
                -known_haps_g ${study_haps} \
                -h ${ref_hapFile}  \
                -l ${ref_legendFile}  \
                -m ${ref_mapFile}  \
                -int ${chunkStart} ${chunkEnd} \
                -Ne 15000 \
                -buffer 250 \
                -align_by_maf_g \
                -o_gz \
                -o ${outfile}.imputed || true

            ## Sometimes there are no (type2) SNP's in a region
            if [ ! -f "${outfile}.imputed.gz" ]; then
                if grep 'ERROR: There are no type 2 SNPs after applying the command-line settings for this run' ${outfile}.imputed_summary || \
                    grep 'Your current command-line settings imply that there will not be any SNPs in the output file, so IMPUTE2 will not perform any analysis or print output files.' ${outfile}.imputed_summary || \
                    grep 'There are no SNPs in the imputation interval' ${outfile}.imputed_summary;
                    then
                    nblines=\$(zcat ${study_vcf} | grep -v '^#' | wc -l)
                    if (( \$nblines >= 1 )); then
                        bcftools convert -gensample \
                            ${study_vcf} \
                            --chrom \
                            -Oz -o ${file(study_vcf.baseName).baseName}
                        mv ${file(study_vcf.baseName).baseName}.gen.gz ${outfile}.imputed.gz
                        touch ${outfile}.imputed_info
                    else
                        touch ${outfile}.imputed
                        bgzip -f ${outfile}.imputed
                        touch ${outfile}.imputed_info
                    fi
                fi
            fi
            """
    }
}


'''
Combine output
'''
impute_all.into{impute_all; impute__imputeCombine_cha}

// Create a dataflow instance of all impute results
imputeCombine_impute = []
imputeCombine_info = []
impute__imputeCombine_cha_list = impute__imputeCombine_cha.toSortedList().val
impute__imputeCombine_cha_list.each{ chromosome, chunkStart, chunkEnd, impute, info, summary, study_vcf, study_sample ->
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

"""
Combine impute chunks to chromosomes
"""
process imputeCombine {
    tag "impComb_chr${chromosome}"
    memory { 2.GB * task.attempt }
    input:
        set chromosome, file(imputed_files) from imputeCombine_impute_cha
    output:
        set chromosome, file(comb_impute) into imputeCombine
    script:
        comb_impute = "${file(params.bedFile).getBaseName()}_chr${chromosome}.imputed.gz"
        """
        zcat ${imputed_files.join(' ')} | bgzip -c > ${comb_impute}
        """
}


"""
Combine impute info chunks to chromosomes
"""
process infoCombine {
    tag "infoComb_chr${chromo}"
    memory { 2.GB * task.attempt }
    publishDir "${params.impute_result}/INFOS", overwrite: true, mode:'copy'
    input:
        set chromo, file(info_files) from imputeCombine_info_cha
    output:
        set chromo, file(comb_info) into infoCombine
    script:
        comb_info = "${file(params.bedFile).getBaseName()}_chr${chromo}.imputed_info"
        """
        echo "snp_id rs_id position a0 a1 exp_freq_a1 info certainty type info_type0 concord_type0 r2_type0" > ${comb_info}
        tail -q -n +2 ${info_files.join(' ')} >> ${comb_info}
        """
}


"""
Convert chromosome impute files to VCF files
"""
imputeCombine.into { imputeCombine; imputeCombine_1 }
process imputeToVCF {
    tag "toVCF_chr${chromosome}"
    memory { 2.GB * task.attempt }
    publishDir "${params.impute_result}/FINAL_VCFS", overwrite: true, mode:'copy'
    input:
        set chromosome, file(chromosome_imputed_gz) from imputeCombine_1
        file study_sample from sample_from_fam
    output:
        set chromosome, file("${base}.vcf.gz"), file("${base}.vcf.gz.tbi"), file(study_sample) into imputeToVCF
    script:
        base = chromosome_imputed_gz.baseName
        """
        gunzip -c ${chromosome_imputed_gz} > ${base}.haps
        ${params.plink} \
            --gen ${base}.haps \
            --sample ${study_sample} \
            --oxford-single-chr ${chromosome} \
            --hard-call-threshold 0.1 \
            --mind ${params.cut_mind} \
            --geno ${params.cut_geno} \
            --hwe ${params.cut_hwe} \
            --recode vcf --out ${base}.temp || true
        bgzip -f ${base}.temp.vcf
        bcftools sort ${base}.temp.vcf.gz -Oz -o ${base}.vcf.gz
        bcftools index --tbi ${base}.vcf.gz
        rm -f ${base}.haps *temp*
        """
}


"""
Generating report
"""
infoCombine.into { infoCombine; infoCombine_1 }
infoCombine_list = infoCombine_1.toSortedList().val
chrm_infos = []
infoCombine_list.each{ chrm, comb_info ->
    chrm_infos << chrm+"=="+comb_info
}
process filter_info {
    tag "filter_${projectName}_${chrms}"
    memory { 2.GB * task.attempt }
    publishDir "${params.output_dir}/INFOS/${projectName}", overwrite: true, mode:'copy'
    input:
        val(projectName) from params.project_name
    output:
        set val(projectName), file(well_out) into info_Well
        set val(projectName), file(acc_out) into info_Acc
    script:
        chrms = chromosomes[0]+"-"+chromosomes[-1]
        comb_info = "${projectName}_${chrms}.imputed_info"
        well_out = "${comb_info}_well_imputed"
        acc_out = "${comb_info}_accuracy"
        infos = chrm_infos.join(',')
        """
        python2.7 ${params.homedir}/scripts/report.py \
            --infoFiles ${infos} \
            --outWell_imputed ${well_out} \
            --outSNP_acc ${acc_out} \
            --infoCutoff ${params.impute_info_cutoff}
        """
}


"""
Report 1: Report Well imputed (info > 0.8)
"""
info_Well.into{ info_Well; info_Well_1}
process report_well_imputed {
    tag "report_wellImputed_${projectName}_${chrms}"
    memory { 2.GB * task.attempt }
    publishDir "${params.output_dir}/REPORTS/${projectName}", overwrite: true, mode:'copy'
    input:
        set val(projectName), file(well_in) from info_Well_1
    output:
        set val(projectName), file(well_report_out) into report_well_imputed
    script:
        chrms = chromosomes[0]+"-"+chromosomes[-1]
        comb_info = "${projectName}_${chrms}.imputed_info"
        well_report_out = "${comb_info}_report_well_imputed.tsv"
        """
        python2.7 -u ${params.homedir}/scripts/report.py \
            --inWell_imputed ${well_in} \
            --outWell_imputed ${well_report_out} \
            --infoCutoff ${params.impute_info_cutoff}
        """
}


"""
Report 2: Accuracy (Concordance)
"""
info_Acc.into{ info_Acc; info_Acc_2}
process report_SNP_acc {
    tag "report_SNP_acc_${projectName}_${chrms}"
    memory { 2.GB * task.attempt }
    publishDir "${params.output_dir}/REPORTS/${projectName}", overwrite: true, mode:'copy'
    input:
        set val(projectName), file(acc_in) from info_Acc_2
    output:
        set val(projectName), file(SNP_acc_out) into report_SNP_acc
    script:
        chrms = chromosomes[0]+"-"+chromosomes[-1]
        comb_info = "${projectName}_${chrms}.imputed_info"
        SNP_acc_out = "${comb_info}_report_SNP_acc.tsv"
        """
        python2.7 -u ${params.homedir}/scripts/report.py \
            --inSNP_acc ${acc_in} \
            --outSNP_acc ${SNP_acc_out} \
            --infoCutoff ${params.impute_info_cutoff}
        """
}
