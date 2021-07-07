#!/usr/bin/env nextflow
nextflow.preview.dsl=2

/*
========================================================================================
=                                 h3achipimputation                                    =
========================================================================================
 h3achipimputation imputation functions and processes.
----------------------------------------------------------------------------------------
 @Authors

----------------------------------------------------------------------------------------
 @Homepage / @Documentation
  https://github.com/h3abionet/chipimputation
----------------------------------------------------------------------------------------

----------------------------------------------------------------------------------------
*/

// Show help message
if (params.help){
    helpMessage()
    exit 0
}

process impute_minimac4 {
    tag "imp_${target_name}_${chrm}:${chunk_start}-${chunk_end}_${ref_name}_${tagName}"
    label "bigmem"
    input:
        tuple val(chrm), val(chunk_start), val(chunk_end), val(target_name), file(target_phased_vcf), val(ref_name), file(ref_vcf), file(ref_m3vcf), val(tagName)
    output:
        tuple val(chrm), val(chunk_start), val(chunk_end), val(target_name), val(ref_name), file("${base}_imputed.dose.vcf.gz"), file("${base}_imputed.info"), val(tagName)
    shell:
        base = "${file(target_phased_vcf.baseName).baseName}_${tagName}_${chrm}_${chunk_start}-${chunk_end}"
        """
        minimac4 \
            --refHaps ${ref_m3vcf} \
            --haps ${target_phased_vcf} \
            --format GT,DS \
            --allTypedSites \
            --minRatio ${params.minRatio} \
            --chr ${chrm} --start ${chunk_start} --end ${chunk_end} --window ${params.buffer_size} \
            --prefix ${base}_imputed \
            --cpus ${task.cpus}
        """
}

process impute_minimac4_1 {
    tag "imp_${target_name1}_${chrm1}:${chunk_start1}-${chunk_end1}_${ref_name1}_${tagName1}"
    label "bigmem"
    input:
        tuple val(chrm1), val(chunk_start1), val(chunk_end1), val(target_name1), file(target_phased_vcf1), val(ref_name1), file(ref_vcf1), file(ref_m3vcf1), val(tagName1)
    output:
        tuple val(chrm1), val(chunk_start1), val(chunk_end1), val(target_name1), val(ref_name1), file("${base}_imputed.dose.vcf.gz"), file("${base}_imputed.info"), val(tagName1)
    shell:
        base = "${file(target_phased_vcf1.baseName).baseName}"
        """
        nblines=\$(zcat ${target_phased_vcf1} | grep -v "^#" | wc -l)
        if (( \$nblines > 1 ))
        then
            bcftools annotate  --set-id '%CHROM\\_%POS\\_%REF\\_%ALT' ${target_phased_vcf1} | 
            vcftools --vcf - --keep-INFO-all --max-missing 0.05 --hwe 0.00001 --mac 1 --recode --stdout | bgzip > ${base}.id.vcf.gz
            minimac4 \
                --cpus ${task.cpus} \
                --refHaps ${ref_m3vcf1} \
                --haps ${base}.id.vcf.gz \
                --format GT,DS \
                --allTypedSites \
                --minRatio ${params.minRatio} \
                --chr ${chrm1} --start ${chunk_start1} --end ${chunk_end1} --window ${params.buffer_size} \
                --prefix ${base}_imputed
        else
             touch ${base}_imputed.dose.vcf && bgzip ${base}_imputed.dose.vcf
             touch ${base}_imputed.info
        fi
        """
}

"""
Combine impute chunks to chromosomes
"""
process combineImpute {
    tag "impComb_${target_name}_${ref_name}_${chrm}"
    publishDir "${params.outDir}/imputed/${ref_name}", overwrite: true, mode:'symlink', pattern: '*imputed.vcf.gz'
    label "bigmem"
    
    input:
        tuple val(target_name), val(ref_name), val(chrm), val(imputed_files)
    output:
        tuple val(target_name), val(ref_name), val(chrm), file(comb_impute)// into combineImpute
    script:
        comb_impute = "${target_name}_${ref_name}_chr${chrm}.imputed.vcf.gz"
        """
        bcftools concat ${imputed_files.join(' ')} | \
        bcftools +fill-tags -- -t AC,AN,AF,MAF | \
        bcftools sort -T . -Oz -o ${comb_impute}
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
        tuple val(target_name), val(ref_name), val(chrm), file(info_files) //from infoCombine.values()
    output:
        tuple val(target_name), val(ref_name), val(chrm), file(comb_info) //into combineInfo_frq
    script:
        comb_info = "${target_name}_${ref_name}_chr${chrm}.imputed_info"
        """
        head -n1 ${info_files[0]} > ${comb_info}
        tail -q -n +2 ${info_files.join(' ')} >> ${comb_info}
        """
}
