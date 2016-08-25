#!/usr/bin/env nextflow
 
//params.study =  "$baseDir/data/example.chr22.{map,study.gens}"
params.study =  "$baseDir/data/{*}.{map,study.gens}"

params.refdir = "/srv/imputation/refdata/1000GP_Phase3"

// currently we only support working on a single chromosome at a time,
// but we should eventually support more than that
params.chr = "chr22"
// this is the starting position on the chromosome that we should look
// at
params.begin = (long)20e6
// this is the ending position on the chromosome that we should be
// looking at
params.end = (long)21e6
// this is the size of imputation blocks that we should be analyzing
params.range = (long)5e5

params.populationsize = 2000
params.reference_hap = "${params.refdir}/1000GP_Phase3_${params.chr}.hap.gz"
params.reference_legend = "${params.refdir}/1000GP_Phase3_${params.chr}.legend.gz"
params.reference_map = "${params.refdir}/genetic_map_${params.chr}_combined_b37.txt"

// these are the chromosome segments to iterate over
chr_segments = (params.begin .. params.end).step((int)params.range+1).each({[it,it+params.range > params.end?params.end:it+params.range]})

// this is the channel which contains the regions to iterate over
Channel
.from(chr_segments).set { chr_segment }


Channel
    .fromFilePairs( params.study)                              
    .set { input_study }



/*
 * pre-phase the study genotypes
 */
process prePhase {
 
    input:
        set val(chunkname), file(studybase) from input_study
 
    output:
        file "${chunkname}.prephasing.impute2" into prePhased
 
    """
    echo "Running prePhase on..." 
    echo ${chunkname}.map
    echo ${chunkname}.study.gens

    impute2 -prephase_g \
        -m ${chunkname}.map \
        -g ${chunkname}.study.gens \
        -int 20.4e6 20.5e6 \
        -Ne 20000 \
        -o ${chunkname}.prephasing.impute2

    """
}
 
/*
 * impute using the prephased genotypes
 */
process imputeStudyWithPrephased {
 
    input:
        file phasedchunkname from prePhased
     
    output:
        stdout result
 
    """
    echo "Running imputeStudyWithPrephased on..." 
    impute2 \
        -known_haps_g ${INPUT}/gwas_data_chr10_phased.haps \
        -h ${INPUT}/pilot1.jun2010.b36.CEU.chr10.snpfilt.haps \
        -l ${INPUT}/pilot1.jun2010.b36.CEU.chr10.snpfilt.legend \
        -m ${INPUT}/genetic_map_chr10_combined_b36.txt \
        -int 20000000 25000000 \
        -Ne 15000 \
        -buffer 250 \
        -o ${OUTPUT}/gwas_data_chr10_imputed.20-25Mb.gen
    echo ${phasedchunkname}
    echo ${phasedchunkname}

    """

}
 
/*
 * print the channel content
 */
result.subscribe { println it }
