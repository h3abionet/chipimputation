#!/usr/bin/env nextflow
 
//params.study =  "$baseDir/data/example.chr22.{map,study.gens}"
params.study =  "$baseDir/data/{*}.{map,study.gens}"

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
    echo ${phasedchunkname}
    echo ${phasedchunkname}

    """

}
 
/*
 * print the channel content
 */
result.subscribe { println it }