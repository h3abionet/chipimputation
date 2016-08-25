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

    #impute2 -prephase_g \
    #    -m ${chunkname}.map \
    #    -g ${chunkname}.study.gens \
    #    -int 20.4e6 20.5e6 \
    #    -Ne 20000 \
    #    -o ${chunkname}.prephasing.impute2

    touch ${chunkname}.prephasing.impute2

    """
}



/*
 * impute using the prephased genotypes
 */
process imputeStudyWithPrephased {
 
    input:
        file(chunkname) from prePhased
     
    output:
        stdout result
 

    script:
        filename = file(chunkname).name
        basename = filename.replaceAll("\\.prephasing\\.impute2","")
        parentname = file(chunkname).parent
        println "---${chunkname}---\n----filename: ${filename}\n----parentname: ${parentname}" 
        println "-----BASENAME: ${basename}"
    """
    echo "Running imputeStudyWithPrephased on..." 
    echo ${chunkname}
    echo ${filename}
    echo ${parentname}


#impute2 \
# -use_prephased_g \
# -m ${chunkname}.map \

# -h ./refs/example.chr22.1kG.haps \
# -l ./refs/example.chr22.1kG.legend \
# -known_haps_g ${chunkname}.prephasing.im" "g ${chunkname}\\..study.str\\.and \
# -int 20.4e6 20.5e6 \
# -Ne 20000 \
# -o ${chunkname}.one.phae2 \
# -phase

    """

}
 
/*
 * print the channel content
 */
result.subscribe { println it }