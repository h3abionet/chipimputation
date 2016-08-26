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

// This takes the start (params.begin) and stop (params.end), makes a
// range operator, then steps through that range by the window size
// (params.range) to create a list, which is turned into the
// start/stop positions by collect
chr_segments = (params.begin .. params.end).step((int)params.range+1).collect({[it,it+params.range > params.end?params.end:it+params.range]})

// We then take all of the file pairs from params.study and spread
// them over the chromosome segments into the input_study channel

Channel
.fromFilePairs( params.study).spread(chr_segments)
    .set { input_study }

/*
 * pre-phase the study genotypes
 */
process prePhase {
 
    input:
    set val(chunkname), file(map),file(genotypes),val(begin),val(end) from input_study
 
    output:
    set file("prephased.impute2"),val(begin),val(end) into prePhased

 
    """
    echo "Running prePhase on..." 
    echo ${chunkname}.map
    echo ${chunkname}.study.gens

    impute2 -prephase_g \
        -m ${map} \
        -g ${genotypes} \
        -int ${begin} ${end} \
        -Ne ${params.populationsize} \
        -o prephased.impute2

    """
}
 
/*
 * impute using the prephased genotypes
 */
process imputeStudyWithPrephased {
 
    input:
    set file(phasedchunkname),val(begin),val(end) from prePhased
     
    output:
    file('imputed.haps') into imputedHaplotypes
 
    """
    impute2 \
        -known_haps_g ${phasedchunkname} \
        -h ${params.reference_hap} \
        -l ${params.reference_legend} \
        -m ${params.reference_map} \
        -int ${begin} ${end} \
        -Ne 15000 \
        -buffer 250 \
        -o imputed.haps

    """

}

// combine the imputed segments

process imputeCombine {
  publishDir "./results/imputation"
  input:
    file haplotype_files from imputedHaplotypes.toList()
  output:
  file ('imputed_haplotypes') into result

  script:
  """
  rm -f imputed_haplotypes;
  for file in ${haplotype_files}; do
     cat \$file >> imputed_haplotypes
  done;
  """
}


/*
 * print the channel content
 */
