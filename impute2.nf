#!/usr/bin/env nextflow
// -*- mode: groovy -*-

//params.study =  "$baseDir/data/example.chr22.{map,study.gens}"
params.study_map =  "$baseDir/data/example.map"
params.study_ped =  "$baseDir/data/example.ped"

params.refdir = "/srv/imputation/refdata/1000GP_Phase3"

// currently we only support working on a single chromosome at a time,
// but we should eventually support more than that
params.positions = ["chr20":
                    [(long)20e6, // this is the starting position on this chr
                     (long)21e6 // this is the ending position
                    ]
                   ]
// this is the size of imputation blocks that we should be analyzing
params.range = (long)5e5

params.populationsize = 2000
params.reference_hap = {chr -> "${params.refdir}/1000GP_Phase3_${chr}.hap.gz"}
params.reference_legend = {chr -> "${params.refdir}/1000GP_Phase3_${chr}.legend.gz"}
params.reference_map = {chr -> "${params.refdir}/genetic_map_${chr}_combined_b37.txt"}

// these are the chromosome segments to iterate over

// This takes the start (params.begin) and stop (params.end), makes a
// range operator, then steps through that range by the window size
// (params.range) to create a list, which is turned into the
// start/stop positions by collect
chr_segments = params.positions.collect{ entry ->
  (entry.value[0] .. entry.value[1]).step((int)params.range+1).
  collect({[entry.key,it,it+params.range > entry.value[1]?entry.value[1]:it+params.range]})
}

chromosomes = params.positions.collect{entry -> entry.key}

ped_maps_per_chr = Channel.fromPath(params.study_ped).\
merge(Channel.fromPath(params.study_map)){o,e->[o,e]}.spread(chromosomes)

ped_maps_per_chr_debug = Channel.fromPath(params.study_ped).\
merge(Channel.fromPath(params.study_map)){o,e->[o,e]}.spread(chromosomes)

ped_maps_per_chr_debug.subscribe {println "Got: $it"}

//chr_segments.each{println "Item: $it"};

// We then take all of the file pairs from params.study and spread
// them over the chromosome segments into the input_study channel

// input_study = Channel.fromPath( params.study).spread(chr_segments)

// input_study_debug = Channel.fromPath( params.study).spread(chr_segments)
// input_study_debug.subscribe { println "Got: $it" }


/* split ped/map by chromosome */
process splitPedMap {
  input:
  set file('full.ped'),file('full.map'),chr from ped_maps_per_chr

  output:
  set val(chr),file('split.ped'),file('split.map') into split_ped_maps

  """
  plink1 --noweb --file full --chr ${chr.replace('chr','')} --recode --out split
  """
}

/* check the study genotypes
 *
 */
process checkGenotypes {
 
    input:
    set chr,file('split.ped'),file('split.map') from split_ped_maps
 
    output:
    set val(chr),file('split.ped'),file('split.map'),file('split_shapeit_log.snp.strand.exclude') \
    into checked_genotypes

    """

    shapeit -check \
        -P split.ped split.map \
        --input-ref ${params.reference_hap(chr)} \
                    ${params.reference_legend(chr)} \
                    ${params.reference_sample(chr)} \
        --output-log split_shapeit_log
    """
}


/*
 * pre-phase each chromosome
 */
process prePhase {
 
    input:
    set chr,file('split.ped'),file('split.map'),file('split.snp.strand.exclude') \
    from checked_genotypes
 
    output:
    set val(chr),file('prephased_gens'),val(begin),val(end) into prePhased

    """
    echo "Running prePhase on..." 
    echo ${genotypes}.study.gens
    echo ${begin}-${end}


    shapeit -phase \
        -P split.ped split.map \
        --input-ref ${params.reference_hap(chr)} \
                    ${params.reference_legend(chr)} \
                    ${params.reference_sample(chr)} \
        --exclude-snps split.snp.strand.exclude \
        -O prephased_gens

    """
}
 
/*
 * impute using the prephased genotypes
 */
process imputeStudyWithPrephased {
 
    input:
    set file('phasedchunkname'),begin,end from prePhased
     
    output:
    file('imputed.haps') into imputedHaplotypes
 
    """
    impute2 \
        -known_haps_g phasedchunkname \
        -h <(gzip -dcf ${params.reference_hap}) \
        -l <(gzip -dcf ${params.reference_legend}) \
        -m <(gzip -dcf ${params.reference_map}) \
        -int ${begin} ${end} \
        -Ne 15000 \
        -buffer 250 \
        -o imputed.haps

    """

}

/*
 *  combine the imputed segments
 * 
 * This particular combination is pretty naive, and a better
 * implementation is probably necessary
 */
process imputeCombine {
  publishDir "./results/imputation"
  input:
    file 'haplotype_files' from imputedHaplotypes.toList()
  output:
  file ('imputed_haplotypes') into result

  script:
  """
  rm -f imputed_haplotypes;
  for file in haplotype_files*; do
     cat \$file >> imputed_haplotypes
  done;
  """
}


/*
 * print the channel content
 */
