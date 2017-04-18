#!/usr/bin/env nextflow
// -*- mode: groovy -*-

//params.study =  "$baseDir/data/example.chr22.{map,study.gens}"
params.map =  "$baseDir/data/example.map"
params.ped =  "$baseDir/data/example.ped"

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
params.reference_hap = {chr -> "${params.refdir}/HRC.r1-1.EGA.GRCh37.chr${chr}.hap.gz"}
params.reference_legend = {chr -> "${params.refdir}/HRC.r1-1.EGA.GRCh37.chr${chr}.legend.gz"}
params.reference_map = {chr -> "${params.refdir}/genetic_map_chr${chr}_combined_b37.txt"}
params.reference_sample = {chr -> "${params.refdir}/HRC.r1-1.EGA.GRCh37.chr${chr}.samples"}

// plink may also be called plink1
params.plink="plink1"

// these are the chromosome segments to iterate over

// This takes the start (params.begin) and stop (params.end), makes a
// range operator, then steps through that range by the window size
// (params.range) to create a list, which is turned into the
// start/stop positions by collect
chr_segments = params.positions.collect{ entry ->
  (entry.value[0] .. entry.value[1]).step((int)params.range+1).
  collect({[entry.key,it,it+params.range > entry.value[1]?entry.value[1]:it+params.range]})
}

// this channel contains the ped and map file
ped_and_maps = Channel.fromPath(params.ped).\
merge(Channel.fromPath(params.map)){o,e->[o,e]}

ped_maps_per_chr = Channel.create();

/* identify chromosomes and start/stop positions per chromosome */
process identifyChromosomes {
  input:
  set file('full.ped'),file('full.map') from ped_and_maps;

  output:
  set file('full.ped'),file('full.map'),file('chr_min_max') into chromosomes,chromosomes2;

  // Read through the map file, and output the min and max for each chromosome
  shell:
  '''
  perl -n -e 'my ($chr,$rs,$m,$pos) = split /\\s+/;
  if ($chr !~ /^(Y|X|0|23|24)$/) {
   $chrs{$chr}{start} = $pos if not defined $chrs{$chr}{start} or $chrs{$chr}{start} > $pos;
   $chrs{$chr}{end} = $pos if not defined $chrs{$chr}{end} or $chrs{$chr}{end} < $pos;
  }
  END { print map { qq($_\\t$chrs{$_}{start}\\t$chrs{$_}{end}\\n)} keys %chrs; }' < \
    full.map > chr_min_max
  '''
}

process duplicatePedMapByChr {
  input:
  set ped,map,chroms from chromosomes;

  executor 'local'
  
  exec:
  for (chrminmax in chroms.readLines()) {
    def cmt = chrminmax.tokenize("\t")
    def chr = cmt[0];
    def start = cmt[1].toInteger();
    start = start - params.range;
    if (start < 0) {
      start = 0;
    }
    def end = cmt[2].toInteger();
    end = end + params.range;
    ped_maps_per_chr.bind(tuple(cmt[0],start,end,ped,map,
                                file(params.reference_hap(chr)),
                                file(params.reference_legend(chr)),     \
                                file(params.reference_map(chr)),
                                file(params.reference_sample(chr))));
    // ped_maps_per_chr.bind(set (chr,file('full.ped'),file('full.map')));
  }
}

/* split ped/map by chromosome */
process splitPedMap {
  input:
  set chr,val(start),val(end),file('full.ped'),file('full.map'),\
  file('refhap'),file('reflegend'),file('refmap'),file('refsample')\
  from ped_maps_per_chr;

  output:
  set val(chr),val(start),val(end),file('split.ped'),file('split.map'), \
  file('refhap'),file('reflegend'),file('refmap'),file('refsample')\
  into split_ped_maps

  """
  ${params.plink} --noweb \
           --file full --chr ${chr.replace('chr','')} \
           --recode --out split
  """
}

/* check the study genotypes
 *
 */
process checkGenotypes {
 
    input:
    set val(chr),val(start),val(end),file('split.ped'),file('split.map'), \
    file('refhap.gz'),file('reflegend.gz'),file('refmap'),file('refsample')   \
    from split_ped_maps;
 
    output:
    set val(chr),val(start),val(end),file('split.ped'),file('split.map'), \
    file('split_shapeit_log.snp.strand.exclude'),\
    file('refhap.gz'),file('reflegend.gz'),file('refmap'),file('refsample')\
    into checked_genotypes

    """
    shapeit -check \
        -P split.ped split.map \
        --input-ref refhap.gz reflegend.gz refsample \
        --output-log split_shapeit_log || true
    """
}


/*
 * pre-phase each chromosome
 */
process prePhase {
    cpus 8
  
    input:
    set val(chr),val(start),val(end),file('split.ped'),\
    file('split.map'),file('split.snp.strand.exclude'),                 \
    file('refhap.gz'),file('reflegend.gz'),file('refmap'),file('refsample')\
    from checked_genotypes
 
    output:
    set val(chr),val(start),val(end),file('prephased_gens.haps'),\
    file('prephased_gens.sample'),\
    file('refhap.gz'),file('reflegend.gz'),file('refmap'),file('refsample') \
    into prePhased;

    """
    shapeit -phase \
        -P split.ped split.map \
        --input-ref refhap.gz reflegend.gz refsample \
        --exclude-snp split.snp.strand.exclude \
        -T 8 \
        -O prephased_gens

    """
}


split_prephased = Channel.create();
process splitPrephased {
  input:
  set val(chr),val(start),val(end),val(haps), \
  val(sample),val(refhap),val(reflegend),val(refmap),val(refsample) from prePhased;

  executor 'local'

  exec:
  chr_segments = \
  (start .. end).step((int)params.range+1).\
  collect({[it,it+params.range > end?end:it+params.range]})
  for (chr_seg in chr_segments) {
    split_prephased.bind(tuple(chr,chr_seg[0],chr_seg[1],haps,sample,\
                               refhap,reflegend,refmap,refsample))
  }
}


/*
 * impute using the prephased genotypes
 */
process imputeStudyWithPrephased {
 
    input:
    set val(chr),val(begin),val(end),file('phasedhaps'),file('phasedsample'),\
    file('refhap'),file('reflegend'),file('refmap'),file('refsample')   \
    from split_prephased;
     
    output:
    file('imputed.haps') into imputedHaplotypes

    errorStrategy 'finish'
    
    """
    impute2 \
        -known_haps_g phasedhaps \
        -h <(gzip -dcf refhap ) \
        -l <(gzip -dcf reflegend ) \
        -m <(gzip -dcf refmap ) \
        -int ${begin} ${end} \
        -Ne 15000 \
        -buffer 250 \
        -o imputed.haps || true;
    if ! [ -e imputed.haps ] && [ -e imputed.haps_summary ]; then 
         touch imputed.haps
    fi;
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
