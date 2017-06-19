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

println "=================================================="

chromosomes = params.chromosomes.split(',')
chunk_size = params.chunk_size[0..-7]
ref_dir = file(params.ref_dir)

// All POP
def POPS_ALL = []
params.POPS.each { entry->
    POPS_ALL.addAll(entry.value.split(','))
}
println "Project : $workflow.projectDir"
//println "Git info: $workflow.repository - $workflow.revision [$workflow.commitId]"
println "Cmd line: $workflow.commandLine"
println "Populations in use: ${POPS_ALL.join(", ")}"
 println "Chromosomes used: ${chromosomes.join(', ')}"


//// Help functions

def chunk_split(map_file, outfile, chromosome, chunk=params.chunk_size) {
    '''
    Return: chunck files in the output folder
    '''
    res = []
    data = map_file.readLines().collect{[it.split()[0], it.split()[3]]}
    for(dat in data){
        if( dat[0] == chromosome && dat[0] != ''){
            res << dat[1]
        }
    }
    myOutfile = file(outfile)
    myOutfile.text = res.join(' ')
}

def read_chunk(chunk_file) {
    '''
    Return: read a chunk file and return a list of chunks
    '''
    data = new File(chunk_file).readLines()[1..-1]

    return data
}
/////

// create channel for the study data from ped and map files
//ped_data = Channel.fromPath(params.pedFile)
//map_data = Channel.fromPath(params.mapFile)
//ped_map_data = ped_data.merge(map_data){ o,e -> [o,e] }
bed_data = Channel.fromPath(params.bedFile)
fam_data = Channel.fromPath(params.famFile)
bim_data = Channel.fromPath(params.bimFile)
ped_map_data = bed_data
                .merge(fam_data){ o,e -> [o,e] }
                .merge(bim_data){ o,e -> o + [e] }
//                .subscribe{ println it.baseName } [bed, fam, bim]

"""
Create a channel with bed/fam/bim for each chromosome
"""
ped_map_data.into{ped_map_data;ped_map_data_chan}
ped_map_data_all = ped_map_data_chan.flatMap { bed_File, fam_File, bim_File ->
    def results = []
    chromosomes.each { chromosome ->
        results.push( [ chromosome, bed_File, fam_File, bim_File] )
    }
    return results
}


"""
Identify chromosomes and start/stop positions per chromosome
"""
ped_map_data_all.into{ped_map_data_all; ped_map_data__identifyChromosomes; ped_map_data__identifyChromosomes_1}
process generate_chunks {

    echo true
    tag "generate_chunks_${chromosome}"
    cpus { 1 * task.attempt }
    publishDir "${params.work_dir}/data/CHUNKS", overwrite: true, mode:'symlink', pattern: '{*.txt}'

    input:
        set val(chromosome), file(data_bed), file(data_fam), file(bim_data) from ped_map_data__identifyChromosomes

    output:
        set val(chromosome), file(data_bed), file(data_fam), file(bim_data), file("analysis_chunks_${chunk_size}Mb_chr${chromosome}.txt") into identifyChromosomes_all

    // Read through the map file, and output the min and max for each chromosome
    script:
        outfile = "analysis_chunks_${chunk_size}Mb_chr${chromosome}.txt"
        //  chunk_split(bim_data, outfile, chromosome, params.chunk_size)
        """
        #!/usr/bin/python
        POS = [int(dat.strip().split()[3]) for dat in open("${bim_data}").readlines() if dat.strip().split()[3] != '' and dat.strip().split()[0] == str(${chromosome})]
        min_POS = int(str(min(POS))[:-1]+'1')
        max_POS = max(POS)
        chunk = int(${params.chunk_size})
        out = open('${outfile}', "wt")
        out.writelines("start" +" "+"end"+"\\n")
        for pos in list(range(min_POS, max_POS+chunk, chunk)):
            start_ = pos
            end_ = start_ + chunk-1
            out.writelines(str(start_) +" "+str(end_)+"\\n")
        out.close()
        """
}


"""
Plink bed/fam/bim per chromosomes
"""
ped_map_data_all.into{ped_map_data_all; ped_map_data__plink_to_chrm}
process plink_to_chrm {

    tag "plink2chrm_${chromosome}"
    cpus { 2 * task.attempt }
    publishDir "${params.work_dir}/data/qc/${chromosome}", overwrite: true, mode:'symlink'

    input:
        set val(chromosome), file(data_bed), file(data_fam), file(bim_data) from ped_map_data__plink_to_chrm

    output:
        set val(chromosome), file("${data_bed.baseName}.chr${chromosome}_clean.ped"), file("${data_bed.baseName}.chr${chromosome}_clean.map") into plink_to_chrm_all

    script:
        /*
        2- Exclude samples with missing more than 10% of genotype calls
        3- Two versions of you dataset, one with minor allele frequencies (MAFs) greater than 5% and one with MAFs less than 5%
        */
        """
        plink2 --bfile ${data_bed.baseName} \
            --chr ${chromosome} --allow-no-sex --make-bed \
            --out ${data_bed.baseName}.chr${chromosome}
        plink2 --bfile ${data_bed.baseName}.chr${chromosome} \
            --mind 0.10 --allow-no-sex --recode \
            --out ${data_bed.baseName}.chr${chromosome}_clean_mind  
        plink2 --file ${data_bed.baseName}.chr${chromosome}_clean_mind \
            --allow-no-sex --list-duplicate-vars ids-only suppress-first \
            --out ${data_bed.baseName}.chr${chromosome}_clean_mind
        plink2 --file ${data_bed.baseName}.chr${chromosome}_clean_mind \
            --allow-no-sex \
            --exclude ${data_bed.baseName}.chr${chromosome}_clean_mind.dupvar \
            --recode --out ${data_bed.baseName}.chr${chromosome}_clean
        rm -f ${data_bed.baseName}.chr${chromosome}.* \
            ${data_bed.baseName}.chr${chromosome}_clean_mind.* \
            ${data_bed.baseName}.chr${chromosome}_clean_mind.dupvar
    """
    /*
    plink2 --bfile ${data_bed.baseName}.chr${chromosome}_clean_mind \
            --maf 0.05 --allow-no-sex --make-bed \
            --out ${data_bed.baseName}.chr${chromosome}_clean_mind_MAF_greater_5
    plink2 --bfile ${data_bed.baseName}.chr${chromosome}_clean_mind \
            --exclude ${data_bed.baseName}.chr${chromosome}_clean_mind_MAF_greater_5.bim \
            --make-bed --allow-no-sex \
            --out ${data_bed.baseName}.chr${chromosome}_clean_mind_MAF_less_5
    plink2 --bfile ${data_bed.baseName}.chr${chromosome}_clean_mind_MAF_greater_5 \
            --geno 0.05 --allow-no-sex \
            --make-bed --out ${data_bed.baseName}.chr${chromosome}_clean_mind_MAF_greater_5_clean
    plink2 --bfile ${data_bed.baseName}.chr${chromosome}_clean_mind_MAF_less_5 \
            --geno 0.01 --allow-no-sex \
            --make-bed --out ${data_bed.baseName}.chr${chromosome}_clean_mind_MAF_less_5_clean
    plink2 --bfile ${data_bed.baseName}.chr${chromosome}_clean_mind_MAF_greater_5_clean \
            --bmerge ${data_bed.baseName}.chr${chromosome}_clean_mind_MAF_less_5_clean \
            --allow-no-sex --recode --out ${data_bed.baseName}.chr${chromosome}_clean_mind_MAF_greater_less_5_clean
    */
}

/* check the study genotypes
 *
 */
plink_to_chrm_all.into{plink_to_chrm__checkGenotypes}
process checkGenotypes {

    tag "checkGenotypes_${chromosome}"
    cpus { 2 * task.attempt }
    validExitStatus 0,1,2
    errorStrategy 'ignore'
    publishDir "${params.work_dir}/data/checkGenotypes_shapeit/${chromosome}", overwrite: true, mode:'symlink'

    input:
        set val(chromosome), file(pedFile), file(mapFile) from plink_to_chrm__checkGenotypes

    output:
        set val(chromosome), \
        file(pedFile), \
        file(mapFile), \
        file("${pedFile.baseName}.alignments.log"), \
        file("${pedFile.baseName}.alignments.snp.strand.exclude") \
        into checkGenotypes_all

    script:
        ref_hapFile = file( sprintf(params.ref_hapFile, chromosome) )
        ref_legendFile = file( sprintf(params.ref_legendFile, chromosome) )
        ref_sampleFile = file( params.ref_sampleFile )

        """
        shapeit -check \
            --input-ped ${pedFile} ${mapFile} \
            --input-ref ${ref_hapFile} ${ref_legendFile} ${ref_sampleFile} \
            --output-log ${pedFile.baseName}.alignments.log
        """
}


"""
Step 3: Pre-phase each chromosome using shapeit
"""

checkGenotypes_all.into{checkGenotypes_all; checkGenotypes_prephase}

process prephase {

    tag "prephase_${chromosome}"
    cpus { 8 * task.attempt }
    memory { 2.GB * task.cpus }
    publishDir "${params.work_dir}/data/prephased/${chromosome}", overwrite: true, mode:'symlink'

    input:
        set val(chromosome), \
        file(data_ped), \
        file(data_map), \
        file(shapeit_check_log), \
        file(shapeit_check_snp_strand_exclude) \
        from checkGenotypes_prephase

    output:
        set val(chromosome), \
        file("${data_ped.baseName}.prephased.haps"), \
        file("${data_ped.baseName}.prephased.sample") \
        into prephase_all, prephase__chunk, prephase__impute

    script:
        ref_hapFile = file( sprintf(params.ref_hapFile, chromosome) )
        ref_legendFile = file( sprintf(params.ref_legendFile, chromosome) )
        ref_sampleFile = file( params.ref_sampleFile )
        """
        shapeit -phase \
            --input-ped ${data_ped} ${data_map} \
            --input-ref ${ref_hapFile} ${ref_legendFile} ${ref_sampleFile} \
            --exclude-snp ${shapeit_check_snp_strand_exclude} \
            -T ${task.cpus} \
            -O ${data_ped.baseName}.prephased
        """
}
prephase_all.subscribe{
    println "|-- Finished for ${it[1][-1]}, ${it[2][-1]}"
}

chunk_prephased = prephase__chunk.flatMap {
    chromosome      = it[0]
    haps            = it[1]
    sample          = it[2]
    ref_hapFile     = file( sprintf(params.ref_hapFile, chromosome) )
    ref_legendFile  = file( sprintf(params.ref_legendFile, chromosome) )
    ref_mapFile     = file( sprintf(params.ref_mapFile, chromosome) )
    ref_sampleFile  = file( params.ref_sampleFile )
    res = []
    chunk_file_data = new File("${params.work_dir}/data/CHUNKS/analysis_chunks_${chunk_size}Mb_chr${chromosome}.txt").readLines()[1..-1]
    for (chunks in chunk_file_data) {
        chunks = chunks.split()
        res << tuple(chromosome, chunks[0], chunks[1], haps, sample, ref_hapFile, ref_legendFile, ref_mapFile, ref_sampleFile)
    }
    res
}
//chunk_prephased.subscribe{ println it }

/*
 * impute using the prephased genotypes
 */
process impute {

    tag "imp_chr${chromosome}_${chunkStart}-${chunkEnd}"
    cpus { 2 * task.attempt }
    memory { 2.GB * task.cpus }
    publishDir "${params.impute_result}/${chromosome}", overwrite: true, mode:'symlink'
    errorStrategy 'finish'

    input:
        set val(chromosome), \
        val(chunkStart), \
        val(chunkEnd), \
        file(haps), \
        file(sample), \
        file(ref_hapFile), \
        file(ref_legendFile), \
        file(ref_mapFile), \
        file(ref_sampleFile) \
        from chunk_prephased

    output:
        set val(chromosome), \
        file("chr${chromosome}_${chunkStart}-${chunkEnd}.imputed.gz"), \
        file("chr${chromosome}_${chunkStart}-${chunkEnd}.imputed_haps.gz"), \
        file("chr${chromosome}_${chunkStart}-${chunkEnd}.imputed_info"), \
        file("chr${chromosome}_${chunkStart}-${chunkEnd}.imputed_summary") \
        into impute_all, impute__sub

    script:
        """
        impute2 \
            -use_prephased_g \
            -known_haps_g ${haps} \
            -h ${ref_hapFile}  \
            -l ${ref_legendFile}  \
            -m ${ref_mapFile}  \
            -int ${chunkStart} ${chunkEnd} \
            -Ne 15000 \
            -buffer 250 \
            -phase \
            -o_gz \
            -o chr${chromosome}_${chunkStart}-${chunkEnd}.imputed || true
        
        #sometimes there are no SNP's in a region
        if [ ! -f "chr${chromosome}_${chunkStart}-${chunkEnd}.imputed.gz" ]; then
            #|| [ ! -f "chr${chromosome}_${chunkStart}-${chunkEnd}.imputed_haps.gz" ]
            touch "chr${chromosome}_${chunkStart}-${chunkEnd}.imputed"
            bgzip -f "chr${chromosome}_${chunkStart}-${chunkEnd}.imputed"
            touch "chr${chromosome}_${chunkStart}-${chunkEnd}.imputed_haps"
            bgzip -f "chr${chromosome}_${chunkStart}-${chunkEnd}.imputed_haps"
            touch "chr${chromosome}_${chunkStart}-${chunkEnd}.imputed_info"
        fi
        """

}
impute__sub.subscribe{
    println "|-- Finished for ${it[1][-1]}, ${it[2][-1]}"
}

workflow.onComplete {
    def subject = 'My pipeline execution'
    def recipient = 'mypandos@gmail.com'

    ['mail', '-s', subject, recipient].execute() << """

    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    Error report: ${workflow.errorReport ?: '-'}
    """
}

