#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { get_chromosome; fill_tags_vcf; check_chromosome; check_files; check_chromosome_vcf; check_mismatch; no_mismatch ; qc_dupl; split_multi_allelic; filter_min_ac; target_qc as target_qc; target_qc as target_qc1; qc_site_missingness as qc_site_missingness1; qc_site_missingness as qc_site_missingness2; sites_only ; combine_vcfs ; combine_infos; combine_csvs as combine_freqs; combine_vcfs_chrm; } from './modules/qc' 
include { vcf_map_simple; extract_site_from_vcf; generate_chunks_vcf; split_target_to_chunk; vcf_map; vcf_freq; info_freq; fill_tags_VCF; sort_vcf ; get_vcf_sites; extract_pop } from './modules/subset_vcf'
include { minimac4_phasing_eagle } from './modules/phasing'
include { impute_minimac4; impute_minimac4_1; combineImpute; combineInfo; filter_info_by_target } from './modules/impute'
include { filter_info; report_site_by_maf; plot_freq_comparison; report_well_imputed_by_target; plot_performance_target; 
report_accuracy_target; plot_accuracy_target; generate_frequency; plot_r2_SNPpos; plot_r2_SNPcount; plot_hist_r2_SNPcount; plot_MAF_r2; 
average_r2 } from './modules/report'


// Header log info
def intro(){
log.info ""
log.info """
=======================================================
h3achipimputation v${params.version}"
======================================================= """
    def summary = [:]
    summary['Pipeline Name']    = 'h3achipimputation'
    summary['Pipeline version'] = params.version
    // summary['Run Name']         = custom_runName ?: workflow.runName
    // summary['Target datasets']  = params.target_datasets.values().join(', ')
    // summary['Reference panels']  = params.ref_panels.keySet().join(', ')
    summary['Max Memory']       = params.max_memory
    summary['Max CPUs']         = params.max_cpus
    summary['Max Time']         = params.max_time
    summary['Output dir']       = params.outDir
    summary['Working dir']      = workflow.workDir
    summary['Script dir']       = workflow.projectDir
    summary['Current path']     = "$PWD"
    summary['Git info']         = "${workflow.repository} - ${workflow.revision} [${workflow.commitId}]"
    summary['Command line']     = workflow.commandLine
    if(workflow.containerEngine) {
        summary['Container Engine'] = workflow.containerEngine
        summary['Container'] = workflow.container
        summary['Current home'] = "$HOME"
        summary['Current user'] = "$USER"
        summary['Current path'] = "$PWD"
        summary['Working dir'] = workflow.workDir
        summary['Output dir'] = params.outDir
        summary['Script dir'] = workflow.projectDir
        summary['Config Profile'] = workflow.profile
    }

    if(params.email) summary['E-mail Address'] = params.email
    log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
    log.info "======================================================="
    log.info ""
}

workflow preprocess {
    take: datasets

    main:
        //// check if study genotype files exist
        target_datasets = []
        datasets.each { name, vcf ->
            check_files([vcf])
            target_datasets << [name, file(vcf)]
        }
        target_datasets = Channel.from(target_datasets)

        //// check if eagle map file exists
        if(params.eagle_genetic_map) {
            check_files([params.eagle_genetic_map])
        }

        //// check if fasta reference genome and index
        if(params.reference_genome) {
            check_files([params.reference_genome, "${params.reference_genome}.fai"])
        }

        // //// Check chromosome
        // get_chromosome( target_datasets )
        // chromosomes = get_chromosome.out
        //     .map{ dataset, dataset_vcf, map_file -> check_chromosome_vcf(dataset, dataset_vcf, map_file, params.chromosomes) }
        //     .map{ dataset, dataset_vcf, map_file, chrms -> chrms.unique() }
        
        //// check if reference panel files exist
        params.chromosomes.split(',').each{ chrm ->
            params.ref_panels.each{ ref_name, ref_m3vcf, ref_vcf ->
                vcf = sprintf(ref_vcf, chrm)
                m3vcf = sprintf(ref_m3vcf, chrm)
                if(vcf.endsWith("vcf.gz")){
                    vcf_idx = "${vcf}.tbi" 
                }
                else if(vcf.endsWith("bcf")){
                    vcf_idx = "${vcf}.csi" 
                }
                check_files([ m3vcf, vcf, vcf_idx ])
            }
        }

        ////// QC
        check_mismatch(target_datasets.map{ dataset, dataset_vcf -> [ dataset, '', '', '', file(dataset_vcf), file(params.reference_genome) ] } )
        qc_dupl( target_datasets.map{ dataset, dataset_vcf -> [ dataset, '', '', '', file(dataset_vcf) ] } )
        split_multi_allelic(qc_dupl.out)
        fill_tags_vcf(split_multi_allelic.out)
        filter_min_ac(fill_tags_vcf.out.map{ dataset, chrm, start, end, vcf -> [ dataset, chrm, start, end, file(vcf), " --min-ac ${params.min_ac} --max-alleles ${params.max_alleles} --min-alleles ${params.min_alleles} -v snps "  ] })

    emit:
        // chrms = chromosomes
        dataset_qc = filter_min_ac.out
        // get_chromosome = get_chromosome.out
}

workflow subset{
    take: data
    
    main:
        get_chromosome( data.map{ dataset, chrm, start, end, vcf, map_file -> [ dataset, file(vcf) ] } )
        data_chrms = get_chromosome.out
            .map{ dataset, dataset_vcf, map_file -> check_chromosome_vcf(dataset, dataset_vcf, map_file, params.chromosomes) }
            .map{ dataset, dataset_vcf, map_file, chrms -> [ dataset, file(dataset_vcf), file(map_file), chrms.unique().join(',') ] }
        
        generate_chunks_vcf(data_chrms.map{ dataset, vcf, map_file, chrms -> [ dataset, file(vcf), file(map_file), chrms, params.chunk_size ] })
        chunks_datas = generate_chunks_vcf.out.flatMap{ dataset, vcf, chunk_file ->
            datas = []
            chunks = file(chunk_file).text.split()
            chunks.each{ chunk_data ->
                data = chunk_data.split(',')
                chrm = data[0]
                chunk_start = data[1]
                chunk_end = data[2]
                datas << [dataset, chrm, chunk_start, chunk_end, dataset, file(vcf)]
            }
            return datas
        }
        split_target_to_chunk(chunks_datas)

    emit: 
        chunks = split_target_to_chunk.out
}

workflow phasing{
    take: data
    
    main:
        minimac4_phasing_eagle(data)

    emit: 
        chunks_phased = minimac4_phasing_eagle.out
}

workflow impute{
    take: data
    
    main:
        impute_minimac4(data)

    emit: 
        data
        chunks_imputed = impute_minimac4.out
}

workflow report_by_ref{
    take: data

    main:
        /// /// By Refecene panel
        // combine by chrom, dataset, refpanel
        imputeCombine_ref = data
                .groupTuple( by:[1] )
                .map{ datasets, refpanel, vcfs, imputed_vcfs, imputed_infos -> [ refpanel, datasets.join(','), '', imputed_infos.join(',') ] }
        filter_info_by_target( imputeCombine_ref )

        /// change to group_by_maf
        report_well_imputed_by_target( filter_info_by_target.out.map{ target_name, ref_panels, wellInfo, accInfo -> [ target_name, ref_panels, file(wellInfo) ]} )
        
        //// Plot performance all targets by maf for a reference panel
        plot_performance_target( report_well_imputed_by_target.out.map{ target_name, ref_panels, wellInfo, wellInfo_summary -> [ target_name, ref_panels, file(wellInfo), file(wellInfo_summary), 'DATASETS' ]} )

        //// Accuracy/Concordance
        report_accuracy_target( filter_info_by_target.out.map{ target_name, ref_panels, wellInfo, accInfo -> [ target_name, ref_panels, file(accInfo), 'DATASETS' ]} )
        plot_accuracy_target ( report_accuracy_target.out )
    emit:
        data
}

workflow report_by_dataset{
    take: data

    main:
        /// /// By Dataset
        // combine by chrom, dataset, refpanel
        
        imputeCombine_ref = data
                .groupTuple( by:[0] )
                .map{ dataset, refpanels, vcfs, imputed_vcfs, imputed_infos -> [ dataset, refpanels.join(','), '', imputed_infos.join(',') ] }
        filter_info_by_target( imputeCombine_ref )

        ///// Number of well imputed snps
        /// change to group_by_maf
        report_well_imputed_by_target( filter_info_by_target.out.map{ target_name, ref_panels, wellInfo, accInfo -> [ target_name, ref_panels, file(wellInfo) ]} )
        /// Plot performance all targets by maf for a reference panel
        plot_performance_target( report_well_imputed_by_target.out.map{ target_name, ref_panels, wellInfo, wellInfo_summary -> [ target_name, ref_panels, file(wellInfo), file(wellInfo_summary), 'REFERENCE_PANELS' ]} )

        //// Accuracy/Concordance
        report_accuracy_target( filter_info_by_target.out.map{ target_name, ref_panels, wellInfo, accInfo -> [ target_name, ref_panels, file(accInfo), 'REFERENCE_PANELS' ]} )
        plot_accuracy_target ( report_accuracy_target.out )

        // Plot number of imputed SNPs over the mean r2 for all reference panels
        input = imputeCombine_ref
        .map{ dataset, refpanels, chrm, infos -> [dataset, refpanels, infos]}

        // Plot number of imputed SNPs over the mean r2 for all reference panels
        plot_r2_SNPcount(input)

        // Plot histograms of number of imputed SNPs over the mean r2 for all reference panels
        plot_hist_r2_SNPcount(input)

        // Plot MAF of imputed SNPs over r2 for all references
        plot_MAF_r2(input)

    emit:
        data
}

workflow {

    intro()

    //// Data preparation
    preprocess(params.target_datasets)
    
    // //// Data chunking
    subset(preprocess.out.dataset_qc)

    //// Phasing
    phasing_data = Channel.from(params.ref_panels)
        .combine(subset.out.chunks)
        .flatMap{ ref_name, ref_m3vcf, ref_vcf, dataset, chrm, start, end, dataset_, dataset_vcf ->
            vcf = sprintf(ref_vcf, chrm)
            m3vcf = sprintf(ref_m3vcf, chrm)
            if(vcf.endsWith("vcf.gz")){
                vcf_idx = "${vcf}.tbi" 
            }
            else if(vcf.endsWith("bcf")){
                vcf_idx = "${vcf}.csi" 
            }
            return [ [ chrm, ref_name, file(m3vcf), file(vcf), file(vcf_idx), file(params.eagle_genetic_map), start, end, dataset, dataset, file(dataset_vcf)] ]
        }
    phasing(phasing_data)

    //// Imputation
    // phasing.out.chunks_phased.view()
    impute(phasing.out.chunks_phased)
    
    // Reporting
    // impute.out.chunks_imputed.view()
    impute_data = impute.out.chunks_imputed
                .map{chr, fwd, rev, test_data, ref, imputed_vcf, 
                imputed_info, tst_data -> [test_data, ref, imputed_vcf, imputed_info]}
                .combine(params.target_datasets, by:0)
                .map {test_data, ref, imputed_vcf, imputed_info, orig_vcf 
                -> [test_data, ref, orig_vcf, imputed_vcf, imputed_info]}

    // //// Report by Reference
    report_by_ref( impute_data )

    // //// Report by datasets
    report_by_dataset( impute_data )

    // // Generate dataset frequencies
    inp = Channel.fromList(params.ref_panels).map{ref, m3vcf, vcf -> [ref, vcf]}
    input = impute_data
    .map{ target_name, ref_name, vcf, impute_vcf, info ->[  ref_name, target_name, file(impute_vcf)]}
    .combine(inp, by:0)
    .map{ ref_name, target_name, impute_vcf, ref_vcf -> [target_name, ref_name, file(impute_vcf), file(ref_vcf)]}
    generate_frequency(input)

    // // Plot frequency Comparison
    freq_comp = impute_data.map {target_name, ref_name, vcf, impute_vcf, info -> 
    [target_name, ref_name, info]}
    .combine(generate_frequency.out, by:[0,1])
    plot_freq_comparison(freq_comp)

    // // Plot number of imputed SNPs over the mean r2 for all reference panels
    combineInfo_frq = impute_data.map{ target_name, ref_name, vcf, impute_vcf, info ->[ target_name, ref_name, info, params.maf_thresh]}
    .combine(generate_frequency.out, by:[0,1])
    .map { target_name, ref_name, info, maf_thresh, target_frq, ref_frq -> 
    [target_name, ref_name, info, maf_thresh, target_frq]}
    plot_r2_SNPpos(combineInfo_frq)

    // // compute for average rsquared values
    rsquared_input = impute_data.map{ target_name, ref_name, vcf, impute_vcf, info ->[ target_name, ref_name, info]}
    average_r2(rsquared_input)
}