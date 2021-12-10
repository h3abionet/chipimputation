#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { check_files } from './modules/qc' 
include { filter_info_by_target } from './modules/impute'
include { filter_info; report_site_by_maf; plot_freq_comparison; 
report_well_imputed_by_target; plot_performance_target; 
report_accuracy_target; plot_accuracy_target; generate_frequency;
plot_r2_SNPpos; plot_r2_SNPcount; plot_hist_r2_SNPcount; plot_MAF_r2; 
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

workflow{

    intro()
    target_datasets = []
        params.datasets.each { dataset, refpanel, vcf, imputed, info ->
            check_files([vcf, imputed, info])
            target_datasets << [dataset, refpanel, file(vcf), file(imputed), file(info) ]
        }
    target_datasets = Channel.from(target_datasets)
    //// Report by Ref
    report_by_ref( target_datasets )

    //// Report by datasets
    report_by_dataset( target_datasets )

    // Generate dataset frequencies
    inp = Channel.fromList(params.ref)
    input = target_datasets
    .map{ target_name, ref_name, vcf, impute_vcf, info ->[  ref_name, target_name, file(impute_vcf)]}
    .combine(inp, by:0)
    .map{ ref_name, target_name, impute_vcf, ref_vcf -> [target_name, ref_name, file(impute_vcf), file(ref_vcf)]}
    generate_frequency(input)

    // Plot frequency Comparison
    freq_comp = target_datasets.map {target_name, ref_name, vcf, impute_vcf, info -> 
    [target_name, ref_name, info]}
    .combine(generate_frequency.out, by:[0,1])
    plot_freq_comparison(freq_comp)

    // Plot number of imputed SNPs over the mean r2 for all reference panels
    combineInfo_frq = target_datasets.map{ target_name, ref_name, vcf, impute_vcf, info ->[ target_name, ref_name, info, params.maf_thresh]}
    .combine(generate_frequency.out, by:[0,1])
    .map { target_name, ref_name, info, maf_thresh, target_frq, ref_frq -> 
    [target_name, ref_name, info, maf_thresh, target_frq]}
    // combineInfo_frq.view()
    plot_r2_SNPpos(combineInfo_frq)

    // compute for average rsquared values

    rsquared_input = target_datasets.map{ target_name, ref_name, vcf, impute_vcf, info ->[ target_name, ref_name, info]}
    average_r2(rsquared_input)
}