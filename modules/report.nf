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

process filter_info {
    tag "filter_${dataset_name}_${tagName}_${ref_panels.join('-')}"
    label "bigmem"

    input:
        tuple val(dataset_name), val(ref_panels), val(ref_infos), val(tagName), val(info_cutoff)
    output:
        tuple val(dataset_name), val(ref_panels), file("${well_out}_${info_cutoff}.tsv"), file("${acc_out}_${info_cutoff}.tsv"), val(tagName)
    script:
        comb_info = "${dataset_name}_${tagName}_${ref_panels.join('-')}.imputed_info"
        well_out = "${comb_info}_well_imputed"
        acc_out = "${comb_info}_accuracy"
        infos = ref_infos.join(',')
        impute_info_cutoff = params.impute_info_cutoff
        template "filter_info_minimac.py"
}


process report_site_by_maf {
    tag "site_by_maf_${dataset_name}"
    label "bigmem"

    input:
        tuple val(dataset_name), file(sites)
    output:
        tuple val(dataset_name), file("${site_by_maf}_number_of_imputed.tsv"), file("${site_by_maf}_number_of_imputed_summary.tsv"), file("${site_by_maf}_avg_info.tsv")
    script:
        site_by_maf = "${sites.baseName}_report_by_maf"
        group = 'REF_PANEL'
        template "report_well_imputed.py"
}


process plot_freq_comparison {
    tag "plot_freq_comparison_${dataset_name}_${ref_name}"
    // publishDir "${params.outDir}/plots/${ref_name}/freq_comparison", overwrite: true, mode:'copy'
    label "bigmem"
    input:
        tuple val(dataset_name), val(ref_name), file(dataset_info), file(dataset_frq), file(ref_frq)
    output:
        tuple val(dataset_name), val(ref_name), file(outputcolor)
    script:
        info = dataset_info
        target = dataset_frq
        frq = ref_frq
        //output = "${dataset_name}_${ref_name}_${chrm}_freq_comparison.png"
        outputcolor = "${dataset_name}_${ref_name}_freq_comparison_color.png"
        template "AF_comparison.R"
}

