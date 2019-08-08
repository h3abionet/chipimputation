# Configuration file


A basic configuration comes with the pipeline, which runs by default (the `standard` config profile - see [`conf/base.config`](../../conf/base.config)). This means that you only need to configure the specifics for your system and overwrite any defaults that you want to change.  

> If you think that there are other people using the pipeline who would benefit from your configuration (eg. other common cluster setups), please let us know. We can add a new configuration and profile which can used by specifying `profile <name>` when running the pipeline.

## Creating own config file

To run the pipeline using your own dataset, you will need to create your config file as `your_project.config` and it will be applied every time you run Nextflow. 
While running the pipeline with the `test` profile, the test configuration file `test.config` will be copied into your output folder (`./output`). Simply update the `test.config` with parameters pertaining to your data and save the file anywhere, and reference it when running the pipeline with `-c path/to/test.config` (see the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more).

## Main Arguments

### `-profile`
Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments. Note that multiple profiles can be loaded, for example: `-profile standard,docker` - the order of arguments is important!

* `standard`
    * The default profile, used if `-profile` is not specified at all.
    * Runs locally and expects all software to be installed and available on the `PATH`.
* `singularity`
    * A generic configuration profile to be used with [Singularity](http://singularity.lbl.gov/)
    * Pulls container quay.io/h3abionet_org/imputation_tools from http://quay.io/h3abionet_org
* `docker`
    * A generic configuration profile to be used with [Docker](http://docker.com/)
    * Pulls container quay.io/h3abionet_org/imputation_tools from http://quay.io/h3abionet_org
* `conda`
    * A generic configuration profile to be used with [conda](https://conda.io/docs/)
    * Pulls most software from [Bioconda](https://bioconda.github.io/)
    * Note that not all tools needed for this pipeline are available through conda. This is the case of `Minimac4` for instance.
* `test`
    * A profile with a complete configuration for automated testing
    * Includes links to test data so needs no other parameters
    * This will copy the test configuration file into your current directory

### Target/Study dataset`--target_datasets`
Use this to specify the location of your input target dataset files in VCF format.  
Multiple target datasets can be specified in `target_datasets` of format `name = dataset`, however each target dataset will be used separately.  

The syntax for this : 
```bash
target_datasets {
    Study_name1 = "https://github.com/h3abionet/chipimputation_test_data/raw/master/testdata_imputation/target_testdata.vcf.gz"
}
```

A test data is provided in https://github.com/h3abionet/chipimputation_test_data/raw/master/testdata_imputation/target_testdata.vcf.gz, which can be used for testing only. 

Please note the following requirements:
1. This is required by the pipeline
1. The path must be enclosed in quotes and must exist otherwise the pipeline will stop.
2. The VCF file can contain a single or multiple chromosomes

### Reference panels `--ref_panels`
The pipeline expects uses minimac4 to imputed genotypes. Therefore, minimac3 reference format [m3vcf](https://genome.sph.umich.edu/wiki/M3VCF_Files) generated used [minimac3](https://genome.sph.umich.edu/wiki/Minimac3) is expected to be used.  
You need to specify both `VCF` and `M3VCF` files for `vcfFile` and `m3vcfFile` respectively in the configuration file before you launch the pipeline. 
A normal glob pattern, enclosed in quotation marks, can then be used. 

The syntax for this :
```bash
ref_panels {
    RefPanel_name1 {
      m3vcfFile   = "refPanel_testdata_22_phased.m3vcf.gz"
      vcfFile     = "refPanel_testdata_22_phased.vcf.gz"
    }
  }
```

A test data is provided in https://github.com/h3abionet/chipimputation_test_data repo:  
- `M3VCF`: https://github.com/h3abionet/chipimputation_test_data/raw/master/testdata_imputation/refPanel_testdata_%s_phased.m3vcf.gz  
- `VCF`: https://github.com/h3abionet/chipimputation_test_data/raw/master/testdata_imputation/refPanel_testdata_%s_phased.vcf.gz  
    
Please note the following requirements:
1. Both `VCF` and `M3VCF` files must be in chromosomes. String extrapolation of `%s` will be used to replace the chromosome
2. The `VCF` files will be used during phasing by `eagle2` and allele frequency comparison by `bcftools` steps
3. The `M3VCF` files will be used during imputation step by `minimac4`
4. The path must be enclosed in quotes and must exist otherwise the pipeline will stop.
5. Must be of the same build as the target dataset.

Some commonly used reference panels are available for download here from [minimac3 website](https://genome.sph.umich.edu/wiki/Minimac3#Reference_Panels_for_Download) including  `1000 Genomes Phase 1` (version 3) and  `1000 Genomes Phase 3` (version 5).  
To generate your own `M3VCF` files from `VCF` files using `minimac3`, please follow the instructions below as described https://genome.sph.umich.edu/wiki/Minimac3_Examples
```bash
Minimac3 --refHaps refPanel.vcf \ 
                --processReference \ 
                --rounds 0 \ 
                --prefix testRun
```

## Reference Genomes --reference_genome

A human reference genome in `fasta` format of the same build as the target dataset is required by the pipeline during the QC step to check the REF mismatch between in the target dataseet.
This can be downloaded from the [AWS-iGenomes](https://ewels.github.io/AWS-iGenomes/) resource.  
An test fasta file that can be used with the test dataset is provided on https://github.com/h3abionet/chipimputation_test_data/raw/master/testdata_imputation/hg19_testdata.fasta.gz

The syntax for this :
```bash
reference_genome = hg19_testdata.fasta.gz
```

## Genetic map for eagle2 --eagle_genetic_map
A genetic map file is required during phasing phase. A full 

## Job Resources
### Automatic resubmission
Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

### Custom resource requests
Wherever process-specific requirements are set in the pipeline, the default value can be changed by creating a custom config file. See the files in [`conf`](../conf) for examples.

Please make sure to also set the `-w/--work-dir` and `--outdir` parameters to a S3 storage bucket of your choice - you'll get an error message notifying you if you didn't.

## Other command line parameters

### `--outdir`
The output directory where the results will be saved.

### `--email`
Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to speicfy this on the command line for every run.

### `-name`
Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

This is used in the MultiQC report (if not default) and in the summary HTML / e-mail (always).

**NB:** Single hyphen (core Nextflow option)

### `-resume`
Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

**NB:** Single hyphen (core Nextflow option)

### `-c`
Specify the path to a specific config file (this is a core NextFlow command).

**NB:** Single hyphen (core Nextflow option)

Note - you can use this to override defaults. For example, you can specify a config file using `-c` that contains the following:

```nextflow
process.$multiqc.module = []
```

### `--max_memory`
Use to set a top-limit for the default memory requirement for each process.
Should be a string in the format integer-unit. eg. `--max_memory '8.GB'``

### `--max_time`
Use to set a top-limit for the default time requirement for each process.
Should be a string in the format integer-unit. eg. `--max_time '2.h'`

### `--max_cpus`
Use to set a top-limit for the default CPU requirement for each process.
Should be a string in the format integer-unit. eg. `--max_cpus 1`

### `--plaintext_email`
Set to receive plain-text e-mails instead of HTML formatted.

### `--multiqc_config`
Specify a path to a custom MultiQC configuration file.
