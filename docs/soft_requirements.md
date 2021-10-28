## Software Requirements
To run the pipeline, several software packages are required. How you satisfy these requirements is essentially up to you and depends on your system.  
If possible, we _highly_ recommend using either Docker or Singularity.

Please see the [installation documentation](docs/installation.md) for how to run using the below as a one-off. These instructions are about configuring a config file for repeated use.

### Docker
Docker is a great way to run h3abionet/chipimputation, as it manages all software installations and allows the pipeline to be run in an identical software environment across a range of systems.

Nextflow has [excellent integration](https://www.nextflow.io/docs/latest/docker.html) with Docker, and beyond installing the two tools, not much else is required - nextflow will automatically fetch the [nfcore/imp](https://quay.io/repository/h3abionet_org/imputation_tools?tag=latest&tab=tags) image that we have created and is hosted at dockerhub at run time.

To add docker support to your own config file, add the following:

```nextflow
docker.enabled = true
process.container = "nfcore/imp"
```
Note that the dockerhub organisation name annoyingly can't have a hyphen, so is `nfcore` and not `nf-core`.

You can now run the pipeline using the command:

```bash
 nextflow run h3abionet/chipimputation -profile docker
```
### Singularity image
Many HPC environments are not able to run Docker due to security issues.
[Singularity](http://singularity.lbl.gov/) is a tool designed to run on such HPC systems which is very similar to Docker.

To specify singularity usage in your pipeline config file, add the following:

```nextflow
singularity.enabled = true
process.container = "shub://h3abionet/chipimputation"
```

If you intend to run the pipeline offline, nextflow will not be able to automatically download the singularity image for you.
Instead, you'll have to do this yourself manually first, transfer the image file and then point to that.

First, pull the image file where you have an internet connection:

```bash
singularity pull --name nf-core-imp.simg shub://h3abionet/chipimputation
```

Then transfer this file and point the config file to the image:

```nextflow
singularity.enabled = true
process.container = "/path/to/nf-core-imp.simg"
```

The pipeline can then be run using the command

```bash
 nextflow run h3abionet/chipimputation -profile singularity
```
### Conda
If you're not able to use Docker or Singularity, you can instead use conda to manage the software requirements.
To use conda in your own config file, add the following:

```nextflow
process.conda = "$baseDir/environment.yml"
```
