# Deploying Workflow on OpenStack

We have tested and deployed this workflow using NCSA's Nebula, an
OpenStack based configuration.

The general layout is to run 

## Base images and deployment

We're currently using Ubuntu 16.04-based images. To deploy them, do
the following steps:

1. Create an ssh key pair (or use an existing one) which you will
   use to log into the images. `ssh-keygen -t rsa -f nebula_key` will
   generate an RSA key, and you can then add the contents of
   `nebula_key.pub` to the "key pairs" pane of the "Access and
   Security" section of the "Compute menu" in OpenStack.
2. Get the OpenStack RC file from the "API Access" pane of the "Access
   and Security" section mentioned above, source it in your shell, and
   enter your password.
3. Write your `ssh_keys_config` file which contains the ssh keys that
   should be allowed to log into the system, and also the `ca.pem`
   which will be trusted by the swarm (once we actually have TLS
   working).
4. Run `./generate_openstack` in this directory and wait for the VMs
   to be generated
5. When they have completed, you can log into one of the two manager
   nodes to start jobs on the docker swarm. (See Running jobs using
   docker below).

# Installing nextflow


# Running jobs using docker

To run a job using docker, ssh into one of the manager nodes. You can
do this using:

`nova ssh -i ssh_identity manager0 --login ubuntu
--extra-opts='-o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no'`

which will log into the manager node manager0.

You can then run docker containers like the following, which starts a
busybox shell, while mounting the /srv/imputation directory which is
shared between all nodes:

`docker -H :4000 run -i -u 1000 -P -v
/srv/imputation:/srv/imputation:Z busybox sh`

## Running the imputation workflow

You can also run our imputation workflow using nextflow:

     NXF_EXECUTOR_CPUS=100 /srv/imputation/nextflow/nextflow run \
     -qs 1000 -ps 1000 run.nf -with-docker 'quay.io/dondelelcaro/impute2:latest' \
     -c ../docker_nextflow

The environmental variable `NXF_EXECUTOR_CPUS=100` tells nextflow that
we have 100 available CPUs so that the local executor will schedule
more than a single job on the available docker swarm. This value may
need to be tuned for a particular workflow depending on how many
swarms you have running.


