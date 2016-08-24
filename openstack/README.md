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
3. Write your ssh-key *TO BE FINISHED*
4. Run `./generate_openstack` in this directory and wait for the VMs
   to be generated
5. When they have completed, you can log into one of the two manager
   nodes to start jobs on the docker swarm.
