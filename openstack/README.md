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
   Security" section of the compute menu in OpenStack.
2. Generate two instances using the "Ubuntu 16.04" base image which
   are at least m1.medium in size. These will be the swarm primary and
   secondary manager servers
3. Generate one instance using same base image which is at least
   m1.medium in size; this will be the discovery backend.
