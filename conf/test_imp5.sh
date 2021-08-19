#!/usr/bin/env bash

cd /cbio/users/mamana

nextflow \
  run /users/mamana/chipimputation/main-dsl2_imp5.nf \
  -c /users/mamana/chipimputation/conf/test_imp5.config \
  -resume -profile slurm,singularity
