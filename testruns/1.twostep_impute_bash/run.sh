#!/bin/bash

set -e

impute2 \
 -prephase_g \
 -m ./input/example.chr22.map \
 -g ./input/example.chr22.study.gens \
 -int 20.4e6 20.5e6 \
 -Ne 20000 \
 -o ./temp/example.chr22.prephasing.impute2


impute2 \
 -use_prephased_g \
 -m ./input/example.chr22.map \
 -h ./refs/example.chr22.1kG.haps \
 -l ./refs/example.chr22.1kG.legend \
 -known_haps_g ./temp/example.chr22.prephasing.impute2_haps \
 -strand_g ./input/example.chr22.study.strand \
 -int 20.4e6 20.5e6 \
 -Ne 20000 \
 -o ./output/example.chr22.one.phased.impute2 \
 -phase
