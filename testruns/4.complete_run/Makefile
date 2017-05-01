#!/usr/bin/make -f

## Use this makefile to generate some sweet test data. You can run it
## with `make;`

## if you need to do something more complicated, like have your
## reference data somewhere else, use make
## IMPUTATION_REFDATA=/path/to/refdata/

## for the testing data, only bother to make chromosomes 19 and 20
CHRS=19 20

SAMPLE_FILES=$(foreach chr,$(CHRS),chr$(chr).sample)
GENS_FILES=$(foreach chr,$(CHRS),chr$(chr).gens)

IMPUTATION_REFDATA=/srv/imputation/refdata/1000GP_Phase3/

HAP_FILES=$(foreach chr,$(CHRS),$(IMPUTATION_REFDATA)1000GP_Phase3_chr$(chr).hap.gz)
LEGEND_FILES=$(foreach chr,$(CHRS),$(IMPUTATION_REFDATA)1000GP_Phase3_chr$(chr).legend.gz)


all: testdata.map testdata.ped

testdata.map testdata.ped: $(HAP_FILES) $(LEGEND_FILES) \
$(IMPUTATION_REFDATA)1000GP_Phase3.sample
	../generate_testdata_ped_map.pl --skip 50 $(foreach				\
	hap,$(HAP_FILES),--hap $(hap)) $(foreach							\
	legend,$(LEGEND_FILES),--legend $(legend)) --sample $(word 3,$^)	\
	--ped-out testdata.ped --map-out testdata.map
