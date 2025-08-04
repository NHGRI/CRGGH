#!/bin/bash

module load qtltools

QTLtools cis --vcf Genotypes.vcf.gz \
--bed Mvalues.bed.gz \
--cov Covariates.txt \
--nominal 1 \
--std-err \
--window 2000000 \
--out mQTLs_out.txt  \
--log mQTLs_log.txt