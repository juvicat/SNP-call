#!/bin/bash
set -eo pipefail

prefix=$1 

# Generate a unique name for each SNP
bcftools annotate --set-id +"%CHROM:%POS" "$prefix".snp.filtered.nocall.vcf > "$prefix".snp.filtered.nodot.vcf

# Filtering
# Missing data and linkage disequilibrium:
plink --vcf-filter --vcf "$prefix".snp.filtered.nodot.vcf --allow-extra-chr --recode --make-bed --geno 0 --const-fid --out "$prefix"
plink --indep 50 5 2 --file "$prefix" --allow-extra-chr --out "$prefix"
plink --extract "$prefix".prune.in --out "$prefix"_pruned --file "$prefix" --make-bed --allow-extra-chr --recode

# Generate individual values for PCA with 20 axes:
plink --pca 20 --file "$prefix"_pruned --allow-extra-chr --out "$prefix"

# Generate statistics:
plink --freq --het 'small-sample' --ibc --file "$prefix"_pruned --allow-extra-chr -out "$prefix"_pruned
