#!/usr/bin/env bash
set -euo pipefail

# Diretório de saída
mkdir -p vep_output
chmod 777 vep_output

docker run --rm \
  -v "$(pwd)":/data \
  -w /data \
  ensemblorg/ensembl-vep \
  vep \
  -i filtered.vcf.gz \
  -o vep_output/filtered.vep.tsv \
  --database --assembly GRCh37 --refseq \
  --pick --pick_allele --force_overwrite --tab --symbol --check_existing \
  --fields "Location,SYMBOL,Consequence,Feature,Amino_acids,CLIN_SIG" \
  --fasta chr9.fa \
  --individual all
