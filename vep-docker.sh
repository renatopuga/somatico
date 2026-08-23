#!/usr/bin/env bash
set -Eeuo pipefail

input=${1:-filtered.vcf.gz}
output=${2:-vep_output/filtered.vep.tsv}

VEP_IMAGE=${VEP_IMAGE:-ensemblorg/ensembl-vep:release_116.0}
REFERENCE_FASTA=${REFERENCE_FASTA:-chr9.fa}

for file in "$input" "$REFERENCE_FASTA"; do
  [[ -f "$file" ]] || {
    printf 'Erro: arquivo não encontrado: %s\n' "$file" >&2
    exit 1
  }
done

command -v docker >/dev/null 2>&1 || {
  printf 'Erro: Docker não está instalado ou não está no PATH.\n' >&2
  exit 1
}

mkdir -p "$(dirname "$output")"

docker run --rm \
  --user "$(id -u):$(id -g)" \
  -v "$(pwd)":/data \
  -w /data \
  "$VEP_IMAGE" \
  vep \
  --input_file "$input" \
  --output_file "$output" \
  --database \
  --assembly GRCh37 \
  --refseq \
  --fasta "$REFERENCE_FASTA" \
  --pick \
  --pick_allele \
  --tab \
  --symbol \
  --check_existing \
  --force_overwrite \
  --fields "Uploaded_variation,Location,Allele,SYMBOL,Consequence,Feature,HGVSc,HGVSp,Amino_acids,Existing_variation,CLIN_SIG"

printf 'Concluído. Anotação: %s\n' "$output"
