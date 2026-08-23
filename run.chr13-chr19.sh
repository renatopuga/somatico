#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
  cat <<'EOF'
Uso:
  bash run.chr13-chr19.sh TUMOR.bam NORMAL.bam TUMOR_ID NORMAL_ID

Variáveis opcionais:
  GATK_CMD   Caminho ou comando do GATK (padrão: ./gatk-4.6.2.0/gatk)
  GENOME     FASTA de referência (padrão: hg19.fa)
  GNOMAD     VCF populacional (padrão: af-only-gnomad-chr13-chr19.vcf.gz)
  INTERVALS  Intervalos do exercício (padrão: hg19.interval_list)
  OUTPUT_DIR Diretório de saída (padrão: results)
EOF
}

if [[ $# -ne 4 ]]; then
  usage >&2
  exit 2
fi

tumor=$1
normal=$2
tumor_id=$3
normal_id=$4

GATK_CMD=${GATK_CMD:-./gatk-4.6.2.0/gatk}
GENOME=${GENOME:-hg19.fa}
GNOMAD=${GNOMAD:-af-only-gnomad-chr13-chr19.vcf.gz}
INTERVALS=${INTERVALS:-hg19.interval_list}
OUTPUT_DIR=${OUTPUT_DIR:-results}

require_file() {
  [[ -f "$1" ]] || {
    printf 'Erro: arquivo não encontrado: %s\n' "$1" >&2
    exit 1
  }
}

if [[ "$GATK_CMD" == */* ]]; then
  [[ -x "$GATK_CMD" ]] || {
    printf 'Erro: GATK não encontrado ou sem permissão de execução: %s\n' "$GATK_CMD" >&2
    exit 1
  }
else
  command -v "$GATK_CMD" >/dev/null 2>&1 || {
    printf 'Erro: comando GATK não encontrado: %s\n' "$GATK_CMD" >&2
    exit 1
  }
fi

for file in "$tumor" "$normal" "$GENOME" "$GNOMAD" "$INTERVALS"; do
  require_file "$file"
done

mkdir -p "$OUTPUT_DIR"

somatic_vcf="$OUTPUT_DIR/${tumor_id}.somatic.vcf.gz"
tumor_pileups="$OUTPUT_DIR/${tumor_id}.pileups.table"
normal_pileups="$OUTPUT_DIR/${normal_id}.pileups.table"
contamination="$OUTPUT_DIR/${tumor_id}.contamination.table"
filtered_vcf="$OUTPUT_DIR/${tumor_id}.filtered.vcf.gz"

printf '[1/5] Mutect2: %s + %s\n' "$tumor_id" "$normal_id"
"$GATK_CMD" Mutect2 \
  -R "$GENOME" \
  -I "$tumor" \
  -I "$normal" \
  -normal "$normal_id" \
  --germline-resource "$GNOMAD" \
  -L "$INTERVALS" \
  -O "$somatic_vcf"

printf '[2/5] GetPileupSummaries: %s\n' "$tumor_id"
"$GATK_CMD" GetPileupSummaries \
  -R "$GENOME" \
  -I "$tumor" \
  -V "$GNOMAD" \
  -L "$INTERVALS" \
  -O "$tumor_pileups"

printf '[3/5] GetPileupSummaries: %s\n' "$normal_id"
"$GATK_CMD" GetPileupSummaries \
  -R "$GENOME" \
  -I "$normal" \
  -V "$GNOMAD" \
  -L "$INTERVALS" \
  -O "$normal_pileups"

printf '[4/5] CalculateContamination\n'
"$GATK_CMD" CalculateContamination \
  -I "$tumor_pileups" \
  -matched "$normal_pileups" \
  -O "$contamination"

printf '[5/5] FilterMutectCalls\n'
"$GATK_CMD" FilterMutectCalls \
  -R "$GENOME" \
  -V "$somatic_vcf" \
  --contamination-table "$contamination" \
  -O "$filtered_vcf"

printf 'Concluído. VCF filtrado: %s\n' "$filtered_vcf"
