#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
  cat <<'EOF'
Uso:
  bash run.chr13-chr19.pon.sh TUMOR.bam TUMOR_ID

Variáveis opcionais:
  GATK_CMD   Caminho ou comando do GATK (padrão: ./gatk-4.6.2.0/gatk)
  GENOME     FASTA de referência (padrão: hg19.fa)
  GNOMAD     VCF populacional (padrão: af-only-gnomad-chr13-chr19.vcf.gz)
  INTERVALS  Intervalos do exercício (padrão: hg19.interval_list)
  PON        Panel of Normals (padrão: Mutect2-exome-panel.vcf)
  OUTPUT_DIR Diretório de saída (padrão: results)
EOF
}

if [[ $# -ne 2 ]]; then
  usage >&2
  exit 2
fi

tumor=$1
tumor_id=$2

GATK_CMD=${GATK_CMD:-./gatk-4.6.2.0/gatk}
GENOME=${GENOME:-hg19.fa}
GNOMAD=${GNOMAD:-af-only-gnomad-chr13-chr19.vcf.gz}
INTERVALS=${INTERVALS:-hg19.interval_list}
PON=${PON:-Mutect2-exome-panel.vcf}
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

for file in "$tumor" "$GENOME" "$GNOMAD" "$INTERVALS" "$PON"; do
  require_file "$file"
done

mkdir -p "$OUTPUT_DIR"

somatic_vcf="$OUTPUT_DIR/${tumor_id}.somatic.pon.vcf.gz"
pileups="$OUTPUT_DIR/${tumor_id}.pon.pileups.table"
contamination="$OUTPUT_DIR/${tumor_id}.contamination.pon.table"
filtered_vcf="$OUTPUT_DIR/${tumor_id}.filtered.pon.vcf.gz"

printf '[1/4] Mutect2 tumor-only com PoN: %s\n' "$tumor_id"
"$GATK_CMD" Mutect2 \
  -R "$GENOME" \
  -I "$tumor" \
  --germline-resource "$GNOMAD" \
  --panel-of-normals "$PON" \
  -L "$INTERVALS" \
  -O "$somatic_vcf"

printf '[2/4] GetPileupSummaries\n'
"$GATK_CMD" GetPileupSummaries \
  -R "$GENOME" \
  -I "$tumor" \
  -V "$GNOMAD" \
  -L "$INTERVALS" \
  -O "$pileups"

printf '[3/4] CalculateContamination\n'
"$GATK_CMD" CalculateContamination \
  -I "$pileups" \
  -O "$contamination"

printf '[4/4] FilterMutectCalls\n'
"$GATK_CMD" FilterMutectCalls \
  -R "$GENOME" \
  -V "$somatic_vcf" \
  --contamination-table "$contamination" \
  -O "$filtered_vcf"

printf 'Concluído. VCF filtrado: %s\n' "$filtered_vcf"
