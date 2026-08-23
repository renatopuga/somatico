# Análise de variantes somáticas com GATK Mutect2

Tutorial prático para compreender uma chamada de variantes somáticas com dados tumor–normal, estimativa de contaminação, filtragem e anotação funcional.

> **Finalidade educacional.** Os arquivos são pequenos recortes genômicos preparados para aula. Este repositório não constitui um pipeline clínico validado e seus resultados não devem ser usados para diagnóstico.

## O que você aprenderá

Ao concluir o exercício, você será capaz de:

- reconhecer os arquivos necessários para uma análise tumor–normal;
- verificar a compatibilidade entre os contigs do BAM, FASTA e VCF;
- executar `Mutect2`, `GetPileupSummaries`, `CalculateContamination` e `FilterMutectCalls`;
- comparar uma análise pareada com uma análise tumor-only usando Panel of Normals (PoN);
- identificar variantes que passaram pelos filtros e preparar o VCF para anotação com VEP.

## Escopo e versões

| Componente | Configuração deste exercício | Observação |
|---|---|---|
| Referência | GRCh37/hg19, somente os contigs usados nos exemplos | Mantida por compatibilidade com os BAMs e VCFs didáticos |
| Nomes dos contigs | `9`, `13` e `19`, sem o prefixo `chr` | Todos os arquivos de entrada precisam usar a mesma convenção |
| GATK | 4.6.2.0 | Versão fixada para reprodutibilidade |
| VEP | imagem Docker `release_116.0` | Executado em modo database; requer internet |

Para projetos novos, escolha a montagem de referência e os recursos populacionais de forma explícita. Não misture GRCh37/hg19 com GRCh38, nem arquivos com e sem o prefixo `chr`.

## Visão geral do fluxo

1. Preparar e indexar a referência FASTA.
2. Confirmar o identificador da amostra normal no read group do BAM.
3. Chamar SNVs e indels somáticos com Mutect2.
4. Estimar contaminação com os sítios populacionais do gnomAD.
5. Filtrar as chamadas.
6. Inspecionar variantes `PASS` e, opcionalmente, anotar com VEP.

## Arquivos incluídos

| Arquivo | Papel no exercício |
|---|---|
| `tumor_JAK2.bam` / `.bai` | Recorte tumoral da região de `JAK2` no cromossomo 9 |
| `normal_JAK2.bam` / `.bai` | Normal pareado do exemplo `JAK2` |
| `tumor_wp017.bam` e `tumor_wp018.bam` | Par tumor–normal adicional |
| `tumor_wp190.bam` e `tumor_wp191.bam` | Par tumor–normal adicional |
| `af-only-gnomad-*.vcf.gz` / `.tbi` | Sítios populacionais usados pelo Mutect2 e pela estimativa de contaminação |
| `somatic.Pon.vcf` | Pequeno VCF de exemplo; não substitui o PoN oficial usado na atividade |
| `run.chr13-chr19.sh` | Pipeline tumor–normal para os exemplos dos cromossomos 13 e 19 |
| `run.chr13-chr19.pon.sh` | Pipeline tumor-only com PoN |
| `vep-docker.sh` | Anotação do VCF filtrado com VEP em Docker |
| `somatico_google_colab.ipynb` | Aula guiada para Google Colab |

## Início rápido: exemplo JAK2

### 1. Pré-requisitos

Para a execução local, instale:

- Bash;
- `wget`, `gzip` e `unzip`;
- [samtools](https://www.htslib.org/);
- Java compatível com a versão do GATK;
- Docker, apenas para a etapa opcional com VEP.

Clone o repositório:

```bash
git clone https://github.com/renatopuga/somatico.git
cd somatico
```

Alternativamente, abra o notebook no Colab:

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/renatopuga/somatico/blob/main/somatico_google_colab.ipynb)

### 2. Baixar o GATK

```bash
wget -c https://github.com/broadinstitute/gatk/releases/download/4.6.2.0/gatk-4.6.2.0.zip
unzip gatk-4.6.2.0.zip
./gatk-4.6.2.0/gatk --version
```

### 3. Preparar a referência do cromossomo 9

Os BAMs do exercício usam o contig `9`, enquanto o FASTA da UCSC usa `chr9`. O comando abaixo remove o prefixo somente do cabeçalho FASTA.

```bash
wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr9.fa.gz
zcat chr9.fa.gz | sed '/^>/ s/^>chr/>/' > chr9.fa

samtools faidx chr9.fa
./gatk-4.6.2.0/gatk CreateSequenceDictionary \
  -R chr9.fa \
  -O chr9.dict

./gatk-4.6.2.0/gatk ScatterIntervalsByNs \
  -R chr9.fa \
  -O chr9.interval_list \
  -OT ACGT
```

Verifique se os três recursos usam o mesmo contig:

```bash
samtools view -H tumor_JAK2.bam | grep '^@SQ' | head
grep '^>' chr9.fa | head
bcftools view -h af-only-gnomad-chr9.vcf.gz | grep '^##contig' | head
```

### 4. Confirmar o nome da amostra normal

O valor informado em `-normal` deve ser exatamente igual ao campo `SM` do read group:

```bash
samtools view -H normal_JAK2.bam \
  | awk -F '\t' '$1=="@RG" {for (i=1; i<=NF; i++) if ($i ~ /^SM:/) {sub(/^SM:/,"",$i); print $i}}'
```

Para os arquivos deste exemplo, o identificador esperado é `WP044`.

### 5. Executar Mutect2 e estimar contaminação

```bash
./gatk-4.6.2.0/gatk Mutect2 \
  -R chr9.fa \
  -I tumor_JAK2.bam \
  -I normal_JAK2.bam \
  -normal WP044 \
  --germline-resource af-only-gnomad-chr9.vcf.gz \
  -L chr9.interval_list \
  -O somatic.vcf.gz

./gatk-4.6.2.0/gatk GetPileupSummaries \
  -R chr9.fa \
  -I tumor_JAK2.bam \
  -V af-only-gnomad-chr9.vcf.gz \
  -L chr9.interval_list \
  -O tumor_JAK2.table

./gatk-4.6.2.0/gatk GetPileupSummaries \
  -R chr9.fa \
  -I normal_JAK2.bam \
  -V af-only-gnomad-chr9.vcf.gz \
  -L chr9.interval_list \
  -O normal_JAK2.table

./gatk-4.6.2.0/gatk CalculateContamination \
  -I tumor_JAK2.table \
  -matched normal_JAK2.table \
  -O contamination.table
```

### 6. Filtrar e inspecionar as chamadas

```bash
./gatk-4.6.2.0/gatk FilterMutectCalls \
  -R chr9.fa \
  -V somatic.vcf.gz \
  --contamination-table contamination.table \
  -O filtered.vcf.gz

bcftools view -f PASS filtered.vcf.gz
```

O VCF bruto preserva evidências candidatas; o VCF filtrado adiciona decisões no campo `FILTER`. Em uma análise real, revise também profundidade, VAF, qualidade, orientação das reads, artefatos técnicos e evidência clínica.

## Exemplos dos cromossomos 13 e 19

Prepare uma referência combinada com nomes de contig compatíveis:

```bash
wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr13.fa.gz
wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr19.fa.gz
zcat chr13.fa.gz chr19.fa.gz | sed '/^>/ s/^>chr/>/' > hg19.fa

samtools faidx hg19.fa
./gatk-4.6.2.0/gatk CreateSequenceDictionary -R hg19.fa -O hg19.dict
./gatk-4.6.2.0/gatk ScatterIntervalsByNs -R hg19.fa -O hg19.interval_list -OT ACGT
```

Pipeline tumor–normal:

```bash
bash run.chr13-chr19.sh tumor_wp017.bam tumor_wp018.bam WP017 WP018
bash run.chr13-chr19.sh tumor_wp190.bam tumor_wp191.bam WP190 WP191
```

Os resultados são gravados em `results/`. Você pode alterar os caminhos por variáveis de ambiente:

```bash
GATK_CMD=/caminho/para/gatk \
GENOME=/referencias/hg19.fa \
GNOMAD=/referencias/af-only-gnomad-chr13-chr19.vcf.gz \
INTERVALS=/referencias/hg19.interval_list \
OUTPUT_DIR=results \
bash run.chr13-chr19.sh tumor.bam normal.bam TUMOR_ID NORMAL_ID
```

## Tumor-only com Panel of Normals

Baixe o PoN compatível com GRCh37/b37:

```bash
wget -c https://storage.googleapis.com/gatk-best-practices/somatic-b37/Mutect2-exome-panel.vcf
wget -c https://storage.googleapis.com/gatk-best-practices/somatic-b37/Mutect2-exome-panel.vcf.idx
```

Execute informando apenas o BAM tumoral e seu identificador:

```bash
bash run.chr13-chr19.pon.sh tumor_wp190.bam WP190
```

O PoN ajuda a remover artefatos recorrentes, mas não substitui um normal pareado. Para uso real, o PoN deve ser construído e validado com amostras normais processadas de modo comparável à coorte analisada.

## Anotação opcional com VEP

Com Docker instalado, o script aceita o VCF de entrada e o TSV de saída:

```bash
bash vep-docker.sh filtered.vcf.gz vep_output/filtered.vep.tsv
```

Por padrão, o script usa GRCh37, RefSeq e a imagem `ensemblorg/ensembl-vep:release_116.0`. Como o modo `database` consulta serviços externos, a execução requer internet e pode variar conforme a disponibilidade do Ensembl.

## Saídas principais

| Saída | Conteúdo |
|---|---|
| `*.somatic.vcf.gz` | Chamadas candidatas produzidas pelo Mutect2 |
| `*.pileups.table` | Contagens usadas para estimar contaminação |
| `*.contamination.table` | Fração estimada de contaminação |
| `*.filtered.vcf.gz` | Chamadas com os filtros do Mutect2 |
| `*.vep.tsv` | Consequências anotadas pelo VEP |

## Problemas frequentes

- **Contigs incompatíveis:** compare os cabeçalhos do BAM, FASTA e VCF. `9` e `chr9` são nomes diferentes.
- **Normal sample not found:** confirme o campo `SM` com `samtools view -H` e use exatamente esse valor em `-normal`.
- **Arquivo sem índice:** BAM requer `.bai`; FASTA requer `.fai` e `.dict`; VCF compactado requer `.tbi` ou `.csi`.
- **Recurso de montagem incorreta:** todos os arquivos devem pertencer à mesma montagem e à mesma convenção de contigs.
- **VEP sem conexão:** repita mais tarde ou configure cache/FASTA local para uma execução offline reprodutível.

## Referências

- [GATK: somatic short variant discovery](https://gatk.broadinstitute.org/hc/en-us/articles/360035894731-Somatic-short-variant-discovery-SNVs-Indels-)
- [GATK Mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2)
- [Diferenças entre GRCh37, hg19 e b37](https://gatk.broadinstitute.org/hc/en-us/articles/360035890711-GRCh37-hg19-b37-humanG1Kv37-Human-Reference-Discrepancies)
- [Ensembl Variant Effect Predictor](https://www.ensembl.org/info/docs/tools/vep/index.html)

## Licença

Distribuído sob a licença MIT. Consulte `LICENSE`.
