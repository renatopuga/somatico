[![Gitpod ready-to-code](https://gitpod.io/button/open-in-gitpod.svg)](https://gitpod.io/#https://github.com/renatopuga/somatico)

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1r4LDUiQqirUFQT2nIrhVW_AHLaBhmS7y?usp=sharing)


# Análise Somática
GATK 4 Mutect2 Somático

# T12020 - Aula Prática

![image](https://user-images.githubusercontent.com/8321336/130251648-7ad77cae-435f-49be-950f-b7af5f59fd7a.png)

## Pipeline

* GATK4 - Mutect2
* Gene JAK2
* Referência chr9
  * Sobre as versões do Genoma Humano: https://gatk.broadinstitute.org/hc/en-us/articles/360035890711-GRCh37-hg19-b37-humanG1Kv37-Human-Reference-Discrepancies#grch37
* Amostras: 
  * WP043 (tumor)
  * WP044 (normal)
* https://gatk.broadinstitute.org/hc/en-us/articles/360035894731-Somatic-short-variant-discovery-SNVs-Indels-


## Amostras Extras

- WP190 (tumor) e WP191 (normal)
- WP017 (tumor) e WP018 (normal)

**Nota 1:** Utilizar o af-gnomad chr13 e chr19. 

**Nota 2:** Será preciso baixar o chr13 e chr19 da UCSC.

**Download chr19**
```bash
wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr19.fa.gz
```

**Download chr13**
```bash
wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr13.fa.gz
```



**Concatenar os arquivos .fa.gz**
> Dica 1: zcat lê arquivos compactados .gz e zip

```bash
zcat chr13.fa.gz chr19.fa.gz | sed -e "s/chr//g" > hg19.fa
```

**Gerar o index do arquivo hg19.fa**
```bash
samtools faidx hg19.fa
```

 ____________________ 
< Primeiro é o chr9 >
 -------------------- 
        \   ^__^
         \  (oo)\_______
            (__)\       )\/\
                ||----w |
                ||     ||



## Workflow
Os arquivos BAM (tumor e normal) já foram gerados e estão prontos para a chamada de variates (ver Anexo 1). Então, agora vamos fazer download da referência `chr9` e gerar o index com o programa `samtools`.

### Download da Referência - chr9

* Download 

```bash
wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr9.fa.gz
```

* Alterar nome do header: DE: >chr9 para >9
> Essa alteração é necessária pois no BAM a referência não tinha `>chr` era apenas `>9`.

```bash
zcat chr9.fa.gz | sed -e "s/chr//g" > chr9.fa
```

* Verificar se alteração foi feita com o comando `head`

```bash
head chr9.fa
```
> O comando `head` lê as 10 primeiras linha de um arquivo texto


## samtools

Samtools is a suite of programs for interacting with high-throughput sequencing data. It consists of three separate repositories:



### samtools install

* samtools install (Mac)

```bash
brew install samtools 
```

* samtools install (Ubuntu)

```bash
sudo apt-get install samtools 
```

* samtools install (Docker)

```bash
docker pull biocontainers/samtools
```



### samtools faidx e index

* samtools faidx

```bash
samtools faidx chr9.fa
```


## GATK4

> Version: 4.2.2.0

Genome Analysis Toolkit - Variant Discovery in High-Throughput Sequencing Data. https://gatk.broadinstitute.org/



### GATK4 install

GATK4 install Docker

```bash
docker pull broadinstitute/gatk:4.2.2.0
```



GATK4 install Mac e Linux

* Download

```bash
wget -c https://github.com/broadinstitute/gatk/releases/download/4.2.2.0/gatk-4.2.2.0.zip
```

* Descompactar

```bash
unzip gatk-4.2.2.0.zip 
```

* Testando gatk

```bash
./gatk-4.2.2.0/gatk
```



### GATK4 .dict

```bash
./gatk-4.2.2.0/gatk CreateSequenceDictionary -R chr9.fa -O chr9.dict
```



### GATK4 intervals

```bash
./gatk-4.2.2.0/gatk ScatterIntervalsByNs -R chr9.fa -O chr9.interval_list -OT ACGT
```



## Mutect2

Call somatic SNVs and indels via local assembly of haplotypes



### Mutect2 Tumor e Normal
> O comando: `samtools view -H normal_JAK2.bam` você consegue pegar o campo SM: que contém o ID da amostra normal.

```bash
./gatk-4.2.2.0/gatk Mutect2 \
	-R chr9.fa \
	-I tumor_JAK2.bam \
	-I normal_JAK2.bam \
	-normal WP044 \
	--germline-resource af-only-gnomad-chr9.vcf.gz \
	-O somatic.vcf.gz \
	-L chr9.interval_list
```



## Calcular Contaminação



### GetPileupSummaries

Tabulates pileup metrics for inferring contamination

* GetPileupSummaries Tumor

```bash
./gatk-4.2.2.0/gatk GetPileupSummaries \
	-I tumor_JAK2.bam \
	-V af-only-gnomad-chr9.vcf.gz \
	-L chr9.interval_list \
	-O tumor_JAK2.table
```

* GetPileupSummaries Normal

```bash
./gatk-4.2.2.0/gatk GetPileupSummaries \
	-I normal_JAK2.bam \
	-V af-only-gnomad-chr9.vcf.gz \
	-L chr9.interval_list \
	-O normal_JAK2.table
```



### CalculateContamination

Calculate the fraction of reads coming from cross-sample contamination

```bash
./gatk-4.2.2.0/gatk CalculateContamination \
	-I tumor_JAK2.table \
	-matched normal_JAK2.table \
	-O contamination.table
```



### FilterMutectCalls

Filter somatic SNVs and indels called by Mutect2

```bash
./gatk-4.2.2.0/gatk FilterMutectCalls \
	-R chr9.fa \
	-V somatic.vcf.gz \
	--contamination-table contamination.table \
	-O filtered.vcf.gz
```

## VEP ensembl - Anotação

### VEP install

```bash
docker pull ensemblorg/ensembl-vep
```

* Criar diretório de output do VEP com permissão total (aplicado apenas no gitpod)

Voltar para a casa
```
cd
```

Verificar se o diretório é o /home/gitpod
```
pwd
```
> ex.: /home/gitpod

Copiar o arquivo filtered.vcf.gz
```
cp /workspace/somatico/filtered.vcf.gz .
```

Copiar o arquivo chr9.fa
```
cp /workspace/somatico/chr9.fa .
```

Copiar o arquivo chr9.fa.fai
```
cp /workspace/somatico/chr9.fa.fai .
```


Criar o diretorio vep_output
```
mkdir -p vep_output
```

Modificar a permissao do diretorio vep_output
```
chmod 777 vep_output
```

* Aplicar apenas no Google Colab

```bash
mkdir -p vep_output
chmod 777 vep_output
```

# rodar o vep

```bash
docker run -it --rm  -v $(pwd):/data ensemblorg/ensembl-vep ./vep  \
	-i /data/filtered.vcf.gz  \
	-o /data/vep_output/filtered.vep.tsv \
	--database --assembly GRCh37 --refseq  \
	--pick --pick_allele --force_overwrite --tab --symbol --check_existing \
	--fields "Location,SYMBOL,Consequence,Feature,Amino_acids,CLIN_SIG" \
	--fasta /data/chr9.fa \
	--individual all 
```



## Panel of Normal (PoN)


# Panel of Normal (PoN)

GATK Best Practices - Exome PoN


* vcf

```bash
wget -c https://storage.googleapis.com/gatk-best-practices/somatic-b37/Mutect2-exome-panel.vcf
```

* vcf.idx

```bash
wget -c https://storage.googleapis.com/gatk-best-practices/somatic-b37/Mutect2-exome-panel.vcf.idx
```

* Mutect2

```bash
./gatk-4.2.2.0/gatk Mutect2 \
  -R hg19.fa \
  -I tumor_wp190.bam \
  --germline-resource af-only-gnomad-chr13-chr19.vcf.gz \
  --panel-of-normals Mutect2-exome-panel.vcf \
  -L hg19.interval_list \
  -O WP190.somatic.pon.vcf.gz
  
```

* CalculateContamination somente com o table do tumor (ex.: wp190)

```bash
./gatk-4.2.2.0/gatk CalculateContamination \
	-I tumor_wp190.table \
	-O WP190.contamination.pon.table
```

* FilterMutectCalls

```bash
./gatk-4.2.2.0/gatk FilterMutectCalls \
	-R hg19.fa \
	-V WP190.somatic.pon.vcf.gz \
	--contamination-table WP190.contamination.pon.table \
	-O WP190.filtered.pon.vcf.gz
```





# Anexo 1

### Converter BAM para FASTQ

```bash
samtools view -h -b /Volumes/Seagate\ Expansion\ Drive/data-lpfap10/projects/proadi/exome/bam/WP043/WP043.bam 9:5021937-5126899 | samtools bam2fq -1 tumor_R1.fq -2 tumor_R2.fq - 
```

```bash
samtools view -h -b /Volumes/Seagate\ Expansion\ Drive/data-lpfap10/projects/proadi/exome/bam/WP044/WP044.bam 9:5021937-5126899 | samtools bam2fq -1 normal_R1.fq -2 normal_R2.fq - 
```



### Converter BAM para JAK2 BAM

Aqui estão os comandos que foram utilizados para gerar os BAMs intermediários e como gerar arquivos FASTQs de regiões específicas do seu arquivo BAM completo.
> Essas etapas não precisam ser executadas nesse pipeline

* BAM para BAM

```bash
samtools view -h -b /Volumes/Seagate\ Expansion\ Drive/data-lpfap10/projects/proadi/exome/bam/WP043/WP043.bam 9:5021937-5126899 > tumor_JAK2.bam
```

```bash
samtools view -h -b /Volumes/Seagate\ Expansion\ Drive/data-lpfap10/projects/proadi/exome/bam/WP044/WP044.bam 9:5021937-5126899 > normal_JAK2.bam
```

* Gerar index do BAM (.BAI)

```bash
samtools index tumor_JAK2.bam 
```

```bash
samtools index normal_JAK2.bam 
```



### af-only-gnomad.vcf.gz (apenas região JAK2)
> [Google Clou af-only-gnomad.vcf.gz](https://console.cloud.google.com/storage/browser/gatk-best-practices/somatic-b37;tab=objects?project=broad-dsde-outreach&pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))&forceOnObjectsSortingFiltering=false&pli=1)

* Header do VCF

```bash
zgrep -w "\#" af-only-gnomad.raw.sites.chr.vcf.gz > header
```

* Apenas Região do Gene JAK2

```bash
zgrep -w "^chr9" af-only-gnomad.raw.sites.chr.vcf.gz  | awk '$2>=5021937 && $2<=5126899' > JAK2.region
```

* Concatenar header + vcf

```bash
cat header JAK2.region > af-only-gnomad-chr9.vcf
```

* Compactar

```bash
bgzip af-only-gnomad-chr9.vcf
```

* Index do VCF

```bash
tabix -p vcf af-only-gnomad-chr9.vcf.gz 
```



