https://gitpod.io/#https://github.com/renatopuga/somatico

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



## Workflow

### Converter BAM para FASTQ

```bash
samtools view -h -b /Volumes/Seagate\ Expansion\ Drive/data-lpfap10/projects/proadi/exome/bam/WP043/WP043.bam 9:5021937-5126899 | samtools bam2fq -1 tumor_R1.fq -2 tumor_R2.fq - 
```

```bash
samtools view -h -b /Volumes/Seagate\ Expansion\ Drive/data-lpfap10/projects/proadi/exome/bam/WP044/WP044.bam 9:5021937-5126899 | samtools bam2fq -1 normal_R1.fq -2 normal_R2.fq - 
```



### Converter BAM para JAK2 BAM

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



### Download da Referência - chr9

* Download 

```bash
wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr9.fa.gz
```

* Descompactar

```bash
gunzip chr9.fa.gz 
```



## BWA

Algoritmo de alinhamento.



### BWA install

* BWA install (Mac)

```bash
brew install bwa
```

* BWA install (Ubuntu)

```bash
sudo apt-get install bwa
```

* BWA install (Docker)

```bash
docker pull comics/bwa
```



### BWA index

* Index chr9

```bash
bwa index chr9.fa 
```



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

* samtools index

```bash
samtools index tumor_JAK2.bam
```

```bash
samtools index normal_JAK2.bam
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

```bash
./gatk-4.2.2.0/gatk Mutect2 -R chr9.fa -I tumor_JAK2.bam -I normal_JAK2.bam -normal WP044 --germline-resource af-only-gnomad-chr9.vcf.gz -O somatic.vcf.gz -L chr9.interval_list
```



## Calcular Contaminação



### GetPileupSummaries

Tabulates pileup metrics for inferring contamination

* GetPileupSummaries Tumor

```bash
./gatk-4.2.2.0/gatk GetPileupSummaries -I tumor_JAK2.bam -V af-only-gnomad-chr9.vcf.gz -L chr9.interval_list  -O tumor_JAK2.table
```

* GetPileupSummaries Normal

```bash
./gatk-4.2.2.0/gatk GetPileupSummaries -I normal_JAK2.bam -V af-only-gnomad-chr9.vcf.gz -L chr9.interval_list  -O normal_JAK2.table
```



### CalculateContamination

Calculate the fraction of reads coming from cross-sample contamination

```bash
./gatk-4.2.2.0/gatk CalculateContamination -I tumor_JAK2.table -matched normal_JAK2.table -O contamination.table
```



### FilterMutectCalls

Filter somatic SNVs and indels called by Mutect2

```bash
./gatk-4.2.2.0/gatk FilterMutectCalls -R chr9.fa -V somatic.vcf.gz --contamination-table contamination.table -O filtered.vcf.gz
```



## Funcotator

Functional Annotator

* 30G Source Somatic (1s)
  *  ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/funcotator/funcotator_dataSources.v1.7.20200521s.tar.gz

* Download Funcotator

```bash
wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/funcotator/funcotator_dataSources.v1.7.20200521s.tar.gz
```

* Descompactar

```bash
tar -zxvf funcotator_dataSources.v1.7.20200521s.tar.gz
```

* Anotar

```bash
./gatk-4.2.2.0/gatk Funcotator --data-sources-path funcotator_dataSources.v1.7.20200521s -O funcotator.maf --output-file-format MAF --ref-version hg19 --reference chr9.fa --variant filtered.vcf.gz
```
