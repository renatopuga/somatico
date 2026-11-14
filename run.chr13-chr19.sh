#####
# rodando o par wp190 (tumor) wp191 (normal)
#####

# variaveis fixas
gnomad="af-only-gnomad-chr13-chr19.vcf.gz"
interval="hg19.interval_list"
genome="hg19.fa"

# bam e ids 
tumor=$1
normal=$2
id_tumor=$3
id_normal=$4

# rodando mutect2
./gatk-4.2.2.0/gatk Mutect2 \
	-R $genome \
	-I $tumor \
	-I $normal \
	-normal $id_normal \
	--germline-resource $gnomad \
	-O $id_tumor.somatic.vcf.gz \
	-L $interval

# getPileup tumor
./gatk-4.2.2.0/gatk GetPileupSummaries \
	-I $tumor \
	-V $gnomad \
	-L $interval \
	-O $id_tumor.table

# getPileup normal
./gatk-4.2.2.0/gatk GetPileupSummaries \
	-I $normal \
	-V $gnomad \
	-L $interval \
	-O $id_normal.table

# CalculateContamination
./gatk-4.2.2.0/gatk CalculateContamination \
	-I $id_tumor.table \
	-matched $id_normal.table \
	-O $id_tumor.contamination.table

# FilterMutectCalls
./gatk-4.2.2.0/gatk FilterMutectCalls \
	-R $genome \
	-V $id_tumor.somatic.vcf.gz \
	--contamination-table $id_tumor.contamination.table \
	-O $id_tumor.filtered.vcf.gz



