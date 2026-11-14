#####
# rodando o par wp190 (tumor) wp191 (normal)
#####

# variaveis fixas
gnomad="af-only-gnomad-chr13-chr19.vcf.gz"
interval="hg19.interval_list"
genome="hg19.fa"
pon="Mutect2-exome-panel.vcf"

# bam e ids 
tumor=$1
id_tumor=$3

# rodando mutect2
./gatk-4.2.2.0/gatk Mutect2 \
	-R $genome \
	-I $tumor \
	--germline-resource $gnomad \
	--panel-of-normals $pon \
	-O $id_tumor.somatic.pon.vcf.gz \
	-L $interval

# getPileup tumor
./gatk-4.2.2.0/gatk GetPileupSummaries \
	-I $tumor \
	-V $gnomad \
	-L $interval \
	-O $id_tumor.pon.table

# CalculateContamination
./gatk-4.2.2.0/gatk CalculateContamination \
	-I $id_tumor.pon.table \
	-O $id_tumor.contamination.pon.table

# FilterMutectCalls
./gatk-4.2.2.0/gatk FilterMutectCalls \
	-R $genome \
	-V $id_tumor.somatic.pon.vcf.gz \
	--contamination-table $id_tumor.contamination.pon.table \
	-O $id_tumor.filtered.pon.vcf.gz



