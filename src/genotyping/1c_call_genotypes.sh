#!/bin/bash

chromosome=$1

echo -e "\n\n******************************\nChromosome: ${chromosome}\n******************************"

input_dir=/home/ramstein/sieve/M2/data/

work_dir=/scratch/${SLURM_JOBID}/
if [ ! -d "${work_dir}" ]; then mkdir ${work_dir}; fi

db_dir=${input_dir}db2/
db=${work_dir}${chromosome}
cp -r ${db_dir}${chromosome} ${db}

vcf_dir=${input_dir}vcf/
if [ ! -d "${vcf_dir}" ]; then mkdir ${vcf_dir}; fi

ref_file=/home/ramstein/sieve/M2/data/reference/BdistachyonBd21_3_537_v1.0.fa

variant_file=${vcf_dir}${chromosome}.vcf.gz
snp_file=${vcf_dir}${chromosome}.snps.vcf.gz
filtered_file=${vcf_dir}${chromosome}.snps.filtered.vcf.gz

###################################################################
# 1. Call all variants
###################################################################
gatk --java-options "-Xmx130g -Xms130g" GenotypeGVCFs \
-R ${ref_file} \
-V gendb://${db} \
-O ${variant_file}

###################################################################
# 2. Call SNP genotypes
###################################################################
gatk --java-options "-Xmx130g -Xms130g" SelectVariants \
-V ${variant_file} \
--select-type-to-include SNP  \
-O ${snp_file}

###################################################################
# 3. Filter SNPs
###################################################################
gatk --java-options "-Xmx130g -Xms130g" VariantFiltration \
-V ${snp_file} \
-filter "QD < 2.0" --filter-name "QD2" \
-filter "QUAL < 30.0" --filter-name "QUAL30" \
-filter "SOR > 3.0" --filter-name "SOR3" \
-filter "FS > 60.0" --filter-name "FS60" \
-filter "MQ < 40.0" --filter-name "MQ40" \
-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
-O ${filtered_file}
