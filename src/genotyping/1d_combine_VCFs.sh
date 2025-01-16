#!/bin/bash

input_dir=/home/ramstein/sieve/M2/data/

tmp_dir=${input_dir}tmp/
if [ ! -d "${tmp_dir}" ]; then mkdir ${tmp_dir}; fi

vcf_dir=${input_dir}vcf/
if [ ! -d "${vcf_dir}" ]; then mkdir ${vcf_dir}; fi

vcf_files=(
    ${vcf_dir}Bd1.snps.filtered.vcf.gz
    ${vcf_dir}Bd2.snps.filtered.vcf.gz
    ${vcf_dir}Bd3.snps.filtered.vcf.gz
    ${vcf_dir}Bd4.snps.filtered.vcf.gz
    ${vcf_dir}Bd5.snps.filtered.vcf.gz
    )

list_file=${tmp_dir}/vcf_files.txt
if [ -f "${list_file}" ]; then rm ${list_file}; fi

for vcf_file in ${vcf_files[*]}
do
	echo "${vcf_file}" >> ${list_file}
done

output_file=${vcf_dir}snps.combined.vcf.gz

###################################################################
# 1. Combine and combine filtered variants
###################################################################
bcftools concat -f ${list_file} -Ou | bcftools view -f 'PASS,.' -Oz -o ${output_file}

tabix -f -p vcf ${output_file}
