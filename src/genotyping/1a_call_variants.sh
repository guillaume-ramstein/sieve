#!/bin/bash

batch=$1
sample=$2

input_dir=/home/ramstein/sieve/M2/data/${batch}/

fq_dir=${input_dir}soapnuke/clean/

gvcf_dir=${input_dir}gvcf/
if [ ! -d "${gvcf_dir}" ]; then mkdir ${gvcf_dir}; fi

ref_file=/home/ramstein/sieve/M2/data/reference/BdistachyonBd21_3_537_v1.0.fa

###################################################################
# Align sequence reads
###################################################################
echo -e "\n\n**************************************************\n${sample}\n**************************************************"

cd ${fq_dir}${sample}

gvcf_file=${gvcf_dir}$(echo ${sample} | sed 's/M2_//').g.vcf.gz

if [ ! -f "${gvcf_file}" ]; then

    # Aligning sequence reads
    if [ ! -f "${sample}.bam" ]; then

        if [ $(ls *.fq.gz | wc -l) -gt 2 ]; then

            if [ ! -d "merged" ]; then mkdir merged; fi

            R1=merged/R1.fq.gz
            R2=merged/R2.fq.gz

            if [ ! -f "${R1}" ]; then cat $(ls *_1.fq.gz) > ${R1}; fi
            if [ ! -f "${R2}" ]; then cat $(ls *_2.fq.gz) > ${R2}; fi

        else

            R1=$(ls *_1.fq.gz)
            R2=$(ls *_2.fq.gz)

        fi

        bwa mem -t 1 ${ref_file} ${R1} ${R2} | samtools view -b -q20 -F0x904 | samtools sort > ${sample}.bam

    fi

    # Processing alignment
    picard MarkDuplicates INPUT=${sample}.bam OUTPUT=dedup_${sample}.bam M=${sample}_metrics.txt
    picard AddOrReplaceReadGroups I=dedup_${sample}.bam O=RG_dedup_${sample}.bam RGID=Bd21_3 RGLB=${sample} RGPL=DNBseq SORT_ORDER=coordinate RGPU=${sample} RGSM=${sample}

    samtools index RG_dedup_${sample}.bam
    samtools flagstat RG_dedup_${sample}.bam > RG_dedup_${sample}.bam.flagstat

    # Calling haplotypes
    gatk --java-options -Xmx16g HaplotypeCaller -R ${ref_file} -I RG_dedup_${sample}.bam -O ${gvcf_file} -ERC GVCF

    if [ -f "${gvcf_file}" ]; then
        rm ${sample}.bam
        rm dedup_${sample}.bam
        rm RG_dedup_${sample}.bam
        rm RG_dedup_${sample}.bam.bai
    fi

fi
