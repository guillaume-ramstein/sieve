#!/bin/bash

chromosome=$1

echo -e "\n\n******************************\nChromosome: ${chromosome}\n******************************"

gvcf_dirs=(
    /home/ramstein/sieve/M2/data/batch_1/gvcf/
    /home/ramstein/sieve/M2/data/batch_2/gvcf/
    /home/ramstein/sieve/M2/data/batch_3/gvcf/
    /home/ramstein/sieve/M2/data/batch_4/gvcf/
    /home/ramstein/sieve/M2/data/batch_5/gvcf/
    /home/ramstein/sieve/M2/data/batch_6/gvcf/
    /home/ramstein/sieve/M2/data/batch_7/gvcf/
    /home/ramstein/sieve/M2/data/batch_8/gvcf/
    /home/ramstein/sieve/M2/data/batch_9/gvcf/
    )

work_dir=/scratch/${SLURM_JOBID}/
if [ ! -d "${work_dir}" ]; then mkdir ${work_dir}; fi

tmp_dir=${work_dir}tmp/
if [ ! -d "${tmp_dir}" ]; then mkdir ${tmp_dir}; fi

map_file=${tmp_dir}${chromosome}.sample_map
if [ -f "${map_file}" ]; then rm ${map_file}; fi

out_dir=/home/ramstein/sieve/M2/data/db2/
if [ ! -d "${out_dir}" ]; then mkdir ${out_dir}; fi

ref_file=/home/ramstein/sieve/M2/data/reference/BdistachyonBd21_3_537_v1.0.fa

###################################################################
# 1. Sample map
###################################################################
# Make sample map
for gvcf_dir in ${gvcf_dirs[*]}
do

    cd ${gvcf_dir}

    for gvcf_file in $(ls *g.vcf.gz)
    do

        sample=$(echo "${gvcf_file}" | sed 's/[.]g[.]vcf[.]gz//g' | sed 's/[AB]//g')

        echo -e "${sample}\t${gvcf_dir}${gvcf_file}\t${gvcf_dir}${gvcf_file}.tbi" >> ${map_file}

    done

done

###################################################################
# 2. Combining genotype files
###################################################################
cd ${work_dir}

db=${work_dir}${chromosome}

# DB import
if [ ! -d "${db}" ]; then

    gatk --java-options "-Xmx265g -Xms265g" GenomicsDBImport \
    -R ${ref_file} \
    --sample-name-map ${map_file} \
    --genomicsdb-workspace-path ${db} \
    --genomicsdb-shared-posixfs-optimizations true \
    -L ${chromosome} \
    --batch-size 250 \
    --tmp-dir ${tmp_dir} \
    --reader-threads 4

else

    gatk --java-options "-Xmx265g -Xms265g" GenomicsDBImport \
    -R ${ref_file} \
    --sample-name-map ${map_file} \
    --genomicsdb-update-workspace-path ${db} \
    --genomicsdb-shared-posixfs-optimizations true \
    --batch-size 250 \
    --tmp-dir ${tmp_dir} \
    --reader-threads 4

fi

# Copying DB from working directory to home directory
cp -r ${db} ${out_dir}${chromosome}
