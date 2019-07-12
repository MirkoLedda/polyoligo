#!/bin/bash

N_CPU=4

mkdir footprinting_data

for i in $(seq 1 $N_CPU)
do
    echo $i
    mprof run -C -o footprinting_data/mprof_${i}.dat polyoligo-kasp ../sample_data/markers.txt footprinting_data/out ~/Desktop/temp/data/Camarosa_genome_v1.2/F_ana_Camarosa_6-28-17.masked.fasta --vcf ~/Desktop/temp/data/wgs/wgs.vcf.gz -nt $i --silent
done
