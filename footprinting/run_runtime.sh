#!/bin/bash

N_CPU=(1 2 4 10 20 30 40 50 64)
N_MARKER=(1 10 100 1000)
MARKERS="~/Projects/KASP/strawberry/markers_chr_subsets/Fvb1-1.txt"
OUTPUT="footprinting_data"

mkdir $OUTPUT

for n_marker in "${N_MARKER[@]}"
do
    echo $OUTPUT/${n_marker}.txt
    #shuf -n $n_marker $MARKERS > $OUTPUT/${n_marker}.txt

    for n_cpu in "${N_CPU[@]}"
    do
        polyoligo-kasp $OUTPUT/${n_marker}.txt $OUTPUT/${n_marker}_${n_cpu} ~/data/Camarosa_genome_v1.2/F_ana_Camarosa_6-28-17.masked.fasta -nt $n_cpu --fast --vcf ~/data/wgs/wgs.vcf.gz --vcf_include ~/data/wgs/vcf_subsets/F.xananassa_VCFinclude.txt --silent
    done
done
