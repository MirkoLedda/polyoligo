#!/bin/bash

N_CPU=(1 10 25 64)
N_MARKER=(1 10 100 1000)
N_REP=5
MARKERS="$HOME/Projects/KASP/strawberry/850K_array/markers_chr_subsets/Fvb1-1.txt"
OUTPUT="footprinting_data"

mkdir $OUTPUT

for n_marker in "${N_MARKER[@]}"
do
    for i in $(seq 1 $N_REP)
    do
    	shuf -n $n_marker $MARKERS > $OUTPUT/${n_marker}_${i}.markers

	for n_cpu in "${N_CPU[@]}"
    	do
		echo ${n_marker} $i ${n_cpu}
        	polyoligo-kasp $OUTPUT/${n_marker}_${i}.markers $OUTPUT/${n_marker}_${i}_${n_cpu} ~/data/Camarosa_genome_v1.2/F_ana_Camarosa_6-28-17.masked.fasta -nt $n_cpu --vcf ~/data/wgs/wgs.vcf.gz --vcf_include ~/data/wgs/vcf_subsets/F.xananassa_VCFinclude.txt --silent
	done
    done
done
