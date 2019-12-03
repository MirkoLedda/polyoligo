#!/bin/bash

CONTIGS="contig_list.txt"

echo "" > $CONTIGS

for file in blastdb/*.fa; do
    echo $file >> $CONTIGS
    grep ">" ${file} >> $CONTIGS
    echo "" >> $CONTIGS
done