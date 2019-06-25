#!/usr/bin/env bash


rm $3

while read roi; do

echo ${roi}
tabix $1 ${roi} | cut -d$'\t' -f1-5 | awk -F $'\t' 'length($4) == 1 && length($5) == 1' >> $3

done < $2
