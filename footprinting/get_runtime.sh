#!/bin/bash

rm runtime.txt

for file; do
    grep "Number of target markers" ${file}
    grep "Searching for KASP candidates using" ${file}
    grep "Total time elapsed" ${file}
done > runtime.txt

cut -f 6 -d "-" runtime.txt > runtime.temp
#grep -v ".log" runtime.temp > runtime.txt
mv runtime.temp runtime.txt
#rm runtime.temp
