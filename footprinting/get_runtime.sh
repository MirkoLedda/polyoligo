#!/bin/bash

DATAPATH="footprinting_data"

rm runtime.txt

for file; do
    grep "Number of target markers" $DATAPATH/${file}
    grep "Searching for KASP candidates using" $DATAPATH/${file}
    grep "Total time elapsed" $DATAPATH/${file}
done > runtime.txt

cut -f 6 -d "-" runtime.txt > runtime.temp
#grep -v ".log" runtime.temp > runtime.txt
mv runtime.temp runtime.txt
#rm runtime.temp
