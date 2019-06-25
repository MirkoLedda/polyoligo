#!/usr/bin/env bash

awk -F $'\t' '{print $1" "$2" "$5" "$3" "$4}' $1 > .temp.txt
cat .temp.txt > $1
rm .temp.txt
