#!/bin/bash

BASEDIR="$HOME/polyoligo/webapp/uploads"
TRASH="$HOME/polyoligo/webapp/trash"
LOGDIR="$HOME/polyoligo/webapp/log.txt"

mkdir $TRASH
find $BASEDIR -mtime +1 -type d | xargs -I {} mv {} $TRASH

for file in $TRASH/*/*.log; do
    head -n 1 ${file} >> $LOGDIR
    grep "- Number of" ${file} >> $LOGDIR
done

rm -r $TRASH
