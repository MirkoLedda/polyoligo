#!/bin/bash

BASEDIR="$HOME/polyoligo/webapp/uploads"
TRASH="$HOME/polyoligo/webapp/trash"
LOGDIR="$HOME/polyoligo/webapp/log.txt

mkdir $TRASH
find $BASEDIR -mtime +1 | xargs -I {} mv {} $TRASH

for file in $TRASH/*/*.log; do
    grep "Number of target markers" ${file}
    grep "Searching for KASP candidates using" ${file}
    grep "Total time elapsed" ${file}
done >> $LOGDIR

rm -r $TRASH
