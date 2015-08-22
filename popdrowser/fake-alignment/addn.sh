#!/bin/sh
for i in `ls *.fasta`
do
        sed '/^$/d' $i > $i.temp && mv $i.temp $i
done
exit 0
