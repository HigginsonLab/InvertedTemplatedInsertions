#!/usr/bin/env bash

#For each vcf file run TINS search
for vcf in *.vcf.gz; do
#Extract sample name from full filename
    NAME="$(basename $vcf .consensus.20161006.somatic.indel.vcf.gz)"
    NAME="${NAME:7}"
    #Only look at the insertions in the vcf file
    cat $vcf | perl -ane '$x=0;for $y (split(",",$F[4])){$x=1 if length($y)>length($F[3])}print if /^#/||$x' | python tins_from_vcf.py --genome genome.fa --vcf $vcf > $NAME.vcf.txt &
    wait
    #Extract the indel count from stats file (vcf-stats)
    sed -n -e 's/^.*number of indels://p' $NAME.stats >> $NAME.vcf.txt
done
