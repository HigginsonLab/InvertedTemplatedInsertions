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

#Format for resulting txt files
format="%s\t%i\t%i\t%i\n"
printf "Sample_ID\tTINS_Raw\tTotal_Indels\tLength\n" > iTINS_size.txt
for v in *.vcf.txt; do
#Extract sample name, raw iTINS count, and total indels
    NAME="$(basename $v .vcf.txt)"
    TINS_Raw="$(tail -3 $v | head -1)"
    Total_Indels="$(tail -1 $v)"
    #If iTINS found, extract insertion length, else return 0
    if grep "iTINS Insertion Length" $v; then
        cat $v | grep "iTINS Insertion Length" | while read line; do
            Length="$(cut -d$':' -f2 <<< $line)"
            printf "$format" $NAME $TINS_Raw $Total_Indels $Length >> iTINS_size.txt
        done
    else
        Length=0
        printf "$format" $NAME $TINS_Raw $Total_Indels $Length >> iTINS_size.txt
    fi
done

#Format for resulting txt files
format="%s\t%i\t%i\t%i\n"
printf "Sample_ID\tTINS_Raw\tTotal_Indels\tLength\n" > uTINS_size.txt
for v in *.vcf.txt; do
#Extract sample name, raw uTINS count, and total indels
    NAME="$(basename $v .vcf.txt)"
    TINS_Raw="$(tail -4 $v | head -1)"
    Total_Indels="$(tail -1 $v)"
    #If uTINS found, extract insertion length, else return 0
    if grep "uTINS Insertion Length" $v; then
        cat $v | grep "uTINS Insertion Length" | while read line; do
            Length="$(cut -d$':' -f2 <<< $line)"
            printf "$format" $NAME $TINS_Raw $Total_Indels $Length >> uTINS_size.txt
        done
    else
        Length=0
        printf "$format" $NAME $TINS_Raw $Total_Indels $Length >> uTINS_size.txt
    fi
done

#Format for resulting txt files
format="%s\t%i\t%i\t%i\n"
printf "Sample_ID\tTINS_Raw\tTotal_Indels\tPosition\n" > iTINS_Position_Mid.txt
for v in *.vcf.txt; do
#Extract sample name, raw iTINS count, and total indels
    NAME="$(basename $v .vcf.txt)"
    TINS_Raw="$(tail -3 $v | head -1)"
    Total_Indels="$(tail -1 $v)"
    #If iTINS found, extract insertion position, else return 0
    if grep "iTINS Position" $v; then
        cat $v | grep "iTINS Position" | while read line; do
            Position="$(cut -d$':' -f2 <<< $line)"
            printf "$format" $NAME $TINS_Raw $Total_Indels $Position >> iTINS_Position_Mid.txt
        done
    else
        Position=-1
        printf "$format" $NAME $TINS_Raw $Total_Indels $Position >> iTINS_Position_Mid.txt
    fi
done

#Format for resulting txt files
format="%s\t%i\t%i\t%i\n"
printf "Sample_ID\tTINS_Raw\tTotal_Indels\tPosition\n" > uTINS_Position_Mid.txt
for v in *.vcf.txt; do
#Extract sample name, raw uTINS count, and total indels
    NAME="$(basename $v .vcf.txt)"
    TINS_Raw="$(tail -4 $v | head -1)"
    Total_Indels="$(tail -1 $v)"
    #If uTINS found, extract insertion position, else return 0
    if grep "uTINS Position" $v; then
        cat $v | grep "uTINS Position" | while read line; do
            Position="$(cut -d$':' -f2 <<< $line)"
            printf "$format" $NAME $TINS_Raw $Total_Indels $Position >> uTINS_Position_Mid.txt
        done
    else
        Position=-1
        printf "$format" $NAME $TINS_Raw $Total_Indels $Position >> uTINS_Position_Mid.txt
    fi
done

