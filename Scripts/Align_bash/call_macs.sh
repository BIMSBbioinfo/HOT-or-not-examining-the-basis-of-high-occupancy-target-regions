#!/bin/bash
#$ -pe smp 5
#$ -l h_vmem=5G

signal=$1 #run id
control=$2
out=$3
log=$4


if [ $control = "/data/akalin/Base/KO_TranscriptonFactors//CONTROLS/NA/MAPPED/NA.bed" ]; then
              # without control
              ~/.guix-profile/bin/macs2 callpeak  --nomodel --extsize 200 -g 1.87e9 -t $signal --outdir $out 2> $log 
            else
              # with control
              ~/.guix-profile/bin/macs2 callpeak --nomodel --extsize 200 -g 1.87e9 -t $signal -c $control --outdir $out 2> $log 

            fi




























