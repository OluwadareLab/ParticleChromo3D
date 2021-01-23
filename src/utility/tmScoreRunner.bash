#!/bin/bash

chrNum=1
infile="../../Results/gm12878/consistency/500kb/chr1/chr1-"
#infile=../../input-and-models/Input/GM12878_input/KR_1mb/chr${chrNum}_matrix.txt

rm tmscore.txt

for i in {1..30}
do
   for j in {1..30}
   do
     if [[ $i -eq $j ]]; then
       break
     fi

     ./TMScore ${infile}${i}.pdb ${infile}${j}.pdb | grep "TM-score    = " | cut -c15-20 >> tmscore.txt

   done
done
