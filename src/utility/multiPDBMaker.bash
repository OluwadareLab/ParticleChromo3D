#!/bin/bash

chrNum=1
infile=../../input-and-models/Input/GM12878_input/KR_500kb/chr${chrNum}_matrix.txt

for i in {1..30}
do
   python3 ../Ps.py ${infile} -o outputfolder/chr${chrNum}-${i}
done
