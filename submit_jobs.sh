#!/bin/bash

for xtal in COF-102.cif COF-103.cif
do
    echo "submitting job for $xtal"
    qsub -N $xtal -o "$xtal.o" -e "$xtal.e" -v xtal="$xtal" gcmc_submit.sh
done
