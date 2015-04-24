#!/bin/bash

if [ $# -ne 3 ]
then
  echo "syntax: $0 [redshift] [annulus file] [N_cores]"
  exit 1
fi

zlong=`printf "%1.7f\n" $1`

n=0

mkdir -p model/conc

for (( m=1300; m<=1600; m++ ))
do
if [ ! -f model/conc/conc_m${m}_${zlong}.fits ]
then
  ./src/resample_conc $1 $2 $m model/conc/conc_m${m}_${zlong}.fits &
fi

n=`expr $n + 1`
if [ $n -ge $3 ]
then
 wait
 n=0
fi
done
