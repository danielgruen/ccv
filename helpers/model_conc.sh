#!/bin/bash

if [ $# -lt 3 ]
then
  echo "syntax: $0 [redshift] [annulus file] [N_cores]"
  exit 1
fi


mkdir -p model/conc
mkdir -p model/gamma/conc
n=0

annuli="${@:(-2):1}"
cores=${@: -1}

while test ${#} -gt 2
do

  zlong=`printf "%1.7f\n" $1`

  for (( m=1300; m<=1600; m++ ))
  do
    if [ ! -s model/conc/conc_m${m}_${zlong}.fits ]
    then
      rm -f model/conc/conc_m${m}_${zlong}.fits
      echo ./src/resample_conc $1 $annuli $m model/conc/conc_m${m}_${zlong}.fits 
      ./src/resample_conc $1 $annuli $m model/conc/conc_m${m}_${zlong}.fits &
      n=`expr $n + 1`
    fi

    if [ ! -s model/gamma/conc/conc_m${m}_${zlong}.fits ]
    then
      rm -f model/gamma/conc/conc_m${m}_${zlong}.fits
      echo ./src/resample_conc_g $1 $annuli $m model/gamma/conc/conc_m${m}_${zlong}.fits
      ./src/resample_conc_g $1 $annuli $m model/gamma/conc/conc_m${m}_${zlong}.fits &
      n=`expr $n + 1`
    fi

    if [ $n -ge $cores ]
    then
      wait
      n=0
    fi
  done

  shift
done

wait
