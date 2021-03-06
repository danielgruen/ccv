#!/bin/bash

if [ $# -lt 3 ]
then
  echo "syntax: $0 ([redshift] [annulus file])s [N_cores]"
  exit 1
fi

mkdir -p model/ell
mkdir -p model/gamma/ell
n=0

cores=${@: -1}

while test ${#} -gt 2
do

  zlong=`printf "%1.7f\n" $1`
  shift
  annuli=$1

  for (( m=1300; m<=1600; m++ ))
  do
    if [ ! -s model/ell/ell_${annuli}_m${m}_${zlong}.fits ]
    then
      rm -f model/ell/ell_${annuli}_m${m}_${zlong}.fits
      ./src/resample_ell $zlong $annuli.tab $m model/ell/ell_${annuli}_m${m}_${zlong}.fits &
      n=`expr $n + 1`
    fi
    if [ ! -s model/gamma/ell/ell_${annuli}_m${m}_${zlong}.fits ]
    then
      rm -f model/gamma/ell/ell_${annuli}_m${m}_${zlong}.fits
      ./src/resample_ell_g $zlong $annuli.tab $m model/gamma/ell/ell_${annuli}_m${m}_${zlong}.fits &
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

