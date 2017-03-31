#!/bin/bash

if [ $# -lt 3 ]
then
  echo "syntax: $0 ([redshift] [annulus file])s [N_cores]"
  exit 1
fi


mkdir -p model/gamma/conc_g
n=0

cores=${@: -1}

while test ${#} -gt 2
do

  zlong=`printf "%1.7f\n" $1`
  shift
  annuli=$1

  for (( m=1300; m<=1600; m++ ))
  do
    if [ ! -s model/gamma/conc_g/conc_${annuli}_m${m}_${zlong}.fits ]
    then
      rm -f model/gamma/conc_g/conc_${annuli}_m${m}_${zlong}.fits
      echo ./src/resample_conc_g_test $zlong $annuli.tab $m model/gamma/conc_g/conc_${annuli}_m${m}_${zlong}.fits
      ./src/resample_conc_g_test $zlong $annuli.tab $m model/gamma/conc_g/conc_${annuli}_m${m}_${zlong}.fits &
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