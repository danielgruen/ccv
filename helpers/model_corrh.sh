#!/bin/bash

if [ $# -lt 2 ]
then
  echo "syntax: $0 [redshift(s)] [annulus file]"
  exit 1
fi

mkdir -p model
mkdir -p model/gamma

annuli=${@: -1}

while test ${#} -gt 1
do
  zlong=`printf "%1.7f\n" $1`

  if [ ! -s model/corrh_${zlong}.fits ]
  then
    rm -f model/corrh_${zlong}.fits
    ./src/resample_corrh $1 $annuli templates/corrh_${zlong}.fits model/corrh_${zlong}.fits
  fi
  if [ ! -s model/gamma/corrh_${zlong}.fits ]
  then
    rm -f model/gamma/corrh_${zlong}.fits
    ./src/resample_corrh_g $1 $annuli templates/corrh_${zlong}.fits model/gamma/corrh_${zlong}.fits
  fi
 
  shift
done

