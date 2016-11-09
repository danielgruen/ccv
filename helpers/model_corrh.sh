#!/bin/bash

if [ $# -lt 2 ]
then
  echo "syntax: $0 [p(z), annulus file](s)"
  exit 1
fi

mkdir -p model
mkdir -p model/gamma


while test ${#} -gt 1
do
  zlong=`printf "%1.7f\n" $1`
  shift
  annuli=$1

  echo "resampling for z=$zlong annuli=$annuli"

  if [ ! -s model/corrh_${annuli}_${zlong}.fits ]
  then
    rm -f model/corrh_${annuli}_${zlong}.fits
    ./src/resample_corrh ${annuli}.tab templates/corrh_${zlong}.fits model/corrh_${annuli}_${zlong}.fits
  fi
  if [ ! -s model/gamma/corrh_${annuli}_${zlong}.fits ]
  then
    rm -f model/gamma/corrh_${annuli}_${zlong}.fits
    ./src/resample_corrh_g ${annuli}.tab templates/corrh_${zlong}.fits model/gamma/corrh_${annuli}_${zlong}.fits
  fi
 
  shift
done

