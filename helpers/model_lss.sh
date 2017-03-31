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
  pz=$1
  shift
  annuli=$1

  echo "resampling for p(z)=$pz annuli=$annuli"

  if [ ! -s model/lss_${pz}_${annuli}.fits ]
  then
    rm -f model/lss_${pz}_${annuli}.fits
    echo ./src/resample_corrh ${annuli}.tab templates/lss_${pz}.fits model/lss_${pz}_${annuli}.fits
    ./src/resample_corrh ${annuli}.tab templates/lss_${pz}.fits model/lss_${pz}_${annuli}.fits
  fi
  if [ ! -s model/gamma/lss_${pz}_${annuli}.fits ]
  then
    rm -f model/gamma/lss_${pz}_${annuli}.fits
    echo ./src/resample_corrh_g ${annuli}.tab templates/lss_${pz}.fits model/gamma/lss_${pz}_${annuli}.fits
    ./src/resample_corrh_g ${annuli}.tab templates/lss_${pz}.fits model/gamma/lss_${pz}_${annuli}.fits
  fi
 
  shift
done

