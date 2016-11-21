#!/bin/bash

if [ $# -lt 1 ]
then
  echo "syntax: $0 [pz_prefix(es)]"
  exit 1
fi

for PZFILE in "$@"
do
  if [ ! -f ./templates/lss_$PZFILE.fits ]
  then
    echo ./src/template_lss lut/Pkappa_$PZFILE.tab templates/lss_$PZFILE.fits
    ./src/template_lss lut/Pkappa_$PZFILE.tab templates/lss_$PZFILE.fits
  else
    echo "./templates/lss_$PZFILE.fits exists"
  fi
done

