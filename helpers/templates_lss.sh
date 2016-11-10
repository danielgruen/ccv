#!/bin/bash

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

