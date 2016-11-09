#!/bin/bash

for PZFILE in "$@"
do
  if [ ! -f ./templates/lss_$PZFILE.fits ]
  then
    ./src/template_lss lut/Pkappa_$PZFILE.tab templates/lss_$PZFILE.fits
  fi
done

