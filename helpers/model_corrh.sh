#!/bin/bash

if [ $# -ne 2 ]
then
  echo "syntax: $0 [redshift] [annulus file]"
  exit 1
fi

zlong=`printf "%1.7f\n" $1`

./src/resample_corrh $1 $2 templates/corrh_${zlong}.fits model/corrh_${zlong}.fits
