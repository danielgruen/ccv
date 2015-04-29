#!/bin/bash

if [ $# -lt 1 ]
then
  echo "syntax: $0 [redshift(s)]"
  exit 1
fi

mkdir -p lut/

# iterate over redshifts
while test ${#} -gt 0
do

  zlong=`printf "%1.7f\n" $1`

  if [ ! -s lut/W${zlong}.tab ]
  then
    ./src/calc_W $zlong > lut/W${zlong}.tab
  fi

  shift
done

