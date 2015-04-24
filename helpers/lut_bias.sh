#!/bin/bash

if [ $# -ne 1 ]
then
  echo "syntax: $0 [redshift]"
  exit 1
fi

zlong=`printf "%1.7f\n" $1`

n=0

mkdir -p lut/

if [ ! -f lut/bias${zlong}.tab ]
then
  ./src/calc_bias $zlong > lut/bias${zlong}.tab
fi
