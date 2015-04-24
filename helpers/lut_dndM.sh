#!/bin/bash

if [ $# -ne 1 ]
then
  echo "syntax: $0 [redshift]"
  exit 1
fi

zlong=`printf "%1.7f\n" $1`

n=0

mkdir -p lut/

if [ ! -f lut/dndM${zlong}.tab ]
then
  ./src/mktinkerconf $zlong > src/tinker/dndm.conf
  ./src/tinker/massfunction.x src/tinker/dndm.conf
  rm -f hod-usedvalues 
  mv tmp.dndM lut/dndM${zlong}.tab
fi
