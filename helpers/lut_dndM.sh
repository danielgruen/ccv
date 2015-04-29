#!/bin/bash

if [ $# -lt 1 ]
then
  echo "syntax: $0 [redshift(s)]"
  exit 1
fi

mkdir -p lut/

while test ${#} -gt 0
do
  zlong=`printf "%1.7f\n" $1`

  if [ ! -s lut/dndM${zlong}.tab ]
  then
    ./src/mktinkerconf $zlong > src/tinker/dndm.conf
    ./src/tinker/massfunction.x src/tinker/dndm.conf
    rm -f hod-usedvalues 
    mv tmp.dndM lut/dndM${zlong}.tab
  fi

  shift
done

