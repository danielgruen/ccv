#!/bin/bash

if [ $# -lt 1 ]
then
  echo "syntax: $0 [pz_prefix(es)]"
  exit 1
fi

mkdir -p lut

for PZFILE in "$@"
do
  if [ ! -s lut/Pkappa_$PZFILE.tab ]
  then
    echo "calculating Pkappa for $PZFILE"
    ./src/mknicaeaconf $PZFILE.tab
    cp src/nicaea_2.5/Demo/cosmo_lens.par .
    ./src/nicaea_2.5/Demo/pkappa
    grep -v "#" P_kappa > lut/Pkappa_$PZFILE.tab
    rm -f pkappa lensingdemo cosmo_lens.par cosmo.par nofz.par P_kappa
  fi
done

