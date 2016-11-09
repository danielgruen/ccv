#!/bin/bash

for PZFILE in "$@"
do
  if [ ! -f lut/Pkappa_$PZFILE.tab ]
  then
    ./src/mknicaeaconf $PZFILE.tab
    cp src/nicaea_2.5/Demo/cosmo_lens.par .
    ./src/nicaea_2.5/Demo/pkappa
    grep -v "#" P_kappa > lut/Pkappa_$PZFILE.tab
    rm -f pkappa lensingdemo cosmo_lens.par cosmo.par nofz.par P_kappa
  fi
done

