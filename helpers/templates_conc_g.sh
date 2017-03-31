#!/bin/bash

if [ $# -ne 1 ]
then
  echo "syntax: $0 [N_cores]"
  exit 1
fi

mkdir -p pub
mkdir -p templates/conc_g

i=0

for (( ic=20; ic<=100; ic+=1 ))
do
  if [ ! -s templates/conc_g/cov_${ic}.fits ]
  then
    rm -f templates/conc_g/cov_${ic}.fits
    c=`calc.sh $ic*0.1`
    ./src/template_conc_g $c templates/conc_g/cov_${ic}.fits &
    i=`expr $i + 1`
    if [ $i -eq $1 ]
    then
      wait
      i=0
    fi
  fi
done

wait
