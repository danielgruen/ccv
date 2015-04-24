#!/bin/bash

if [ $# -ne 1 ]
then
  echo "syntax: $0 [N_cores]"
  exit 1
fi

i=0

for (( ic=20; ic<=100; ic+=1 ))
do
  if [ ! -f templates/ell/cov_${ic}.fits ]
  then
    c=`calc.sh $ic*0.1`
    ./src/template_ell $c templates/ell/cov_${ic}.fits &
    i=`expr $i + 1`
    if [ $i -eq $1 ]
    then
      wait
      i=0
    fi
  fi
done

wait
