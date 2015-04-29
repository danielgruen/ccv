#!/bin/bash

if [ $# -lt 2 ]
then
  echo "syntax: $0 [redshift(s)] [N_cores]"
  exit 1
fi

n=0
ncores=${@: -1}

while test ${#} -gt 1
do

  zlong=`printf "%1.7f\n" $1`

  for (( m=80; m<=159; m++ ))
  do
    if [ ! -s templates/corrh/cov_${zlong}_$m.fits ]
    then
      rm -f templates/corrh/cov_${zlong}_$m.fits
      ./src/template_corrh $1 `bash helpers/calc.sh $m/10` templates/corrh/cov_${zlong}_$m.fits &
      n=`expr $n + 1`
    fi

    if [ $n -ge $ncores ]
    then
      echo "all $ncores cores busy; waiting" 
      echo
      wait
      n=0
    fi
  done

  shift
done

wait
