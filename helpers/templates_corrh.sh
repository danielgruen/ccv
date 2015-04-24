#!/bin/bash

if [ $# -ne 2 ]
then
  echo "syntax: $0 [redshift] [N_cores]"
  exit 1
fi

zlong=`printf "%1.7f\n" $1`

n=0

for (( m=80; m<=159; m++ ))
do
if [ ! -f templates/corrh/cov_${zlong}_$m.fits ]
then
  ./src/template_corrh $1 `bash helpers/calc.sh $m/10` templates/corrh/cov_${zlong}_$m.fits &
fi

n=`expr $n + 1`
if [ $n -ge $2 ]
then
 wait
 n=0
fi
done
