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

  if [ -s templates/corrh_${zlong}.fits ]
  then
    continue
  fi 

  URL=http://www.usm.uni-muenchen.de/~dgruen/code/templates/corrh_${zlong}.fits

  curl -s --head $URL | head -n 1 | grep "HTTP/1.[01] [23].." > /dev/null
  if [ $? -eq 0 ]
  then    
    echo "I am a lazy script and will download corrh template rather than calculating it"
    curl $URL > templates/corrh_${zlong}.fits
    continue
  fi


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
