#!/bin/bash

if [ $# -ne 1 ]
then
  echo "syntax: $0 [N_cores]"
  exit 1
fi

mkdir -p pub
mkdir -p templates/off

#if [ ! -s templates/off/cov_100.fits ]
#then

#  URL=http://www.usm.uni-muenchen.de/~dgruen/code/templates/templates_ell.tar.gz

#  curl -s --head $URL | head -n 1 | grep "HTTP/1.[01] [23].." > /dev/null
#  if [ $? -eq 0 ]
#  then    
#    echo "I am a lazy script and will download ell templates rather than calculating them"
#    curl $URL > pub/templates_ell.tar.gz
#    tar xzf pub/templates_ell.tar.gz
#    exit 0
#  fi
#fi



i=0

for (( ic=1; ic<=200; ic+=1 ))
do
  if [ ! -s templates/off/cov_${ic}.fits ]
  then
    rm -f templates/off/cov_${ic}.fits
    c=`calc.sh $ic*0.01`
    echo ./src/template_off $c templates/off/cov_${ic}.fits
    ./src/template_off $c templates/off/cov_${ic}.fits &
    i=`expr $i + 1`
    if [ $i -eq $1 ]
    then
      wait
      i=0
    fi
  fi
done

wait
