#!/bin/bash

if [ $# -ne 1 ]
then
  echo "syntax: $0 [N_cores]"
  exit 1
fi


if [ ! -s templates/conc/cov_100.fits ]
then

  URL=http://www.usm.uni-muenchen.de/~dgruen/code/templates/templates_conc.tar.gz

  curl -s --head $URL | head -n 1 | grep "HTTP/1.[01] [23].." > /dev/null
  if [ $? -eq 0 ]
  then
    echo "I am a lazy script and will download conc templates rather than calculating them"
    curl $URL > pub/templates_conc.tar.gz
    tar xzf pub/templates_conc.tar.gz
    exit 0
  fi
fi


i=0

for (( ic=20; ic<=100; ic+=1 ))
do
  if [ ! -s templates/conc/cov_${ic}.fits ]
  then
    rm -f templates/conc/cov_${ic}.fits
    c=`calc.sh $ic*0.1`
    ./src/template_conc $c templates/conc/cov_${ic}.fits &
    i=`expr $i + 1`
    if [ $i -eq $1 ]
    then
      wait
      i=0
    fi
  fi
done

wait
