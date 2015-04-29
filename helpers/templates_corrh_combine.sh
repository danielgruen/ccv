#!/bin/bash

if [ $# -lt 1 ]
then
  echo "syntax: $0 [redshift(s)]"
  exit 1
fi


while test ${#} -gt 0
do
  zlong=`printf "%1.7f\n" $1`

  if [ ! -s templates/corrh_${zlong}.fits ]
  then
    rm -f templates/corrh_${zlong}.fits
    ./src/template_corrh_combine $1 templates/corrh_${zlong}.fits templates/corrh/cov_${zlong}_
  fi

  shift
done

