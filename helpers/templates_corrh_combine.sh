#!/bin/bash

if [ $# -ne 1 ]
then
  echo "syntax: $0 [redshift]"
  exit 1
fi

zlong=`printf "%1.7f\n" $1`

if [ ! -f templates/corrh_${zlong}.fits ]
then
  ./src/template_corrh_combine $1 templates/corrh_${zlong}.fits templates/corrh/cov_${zlong}_
fi

