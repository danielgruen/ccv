#!/bin/bash

if [ $# -lt 1 ]
then
  echo "syntax: $0 [redshift(s)]"
  exit 1
fi

mkdir -p lut/

# iterate over redshifts
while test ${#} -gt 0
do

  zlong=`printf "%1.7f\n" $1`

  if [ ! -s lut/W${zlong}.tab ]
  then

    URL=http://www.usm.uni-muenchen.de/~dgruen/code/templates/W${zlong}.tab

    curl -s --head $URL | head -n 1 | grep "HTTP/1.[01] [23].." > /dev/null
    if [ $? -eq 0 ]
    then
      echo "I am a lazy script and will download W lut rather than calculating it"
      curl $URL > lut/W${zlong}.tab
      continue
    fi

    ./src/calc_W $zlong > lut/W${zlong}.tab
  fi
  
  shift

done

