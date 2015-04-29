#!/bin/bash

if [ $# -lt 1 ]
then
  echo "syntax: $0 [redshift(s)]"
  exit 1
fi

mkdir -p lut/

while test ${#} -gt 0
do
  zlong=`printf "%1.7f\n" $1`

  if [ ! -s lut/dndM${zlong}.tab ]
  then

    URL=http://www.usm.uni-muenchen.de/~dgruen/code/templates/dndM${zlong}.tab

    curl -s --head $URL | head -n 1 | grep "HTTP/1.[01] [23].." > /dev/null
    if [ $? -eq 0 ]
    then
      echo "I am a lazy script and will download dndM lut rather than calculating it"
      curl $URL > lut/dndM${zlong}.tab
      continue
    fi

    ./src/mktinkerconf $zlong > src/tinker/dndm.conf
    ./src/tinker/massfunction.x src/tinker/dndm.conf
    rm -f hod-usedvalues 
    mv tmp.dndM lut/dndM${zlong}.tab
  fi

  shift
done

