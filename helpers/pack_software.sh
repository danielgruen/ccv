#!/bin/bash

if [ $# -ne 1 ]
then
  echo "syntax: $0 [version number]"
  exit 1
fi

git commit -a
dir=`pwd`

mkdir -p ~/tmp
cd ~/tmp
git clone $dir
rm -rf ccv/.git
tar czf ccv-${1}.tar.gz ccv
mv ccv-${1}.tar.gz $dir/pub
rm -rf ccv
