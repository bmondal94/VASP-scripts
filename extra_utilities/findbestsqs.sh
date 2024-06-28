#!/bin/sh
stem=$(grep -n Objective bestcorr?*.out | sort -n -k 2 | head -1 | sed 's/:.*$//g')
if [ "${stem}" == "" ]; then
  echo No output yet.
  exit 1
fi
cp $stem bestcorr.out
cp $(echo $stem | sed 's/corr/sqs/g') bestsqs.out
#cat bestcorr.out
