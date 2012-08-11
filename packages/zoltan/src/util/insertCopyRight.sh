#!/bin/sh

startDir="/Users/kddevin/code/zoltan/"
theFiles=`find ${startDir} -name "*.[ch]" -print`

for origf in $theFiles; do
  # pull out file name, removing all of the path ahead of it
  tt=`echo ${origf} | sed "s/\// /g" | awk '{ print $NF }'`
  echo ${tt}
  tmpfile="/tmp/${tt}_tmp"
  rm -f ${tmpfile}
  cat ${startDir}/COPYRIGHT_AND_LICENSE $origf > $tmpfile
  mv -f ${tmpfile} ${origf}
  rm -f ${tmpfile}
done
