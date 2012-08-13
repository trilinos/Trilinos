#!/bin/sh

startDir="/Users/kddevin/code/Trilinos/preCopyrightTrilinos/zoltan2"
theFiles=`find ${startDir} -name "*.[ch]pp" -print`

for origf in $theFiles; do
  # pull out file name, removing all of the path ahead of it
  tt=`echo ${origf} | sed "s/\// /g" | awk '{ print $NF }'`
  tmpfile="/tmp/${tt}_tmp"
  rm -f ${tmpfile}
  cat ${startDir}/doc/COPYRIGHT_AND_LICENSE $origf > $tmpfile
  mv -f ${tmpfile} ${origf}
  rm -f ${tmpfile}
done
