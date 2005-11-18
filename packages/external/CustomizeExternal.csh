#!/bin/csh
# see README in this directory for details on how to run this script
#
if( "$1" == "" ) then
  echo "External package name must be passed as argument.  Please see README."
else
  set external_package=$1
  sed "s/sample/$external_package/g" configure.txt > configure
  sed "s/sample/$external_package/g" Makefile.txt > Makefile.in

endif
