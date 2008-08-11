#!/bin/csh
# see README in this directory for details on how to run this script
#
if( "$1" == "" ) then
  echo "External package name must be passed as argument.  Please see README."
else
  if ( "$2" == "" ) then

    set external_package=$1
    sed "s/sample/$external_package/g" backup/configure > configure
    sed "s/sample/$external_package/g" backup/Makefile.in > Makefile.in

  else

    sed "s/sample/$argv[*]/g" backup/configure.ac > configure.ac
    sed "s/sample/$argv[*]/g" backup/Makefile.am > Makefile.am
    ./bootstrap

  endif

endif
