#!/bin/sh
#
# Syntax:
#
#    toMatlab.sh inputFile outputfile.m
#
# Purpose:
#
#    Prepare output for analysis in Matlab.  Input is raw output from Scaling
#    example.  Output is a cell array containing two columns of data.
#    First column is labels, second column is average times, e.g.,
#        'SaPFactory:APtent1'           [0.0335181]
#        'SaPFactory:Dinv_APtent_1'     [0.004318]
#        'SaPFactory:eigen_estimate_1'  [0.00355911]
#
# Note: This script was written with the assumption that timing data lines
#       are of the form
#
#          &&&TimerName max=1.23 min=1.23 avg=1.23
#

if [ $# -ne 2 ]
then
  echo "usage: toMatlab.sh inputfile outputfile"
else
  if [ -f $2 ]
  then
    echo "output file '$2' exists"
  else
    echo "data = {..." > $2
  
    # stick single quotes around first field, end line with ";..."
    cat $1 | egrep "^&&&" | egrep -v "MemFree" | sed "s/^&&&//" | sed "s/avg=//" | awk -v qt="'" '{printf "%s%s%s %s; ...\n",qt,$1,qt,$4}' >> $2
  
    echo "};" >> $2
  fi
fi
