#!/bin/bash

froot="somefile"    ## root without exodus extension
fext1="e"           ## exodus extension
fext2="txt"         ## new extension
fnum=64             ## number of files

digits="${#fnum}"
ct=0
while (($ct < $fnum))
do
  printf -v ctpad "%0${digits}d" $ct
  finput=$froot.$fext1.$fnum.$ctpad   ## integers assumed padded with zeros
  foutput=$froot.$fext2.$fnum.$ct     ## integers not padded with zeros
  echo "Running command:  ncdump" $finput ">" $foutput
  ncdump $finput > $foutput
  ct=$(( ct+1 ))
done

