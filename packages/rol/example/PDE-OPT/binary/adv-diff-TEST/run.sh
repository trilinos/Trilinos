#!/bin/bash

size="$1x$2"
file="output_$1x$2.txt"

sed -i '11s/.*/    <Parameter name="Number of X-Cells"           type="int"     value="'$1'"        \/>/' input_ex01.xml
sed -i '12s/.*/    <Parameter name="Number of Y-Cells"           type="int"     value="'$2'"        \/>/' input_ex01.xml

./ROL_example_PDE-OPT_binary_adv-diff-TEST_RelaxedBinary.exe --maxWallMinutes=600 --relTolerance=1e-6 $size |tee $file

plot="${3:-1000}"
if [ $plot -eq "1" ]
then
  matlab -nodisplay -r "plotPEBBL($1,$2); quit"; stty sane
  epstopdf 'control_'$1'x'$2'.eps'
  xdg-open 'control_'$1'x'$2'.pdf' &
fi
