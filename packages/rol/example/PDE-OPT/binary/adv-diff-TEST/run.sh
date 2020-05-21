#!/bin/bash

size="$1x$2"
file="output_$1x$2.txt"

sed -i '37s/.*/    <Parameter name="NX"                         type="int"     value="'$1'"    \/>/' input_ex01.xml
sed -i '38s/.*/    <Parameter name="NY"                         type="int"     value="'$2'"    \/>/' input_ex01.xml

./ROL_example_PDE-OPT_binary_adv-diff-TEST_Binary.exe --maxWallMinutes=600 --relTolerance=1e-6 $size |tee $file

mv control.txt control_$1x$2.txt
mv map_control.txt map_control_$1x$2.txt

mv state.txt state_$1x$2.txt
mv map_state.txt map_state_$1x$2.txt

mv nodes.txt nodes_$1x$2.txt
mv cell_to_node_quad.txt cell_to_node_quad_$1x$2.txt

plot="${3:-1000}"
if [ $plot -eq "1" ]
then
  matlab -nodisplay -r "plotresultsC0($1,$2); quit"; stty sane
  epstopdf 'control_'$1'x'$2'.eps'
  xdg-open 'control_'$1'x'$2'.pdf' &
fi
