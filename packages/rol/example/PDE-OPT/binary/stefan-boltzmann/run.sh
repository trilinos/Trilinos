#!/bin/bash

size="$1x$1"
file="output_$1x$1.txt"

sed -i '37s/.*/    <Parameter name="NX"                         type="int"     value="'$1'"    \/>/' input_ex01.xml
sed -i '38s/.*/    <Parameter name="NY"                         type="int"     value="'$1'"    \/>/' input_ex01.xml

./ROL_example_PDE-OPT_binary_adv-diff-TEST_Binary.exe --maxWallMinutes=600 --relTolerance=1e-6 $size |tee $file

mv control.txt control_$size.txt
mv map_control.txt map_control_$size.txt

mv state.txt state_$size.txt
mv map_state.txt map_state_$size.txt

mv nodes.txt nodes_$size.txt
mv cell_to_node_quad.txt cell_to_node_quad_$size.txt

plot="${2:-1000}"
if [ $plot -eq "1" ]
then
  matlab -nodisplay -r "plotresultsC0($1); quit"; stty sane
  epstopdf 'control_'$1'x'$1'.eps'
  xdg-open 'control_'$1'x'$1'.pdf' &
fi
