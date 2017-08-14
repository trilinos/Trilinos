#!/bin/bash

NN=$(cat nodes.txt | wc -l)
NN0=$((NN * 4 + 2))
NN1=$((NN0 + 1))

sed -n "3~4p;$NN0 q" state.txt > state-velx.txt
sed -n "4~4p;$NN0 q" state.txt > state-vely.txt
sed -n "5~4p;$NN0 q" state.txt > state-velz.txt
sed -n "6~4p;$NN0 q" state.txt > state-pres.txt
sed -n "7~4p;$NN1 q" state.txt > state-thrm.txt

cat nodes.txt | tr -s '[:space:]' > nodes-mod.txt

paste -d'\0    ' nodes-mod.txt state-velx.txt state-vely.txt state-velz.txt state-pres.txt state-thrm.txt > results-state.txt

sed -n "3~4p;$NN0 q" control.txt > control-velx.txt
sed -n "4~4p;$NN0 q" control.txt > control-vely.txt
sed -n "5~4p;$NN0 q" control.txt > control-velz.txt
sed -n "6~4p;$NN0 q" control.txt > control-pres.txt
sed -n "7~4p;$NN1 q" control.txt > control-thrm.txt

cat nodes.txt | tr -s '[:space:]' > nodes-mod.txt

paste -d'\0    ' nodes-mod.txt control-velx.txt control-vely.txt control-velz.txt control-pres.txt control-thrm.txt > results-control.txt

rm state-velx.txt state-vely.txt state-velz.txt state-pres.txt state-thrm.txt
rm control-velx.txt control-vely.txt control-velz.txt control-pres.txt control-thrm.txt
rm nodes-mod.txt
