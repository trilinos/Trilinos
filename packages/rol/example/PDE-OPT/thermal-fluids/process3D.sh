#!/bin/bash

# Grab number of nodes
NN=$(cat nodes.txt | wc -l)
NN0=$((NN * 4 + 2))
NN1=$((NN0 + 1))

# Replace double space with single space in nodes
cat nodes.txt | tr -s '[:space:]' > nodes-mod.txt

# Extract nodal state data
sed -n "3~4p;$NN0 q" state.txt > state-velx.txt
sed -n "4~4p;$NN0 q" state.txt > state-vely.txt
sed -n "5~4p;$NN0 q" state.txt > state-velz.txt
sed -n "6~4p;$NN0 q" state.txt > state-pres.txt
sed -n "7~4p;$NN1 q" state.txt > state-thrm.txt

# Add nodal state data to nodes
paste -d'\0    ' nodes-mod.txt state-velx.txt state-vely.txt state-velz.txt state-pres.txt state-thrm.txt > results-state.txt

# Extract nodal control data
sed -n "3~4p;$NN0 q" control.txt > control-velx.txt
sed -n "4~4p;$NN0 q" control.txt > control-vely.txt
sed -n "5~4p;$NN0 q" control.txt > control-velz.txt
sed -n "6~4p;$NN0 q" control.txt > control-pres.txt
sed -n "7~4p;$NN1 q" control.txt > control-thrm.txt

# Add nodal control data to nodes
paste -d'\0    ' nodes-mod.txt control-velx.txt control-vely.txt control-velz.txt control-pres.txt control-thrm.txt > results-control.txt

# Remove temporary files
rm state-velx.txt state-vely.txt state-velz.txt state-pres.txt state-thrm.txt
rm control-velx.txt control-vely.txt control-velz.txt control-pres.txt control-thrm.txt
rm nodes-mod.txt
