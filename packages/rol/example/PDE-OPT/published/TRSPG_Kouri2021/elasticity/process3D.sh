#!/bin/bash

# Grab number of nodes
NN=$(cat nodes.txt | wc -l)
NN0=$((NN * 3 + 2))

# Replace double space with single space in nodes
cat nodes.txt | tr -s '[:space:]' > nodes-mod.txt

# Extract nodal state data
sed -n "3~3p;$NN0 q" state.txt > state-dispx.txt
sed -n "4~3p;$NN0 q" state.txt > state-dispy.txt
sed -n "5~3p;$NN0 q" state.txt > state-dispz.txt

# Add nodal state data to nodes
paste -d'\0    ' nodes-mod.txt state-dispx.txt state-dispy.txt state-dispz.txt > results-state-0.txt

# Extract nodal control data
sed -n "3~3p;$NN0 q" density.txt > density-0.txt

# Add nodal state data to nodes
paste -d'\0    ' nodes-mod.txt density-0.txt > results-density-0.txt

# Add header to state and control results files
echo 'X                      Y                      Z                      DX                      DY                      DZ' > header-state.txt
cat header-state.txt results-state-0.txt > results-state.txt
echo 'X                      Y                      Z                      D' > header-density.txt
cat header-density.txt results-density-0.txt > results-density.txt

# Remove temporary files
rm state-dispx.txt state-dispy.txt state-dispz.txt density-0.txt
rm results-state-0.txt results-density-0.txt nodes-mod.txt
rm header-state.txt header-density.txt
