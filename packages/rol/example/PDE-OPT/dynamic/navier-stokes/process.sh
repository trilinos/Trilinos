#!/bin/bash

# Read user input file name
NAME="$1";

# Get total number of files with name $NAME.*.txt
NUMFILES=$(ls -dq $NAME.*.txt | wc -l);
MAX=$((NUMFILES-1));

# Concatenate all files horizontally
eval paste "$NAME".{0..$MAX}.txt > "$NAME".txt;

# Remove first two lines of file
sed -i '1,2d' "$NAME".txt;

# Copy map file
cat map_"$NAME".0.txt > map_"$NAME".txt;
sed -i '1,9d' map_"$NAME".txt;

# Clean up
rm $NAME.*.txt;
rm map_"$NAME".*.txt
