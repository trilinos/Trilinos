#!/bin/bash

# Read user input file name
NAME="$1";

# Get total number of files with name $NAME.*.txt
NUMFILES=$(ls -dq $NAME.*.txt | wc -l);
MAX=$((NUMFILES-1));

# Concatenate all files horizontally
eval paste "$NAME".{0..$MAX}.txt > "$NAME".txt;

# Clean up
rm $NAME.*.txt;
