#!/bin/bash

rm example.txt

gnuplot <<_EOF_
set dgrid3d 10,10
set style data lines
splot "example0.txt" using 3:4:5
_EOF_