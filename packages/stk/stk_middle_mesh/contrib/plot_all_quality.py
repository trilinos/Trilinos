#!/bin/bash

pth="$HOME"/build/triangulator/
for i in `seq 1 17`
do
  "$pth"/contrib/plot_quality.py mesh"$i"_vert_qualities.txt &
done

wait
