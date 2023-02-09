#!/bin/bash

nsteps=4

pth="$HOME"/build/triangulator/

$pth/contrib/plot_c.py mesh2_initial
for i in `seq 0 4`
do
  $pth/contrib/plot_c.py mesh2_after_snap"$i" &
  $pth/contrib/plot_c.py mesh2_after_smooth"$i" &
done

wait
