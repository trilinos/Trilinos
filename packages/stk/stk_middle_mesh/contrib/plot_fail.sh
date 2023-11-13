#!/bin/bash

../contrib/plot_c.py "mesh1_initial" &
../contrib/plot_c.py "mesh2_initial" &

../contrib/plot_c.py "mesh1_snapped" &
../contrib/plot_c.py "mesh2_snapped" &

../contrib/plot_c.py "mesh1_fixed" &
../contrib/plot_c.py "mesh2_fixed" &

for i in `seq 0 99`;
do
  ../contrib/plot_c.py "mesh_invalid_laplacian$i" &
  ../contrib/plot_c.py "mesh_invalid_quality$i" &
done

../contrib/plot_c2.py "mesh1_initial" "mesh2_initial" &
../contrib/plot_c2.py "mesh1_snapped" "mesh2_snapped" &
../contrib/plot_c2.py "mesh1_fixed" "mesh2_fixed" &

wait
