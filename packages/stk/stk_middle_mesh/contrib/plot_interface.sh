#!/bin/bash

#../contrib/plot_c.py mesh1_element &
#../contrib/plot_c.py mesh2_elements &
#../contrib/plot_c.py mesh_in &
#../contrib/plot_c2.py mesh1_element mesh2_elements &

../contrib/plot_c.py mesh_fixed &
../contrib/plot_c.py mesh_in_projected &
../contrib/plot_c_3d.py mesh_in3d &
../contrib/plot_c_3d_both.py mesh1_element3d mesh2_elements3d &
../contrib/plot_c_3d.py mesh1_snapped &
../contrib/plot_c_3d.py mesh2_snapped

wait
