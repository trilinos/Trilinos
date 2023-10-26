#!/bin/bash

# $1 can (opionally) be a suffix for the file names

pth="$HOME"/build/triangulator/
$pth/contrib/plot.py "mesh$1" &
$pth/contrib/plot.py "mesh_fixed$1" &
$pth/contrib/plot_c.py "mesh_constraints$1" &
$pth/contrib/plot_c2.py "mesh1_element$1" "mesh2_elements$1" &
