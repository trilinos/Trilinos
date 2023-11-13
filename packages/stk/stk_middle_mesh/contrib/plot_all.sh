#!/bin/bash

# $1 can (optionally) be a suffix for the first mesh file name
# $2 can (optionallY) be a suffix for the second mesh file name

pth="$HOME"/build/triangulator/
$pth/contrib/plot_c.py mesh"$1"_initial &
$pth/contrib/plot_c.py mesh"$1"_snapped &
$pth/contrib/plot_c.py mesh"$1"_fixed &

$pth/contrib/plot_c.py mesh"$2"_initial &
$pth/contrib/plot_c.py mesh"$2"_snapped &
$pth/contrib/plot_c.py mesh"$2"_fixed &


$pth/contrib/plot_c2.py mesh"$1"_initial mesh"$2"_initial &
$pth/contrib/plot_c2.py mesh"$1"_snapped mesh"$2"_snapped &
$pth/contrib/plot_c2.py mesh"$1"_fixed mesh"$2"_fixed &

wait
