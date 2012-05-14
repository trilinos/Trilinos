Some of the coordinate files were obtained from the
University of Florida sparse matrix collection.

Those listed below were converted from Zoltan
test input files.

To compare to Zoltan nightly tests:

brack2_3_coord.mtx:
  5 processes

bug_coord.mtx
  3 processes

degenerateAA_coord.mtx
  6 processes

degenerate_coord.mtx
  6 processes
  rectilinear_blocks

ewgt_coord.mtx
  4 processes

drake_coord.mtx
  3 processes

grid20x19_coord.mtx
  4 processes
  average cuts
  rectilinear_blocks

hammond_coord.mtx
  8 processes

hammond2_coord.mtx
  6 processes

nograph_coord.mtx
  4 processes

onedbug_coord.mtx
  3 processes

simple_coord.mtx
  4 processes
  average cuts

vwgt_coord.mtx
  3 processes

vwgt2_coord.mtx
  2 processes

To convert from Zoltan format to UF format:

#!/bin/bash
#
#  The script to convert Chaco coordinates to "mtx" coordinates
#  is in the Zoltan_Tests/converters repository.
#

fnames="degenerateAA ewgt hammond simple3d vwgt bug degenerate grid20x19 nograph simple drake hammond2 onedbug vwgt2"
suffix=_coord.mtx

for fname in $fnames ; do
  echo "chacoToDavisCoords.py $fname.coords  $fname$suffix"
  chacoToDavisCoords.py $fname.coords  $fname$suffix
done
