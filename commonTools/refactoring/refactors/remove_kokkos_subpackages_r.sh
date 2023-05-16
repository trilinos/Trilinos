#!/bin/bash
#
# Run this in a subdir to strip out usage of Kokkos subpackages for TriBITS
# packages downstream from the Kokkos package (4/19/2023).
#
# Run as:
#
#   $ cd <some-base-dir>
#   $ <this-script-dir>/remove_kokkos_subpackages_r.sh
#
# NOTE: Do **NOT** run this script on the packages/kokkos/subdirectory itself
# (or undo the changes if you do by accident).
#
# NOTE: Some manual modifications may be desired after running this script
# like removing duplicate if clauses like:
#
#   if (<Package>_ENABLE_KokkosCore AND <Package>_ENABLE_KokkosContainers ...)
#
# becoming:
#
#   if (<Package>_ENABLE_Kokkos AND <Package>_ENABLE_Kokkos ...)
#
# which should be manually refactored to be
#
#   if (<Package>_ENABLE_Kokkos ...)
#
# NOTE: This script can't be in the packages/kokkos/ directory because the
# next Kokos snapshot would delete it and this refactoring is really
# Trilinos-specific and therefore does not need to be in the external kokkos
# repo.
#

_SCRIPT_DIR=`echo $0 | sed "s/\(.*\)\/.*remove_kokkos_subpackages_r.sh/\1/g"`
#echo $_SCRIPT_DIR

refactDir=${_SCRIPT_DIR}/..

kokkosSubpkgs=(Core Containers Algorithms Simd)

for kokkosSubpkg in ${kokkosSubpkgs[@]}; do

  echo
  echo "Replace Kokkos${kokkosSubpkg} with Kokkos in Dependencies.cmake file ..."
  echo

  find . -name Dependencies.cmake -exec ${refactDir}/token-replace.pl Kokkos${kokkosSubpkg} Kokkos {} {} \;

  echo
  echo "Replace Kokkos${kokkosSubpkg}_ with Kokkos_ in CMakeLists.txt and *.cmake files ..."
  echo

  find . \( -name CMakeLists.txt -or -name "*.cmake" \) -exec ${refactDir}/string-replace.pl Kokkos${kokkosSubpkg}_ Kokkos_ {} {} \;

  echo
  echo "Replace _Kokkos${kokkosSubpkg} with _Kokkos in CMakeLists.txt and *.cmake files ..."
  echo

  find . \( -name CMakeLists.txt -or -name "*.cmake" \) -exec ${refactDir}/string-replace.pl _Kokkos${kokkosSubpkg} _Kokkos {} {} \;

  echo
  echo "Replace _KOKKOS${kokkosSubpkg^^} with _KOKKOS in all files ..."
  echo

  find . -type f -exec ${refactDir}/string-replace.pl _KOKKOS${kokkosSubpkg^^} _KOKKOS {} {} \;

done

echo
echo "Replace duplicates of 'Kokkos' in Dependencies.cmake file ..."
echo

find . -name Dependencies.cmake -exec sed -i 's/Kokkos Kokkos\([ )\n]\)/Kokkos\1/g' {} \; 
find . -name Dependencies.cmake -exec sed -i 's/Kokkos Kokkos\([ )\n]\)/Kokkos\1/g' {} \; 
find . -name Dependencies.cmake -exec sed -i 's/Kokkos Kokkos\([ )\n]\)/Kokkos\1/g' {} \; 
