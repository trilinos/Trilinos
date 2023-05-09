#!/bin/bash
#
# Run this script from the Trilinos/packages/ subdir to remove the usage of Kokkos subpackages.
#
# Usage:
#
#   $ cd Trilinos/packages/
#   $ ../commonTools/refactoring/refactors/remove_kokkos_subpackages_from_trilinos_packages.sh
# 

pkgsToRefactor=$(ls | grep -v "kokkos$" | grep -v framework | grep -v TrilinosInstallTests | grep -v new_package)

for pkg in ${pkgsToRefactor} ; do

  if [[ -d ${pkg} ]] ; then

    echo
    echo "***"
    echo "*** Refactoring ${pkg}"
    echo "***"
    echo

    cd ${pkg}

    time ../../commonTools/refactoring/refactors/remove_kokkos_subpackages_r.sh
    
    cd ..

  fi

done
