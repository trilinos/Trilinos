#!/usr/bin/env bash
SCRIPTFILE=$(realpath ${WORKSPACE:?}/Trilinos/packages/framework/pr_tools/SetKokkosDevelop.sh)
SCRIPTPATH=$(dirname $SCRIPTFILE)
source ${SCRIPTPATH:?}/common.bash
PACKAGESPATH=$(realpath ${WORKSPACE:?}/Trilinos/packages)

# Ensures git is loaded properly
if ! command -v git &> /dev/null; then
    message_std "SetKokkosDevelop> ERROR: git is not available"
    exit 1
fi

cd $PACKAGESPATH
rm -rf kokkos kokkos-kernels
git clone --depth=1 --single-branch --branch=develop --shallow-submodules https://github.com/kokkos/kokkos.git
git clone --depth=1 --single-branch --branch=develop --shallow-submodules https://github.com/kokkos/kokkos-kernels.git
message_std "SetKokkosDevelop> INFO: updated kokkos and kokkos-kernels packages with current develop"

# Returns to previous path from before running this script
cd -
