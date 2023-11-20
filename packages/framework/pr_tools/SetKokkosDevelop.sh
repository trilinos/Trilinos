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

ln -s "$PACKAGESPATH/kokkos" "${WORKSPACE:?}/Trilinos/packages/kokkos"

cmake ${WORKSPACE:?}/Trilinos -DKokkos_SOURCE_DIR_OVERRIDE:STRING=kokkos

ln -s "$PACKAGESPATH/kokkos-kernels" "${WORKSPACE:?}/Trilinos/packages/kokkos-kernels"

cmake ${WORKSPACE:?}/Trilinos -DKokkos_SOURCE_DIR_OVERRIDE:STRING=kokkos-kernels

export CXX=$PACKAGESPATH/kokkos/bin/nvcc_wrapper

export OMPI_CXX=$PACKAGESPATH/kokkos-kernels/bin/nvcc_wrapper