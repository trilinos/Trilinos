#!/usr/bin/env bash
SCRIPTFILE=$(realpath ${WORKSPACE:?}/Trilinos/packages/framework/pr_tools/SetKokkosDevelop.sh)
SCRIPTPATH=$(dirname $SCRIPTFILE)
source ${SCRIPTPATH:?}/common.bash
DIR_CONTAINING_TRILINOS=$(realpath ${WORKSPACE:?})
TRILINOS_SRC=${DIR_CONTAINING_TRILINOS}/Trilinos

# Ensures git is loaded properly
if ! command -v git &> /dev/null; then
    message_std "SetKokkosDevelop> ERROR: git is not available"
    exit 1
fi

cd $DIR_CONTAINING_TRILINOS
git clone --depth=1 --single-branch --branch=develop --shallow-submodules https://github.com/kokkos/kokkos.git kokkos
git clone --depth=1 --single-branch --branch=develop --shallow-submodules https://github.com/kokkos/kokkos-kernels.git kokkos-kernels
message_std "SetKokkosDevelop> INFO: updated kokkos and kokkos-kernels packages with current develop"

ln -s "${DIR_CONTAINING_TRILINOS}/kokkos" "$TRILINOS_SRC/kokkos"
ln -s "${DIR_CONTAINING_TRILINOS}/kokkos-kernels" "$TRILINOS_SRC/kokkos-kernels"

cd -