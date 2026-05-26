#!/bin/bash

BINDER_SOURCE=/binder
BINDER_BUILD=/binder/binder-build
BINDER_INSTALL=/binder-install

git clone https://github.com/RosettaCommons/binder.git ${BINDER_SOURCE}

mkdir -p ${BINDER_BUILD}

cmake -DCMAKE_CXX_COMPILER="$(which clang++)" -B ${BINDER_BUILD} -S ${BINDER_SOURCE} -DCMAKE_INSTALL_PREFIX=${BINDER_INSTALL}
cmake --build ${BINDER_BUILD} --target install

rm -rf ${BINDER_SOURCE}
