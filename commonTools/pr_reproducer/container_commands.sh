#!/bin/bash

clean_trilinos() {
    cmake --build $TRILINOS_BUILD_DIR --target clean
}

get_dependencies() {
    ${TRILINOS_DIR}/packages/framework/get_dependencies.sh
}

configure_trilinos() {
    source ${TRILINOS_DIR}/packages/framework/GenConfig/gen-config.sh --force --yes --ci-mode --cmake-fragment /workspace/pr.cmake ${GENCONFIG_BUILD_ID} ${TRILINOS_DIR}
    rm -rf $TRILINOS_BUILD_DIR/CMakeCache.txt $TRILINOS_BUILD_DIR/CMakeFiles/
    cmake -S $TRILINOS_DIR -B $TRILINOS_BUILD_DIR -C /workspace/pr.cmake -C /workspace/packageEnables.cmake -G Ninja |& tee $TRILINOS_BUILD_DIR/configure.log
}

build_trilinos() {
    cmake --build $TRILINOS_BUILD_DIR $* -- -j 30
}

test_trilinos() {
    ctest --test-dir $TRILINOS_BUILD_DIR $*
}

install_trilinos() {
    rm -rf ${TRILINOS_INSTALL_DIR}/*
    cmake --build $TRILINOS_BUILD_DIR --target install $* -- -j 30
}

reproduce_build() {
    get_dependencies
    configure_trilinos
    build_trilinos
    test_trilinos
}

export -f clean_trilinos
export -f get_dependencies
export -f configure_trilinos
export -f build_trilinos
export -f test_trilinos
export -f install_trilinos

mkdir -p ${TRILINOS_BUILD_DIR}
mkdir -p ${TRILINOS_INSTALL_DIR}

echo ""
echo "The Trilinos source repository is ${TRILINOS_DIR}"
echo "The build directory is            ${TRILINOS_BUILD_DIR}"
echo "The install directory is          ${TRILINOS_INSTALL_DIR}"
echo ""
echo "Available commands:"
echo " clean_trilinos     -- cleans up CMake files and allows to configure from scratch"
echo " get_dependencies   -- installs GenConfig dependencies"
echo " configure_trilinos -- configures Trilinos using CMake"
echo " build_trilinos     -- build Trilinos"
echo " test_trilinos      -- runs tests"
echo " install_trilinos   -- installs Trilinos"
echo " reproduce_build    -- installs dependencies, configures, builds and tests"
echo ""
echo "Normally, this is what you want to do:"
echo " reproduce_build"
echo ""
