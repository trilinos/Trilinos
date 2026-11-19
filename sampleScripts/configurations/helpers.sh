#!/bin/bash

# Assumes that these environment variables are set
#  TRILINOS_DIR
#  TRILINOS_BUILD_DIR
#  TRILINOS_INSTALL_DIR


trilinos_clean() {
    if [ "${TRILINOS_BUILD_DIR}" == "" ]; then
        echo "TRILINOS_BUILD_DIR is not set or is empty."
        exit 1
    fi
    cmake --build ${TRILINOS_BUILD_DIR} --target clean
}

trilinos_configure() {
    # Configure Trilinos
    # All arguments get passed on to CMake

    if [ "${TRILINOS_DIR}" == "" ]; then
        echo "TRILINOS_DIR is not set or is empty."
        exit 1
    fi
    if [ "${TRILINOS_BUILD_DIR}" == "" ]; then
        echo "TRILINOS_BUILD_DIR is not set or is empty."
        exit 1
    fi
    if [ "${TRILINOS_INSTALL_DIR}" == "" ]; then
        echo "TRILINOS_INSTALL_DIR is not set or is empty."
        exit 1
    fi

    sourceDir=${TRILINOS_DIR}
    buildDir=${TRILINOS_BUILD_DIR}
    installDir=${TRILINOS_INSTALL_DIR}

    if [ ! -d ${sourceDir} ] ; then
       echo "Source directory ${sourceDir} does not exist."
       exit 1
    fi

    # make sure build dir exists
    mkdir -p ${buildDir}
    mkdir -p ${installDir}

    # delete CMake files
    rm -rf ${buildDir}/CMakeCache.txt ${buildDir}/CMakeFiles/

    # configure
    cmake -S ${sourceDir} -B ${buildDir} -D CMAKE_INSTALL_PREFIX:PATH=${installDir} $* |& tee ${buildDir}/configure.log

    if [ "$*" == "" ]; then
        echo ""
        echo "No configuration options were set for \"trilinos_configure\"."
        echo "You could try"
        echo " trilinos_configure -C \$TRILINOS_DIR/sampleScripts/configurations/cpu/config.cmake"
    fi
}

NP=${NP:=20}
export NP

trilinos_build() {
    if [ "${TRILINOS_BUILD_DIR}" == "" ]; then
        echo "TRILINOS_BUILD_DIR is not set or is empty."
        exit 1
    fi
    cmake --build ${TRILINOS_BUILD_DIR} $* -- -j $NP
}

trilinos_test() {
    if [ "${TRILINOS_BUILD_DIR}" == "" ]; then
        echo "TRILINOS_BUILD_DIR is not set or is empty."
        exit 1
    fi
    ctest --test-dir ${TRILINOS_BUILD_DIR} $*
}

trilinos_install() {
    if [ "${TRILINOS_BUILD_DIR}" == "" ]; then
        echo "TRILINOS_BUILD_DIR is not set or is empty."
        exit 1
    fi
    if [ "${TRILINOS_INSTALL_DIR}" == "" ]; then
        echo "TRILINOS_INSTALL_DIR is not set or is empty."
        exit 1
    fi
    rm -rf ${TRILINOS_INSTALL_DIR}/*
    cmake --build ${TRILINOS_BUILD_DIR} --target install $* -- -j $NP
}


trilinos_help() {
    echo ""
    echo "The Trilinos source repository is  ${TRILINOS_DIR}"
    echo "The build directory is             ${TRILINOS_BUILD_DIR}"
    echo "The install directory is           ${TRILINOS_INSTALL_DIR}"
    echo ""
    echo "Available commands:"
    echo " trilinos_configure  -- Configure Trilinos"
    echo "                        Examples:"
    echo "                          trilinos_configure -D Trilinos_ENABLE_Tpetra=ON"
    echo "                          trilinos_configure -C \$TRILINOS_DIR/sampleScripts/configurations/cpu/config.cmake"
    echo " trilinos_build      -- Build Trilinos"
    echo "                        Build parallelism is set via the environment variable NP (Current value: ${NP})."
    echo " trilinos_test       -- Test Trilinos"
    echo " trilinos_install    -- Install Trilinos"
    echo " trilinos_clean      -- Clean build"
    echo " trilinos_help       -- Display this help message"
    echo ""
}

export -f trilinos_clean
export -f trilinos_configure
export -f trilinos_build
export -f trilinos_test
export -f trilinos_install
export -f trilinos_help


trilinos_help
