#!/bin/bash


export NP=${BUILD_PARALLELISM:=2}

clean_trilinos() {
    cmake --build $TRILINOS_BUILD_DIR --target clean
}

configure_trilinos() {
    rm -rf $TRILINOS_BUILD_DIR/CMakeCache.txt $TRILINOS_BUILD_DIR/CMakeFiles/
    cmake -S $TRILINOS_DIR -B $TRILINOS_BUILD_DIR -C /scripts/trilinos-build.cmake |& tee $TRILINOS_BUILD_DIR/configure.log

    # fix up compile_commands.json so that it can be used from host system
    if [ -f $TRILINOS_BUILD_DIR/compile_commands.json ] ; then
        rm -f $TRILINOS_BUILD_DIR/compile_commands_host.json
        cp $TRILINOS_BUILD_DIR/compile_commands.json $TRILINOS_BUILD_DIR/compile_commands_host.json
        ESCAPED_TRILINOS_SRC=$(printf '%s\n' "${TRILINOS_DIR}" | sed -e 's/[\/&]/\\&/g')
        ESCAPED_TRILINOS_BUILD=$(printf '%s\n' "${TRILINOS_BUILD_DIR}" | sed -e 's/[\/&]/\\&/g')
        ESCAPED_TRILINOS_SRC_HOST=$(printf '%s\n' "${TRILINOS_SRC_HOST}" | sed -e 's/[\/&]/\\&/g')
        ESCAPED_TRILINOS_BUILD_HOST=$(printf '%s\n' "${TRILINOS_BUILD_HOST}" | sed -e 's/[\/&]/\\&/g')
        sed -i "s/${ESCAPED_TRILINOS_SRC}/${ESCAPED_TRILINOS_SRC_HOST}/g" ${TRILINOS_BUILD_DIR}/compile_commands_host.json
        sed -i "s/${ESCAPED_TRILINOS_BUILD}/${ESCAPED_TRILINOS_BUILD_HOST}/g" ${TRILINOS_BUILD_DIR}/compile_commands_host.json
    fi
}

build_trilinos() {
    cmake --build $TRILINOS_BUILD_DIR $* -- -j $NP
}

test_trilinos() {
    ctest --test-dir $TRILINOS_BUILD_DIR $*
}

install_trilinos() {
    rm -rf ${TRILINOS_INSTALL_DIR}/*
    cmake --build $TRILINOS_BUILD_DIR --target install $* -- -j $NP
}

export -f clean_trilinos
export -f configure_trilinos
export -f build_trilinos
export -f test_trilinos
export -f install_trilinos

export CTEST_OUTPUT_ON_FAILURE=1
export CTEST_PROGRESS_OUTPUT=1
export CTEST_PARALLEL_LEVEL=20

clean_kokkos_tools() {
    cmake --build $KOKKOS_TOOLS_BUILD_DIR --target clean
}


configure_kokkos_tools() {
    rm -rf $KOKKOS_TOOLS_BUILD_DIR/CMakeCache.txt $KOKKOS_TOOLS_BUILD_DIR/CMakeFiles/
    cmake -S ${KOKKOS_TOOLS_DIR} -B ${KOKKOS_TOOLS_BUILD_DIR} -DCMAKE_INSTALL_PREFIX=${KOKKOS_TOOLS_INSTALL_DIR} |& tee ${KOKKOS_TOOLS_BUILD_DIR}/kokkos-tools-configure.log
}

build_kokkos_tools() {
    cmake --build ${KOKKOS_TOOLS_BUILD_DIR} $* -- -j $NP
}

install_kokkos_tools() {
    cmake --build ${KOKKOS_TOOLS_BUILD_DIR} --target install $* -- -j $NP
}

export -f clean_kokkos_tools
export -f configure_kokkos_tools
export -f build_kokkos_tools
export -f install_kokkos_tools


launch_jupyter_notebook() {
    jupyter notebook --port=9999 --no-browser --allow-root --ip=0.0.0.0 \
        --NotebookApp.token='' --NotebookApp.password='' \
        --notebook-dir=/workspace/trilinos/build/packages/PyTrilinos2/examples/notebooks/  > /dev/null 2>&1 &
}
