if ([ -n "$TRILINOS_DIR" ]); then
    echo "TRILINOS_DIR: $TRILINOS_DIR"
else
    echo "TRILINOS_DIR not set! Please set to the directory containing Trilinos."
    return
fi

. $TRILINOS_DIR/cmake/std/atdm/load-env.sh RHEL6-None-gnu-7.2.0-openmp-release-debug

module load sems-yaml_cpp
