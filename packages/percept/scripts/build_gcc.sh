if ([ -n "$TRILINOS_DIR" ]); then
    echo "TRILINOS_DIR: $TRILINOS_DIR"
else
    echo "TRILINOS_DIR not set! Please set to the directory containing Trilinos."
    return
fi

rm -rf CMakeFiles CMakeCache.txt

cmake \
-GNinja \
-DTrilinos_CONFIGURE_OPTIONS_FILE:STRING=cmake/std/atdm/ATDMDevEnv.cmake \
-DTrilinos_ENABLE_Intrepid=ON \
-DTrilinos_ENABLE_Percept=ON \
-DTPL_ENABLE_yaml-cpp=ON \
-D yaml-cpp_INCLUDE_DIRS:FILEPATH="$SEMS_YAML_CPP_INCLUDE_PATH" \
-D yaml-cpp_LIBRARY_DIRS:FILEPATH="$SEMS_YAML_CPP_LIBRARY_PATH" \
$TRILINOS_DIR
