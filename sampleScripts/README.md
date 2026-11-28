This directory contains sample environment scripts and CMake
configuration snippets for building Trilinos on a couple of different
HPC systems. 

The environment can be loaded like 
```
source ${TRILINOS_SRC_DIR}/sampleScripts/environments/frontier/env.sh
```

The CMake snippets can be used like this:

```cmake
cmake \
  -S ${TRILINOS_SRC_DIR} \
  -B ${TRILINOS_BUILD_DIR} \
  -C ${TRILINOS_SRC_DIR}/sampleScripts/configurations/frontier/config_nightly_performance.cmake
```
