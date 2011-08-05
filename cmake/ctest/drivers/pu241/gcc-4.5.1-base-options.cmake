INCLUDE(${TRILINOS_HOME_DIR}/cmake/ctest/drivers/pu241/mkl-12.0.4-options.cmake)
INCLUDE(${TRILINOS_HOME_DIR}/cmake/ctest/drivers/pu241/casl-vri-tpls.cmake)

SET(PVM_LIBRARY_DIRS /opt/gcc-4.5.1/tpls/pvm3/lib/LINUX64                    CACHE FILEPATH "")
SET(PVM_INCLUDE_DIRS /opt/gcc-4.5.1/tpls/pvm3/include                        CACHE FILEPATH "")

SET(TRILINOS_TOOLSET_BASE /opt/gcc-4.5.1/trilinos-toolset)
SET(Trilinos_EXTRA_LINK_FLAGS "-Wl,-rpath,${TRILINOS_TOOLSET_BASE}/lib64" CACHE STRING "")
