#
# checkin-test-fissile4.sh INTEL12_SERIAL_DEBUG build
#

SET(Trilinos_ENABLE_SECONDARY_STABLE_CODE  ON  CACHE BOOL "")

# Include last so that above override these cache variables
INCLUDE(${TRILINOS_HOME_DIR}/cmake/ctest/drivers/pu241/intel-12.0.4-serial-options.cmake)
INCLUDE(${TRILINOS_HOME_DIR}/cmake/ctest/drivers/pu241/intel-debug-options.cmake)
INCLUDE(${TRILINOS_HOME_DIR}/cmake/ctest/drivers/pu241/boost-1.46.1-options.cmake)
