#
# checkin-test-fissile4.sh INTEL11_SERIAL_RELEASE build
#

SET(Trilinos_ENABLE_SECONDARY_STABLE_CODE  ON  CACHE BOOL "")

# Include last so that above override these cache variables
INCLUDE(${TRILINOS_HOME_DIR}/cmake/ctest/drivers/pu241/intel-11.1.064-serial-options.cmake)
INCLUDE(${TRILINOS_HOME_DIR}/cmake/ctest/drivers/pu241/intel-release-options.cmake)
INCLUDE(${TRILINOS_HOME_DIR}/cmake/ctest/drivers/pu241/boost-1.46.1-options.cmake)
