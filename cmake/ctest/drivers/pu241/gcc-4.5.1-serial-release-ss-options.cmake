#
# checkin-test-fissile4.sh SERIAL_RELEASE_SS build
#

SET(Trilinos_ENABLE_SECONDARY_STABLE_CODE  ON  CACHE BOOL "")
SET(TPL_ENABLE_MPI  OFF  CACHE  BOOL "")

# Include last so that above override these cache variables
INCLUDE(${${PROJECT_NAME}_HOME_DIR}/cmake/ctest/drivers/pu241/gcc-4.5.1-serial-ss-options.cmake)
INCLUDE(${${PROJECT_NAME}_HOME_DIR}/cmake/ctest/drivers/pu241/gcc-4.5.1-release-options.cmake)
