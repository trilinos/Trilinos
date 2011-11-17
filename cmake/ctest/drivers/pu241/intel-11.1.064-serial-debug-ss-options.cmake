#
# checkin-test-fissile4.sh INTEL11_SERIAL_DEBUG build
#

SET(Trilinos_ENABLE_SECONDARY_STABLE_CODE  ON  CACHE BOOL "")

# Include last so that above override these cache variables
INCLUDE(${${PROJECT_NAME}_HOME_DIR}/cmake/ctest/drivers/pu241/intel-11.1.064-serial-options.cmake)
INCLUDE(${${PROJECT_NAME}_HOME_DIR}/cmake/ctest/drivers/pu241/intel-debug-options.cmake)
INCLUDE(${${PROJECT_NAME}_HOME_DIR}/cmake/ctest/drivers/pu241/boost-1.46.1-options.cmake)
