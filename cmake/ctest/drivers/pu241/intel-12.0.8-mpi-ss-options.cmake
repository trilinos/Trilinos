#
# checkin-test-fissile4.sh INTEL12_SERIAL_DEBUG build
#

SET(${PROJECT_NAME}_ENABLE_SECONDARY_TESTED_CODE  ON  CACHE BOOL "")

# Include last so that above override these cache variables
INCLUDE(${CMAKE_CURRENT_LIST_DIR}/intel-12.0.8-mpi-options.cmake)
INCLUDE(${CMAKE_CURRENT_LIST_DIR}/boost-1.49.0-options.cmake)
