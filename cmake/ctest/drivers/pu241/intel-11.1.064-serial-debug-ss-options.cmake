#
# checkin-test-fissile4.sh INTEL11_SERIAL_DEBUG build
#

SET(${PROJECT_NAME}_ENABLE_SECONDARY_TESTED_CODE  ON  CACHE BOOL "")

# Include last so that above override these cache variables
INCLUDE(${CMAKE_CURRENT_LIST_DIR}/intel-11.1.064-serial-options.cmake)
INCLUDE(${CMAKE_CURRENT_LIST_DIR}/intel-debug-options.cmake)
INCLUDE(${CMAKE_CURRENT_LIST_DIR}/boost-1.46.1-options.cmake)
