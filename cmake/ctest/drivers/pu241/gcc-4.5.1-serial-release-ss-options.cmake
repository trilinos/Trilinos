#
# checkin-test-fissile4.sh SERIAL_RELEASE_SS build
#

SET(TPL_ENABLE_MPI  OFF  CACHE  BOOL "")

# Include last so that above override these cache variables
INCLUDE(${CMAKE_CURRENT_LIST_DIR}/gcc-4.5.1-serial-ss-options.cmake)
INCLUDE(${CMAKE_CURRENT_LIST_DIR}/gcc-4.5.1-release-options.cmake)
