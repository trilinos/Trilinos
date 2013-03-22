#
# checkin-test-fissile4.sh SERIAL_RELEASE_SS build
#

SET(${PROJECT_NAME}_ENABLE_SECONDARY_STABLE_CODE  ON  CACHE BOOL "")
SET(TPL_ENABLE_MPI  OFF  CACHE  BOOL "")
SET(${PROJECT_NAME}_EXCLUDE_PACKAGES ${${PROJECT_NAME}_EXCLUDE_PACKAGES} Panzer)

# Include last so that above override these cache variables
INCLUDE(${CMAKE_CURRENT_LIST_DIR}/gcc-4.6.1-serial-ss-options.cmake)
# these options are the same for 4.6.1 and 4.5.1
INCLUDE(${CMAKE_CURRENT_LIST_DIR}/gcc-4.5.1-release-options.cmake)
