#
# Secondary Stable Serial build with GCC 4.6.1
#

SET(${PROJECT_NAME}_ENABLE_EXPLICIT_INSTANTIATION  ON  CACHE  BOOL  "")

INCLUDE(${CMAKE_CURRENT_LIST_DIR}/gcc-4.6.1-serial-options.cmake)
INCLUDE(${CMAKE_CURRENT_LIST_DIR}/boost-1.46.1-options.cmake)

