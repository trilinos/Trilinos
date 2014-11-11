#
# checkin-test-fissile4.sh MPI_DEBUG_SS build
#

SET(TPL_ENABLE_MPI  ON  CACHE  BOOL  "")
# SET(${PROJECT_NAME}_ENABLE_CHECKED_STL  ON  CACHE  BOOL  "")
SET(${PROJECT_NAME}_ENABLE_DEBUG_SYMBOLS  ON  CACHE  BOOL  "")
SET(Teuchos_ENABLE_DEFAULT_STACKTRACE  OFF  CACHE  BOOL  "")

# Disable boost libs since not build with checked-stl turned on
SET(STK_ENABLE_BoostLib  OFF  CACHE  BOOL  "")

# Include last so that above override these cache variables
INCLUDE(${CMAKE_CURRENT_LIST_DIR}/gcc-4.6.1-mpi-ss-options.cmake)
# these options are the same for 4.6.1 and 4.5.1
INCLUDE(${CMAKE_CURRENT_LIST_DIR}/gcc-4.5.1-debug-options.cmake)
