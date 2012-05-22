INCLUDE(AddSubdirectories)

# For Tpetra_config.h to get explicit instantiation info
INCLUDE_DIRECTORIES(${PROJECT_BINARY_DIR}/packages/tpetra/src)

# For Zoltan2_TestHelpers.hpp
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/helpers)

# For use of input file readers in Zoltan's zdrive
IF (HAVE_${PACKAGE_NAME}_ZOLTAN)

SET(ZOLTAN_SRC "${PROJECT_SOURCE_DIR}/packages/zoltan/src")
SET(ZOLTAN_BIN "${PROJECT_BINARY_DIR}/packages/zoltan/src")

INCLUDE_DIRECTORIES(${ZOLTAN_BIN})

INCLUDE_DIRECTORIES(${ZOLTAN_SRC}/include)
INCLUDE_DIRECTORIES(${ZOLTAN_SRC}/ch)
INCLUDE_DIRECTORIES(${ZOLTAN_SRC}/driver)

ENDIF()

ADD_SUBDIRECTORIES(
scaling
correctness
order
partition
unit
)
