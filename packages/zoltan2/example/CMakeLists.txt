INCLUDE(AddSubdirectories)

# For Tpetra_config.h to get explicit instantiation info
INCLUDE_DIRECTORIES(${PROJECT_BINARY_DIR}/packages/tpetra/src)

# For Zoltan2_TestHelpers.hpp
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../test/helpers)


ADD_SUBDIRECTORIES(
  block
  geometric
)
