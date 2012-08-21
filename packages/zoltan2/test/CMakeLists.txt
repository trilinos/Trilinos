INCLUDE(AddSubdirectories)

# For Tpetra_config.h to get explicit instantiation info
INCLUDE_DIRECTORIES(${PROJECT_BINARY_DIR}/packages/tpetra/src)

# For Zoltan2_TestHelpers.hpp
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/helpers)

ADD_SUBDIRECTORIES(
temp
scaling
correctness
order
partition
unit
)
