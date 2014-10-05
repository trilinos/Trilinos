
# For Tpetra_config.h to get explicit instantiation info
TRIBITS_INCLUDE_DIRECTORIES(${PROJECT_BINARY_DIR}/packages/tpetra/src)

# For Zoltan2_TestHelpers.hpp
TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/helpers)

ADD_SUBDIRECTORIES(
unit
partition
order
color
scaling
correctness
temp
driver
)
