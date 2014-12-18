
# For Tpetra_config.h to get explicit instantiation info
#
# FIXME (mfh 17 Dec 2014) Tpetra_config.h doesn't exist anymore.  I'll
# do my best by pointing to tpetra/core/src, where one could find the
# automatically generated header file TpetraCore_config.h.  However,
# the best policy is still to include Tpetra_ConfigDefs.hpp.
# Furthermore, if Zoltan2 depends on Tpetra, why does Zoltan2 need to
# add a Tpetra include directory explicitly to its list of include
# directories?
#TRIBITS_INCLUDE_DIRECTORIES(${PROJECT_BINARY_DIR}/packages/tpetra/src)
TRIBITS_INCLUDE_DIRECTORIES(${PROJECT_BINARY_DIR}/packages/tpetra/core/src)

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
