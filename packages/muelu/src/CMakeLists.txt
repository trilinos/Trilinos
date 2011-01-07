
INCLUDE(PackageLibraryMacros)

#
# A) Package-specific configuration options
#

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.hpp)

#
# B) Define the header and source files (and directories)
#

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.hpp
  )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/Headers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/MueCentral)

#TODO: APPEND_GLOB(HEADERS */*.hpp) instead ?
APPEND_SET(HEADERS
  Headers/MueLu_Version.hpp
  Headers/MueLu_Types.hpp
  MueCentral/MueLu_Level.hpp
  MueCentral/MueLu_SaLevel.hpp
  MueCentral/MueLu_Hierarchy.hpp
  MueCentral/MueLu_Needs.hpp
  MueCentral/MueLu_RAPFactory.hpp
  Transfers/MueLu_SaPFactory.hpp
  Transfers/MueLu_TransPFactory.hpp
  Transfers/MueLu_TentativePFactory.hpp
  Smoothers/MueLu_GaussSeidel.hpp
  Smoothers/MueLu_SmootherPrototype.hpp
  Smoothers/MueLu_SmootherBase.hpp
  Smoothers/MueLu_SmootherFactoryBase.hpp
  Smoothers/MueLu_SmootherFactory.hpp
  Smoothers/MueLu_Smoother.hpp
  Smoothers/MueLu_IfpackSmoother.hpp
  Graph/MueLu_AggAlgorithm.hpp
  Graph/MueLu_AggAlgorithm2.hpp
  Graph/MueLu_AggAlgorithm2Comm.hpp
  Graph/MueLu_AggregationOptions.hpp
  Graph/MueLu_Aggregates.hpp
  Graph/MueLu_Graph.hpp
  )

#TODO: APPEND_GLOB(SOURCES */*.cpp) instead ?
APPEND_SET(SOURCES
  mac.cpp
  )

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  muelu
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

