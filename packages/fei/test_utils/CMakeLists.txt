
INCLUDE(PackageLibraryMacros)

#
# A) Package-specific configuration options
#

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# B) Define the header and source files (and directories)
#

#
# src
#

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../support-Trilinos/)

SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

#
# Core headers
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
  AccessPattern.hpp
  BCNodeSet.hpp
  CoefAccess.hpp
  CommNodeSet.hpp
  CRSet.hpp
  DataReader.hpp
  driverData.hpp
  ElemBlock.hpp
  FEData.hpp
  FEI_tester.hpp
  feitester.hpp
  fei_test_utils.hpp
  HexBeamCR.hpp
  HexBeam.hpp
  InputData.hpp
  InputFileReader.hpp
  LibraryFactory.hpp
  PoissonData.hpp
  Poisson_Elem.hpp
  poisson_beam_mains.hpp
  snl_fei_tester.hpp
  SolnCheck.hpp
  test_Algebraic.hpp
  test_AztecWrappers.hpp
  test_Database.hpp
  testData.hpp
  test_EqnBuffer.hpp
  test_EqnCommMgr.hpp
  tester.hpp
  test_Factory.hpp
  test_Factory_helper.hpp
  test_FEI.hpp
  test_FEI_Implementation.hpp
  test_LinearSystem.hpp
  test_MatrixGraph.hpp
  test_Matrix.hpp
  test_misc.hpp
  test_PointBlockMap.hpp
  test_Set.hpp
  test_SNL_FEI_Structure.hpp
  test_FEI_Impl.hpp
  test_Tables.hpp
  test_Utils.hpp
  test_Vector.hpp
  test_VectorSpace.hpp
  test_benchmarks.hpp
  )

APPEND_SET(SOURCES
  BCNodeSet.cpp
  CommNodeSet.cpp
  CRSet.cpp
  beam_main.cpp
  beam_oldfei_main.cpp
  DataReader.cpp
  driverData.cpp
  ElemBlock.cpp
  feiDriver_main.cpp
  fei_test_utils.cpp
  FEData.cpp
  FEI_tester.cpp
  HexBeam.cpp
  HexBeamCR.cpp
  InputFileReader.cpp
  LibraryFactory.cpp
  PoissonData.cpp
  Poisson_Elem.cpp
  poisson3_main.cpp
  poisson_main.cpp
  snl_fei_tester.cpp
  SolnCheck.cpp
  tester.cpp
  test_Algebraic.cpp
  test_AztecWrappers.cpp
  test_Database.cpp
  test_EqnBuffer.cpp
  test_EqnCommMgr.cpp
  test_Factory.cpp
  test_Factory_helper.cpp
  test_FEI.cpp
  test_FEI_Implementation.cpp
  test_LinearSystem.cpp
  test_Matrix.cpp
  test_MatrixGraph.cpp
  test_misc.cpp
  test_PointBlockMap.cpp
  test_Set.cpp
  test_SNL_FEI_Structure.cpp
  test_FEI_Impl.cpp
  test_Tables.cpp
  test_Utils.cpp
  test_Vector.cpp
  test_VectorSpace.cpp
  test_benchmarks.cpp
  )

IF (TPL_ENABLE_MPI)
ENDIF()

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  fei_test_utils
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
