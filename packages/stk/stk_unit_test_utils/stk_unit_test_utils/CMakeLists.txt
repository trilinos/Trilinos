
INCLUDE_DIRECTORIES(${${PARENT_PACKAGE_NAME}_SOURCE_DIR})
INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR})
INCLUDE_DIRECTORIES(${${PARENT_PACKAGE_NAME}_SOURCE_DIR}/stk_ngp_test)
INCLUDE_DIRECTORIES(${${PARENT_PACKAGE_NAME}_SOURCE_DIR}/stk_util)
INCLUDE_DIRECTORIES(${${PARENT_PACKAGE_NAME}_BINARY_DIR}/stk_util)

SET(HEADERS_NEED_IO
  ElemGraphMultipleSharedSidesUtils.hpp
  BulkDataTester.hpp
  FaceCreationTestUtils.hpp
  GeneratedMeshToFile.hpp
  ReadWriteSidesetTester.hpp
  StkMeshFromGeneratedMesh.hpp
  TextMeshFixture.hpp
  TextMeshStkTopologyMapping.hpp
  TextMeshToFile.hpp
  TextMeshUtils.hpp
  TextMesh.hpp
  ioUtils.hpp
  )
SET(SOURCES_NEED_IO
  BulkDataTester.cpp
  ConstructedMesh.cpp
  FaceTestingUtils.cpp
  GenerateALefRAMesh.cpp
  GeneratedMeshToFile.cpp
  ReadWriteSidesetTester.cpp
  TextMeshFixture.cpp
  TextMeshToFile.cpp
  TextMesh.cpp
  ioUtils.cpp
  meshCreationHelpers.cpp
  )

SET(HEADERS_NEED_MESH_BUT_NOT_IO
  BucketTester.hpp
  MeshFixture.hpp
  MeshFileFixture.hpp
  MeshUtilsForBoundingVolumes.hpp
  ParticleUtils.hpp
  WriteSelectedEntities.hpp
  unittestMeshUtils.hpp
  )
SET(SOURCES_NEED_MESH_BUT_NOT_IO
  ParticleUtils.cpp
  WriteSelectedEntities.cpp
  unittestMeshUtils.cpp
  )

SET(HEADERS_NEED_MPI
  ParallelGtestOutput.hpp
  getOption.h
  StkReportRedirector.hpp
  stringAndNumberComparisons.hpp
  )
SET(SOURCES_NEED_MPI
  ParallelGtestOutput.cpp
  StkReportRedirector.cpp
  stringAndNumberComparisons.cpp
  )


SET(HEADERS
  getOption.h
  StkReportRedirector.hpp
  stringAndNumberComparisons.hpp
  )
SET(SOURCES
  StkReportRedirector.cpp
  stringAndNumberComparisons.cpp
  )


IF(TPL_ENABLE_MPI)
LIST(APPEND SOURCES ${SOURCES_NEED_MPI})
LIST(APPEND HEADERS ${HEADERS_NEED_MPI})
ENDIF()

IF(Trilinos_ENABLE_STKMesh)
LIST(APPEND SOURCES ${SOURCES_NEED_MESH_BUT_NOT_IO})
LIST(APPEND HEADERS ${HEADERS_NEED_MESH_BUT_NOT_IO})
ENDIF()

IF(Trilinos_ENABLE_STKIO)
LIST(APPEND SOURCES ${SOURCES_NEED_IO})
LIST(APPEND HEADERS ${HEADERS_NEED_IO})
ENDIF()


IF(Trilinos_ENABLE_STKMesh)
  ADD_SUBDIRECTORY(stk_mesh_fixtures)
ENDIF()

INSTALL(FILES ${HEADERS} DESTINATION
    ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/stk_unit_test_utils/)
    
TRIBITS_ADD_LIBRARY(
  stk_unit_test_utils
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  ${Gtest_NO_INSTALL_LIB_OR_HEADERS_ARG}
  )

IF (NOT DEFINED STK_ENABLE_UnitMain)
  SET(STK_ENABLE_UnitMain ON)
ENDIF()

IF (STK_ENABLE_UnitMain)
  FILE(GLOB UNIT_MAIN_SOURCES unit_main_lib/*.cpp)

  TRIBITS_ADD_LIBRARY(
    stk_unit_main
    TESTONLY
    SOURCES ${UNIT_MAIN_SOURCES}
  )
ENDIF()

