
if(HAVE_STK_Trilinos)
TRIBITS_INCLUDE_DIRECTORIES(${${PARENT_PACKAGE_NAME}_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${${PARENT_PACKAGE_NAME}_SOURCE_DIR}/stk_ngp_test)
TRIBITS_INCLUDE_DIRECTORIES(${${PARENT_PACKAGE_NAME}_SOURCE_DIR}/stk_util)
TRIBITS_INCLUDE_DIRECTORIES(${${PARENT_PACKAGE_NAME}_BINARY_DIR}/stk_util)
endif()

SET(HEADERS_NEED_IO
  ElemGraphMultipleSharedSidesUtils.hpp
  BulkDataTester.hpp
  FaceCreationTestUtils.hpp
  GeneratedMeshToFile.hpp
  ReadWriteSidesetTester.hpp
  StkMeshFromGeneratedMesh.hpp
  TextMeshAdjacencyGraph.hpp
  TextMeshAssembly.hpp
  TextMeshDataTypes.hpp
  TextMeshEntityGroup.hpp
  TextMeshFixture.hpp
  TextMeshFuncs.hpp
  TextMeshNodeset.hpp
  TextMeshSideset.hpp
  TextMeshSidesetSkinner.hpp
  TextMeshSidesetSplitter.hpp
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
  BuildMesh.hpp
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
  AssemblyUtils.hpp
  CommandLineArgs.hpp
  getOption.h
  StkReportRedirector.hpp
  stringAndNumberComparisons.hpp
  )
SET(SOURCES
  AssemblyUtils.cpp
  StkReportRedirector.cpp
  stringAndNumberComparisons.cpp
  )


IF(STK_HAS_MPI)
  LIST(APPEND SOURCES ${SOURCES_NEED_MPI})
  LIST(APPEND HEADERS ${HEADERS_NEED_MPI})
ENDIF()

IF(${PROJECT_NAME}_ENABLE_STKMesh)
  LIST(APPEND SOURCES ${SOURCES_NEED_MESH_BUT_NOT_IO})
  LIST(APPEND HEADERS ${HEADERS_NEED_MESH_BUT_NOT_IO})
ENDIF()

IF(${PROJECT_NAME}_ENABLE_STKIO)
  LIST(APPEND SOURCES ${SOURCES_NEED_IO})
  LIST(APPEND HEADERS ${HEADERS_NEED_IO})
ENDIF()


IF(${PROJECT_NAME}_ENABLE_STKMesh)
  ADD_SUBDIRECTORY(stk_mesh_fixtures)
ENDIF()

IF(${PROJECT_NAME}_ENABLE_STKTransfer AND ${PROJECT_NAME}_ENABLE_STKMiddle_mesh)
  ADD_SUBDIRECTORY(stk_transfer_fixtures)
ENDIF()

INSTALL(FILES ${HEADERS} DESTINATION
    ${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/stk_unit_test_utils/)
    
if(HAVE_STK_Trilinos)
  TRIBITS_ADD_LIBRARY(stk_unit_test_utils
    HEADERS ${HEADERS}
    SOURCES ${SOURCES}
    ${Gtest_NO_INSTALL_LIB_OR_HEADERS_ARG}
  )
else()
  add_library(stk_unit_test_utils ${SOURCES})

  target_link_libraries(stk_unit_test_utils PUBLIC stk_util_parallel)

  find_package(GTest REQUIRED)
  target_link_libraries(stk_unit_test_utils PUBLIC GTest::gtest)
endif()

target_include_directories(stk_unit_test_utils PUBLIC
  $<BUILD_INTERFACE:${STK_TOPLEVEL_SOURCE_DIR}/stk_unit_test_utils>
  $<INSTALL_INTERFACE:include>
)

if (NOT DEFINED STK_ENABLE_UnitMain)
  set(STK_ENABLE_UnitMain ON)
endif()

if (STK_ENABLE_UnitMain)
  FILE(GLOB UNIT_MAIN_SOURCES unit_main_lib/*.cpp)

  if(HAVE_STK_Trilinos)
    TRIBITS_ADD_LIBRARY(stk_unit_main
      TESTONLY
      SOURCES ${UNIT_MAIN_SOURCES}
    )
  else()
    add_library(stk_unit_main ${UNIT_MAIN_SOURCES})
    target_link_libraries(stk_unit_main PUBLIC stk_unit_test_utils)
    target_link_libraries(stk_unit_main PUBLIC stk_ngp_test)
    target_link_libraries(stk_unit_main PUBLIC stk_util_env)
    target_link_libraries(stk_unit_main PUBLIC stk_util_parallel)
    target_link_libraries(stk_unit_main PUBLIC Kokkos::kokkos)
  endif()
endif()

