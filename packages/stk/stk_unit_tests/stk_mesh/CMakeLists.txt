# Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
# Solutions of Sandia, LLC (NTESS). Under the terms of Contract
# DE-NA0003525 with NTESS, the U.S. Government retains certain rights
# in this software.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
# 
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
# 
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
# 
#     * Neither the name of NTESS nor the names of its contributors
#       may be used to endorse or promote products derived from this
#       software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 

FILE(GLOB SOURCES
   *.cpp
   change_parts/*.cpp
   entitySorting/*.cpp
   edge_creation/*.cpp
   face_creation/*.cpp
   face_creation/element_graph/*.cpp
   face_creation/integration/*.cpp
   face_creation/mesh_diagnostics/*.cpp
   face_creation/skin_mesh/*.cpp
   face_creation/user_created/*.cpp
   ngp/*.cpp
)

IF(NOT ${PROJECT_NAME}_ENABLE_STKIO)
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/UnitTestFieldBLAS.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/UnitTestAddNodeSharingWithInternalSideset.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/UnitTestAura.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/UnitTestAuraToShared.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/UnitTestBoundaryAnalysis.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/UnitTestBucket.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/UnitTestBulkData.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/UnitTestBulkDataIdMapper.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/UnitTestBulkDataNotifications.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/UnitTestBulkDataSharing.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/UnitTestCheckOwnedOrphanedSidesOrEdges.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/UnitTestCEO.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/UnitTestCEOME.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/UnitTestCEONodes.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/UnitTestChangeEntityOwnerCommMaps.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/UnitTestCheckUniqueGlobalIds.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/UnitTestChangeParts.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/UnitTestCommInfoObserver.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/UnitTestCommunicateFieldData.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/UnitTestCreateEdges.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/UnitTestCreateFaces.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/UnitTestDeletedEntityCache.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/UnitTestDeleteEntities.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/UnitTestDestroyElements.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/UnitTestDistributedIndexWithBulkData.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/UnitTestElemGraphCoincidentElements.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/UnitTestEntitiesNodesHaveInCommon.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/UnitTestEntityProcMapping.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/UnitTestFEMHelper.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/UnitTestField.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/UnitTestFieldQueryFunctions.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/UnitTestForEachEntityLoopAbstractions.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/UnitTestGenIds.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/UnitTestGetEntities.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/UnitTestGhostingWithModification.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/UnitTestGhostingWithShared.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/UnitTestInducedPart.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/UnitTestKeyhole.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/UnitTestLocalIds.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/UnitTestMeshImplUtils.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/UnitTestMeshModLogObserver.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/UnitTestModificationEnd.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/UnitTestModificationSummary.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/UnitTestNoUpwardConnectivity.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/UnitTestParallelGraphInfo.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/UnitTestPartAfterCommit.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/UnitTestPartAlias.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/UnitTestPartToBucket.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/UnitTestSelfContainedModCycles.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/UnitTestSidePolarity.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/UnitTestSideSet.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/UnitTestStkTextMesh.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/UnitTestTextMesh.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/UnitTestTopology.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/UnitTestVisitAura.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/UnitTest3Tets3Procs.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/change_parts/CustomGhostEntities.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/entitySorting/UnitTestEntitySorting.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/edge_creation/UnitTestEdgeConnection.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/edge_creation/UnitTestEdgesOnFaces.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/face_creation/FaceCreatorFixture.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/face_creation/TwoHexesTwoProcs.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/face_creation/TwoHexesTwoProcsEGWithPermutations.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/face_creation/TwoHexesTwoProcsElemGraph.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/face_creation/TwQuads2DTwoProcsElemGraph.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/face_creation/element_graph/BulkDataElementGraphTester.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/face_creation/element_graph/ElementGraphAddElements.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/face_creation/element_graph/ElementGraphDeleteElements.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/face_creation/element_graph/UnitTestElemGraphMultipleSharedSides.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/face_creation/element_graph/UnitTestElemElemGraph.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/face_creation/element_graph/UnitTestElemElemGraphChangeOwner.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/face_creation/element_graph/UnitTestElemElemGraphDeathIntegration.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/face_creation/element_graph/UnitTestElementDeath.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/face_creation/element_graph/UnitTestSideIdPool.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/face_creation/integration/UnitTestSkinMeshElementDeath.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/face_creation/mesh_diagnostics/UnitTestRuleThreeViolation.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/face_creation/skin_mesh/UnitTestSkinMesh.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/face_creation/skin_mesh/UnitTestSkinMeshCoincidentElements.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/face_creation/skin_mesh/UnitTestSkinMeshDegenerate.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/face_creation/skin_mesh/UnitTestSkinMeshHeterogenous.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/face_creation/skin_mesh/UnitTestSkinMeshExposedBoundary.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/face_creation/skin_mesh/UnitTestSkinMeshMultipleSharedSides.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/face_creation/skin_mesh/UnitTestSkinMeshRefined.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/face_creation/skin_mesh/UnitTestSkinMeshShell.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/face_creation/skin_mesh/UnitTestSkinMeshSkinPart.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/face_creation/skin_mesh/UnitTestSkinMeshUseCaseUtils.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/face_creation/user_created/UnitTestFaceSharingUsingGraph.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/ngp/NgpDebugFieldSync_Fixtures.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/ngp/NgpParallelSumTest.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/ngp/NgpMeshManagerTest.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/ngp/TestNgpMeshUpdate.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/ngp/howToNgp.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/ngp/howToNgpMeshMod.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/ngp/implNgpMeshMod.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/ngp/NgpFieldAsyncTest.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/ngp/ngpFieldBLASTest.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/ngp/ngpFieldTest.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/ngp/UnitTestNgp.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/ngp/UnitTestNgpDebugFieldSync.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/ngp/UnitTestNgpDebugFieldSync_AccessDuringMeshModification.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/ngp/UnitTestNgpDebugFieldSync_MeshModification.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/ngp/UnitTestNgpDebugFieldSync_MultiStateRotation.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/ngp/UnitTestNgpDebugFieldSync_PartialAllocation.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/ngp/UnitTestNgpDebugFieldSync_PartialSync.cpp")
ENDIF()

if(HAVE_STK_Trilinos)
  TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

  TRIBITS_ADD_EXECUTABLE_AND_TEST(stk_mesh_unit_tests
    SOURCES ${SOURCES}
    TESTONLYLIBS stk_unit_main
    ARGS ""
    COMM serial mpi
    NUM_MPI_PROCS 1-4
    NOEXEPREFIX
  )
else()
  add_executable(stk_mesh_unit_tests ${SOURCES})
  target_link_libraries(stk_mesh_unit_tests stk_unit_test_utils)

  target_include_directories(stk_mesh_unit_tests PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
  )

  if(${PROJECT_NAME}_ENABLE_STKIO)
    target_link_libraries(stk_mesh_unit_tests stk_io)
  endif()
  target_link_libraries(stk_mesh_unit_tests stk_mesh_fixtures)
  target_link_libraries(stk_mesh_unit_tests stk_mesh_base)
  target_link_libraries(stk_mesh_unit_tests stk_unit_main)

  add_test(NAME "stk_mesh_unit_tests" COMMAND stk_mesh_unit_tests)
endif()

