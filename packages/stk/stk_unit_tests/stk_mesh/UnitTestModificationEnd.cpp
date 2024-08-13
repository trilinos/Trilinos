// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#include <stddef.h>                     // for size_t
#include <stdlib.h>                     // for exit
#include <exception>                    // for exception
#include <iostream>                     // for ostringstream, etc
#include <iterator>                     // for distance
#include <map>                          // for _Rb_tree_const_iterator, etc
#include <stdexcept>                    // for logic_error, runtime_error
#include <stk_mesh/base/BulkData.hpp>   // for BulkData, etc
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/GetEntities.hpp>  // for count_entities, etc
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine, etc
#include <gtest/gtest.h>
#include <string>                       // for string, basic_string, etc
#include <vector>                       // for vector, etc

#include "stk_mesh/base/Bucket.hpp"     // for Bucket, has_superset
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/EntityKey.hpp"  // for EntityKey
#include "stk_mesh/base/Ghosting.hpp"   // for Ghosting
#include "stk_mesh/base/MetaData.hpp"   // for MetaData, entity_rank_names, etc
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_mesh/base/Relation.hpp"
#include "stk_mesh/base/Selector.hpp"   // for Selector, operator|
#include "stk_mesh/base/Types.hpp"      // for EntityProc, EntityVector, etc
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/util/PairIter.hpp"   // for PairIter

#include "stk_io/StkMeshIoBroker.hpp"
#include <stk_mesh/base/Comm.hpp>

#include "stk_mesh/base/CreateEdges.hpp"
#include "stk_mesh/base/CreateFaces.hpp"
#include "stk_mesh/base/SkinMesh.hpp"
#include "stk_mesh/baseImpl/MeshImplUtils.hpp"
#include <stk_unit_test_utils/BulkDataTester.hpp>
#include "stk_unit_test_utils/TextMesh.hpp"
#include "UnitTestModificationEnd.hpp"

#include <stdio.h> // getline

namespace stk { namespace mesh { namespace unit_test {

//TEST(BulkDataModificationEnd, test_IR_ghosted_modify_delete_with_element_being_marked_as_modified)
//{
//  MPI_Comm communicator = MPI_COMM_WORLD;
//  int numProcs = -1;
//  MPI_Comm_size(communicator, &numProcs);
//
//  if(numProcs == 2)
//  {
//    //============== Load 1x1x4 mesh into BulkData
//
//    const int spatialDim = 3;
//    stk::mesh::MetaData stkMeshMetaData(spatialDim);
//    stk::unit_test_util::BulkDataTester stkMeshBulkData(stkMeshMetaData, communicator);
//
//    std::string exodusFileName = getOption("-i", "generated:1x1x4");
//    populateBulkDataWithFile(exodusFileName, communicator, stkMeshBulkData);
//
//    //============== Get Global Count Of Elements and Nodes
//
//    std::vector<size_t> globalCounts;
//    stk::mesh::comm_mesh_counts(stkMeshBulkData, globalCounts);
//
//    size_t numNodes = 20;
//    size_t numEdges = 0;
//    size_t numFaces = 0;
//    size_t numElements = 4;
//
//    ASSERT_EQ(numNodes, globalCounts[stk::topology::NODE_RANK]);
//    ASSERT_EQ(numEdges, globalCounts[stk::topology::EDGE_RANK]);
//    ASSERT_EQ(numFaces, globalCounts[stk::topology::FACE_RANK]);
//    ASSERT_EQ(numElements, globalCounts[stk::topology::ELEMENT_RANK]);
//
//    checkThatMeshIsParallelConsistent(stkMeshBulkData);
//
//    bool isCheckAfterCallToIRGMD = false;
//
//    checkCommListAndMap(stkMeshBulkData, isCheckAfterCallToIRGMD);
//    stkMeshBulkData.modification_begin();
//    stkMeshBulkData.my_internal_resolve_ghosted_modify_delete();
//    checkCommListAndMap(stkMeshBulkData, isCheckAfterCallToIRGMD);
//
//    std::vector<std::vector<stk::mesh::EntityState> > elementStates(numProcs);
//    for (size_t i=0;i<elementStates.size();i++)
//    {
//      elementStates[i].resize(numElements, stk::mesh::Unchanged);
//    }
//
//    std::vector<std::vector<stk::mesh::EntityState> > nodeStates(numProcs);
//    for (size_t i=0;i<nodeStates.size();i++)
//    {
//      nodeStates[i].resize(numNodes, stk::mesh::Unchanged);
//    }
//
//    bool areNodesValid[2][20] = {
//      { true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, false, false, false, false, },
//      { false, false, false, false, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true }
//    };
//
//    bool areElementsValid[2][4] = {
//      { true, true, true, false, },
//      { false, true, true, true }
//    };
//
//    checkStatesOfEntities(nodeStates, elementStates, areNodesValid, areElementsValid, stkMeshBulkData);
//
//    mark_element3_as_modified(stkMeshBulkData);
//
//    int element3Index = 2;
//    int proc_that_modified_element3 = 1;
//    elementStates[proc_that_modified_element3][element3Index] = stk::mesh::Modified;
//
//    checkStatesOfEntities(nodeStates, elementStates, areNodesValid, areElementsValid, stkMeshBulkData);
//
//    //============== Call IRGMD and check what happens to states of entities and comm maps and lists
//
//    stkMeshBulkData.my_internal_resolve_ghosted_modify_delete();
//
//    int proc_that_no_longer_has_element3 = 0;
//    areElementsValid[proc_that_no_longer_has_element3][element3Index] = false;
//
//    for (size_t i=9;i<=16;i++)
//    {
//      int proc_on_which_nodes_9thru16_are_modified = 0;
//      nodeStates[proc_on_which_nodes_9thru16_are_modified][i-1] = stk::mesh::Modified;
//    }
//
//    int elementNextToModifiedElementOnProc0 = 1;
//    int proc_that_has_element2_modified_as_result_of_element3_modification = 0;
//    elementStates[proc_that_has_element2_modified_as_result_of_element3_modification][elementNextToModifiedElementOnProc0] = stk::mesh::Modified;
//
//    checkStatesOfEntities(nodeStates, elementStates, areNodesValid, areElementsValid, stkMeshBulkData);
//
//    isCheckAfterCallToIRGMD = true;
//    checkCommListAndMap(stkMeshBulkData, isCheckAfterCallToIRGMD);
//  }
//}

//TEST(BulkDataModificationEnd, create_an_edge_and_test_pieces_of_internal_modification_end_for_edges)
//{
//  MPI_Comm communicator = MPI_COMM_WORLD;
//  int numProcs = -1;
//  MPI_Comm_size(communicator, &numProcs);
//
//  if(numProcs == 2)
//  {
//    //============== Load 1x1x4 mesh into Bulk Data
//
//    const int spatialDim = 3;
//    stk::mesh::MetaData stkMeshMetaData(spatialDim);
//    stk::unit_test_util::BulkDataTester stkMeshBulkData(stkMeshMetaData, communicator);
//
//    // Elements 1 and 2 on proc 0, Elements 3 and 4 on proc 1
//    // Elements 2 and 3 are shared because of nodes 9, 10, 11, 12
//    // Element 2 is ghosted onto Proc 1, and Element 0 is ghosted onto Proc 0
//
//    std::string exodusFileName = getOption("-i", "generated:1x1x4");
//    populateBulkDataWithFile(exodusFileName, communicator, stkMeshBulkData);
//
//    //============== Before starting, make sure mesh is parallel consistent
//
//    std::vector<std::string> meshStart;
//    getMeshLineByLine(stkMeshBulkData, meshStart);
//
//    checkThatMeshIsParallelConsistent(stkMeshBulkData);
//
//    int myProcId = stkMeshBulkData.parallel_rank();
//
//    std::vector<stk::mesh::EntityId> edgeIds;
//    std::vector<std::vector<stk::mesh::EntityId> > nodeIdsForEdge;
//    std::vector<std::vector<stk::mesh::EntityId> > elementRelations;
//
//    edgeIds.push_back(100+myProcId);
//    nodeIdsForEdge.resize(edgeIds.size());
//    elementRelations.resize(edgeIds.size());
//    for ( size_t i=0;i<edgeIds.size();i++)
//    {
//      nodeIdsForEdge[i].resize(2);
//      nodeIdsForEdge[i][0] = 9;
//      nodeIdsForEdge[i][1] = 10;
//
//      elementRelations[i].resize(1);
//      elementRelations[i][0] = 2+myProcId;
//    }
//
//    std::vector<stk::mesh::Entity> edgeEntities(edgeIds.size());
//    stk::mesh::Part& edge_part = stkMeshBulkData.mesh_meta_data().declare_part("edge_part", stk::topology::EDGE_RANK);
//
//    stkMeshBulkData.modification_begin();
//
//    create_edges(stkMeshBulkData, edgeIds, nodeIdsForEdge, elementRelations, edgeEntities, edge_part);
//
//    stkMeshBulkData.my_internal_resolve_shared_modify_delete();
//
//    checkResultsOfIRSMD_for_edges(stkMeshBulkData);
//
//    //============== Make sure we have 2 edges
//
//    std::vector<size_t> globalCounts;
//    stk::mesh::comm_mesh_counts(stkMeshBulkData, globalCounts);
//    size_t numEdgesTotal = 2;
//    EXPECT_EQ(numEdgesTotal, globalCounts[stk::topology::EDGE_RANK]);
//
//    //============== Check results of IRGMD
//
//    stkMeshBulkData.my_internal_resolve_ghosted_modify_delete();
//
//    checkResultsOfIRGMD_for_edges(stkMeshBulkData, edgeEntities);
//
//    //============== Before starting, make sure mesh is parallel consistent
//
//    checkThatMeshIsParallelConsistent(stkMeshBulkData);
//
//    //============== Check the result of internal_update_distributed_index
//
//    for (size_t i=0;i<edgeIds.size();i++)
//    {
//      stk::mesh::EntityKey edgeKey(stk::topology::EDGE_RANK, edgeIds[i]);
//      stk::mesh::EntityCommListInfoVector::const_iterator iter = std::lower_bound(stkMeshBulkData.my_internal_comm_list().begin(), stkMeshBulkData.my_internal_comm_list().end(), edgeKey);
//      bool edge_is_not_in_comm_list = iter == stkMeshBulkData.my_internal_comm_list().end() || iter->key != edgeKey;
//      EXPECT_TRUE(edge_is_not_in_comm_list);
//
//      const bool is_entity_ghosted = !stkMeshBulkData.my_internal_entity_comm_map(edgeKey, stkMeshBulkData.aura_ghosting()).empty();
//      EXPECT_FALSE( is_entity_ghosted );
//
//      const bool is_entity_shared = !stkMeshBulkData.my_internal_entity_comm_map(edgeKey, stkMeshBulkData.shared_ghosting()).empty();
//      EXPECT_FALSE( is_entity_shared );
//
//    }
//
//    for (size_t i=0;i<elementRelations.size();i++)
//    {
//      for (size_t j=0;j<elementRelations[i].size();j++)
//      {
//        stk::mesh::EntityKey elementKey(stk::topology::ELEMENT_RANK, elementRelations[i][j]);
//        stk::mesh::EntityCommListInfoVector::const_iterator iter = std::lower_bound(stkMeshBulkData.my_internal_comm_list().begin(), stkMeshBulkData.my_internal_comm_list().end(), elementKey);
//        bool element_not_in_comm_list = iter == stkMeshBulkData.my_internal_comm_list().end();
//        EXPECT_TRUE(element_not_in_comm_list);
//
//        const bool is_entity_ghosted = !stkMeshBulkData.my_internal_entity_comm_map(elementKey, stkMeshBulkData.aura_ghosting()).empty();
//        EXPECT_FALSE( is_entity_ghosted );
//
//        const bool is_entity_shared = !stkMeshBulkData.my_internal_entity_comm_map(elementKey, stkMeshBulkData.shared_ghosting()).empty();
//        EXPECT_FALSE( is_entity_shared );
//      }
//    }
//
//    std::vector<stk::mesh::Entity> shared_new;
//    stkMeshBulkData.my_internal_update_distributed_index(stk::topology::EDGE_RANK, shared_new);
//    ASSERT_EQ(1u, shared_new.size());
//    for (size_t i=0;i<edgeIds.size();i++)
//    {
//      stk::mesh::EntityId shared_edge_id = 100;
//      EXPECT_EQ(shared_edge_id, stkMeshBulkData.identifier(shared_new[i]));
//
//      //============== Funky! Edge is marked as created, not modified. So resolved ownership of "modified" doesn't mess with created.
//
//      EXPECT_EQ(myProcId, stkMeshBulkData.parallel_owner_rank(shared_new[i]));
//      EXPECT_EQ(stk::mesh::Created, stkMeshBulkData.state(shared_new[i]));
//    }
//
//    stk::mesh::EntityKey entity_key(stk::topology::NODE_RANK, 9);
//    stk::mesh::Entity node = stkMeshBulkData.get_entity(entity_key);
//
//    EXPECT_TRUE(stkMeshBulkData.bucket(node).member(edge_part));
//    for (size_t i=0;i<edgeEntities.size();i++)
//    {
//      EXPECT_FALSE(stkMeshBulkData.bucket(edgeEntities[i]).shared());
//    }
//
//    stkMeshBulkData.my_resolve_ownership_of_modified_entities(shared_new);
//
//    for (size_t i=0;i<edgeEntities.size();i++)
//    {
//      EXPECT_FALSE(stkMeshBulkData.bucket(edgeEntities[i]).shared());
//
//      EXPECT_EQ(myProcId, stkMeshBulkData.parallel_owner_rank(shared_new[i])) <<
//                                                                                 "Proc " << stkMeshBulkData.parallel_rank() << " failed.";
//
//      EXPECT_FALSE(stkMeshBulkData.bucket(shared_new[i]).shared());
//      EXPECT_TRUE(stkMeshBulkData.bucket(shared_new[i]).owned());
//      EXPECT_FALSE(stkMeshBulkData.bucket(shared_new[i]).member(stkMeshBulkData.ghosting_part(stkMeshBulkData.aura_ghosting())));
//
//      EXPECT_FALSE(stkMeshBulkData.bucket(edgeEntities[i]).shared());
//    }
//    EXPECT_TRUE(stkMeshBulkData.bucket(node).member(edge_part));
//
//    //============== Hopefully this works!
//
//    stkMeshBulkData.my_move_entities_to_proper_part_ownership( shared_new );
//
//    for (size_t i=0;i<edgeEntities.size();i++)
//    {
//      EXPECT_TRUE(stkMeshBulkData.bucket(shared_new[i]).shared());
//      if ( stkMeshBulkData.parallel_rank() == 0 )
//      {
//        EXPECT_TRUE(stkMeshBulkData.bucket(shared_new[i]).owned());
//      }
//      else
//      {
//        EXPECT_FALSE(stkMeshBulkData.bucket(shared_new[i]).owned());
//      }
//      EXPECT_FALSE(stkMeshBulkData.bucket(shared_new[i]).member(stkMeshBulkData.ghosting_part(stkMeshBulkData.aura_ghosting())));
//      EXPECT_EQ(0, stkMeshBulkData.parallel_owner_rank(shared_new[i])) ;
//    }
//
//    stkMeshBulkData.my_add_comm_list_entries_for_entities( shared_new );
//
//    for (size_t i=0;i<edgeEntities.size();i++)
//    {
//      stk::mesh::EntityKey edgeKey = stkMeshBulkData.entity_key(shared_new[i]);
//      stk::mesh::EntityCommListInfoVector::const_iterator  iter = std::lower_bound(stkMeshBulkData.my_internal_comm_list().begin(), stkMeshBulkData.my_internal_comm_list().end(), edgeKey);
//      bool edge_is_in_comm_list = iter != stkMeshBulkData.my_internal_comm_list().end() && iter->key == edgeKey;
//      EXPECT_TRUE(edge_is_in_comm_list);
//
//      const bool is_entity_ghosted = !stkMeshBulkData.my_internal_entity_comm_map(edgeKey, stkMeshBulkData.aura_ghosting()).empty();
//      EXPECT_FALSE( is_entity_ghosted );
//
//      const bool is_entity_shared = !stkMeshBulkData.my_internal_entity_comm_map(edgeKey, stkMeshBulkData.shared_ghosting()).empty();
//      EXPECT_TRUE( is_entity_shared );
//    }
//
//    //============== Not sure what IRSM is supposed to do?
//
//    std::vector<std::string> meshBefore;
//    getMeshLineByLine(stkMeshBulkData, meshBefore);
//
//    stkMeshBulkData.my_internal_resolve_shared_membership();
//
//    std::vector<std::string> meshAfter;
//    getMeshLineByLine(stkMeshBulkData, meshAfter);
//
//    stkMeshBulkData.my_internal_regenerate_aura();
//
//    std::vector<std::string> meshAfterAura;
//    getMeshLineByLine(stkMeshBulkData, meshAfterAura);
//
//    ////        compareMeshes(stkMeshBulkData.parallel_rank(), meshStart, meshAfterAura);
//    ////        writeMesh(myProcId, "Before::", meshStart);
//    //        writeMesh(myProcId, "After::", meshAfterAura);
//
//    checkItAllForThisCase(stkMeshBulkData);
//  }
//}

TEST(BulkDataModificationEnd, create_an_edge_and_test_up_to_IR_parallel_create)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  int numProcs = -1;
  MPI_Comm_size(communicator, &numProcs);

  if(numProcs == 2)
  {
    //============== Load 1x1x4 mesh into Bulk Data

    const int spatialDim = 3;
    stk::mesh::MetaData stkMeshMetaData(spatialDim);
    stk::unit_test_util::BulkDataTester stkMeshBulkData(stkMeshMetaData, communicator);

    // Elements 1 and 2 on proc 0, Elements 3 and 4 on proc 1
    // Elements 2 and 3 are shared because of nodes 9, 10, 11, 12
    // Element 2 is ghosted onto Proc 1, and Element 0 is ghosted onto Proc 0

    std::string exodusFileName = getOption("-i", "generated:1x1x4");
    populateBulkDataWithFile(exodusFileName, communicator, stkMeshBulkData);

    checkThatMeshIsParallelConsistent(stkMeshBulkData);

    int myProcId = stkMeshBulkData.parallel_rank();

    std::vector<stk::mesh::EntityId> edgeIds;
    std::vector<std::vector<stk::mesh::EntityId> > nodeIdsForEdge;
    std::vector<std::vector<stk::mesh::EntityId> > elementRelations;

    edgeIds.push_back(100+myProcId);
    nodeIdsForEdge.resize(edgeIds.size());
    elementRelations.resize(edgeIds.size());
    for ( size_t i=0;i<edgeIds.size();i++)
    {
      nodeIdsForEdge[i].resize(2);
      nodeIdsForEdge[i][0] = 9;
      nodeIdsForEdge[i][1] = 10;

      elementRelations[i].resize(1);
      elementRelations[i][0] = 2+myProcId;
    }

    std::vector<stk::mesh::Entity> edgeEntities(edgeIds.size());
    stk::mesh::Part& edge_part = stkMeshBulkData.mesh_meta_data().declare_part("edge_part", stk::topology::EDGE_RANK);

    stkMeshBulkData.modification_begin();

    create_edges(stkMeshBulkData, edgeIds, nodeIdsForEdge, elementRelations, edgeEntities, edge_part);

    std::vector<size_t> globalCounts;
    stk::mesh::comm_mesh_counts(stkMeshBulkData, globalCounts);
    size_t numEdgesTotal = 2;
    EXPECT_EQ(numEdgesTotal, globalCounts[stk::topology::EDGE_RANK]);

    checkThatMeshIsParallelConsistent(stkMeshBulkData);

    std::vector<stk::mesh::Entity> shared_new;
    stkMeshBulkData.my_internal_update_distributed_index(stk::topology::EDGE_RANK, shared_new);
    ASSERT_EQ(1u, shared_new.size());
    for (size_t i=0;i<edgeIds.size();i++)
    {
      stk::mesh::EntityId shared_edge_id = 100;
      EXPECT_EQ(shared_edge_id, stkMeshBulkData.identifier(shared_new[i]));

      //============== Funky! Edge is marked as created, not modified. So resolved ownership of "modified" doesn't mess with created.

      EXPECT_EQ(myProcId, stkMeshBulkData.parallel_owner_rank(shared_new[i]));
      EXPECT_EQ(stk::mesh::Created, stkMeshBulkData.state(shared_new[i]));
    }

    stk::mesh::EntityKey entity_key(stk::topology::NODE_RANK, 9);
    stk::mesh::Entity node = stkMeshBulkData.get_entity(entity_key);

    EXPECT_TRUE(stkMeshBulkData.bucket(node).member(edge_part));
    for (size_t i=0;i<edgeEntities.size();i++)
    {
      EXPECT_FALSE(stkMeshBulkData.bucket(edgeEntities[i]).shared());
    }

    stkMeshBulkData.my_resolve_ownership_of_modified_entities(shared_new);

    for (size_t i=0;i<edgeEntities.size();i++)
    {
      EXPECT_FALSE(stkMeshBulkData.bucket(edgeEntities[i]).shared());

      EXPECT_EQ(myProcId, stkMeshBulkData.parallel_owner_rank(shared_new[i])) <<
                                                                                 "Proc " << stkMeshBulkData.parallel_rank() << " failed.";

      EXPECT_FALSE(stkMeshBulkData.bucket(shared_new[i]).shared());
      EXPECT_TRUE(stkMeshBulkData.bucket(shared_new[i]).owned());
      EXPECT_FALSE(stkMeshBulkData.bucket(shared_new[i]).member(stkMeshBulkData.ghosting_part(stkMeshBulkData.aura_ghosting())));

      EXPECT_FALSE(stkMeshBulkData.bucket(edgeEntities[i]).shared());
    }
    EXPECT_TRUE(stkMeshBulkData.bucket(node).member(edge_part));

    //============== Hopefully this works!

    stkMeshBulkData.my_move_entities_to_proper_part_ownership( shared_new );

    for (size_t i=0;i<edgeEntities.size();i++)
    {
      EXPECT_TRUE(stkMeshBulkData.bucket(shared_new[i]).shared());
      if ( stkMeshBulkData.parallel_rank() == 0 )
      {
        EXPECT_TRUE(stkMeshBulkData.bucket(shared_new[i]).owned());
      }
      else
      {
        EXPECT_FALSE(stkMeshBulkData.bucket(shared_new[i]).owned());
      }
      EXPECT_FALSE(stkMeshBulkData.bucket(shared_new[i]).member(stkMeshBulkData.ghosting_part(stkMeshBulkData.aura_ghosting())));
      EXPECT_EQ(0, stkMeshBulkData.parallel_owner_rank(shared_new[i])) ;
    }

    stkMeshBulkData.my_add_comm_list_entries_for_entities( shared_new );

    for (size_t i=0;i<edgeEntities.size();i++)
    {
      stk::mesh::EntityKey edgeKey = stkMeshBulkData.entity_key(shared_new[i]);
      stk::mesh::EntityCommListInfoVector::const_iterator  iter = std::lower_bound(stkMeshBulkData.my_internal_comm_list().begin(), stkMeshBulkData.my_internal_comm_list().end(), edgeKey);
      bool edge_is_in_comm_list = iter != stkMeshBulkData.my_internal_comm_list().end() && iter->key == edgeKey;
      EXPECT_TRUE(edge_is_in_comm_list);

      const bool is_entity_ghosted = !stkMeshBulkData.my_internal_entity_comm_map(edgeKey, stkMeshBulkData.aura_ghosting()).empty();
      EXPECT_FALSE( is_entity_ghosted );

      const bool is_entity_shared = !stkMeshBulkData.my_internal_entity_comm_map(edgeKey, stkMeshBulkData.shared_ghosting()).empty();
      EXPECT_TRUE( is_entity_shared );
    }

    checkItAllForThisCase(stkMeshBulkData);
  }
}

TEST(BulkDataModificationEnd, create_a_ghosted_edge_and_test_internal_modification_end_for_edges)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  int numProcs = -1;
  MPI_Comm_size(communicator, &numProcs);

  if(numProcs == 2)
  {
    //============== Load 1x1x4 mesh into Bulk Data

    const int spatialDim = 3;
    stk::mesh::MetaData stkMeshMetaData(spatialDim);
    stk::unit_test_util::BulkDataTester stkMeshBulkData(stkMeshMetaData, communicator);

    // Elements 1 and 2 on proc 0, Elements 3 and 4 on proc 1
    // Elements 2 and 3 are shared because of nodes 9, 10, 11, 12
    // Element 2 is ghosted onto Proc 1, and Element 0 is ghosted onto Proc 0

    std::string exodusFileName = getOption("-i", "generated:1x1x4");
    populateBulkDataWithFile(exodusFileName, communicator, stkMeshBulkData);

    //============== Before starting, make sure mesh is parallel consistent

    std::vector<std::string> meshStart;
    getMeshLineByLine(stkMeshBulkData, meshStart);

    checkThatMeshIsParallelConsistent(stkMeshBulkData);

    int myProcId = stkMeshBulkData.parallel_rank();

    std::vector<stk::mesh::EntityId> edgeIds;
    std::vector<std::vector<stk::mesh::EntityId> > nodeIdsForEdge;
    std::vector<std::vector<stk::mesh::EntityId> > elementRelations;

    edgeIds.push_back(100+myProcId);
    nodeIdsForEdge.resize(edgeIds.size());
    elementRelations.resize(edgeIds.size());
    for ( size_t edge_index=0;edge_index<edgeIds.size();edge_index++)
    {
      nodeIdsForEdge[edge_index].resize(2);
      elementRelations[edge_index].resize(2);
      if ( myProcId == 0 )
      {
        nodeIdsForEdge[edge_index][0] = 5;
        nodeIdsForEdge[edge_index][1] = 6;
        elementRelations[edge_index][0] = 1;
        elementRelations[edge_index][1] = 2;
      }
      else
      {
        nodeIdsForEdge[edge_index][0] = 13;
        nodeIdsForEdge[edge_index][1] = 14;
        elementRelations[edge_index][0] = 3;
        elementRelations[edge_index][1] = 4;
      }
    }

    std::vector<stk::mesh::Entity> edgeEntities(edgeIds.size());
    stk::mesh::Part& edge_part = stkMeshBulkData.mesh_meta_data().declare_part("edge_part", stk::topology::EDGE_RANK);

    stkMeshBulkData.modification_begin();

    create_edges(stkMeshBulkData, edgeIds, nodeIdsForEdge, elementRelations, edgeEntities, edge_part);

    stkMeshBulkData.my_internal_resolve_shared_modify_delete();

    //============== Make sure we have 2 edges

    std::vector<size_t> globalCounts;
    stk::mesh::comm_mesh_counts(stkMeshBulkData, globalCounts);
    size_t numEdgesTotal = 2;
    EXPECT_EQ(numEdgesTotal, globalCounts[stk::topology::EDGE_RANK]);

    //============== Check results of IRGMD

    stkMeshBulkData.my_internal_resolve_ghosted_modify_delete();

    stkMeshBulkData.my_update_comm_list_based_on_changes_in_comm_map();

    std::vector<stk::mesh::Entity> shared_new;
    //        stkMeshBulkData.my_internal_update_distributed_index(stk::topology::EDGE_RANK, shared_new);
    stkMeshBulkData.my_internal_update_distributed_index(shared_new);
    //        ASSERT_EQ(0u, shared_new.size());

    stkMeshBulkData.my_resolve_ownership_of_modified_entities(shared_new);

    stkMeshBulkData.my_move_entities_to_proper_part_ownership( shared_new );

    stkMeshBulkData.my_add_comm_list_entries_for_entities( shared_new );

    stkMeshBulkData.my_internal_resolve_shared_membership();

    stkMeshBulkData.my_internal_regenerate_aura();

    //        std::vector<std::string> meshAfterAura;
    //        getMeshLineByLine(stkMeshBulkData, meshAfterAura);
    //        writeMesh(myProcId, "After::", meshAfterAura);

    checkItAllForThisGhostedCase(stkMeshBulkData);
  }
}

//BeginDocTest1
TEST(BulkDataModificationEnd, create_a_ghosted_edge_using_only_needed_pieces)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  int numProcs = -1;
  MPI_Comm_size(communicator, &numProcs);

  if(numProcs == 2)
  {
    //============== Load 1x1x4 mesh into Bulk Data

    const int spatialDim = 3;
    stk::mesh::MetaData stkMeshMetaData(spatialDim);
    stk::unit_test_util::BulkDataTester stkMeshBulkData(stkMeshMetaData, communicator);

    // Elements 1 and 2 on proc 0, Elements 3 and 4 on proc 1
    // Elements 2 and 3 are shared because of nodes 9, 10, 11, 12
    // Element 2 is ghosted onto Proc 1, and Element 0 is ghosted onto Proc 0

    std::string exodusFileName = getOption("-i", "generated:1x1x4");
    populateBulkDataWithFile(exodusFileName, communicator, stkMeshBulkData);

    checkThatMeshIsParallelConsistent(stkMeshBulkData);

    int myProcId = stkMeshBulkData.parallel_rank();

    std::vector<stk::mesh::EntityId> edgeIds;
    std::vector<std::vector<stk::mesh::EntityId> > nodeIdsForEdge;
    std::vector<std::vector<stk::mesh::EntityId> > elementRelations;

    edgeIds.push_back(100+myProcId);
    nodeIdsForEdge.resize(edgeIds.size());
    elementRelations.resize(edgeIds.size());

    // Proc 0 creates an edge with id=100 on nodes 5 and 6 related to elements 1 and 2
    // Proc 1 creates an edge with id=101 on nodes 13 and 14 related to elements 3 and 4
    for ( size_t i=0;i<edgeIds.size();i++)
    {
      nodeIdsForEdge[i].resize(2);
      elementRelations[i].resize(2);
      if ( myProcId == 0 )
      {
        nodeIdsForEdge[i][0] = 5;
        nodeIdsForEdge[i][1] = 6;
        elementRelations[i][0] = 1;
        elementRelations[i][1] = 2;
      }
      else
      {
        nodeIdsForEdge[i][0] = 13;
        nodeIdsForEdge[i][1] = 14;
        elementRelations[i][0] = 3;
        elementRelations[i][1] = 4;
      }
    }

    std::vector<stk::mesh::Entity> edgeEntities(edgeIds.size());
    stk::mesh::Part& edge_part = stkMeshBulkData.mesh_meta_data().declare_part("edge_part", stk::topology::EDGE_RANK);

    stkMeshBulkData.modification_begin();
    create_edges(stkMeshBulkData, edgeIds, nodeIdsForEdge, elementRelations, edgeEntities, edge_part);

    std::vector<size_t> globalCounts;
    stk::mesh::comm_mesh_counts(stkMeshBulkData, globalCounts);
    size_t numEdgesTotal = 2;
    EXPECT_EQ(numEdgesTotal, globalCounts[stk::topology::EDGE_RANK]);

    stkMeshBulkData.my_modification_end_for_entity_creation({stk::topology::EDGE_RANK});

    checkItAllForThisGhostedCase(stkMeshBulkData);
  }
}
//EndDocTest1


TEST(BulkDataModificationEnd, create_edges)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  int numProcs = -1;
  MPI_Comm_size(communicator, &numProcs);

  if(numProcs == 2)
  {
    //============== Load 1x1x4 mesh into Bulk Data

    const int spatialDim = 3;
    stk::mesh::MetaData stkMeshMetaData(spatialDim);
    stk::unit_test_util::BulkDataTester stkMeshBulkData(stkMeshMetaData, communicator);

    std::string exodusFileName = getOption("-i", "generated:1x1x4");
    populateBulkDataWithFile(exodusFileName, communicator, stkMeshBulkData);

    checkThatMeshIsParallelConsistent(stkMeshBulkData);

    stk::mesh::Part& edge_part = stkMeshBulkData.mesh_meta_data().declare_part("edge_part", stk::topology::EDGE_RANK);
    stk::mesh::create_edges(stkMeshBulkData, stkMeshBulkData.mesh_meta_data().universal_part(), &edge_part);

    checkThatMeshIsParallelConsistent(stkMeshBulkData);
  }
}

TEST(BulkDataModificationEnd, test_invalid_add_node_sharing)
{
  // this is to reproduce error seen by Steve Kennon, ticket #12829
  int numProcs = -1;
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  if ( numProcs == 3 )
  {
    int myProcId = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &myProcId);

    const unsigned spatial_dim = 3;
    stk::mesh::MetaData meta_data(spatial_dim);
    stk::mesh::Part &node_part = meta_data.get_topology_root_part(stk::topology::NODE);
    meta_data.commit();
    stk::unit_test_util::BulkDataTester mesh(meta_data, MPI_COMM_WORLD);
    mesh.modification_begin();

    stk::mesh::Entity node1 = mesh.declare_node(1, stk::mesh::ConstPartVector{&node_part});

    if ( myProcId == 2 )
    {
      mesh.add_node_sharing(node1, 0);
      mesh.add_node_sharing(node1, 1);
    }
    else if ( myProcId == 1)
    {
      mesh.add_node_sharing(node1, 2);
    }
    else if ( myProcId == 0 )
    {
      mesh.add_node_sharing(node1, 2);
    }

    mesh.my_internal_resolve_shared_modify_delete();
    mesh.my_internal_resolve_ghosted_modify_delete();
    mesh.my_update_comm_list_based_on_changes_in_comm_map();
    mesh.my_internal_resolve_parallel_create();
    EXPECT_THROW(mesh.check_sharing_comm_maps(), std::logic_error);
  }
}

TEST(BulkDataModificationEnd, create_edges_with_min_map)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  int numProcs = -1;
  MPI_Comm_size(communicator, &numProcs);

  if(numProcs == 2)
  {
    //============== Load 1x1x4 mesh into Bulk Data

    const int spatialDim = 3;
    stk::mesh::MetaData stkMeshMetaData(spatialDim);
    stk::unit_test_util::BulkDataTester stkMeshBulkData(stkMeshMetaData, communicator);

    // Elements 1 and 2 on proc 0, Elements 3 and 4 on proc 1
    // Elements 2 and 3 are shared because of nodes 9, 10, 11, 12
    // Element 2 is ghosted onto Proc 1, and Element 0 is ghosted onto Proc 0

    std::string exodusFileName = getOption("-i", "generated:1x1x4");
    populateBulkDataWithFile(exodusFileName, communicator, stkMeshBulkData);

    checkThatMeshIsParallelConsistent(stkMeshBulkData);

    stk::mesh::Part& face_part = stkMeshBulkData.mesh_meta_data().declare_part("face_part", stk::topology::FACE_RANK);
    stk::mesh::PartVector partVec;
    partVec.push_back(&face_part);
    stk::mesh::skin_mesh(stkMeshBulkData, stkMeshBulkData.mesh_meta_data().universal_part(), partVec);

    stk::mesh::Part& edge_part = stkMeshBulkData.mesh_meta_data().declare_part("edge_part", stk::topology::EDGE_RANK);
    stk::mesh::create_edges(stkMeshBulkData, stkMeshBulkData.mesh_meta_data().universal_part(), &edge_part);

    checkThatMeshIsParallelConsistent(stkMeshBulkData);
  }
}


void build_one_hex_on_p0(stk::unit_test_util::BulkDataTester& bulkData, const stk::mesh::MetaData& metaData)
{
  bulkData.modification_begin();
  if (bulkData.parallel_rank() == 0) {
    stk::mesh::Part& elem_part = metaData.get_topology_root_part(stk::topology::HEX_8);
    const int elemID = 1;
    stk::mesh::EntityIdVector nodeIDs = { 1, 2, 3, 4, 5, 6, 7, 8 };
    stk::mesh::declare_element(bulkData, elem_part, elemID, nodeIDs);
  }
  bulkData.modification_end();
}

void ghost_one_hex_to_p1(stk::unit_test_util::BulkDataTester & bulkData)
{
  // Ghost element and downward entities to other proc
  std::vector< std::pair<stk::mesh::Entity, int> > addGhost;
  if (bulkData.parallel_rank() == 0) {
    stk::mesh::Entity elem1 = bulkData.get_entity(stk::topology::ELEM_RANK, 1);
    addGhost.push_back(std::make_pair(elem1, 1));
  }
  bulkData.modification_begin();
  stk::mesh::Ghosting &ghosting = bulkData.create_ghosting("custom ghosting 1");
  bulkData.change_ghosting(ghosting, addGhost);
  bulkData.modification_end();
}

// DISABLED because automatic promotion of ghosted to shared is not working
// properly in pure STK apps.  This will throw with mesh consistency errors
// in debug.  Adding the manual add_node_sharing() calls fixes node sharing
// but leaves edges and faces unshared.
//
TEST(ModEndForEntityCreation, DISABLED_promotion_of_ghosted_to_shared)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  int numProcs = -1;
  MPI_Comm_size(communicator, &numProcs);

  if (numProcs != 2) return;

  const int spatialDim = 3;
  stk::mesh::MetaData metaData(spatialDim);
  stk::unit_test_util::BulkDataTester bulkData(metaData, communicator, stk::mesh::BulkData::NO_AUTO_AURA);
  const int p_rank = bulkData.parallel_rank();

  build_one_hex_on_p0(bulkData, metaData);

  stk::mesh::create_edges(bulkData);
  stk::mesh::create_faces(bulkData);

  ghost_one_hex_to_p1(bulkData);

  std::vector<size_t> localCounts;
  const Selector universalPart = Selector(metaData.universal_part());
  stk::mesh::count_entities(universalPart, bulkData, localCounts);

  // Count all entities on each proc
  EXPECT_EQ( 8u, localCounts[stk::topology::NODE_RANK]);
  EXPECT_EQ(12u, localCounts[stk::topology::EDGE_RANK]);
  EXPECT_EQ( 6u, localCounts[stk::topology::FACE_RANK]);
  EXPECT_EQ( 1u, localCounts[stk::topology::ELEM_RANK]);

  stk::mesh::Part& elem_part = metaData.get_topology_root_part(stk::topology::HEX_8);
  bulkData.modification_begin();
  if (p_rank == 1) {
    stk::mesh::EntityIdVector nodeIDs = {5, 6, 7, 8, 9, 10, 11, 12};
    stk::mesh::declare_element(bulkData, elem_part, 2, nodeIDs);

    // Rather not have to manually call this.  This will make nodes shared
    // but will leave edges and faces not shared.
    //bulkData.add_node_sharing(bulkData.get_entity(stk::mesh::EntityKey(stk::topology::NODE_RANK, 5)), 0);
    //bulkData.add_node_sharing(bulkData.get_entity(stk::mesh::EntityKey(stk::topology::NODE_RANK, 6)), 0);
    //bulkData.add_node_sharing(bulkData.get_entity(stk::mesh::EntityKey(stk::topology::NODE_RANK, 7)), 0);
    //bulkData.add_node_sharing(bulkData.get_entity(stk::mesh::EntityKey(stk::topology::NODE_RANK, 8)), 0);
  }
  if (p_rank == 0) {
    // Rather not have to manually call this.  This will make nodes shared
    // but will leave edges and faces not shared.
    //bulkData.add_node_sharing(bulkData.get_entity(stk::mesh::EntityKey(stk::topology::NODE_RANK, 5)), 1);
    //bulkData.add_node_sharing(bulkData.get_entity(stk::mesh::EntityKey(stk::topology::NODE_RANK, 6)), 1);
    //bulkData.add_node_sharing(bulkData.get_entity(stk::mesh::EntityKey(stk::topology::NODE_RANK, 7)), 1);
    //bulkData.add_node_sharing(bulkData.get_entity(stk::mesh::EntityKey(stk::topology::NODE_RANK, 8)), 1);
  }

  std::vector<stk::mesh::EntityRank> entityRanksToConsider;
  for (stk::mesh::EntityRank rank = stk::mesh::EntityRank::BEGIN_RANK; rank <= stk::mesh::EntityRank::ELEM_RANK; ++rank) {
    entityRanksToConsider.push_back(rank);
  }
  bulkData.my_modification_end_for_entity_creation(entityRanksToConsider);

  const Selector sharedPart = Selector(metaData.globally_shared_part());
  stk::mesh::count_entities(sharedPart, bulkData, localCounts);

  // Count shared entities on each proc
  EXPECT_EQ( 4u, localCounts[stk::topology::NODE_RANK]);
  EXPECT_EQ( 4u, localCounts[stk::topology::EDGE_RANK]);
  EXPECT_EQ( 1u, localCounts[stk::topology::FACE_RANK]);
  EXPECT_EQ( 0u, localCounts[stk::topology::ELEM_RANK]);
}

std::shared_ptr<stk::mesh::BulkData> create_mesh(stk::ParallelMachine comm,
                                                 unsigned spatialDim,
                                                 stk::mesh::BulkData::AutomaticAuraOption auraOption)
{
  stk::mesh::MeshBuilder builder(comm);
  builder.set_spatial_dimension(spatialDim);
  builder.set_aura_option(auraOption);
  std::shared_ptr<stk::mesh::BulkData> bulk = builder.create();
  return bulk;
}

TEST(TestModificationEnd, destroySharedNode_twoSharers_deinducePartMembership)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) return;

  const unsigned spatialDim = 2;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = create_mesh(MPI_COMM_WORLD, spatialDim, stk::mesh::BulkData::NO_AUTO_AURA);
  stk::mesh::BulkData& bulk = *bulkPtr;
  stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  std::string meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
                         "1,2,QUAD_4_2D,2,3,6,5,block_2";
  std::vector<double> coordinates = {
    0,0, 1,0, 2,0,
    0,1, 1,1, 2,1
  };
  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

  bulk.modification_begin();
  if (bulk.parallel_rank() == 1) {
    bulk.destroy_entity(bulk.get_entity(stk::topology::ELEM_RANK, 2));
    bulk.destroy_entity(bulk.get_entity(stk::topology::NODE_RANK, 2));
    bulk.destroy_entity(bulk.get_entity(stk::topology::NODE_RANK, 3));
    bulk.destroy_entity(bulk.get_entity(stk::topology::NODE_RANK, 6));
    bulk.destroy_entity(bulk.get_entity(stk::topology::NODE_RANK, 5));
  }
  bulk.modification_end();

  if (bulk.parallel_rank() == 0) {
    const stk::mesh::Entity node2 = bulk.get_entity(stk::topology::NODE_RANK, 2);
    const stk::mesh::Entity node5 = bulk.get_entity(stk::topology::NODE_RANK, 5);
    const stk::mesh::Part & block1 = *meta.get_part("BLOCK_1");
    const stk::mesh::Part & block2 = *meta.get_part("BLOCK_2");

    EXPECT_TRUE (bulk.bucket(node2).member(block1));
    EXPECT_FALSE(bulk.bucket(node2).member(block2));

    EXPECT_TRUE (bulk.bucket(node5).member(block1));
    EXPECT_FALSE(bulk.bucket(node5).member(block2));
  }
}

TEST(TestModificationEnd, destroySharedNode_threeSharers_deinducePartMembership)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 3) return;

  const unsigned spatialDim = 2;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = create_mesh(MPI_COMM_WORLD, spatialDim, stk::mesh::BulkData::NO_AUTO_AURA);
  stk::mesh::BulkData& bulk = *bulkPtr;
  stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  std::string meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
                         "1,2,QUAD_4_2D,2,3,6,5,block_2\n"
                         "2,3,QUAD_4_2D,4,5,8,7,block_3";
  std::vector<double> coordinates = {
    0,0, 1,0, 2,0,
    0,1, 1,1, 2,1,
    0,2, 1,2
  };
  stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

  bulk.modification_begin();
  if (bulk.parallel_rank() == 2) {
    bulk.destroy_entity(bulk.get_entity(stk::topology::ELEM_RANK, 3));
    bulk.destroy_entity(bulk.get_entity(stk::topology::NODE_RANK, 4));
    bulk.destroy_entity(bulk.get_entity(stk::topology::NODE_RANK, 5));
    bulk.destroy_entity(bulk.get_entity(stk::topology::NODE_RANK, 7));
    bulk.destroy_entity(bulk.get_entity(stk::topology::NODE_RANK, 8));
  }
  bulk.modification_end();

  if (bulk.parallel_rank() == 0) {
    const stk::mesh::Entity node2 = bulk.get_entity(stk::topology::NODE_RANK, 2);
    const stk::mesh::Entity node4 = bulk.get_entity(stk::topology::NODE_RANK, 4);
    const stk::mesh::Entity node5 = bulk.get_entity(stk::topology::NODE_RANK, 5);
    const stk::mesh::Part & block1 = *meta.get_part("BLOCK_1");
    const stk::mesh::Part & block2 = *meta.get_part("BLOCK_2");
    const stk::mesh::Part & block3 = *meta.get_part("BLOCK_3");

    EXPECT_TRUE(bulk.bucket(node2).member(block1));
    EXPECT_TRUE(bulk.bucket(node2).member(block2));

    EXPECT_TRUE (bulk.bucket(node4).member(block1));
    EXPECT_FALSE(bulk.bucket(node4).member(block3));

    EXPECT_TRUE (bulk.bucket(node5).member(block1));
    EXPECT_TRUE (bulk.bucket(node5).member(block2));
    EXPECT_FALSE(bulk.bucket(node5).member(block3));
  }
  else if (bulk.parallel_rank() == 1) {
    const stk::mesh::Entity node2 = bulk.get_entity(stk::topology::NODE_RANK, 2);
    const stk::mesh::Entity node5 = bulk.get_entity(stk::topology::NODE_RANK, 5);
    const stk::mesh::Part & block1 = *meta.get_part("BLOCK_1");
    const stk::mesh::Part & block2 = *meta.get_part("BLOCK_2");
    const stk::mesh::Part & block3 = *meta.get_part("BLOCK_3");

    EXPECT_TRUE(bulk.bucket(node2).member(block1));
    EXPECT_TRUE(bulk.bucket(node2).member(block2));

    EXPECT_TRUE (bulk.bucket(node5).member(block1));
    EXPECT_TRUE (bulk.bucket(node5).member(block2));
    EXPECT_FALSE(bulk.bucket(node5).member(block3));
  }
  else if (bulk.parallel_rank() == 2) {
    const stk::mesh::Entity node4 = bulk.get_entity(stk::topology::NODE_RANK, 4);
    const stk::mesh::Entity node5 = bulk.get_entity(stk::topology::NODE_RANK, 5);
    EXPECT_FALSE(bulk.is_valid(node4));
    EXPECT_FALSE(bulk.is_valid(node5));
  }
}

} } } // namespace stk mesh unit_test
