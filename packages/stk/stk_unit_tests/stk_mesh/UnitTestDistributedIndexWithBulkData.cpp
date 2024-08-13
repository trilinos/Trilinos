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

#include <gtest/gtest.h>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine, etc
#include <stk_util/parallel/DistributedIndex.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/baseImpl/MeshImplUtils.hpp>
#include <stk_unit_test_utils/StkMeshFromGeneratedMesh.hpp>

namespace {

class DistributedIndexTester : public stk::parallel::DistributedIndex
{
public:
  DistributedIndexTester( stk::ParallelMachine comm ,const KeySpanVector & partition_spans ) :
    DistributedIndex(comm, partition_spans)  {}

  KeyProcVector getKeys() const { return m_key_usage ; }
};



void updateDistributedIndexUsingStkMesh(stk::unit_test_util::BulkDataTester &stkMeshBulkData, stk::parallel::DistributedIndex& distributedIndex)
{
  stk::parallel::DistributedIndex::KeyTypeVector local_created_or_modified; // only store locally owned/shared entities

  size_t num_created_or_modified = 0;

  for(stk::mesh::const_entity_iterator i = stkMeshBulkData.begin_entities(stk::topology::NODE_RANK); i != stkMeshBulkData.end_entities(stk::topology::NODE_RANK); ++i)
  {
    stk::mesh::Entity entity = i->second;

    if(stkMeshBulkData.owned_closure(entity))
    {
      ++num_created_or_modified;
    }
  }

  local_created_or_modified.reserve(num_created_or_modified);

  for(stk::mesh::const_entity_iterator i = stkMeshBulkData.begin_entities(stk::topology::NODE_RANK); i != stkMeshBulkData.end_entities(stk::topology::NODE_RANK); ++i)
  {
    stk::mesh::Entity entity = i->second;

    if(stkMeshBulkData.owned_closure(entity))
    {
      local_created_or_modified.push_back(stkMeshBulkData.entity_key(entity));
    }
  }

  stk::parallel::DistributedIndex::KeyTypeVector::const_iterator begin = local_created_or_modified.begin();
  stk::parallel::DistributedIndex::KeyTypeVector::const_iterator end = local_created_or_modified.end();
  distributedIndex.update_keys(begin, end);
}

TEST( UnderstandingDistributedIndex, WithoutStkMeshBulkData)
{
  stk::ParallelMachine communicator = MPI_COMM_WORLD;
  int procCount = stk::parallel_machine_size(communicator);
  int myProc = stk::parallel_machine_rank(communicator);

  if(procCount == 2)
  {
    const std::string generatedMeshSpec = "generated:2x2x2|sideset:xXyYzZ|nodeset:xXyYzZ";
    stk::unit_test_util::StkMeshCreator stkMesh(generatedMeshSpec, communicator);

    stk::mesh::MetaData &stkMeshMetaData = *stkMesh.getMetaData();
    stk::unit_test_util::BulkDataTester &stkMeshBulkData = *stkMesh.getBulkData();

    stk::parallel::DistributedIndex::KeySpanVector spans = stk::mesh::impl::convert_entity_keys_to_spans(stkMeshMetaData);
    DistributedIndexTester distributedIndex(communicator, spans);

    updateDistributedIndexUsingStkMesh(stkMeshBulkData, distributedIndex);

    const unsigned rankCount = stkMeshMetaData.entity_rank_count();
    EXPECT_EQ(4u, rankCount);

    for(size_t i = 0; i < spans.size(); i++)
    {
      stk::parallel::DistributedIndex::KeyType keyMin = spans[i].first;
      stk::parallel::DistributedIndex::KeyType keyMax = spans[i].second;
      stk::mesh::EntityKey key_1(static_cast<stk::mesh::EntityKey::entity_key_t>((keyMin)));
      stk::mesh::EntityKey key_2(static_cast<stk::mesh::EntityKey::entity_key_t>((keyMax)));

      stk::mesh::EntityRank eRank = static_cast<stk::mesh::EntityRank>(i);

      EXPECT_EQ(eRank, key_1.rank());
      EXPECT_EQ(eRank, key_2.rank());
    }

    ///////////////////////////////////////////////////////////////////////////////

    std::vector<stk::parallel::DistributedIndex::KeyTypeVector> requested_key_types;
    std::vector<size_t> requests(rankCount, 0);
    if ( myProc == 0 )
    {
      requests[0] = 4;
    }

    distributedIndex.generate_new_keys(requests, requested_key_types);

    size_t numNodesInMesh = 27;

    for(size_t i = 0; i < requested_key_types.size(); i++)
    {
      stk::mesh::EntityRank eRank = static_cast<stk::mesh::EntityRank>(i);
      if(myProc == 0)
      {
        for(size_t j = 0; j < requested_key_types[i].size(); j++)
        {
          stk::mesh::EntityKey keyType(static_cast<stk::mesh::EntityKey::entity_key_t>(requested_key_types[i][j]));
          EXPECT_EQ(eRank, keyType.rank());
          size_t goldId = numNodesInMesh + j + 1;
          EXPECT_EQ(goldId, keyType.id());
        }
      }
    }

    stk::parallel::DistributedIndex::KeyProcVector keys = distributedIndex.getKeys();
    size_t numNewNodes = requests[stk::topology::NODE_RANK];
    size_t numNodesLocalProc0 = 18;
    size_t numNodesLocalProc1 = 18;
    if(myProc == 0)
    {
      EXPECT_EQ(numNodesLocalProc0+numNodesLocalProc1+numNewNodes, keys.size());
    }
    MPI_Barrier(communicator);
  }
}

TEST( UnderstandingDistributedIndex, ViaStkMeshBulkData)
{
  stk::ParallelMachine communicator = MPI_COMM_WORLD;
  int procCount = stk::parallel_machine_size(communicator);
  int myProc = stk::parallel_machine_rank(communicator);

  if(procCount == 2)
  {
    const std::string generatedMeshSpec = "generated:2x2x2|sideset:xXyYzZ|nodeset:xXyYzZ";
    stk::unit_test_util::StkMeshCreator stkMesh(generatedMeshSpec, communicator);

    stk::mesh::MetaData &stkMeshMetaData = *stkMesh.getMetaData();
    stk::mesh::Part &line2_part = stkMeshMetaData.get_topology_root_part(stk::topology::LINE_2);
    stk::mesh::Part &quad4_part = stkMeshMetaData.get_topology_root_part(stk::topology::QUAD_4);
    stk::mesh::Part &hex8_part = stkMeshMetaData.get_topology_root_part(stk::topology::HEX_8);

    stk::mesh::BulkData &stkMeshBulkData = *stkMesh.getBulkData();

    MPI_Barrier(communicator);

    const int rankCount = 4; // nodes, edges, faces, and elements
    std::vector<size_t> requests(rankCount, 0);
    if ( myProc == 0 )
    {
      requests[0] = 4;
      requests[1] = 5;
      requests[2] = 6;
      requests[3] = 7;
    }

    size_t totalCount = 0;
    for(size_t i = 0; i < requests.size(); i++)
    {
      totalCount += requests[i];
    }

    std::vector<stk::mesh::Entity> requested_entities;

    stkMeshBulkData.modification_begin();
    stkMeshBulkData.generate_new_entities(requests, requested_entities);

    stk::mesh::PartVector element_topo_parts, face_topo_parts, edge_topo_parts, empty_parts;
    edge_topo_parts.push_back(&line2_part);
    face_topo_parts.push_back(&quad4_part);
    element_topo_parts.push_back(&hex8_part);

    //each entity in requested_entities needs to be connected to at least one node
    //because stk-mesh requires that.
    size_t idx = requests[0];
    stk::mesh::Entity node = totalCount>0 ? requested_entities[0] : stk::mesh::Entity();
    for(size_t i=1; i<requests.size(); ++i) {
      for(size_t j=0; j<requests[i]; ++j) {
        stk::mesh::Entity entity = requested_entities[idx++];
        stk::topology entity_topo;
        stk::mesh::EntityRank entity_rank = stkMeshBulkData.entity_rank(entity);
        switch (entity_rank)
        {
        case stk::topology::EDGE_RANK:
          stkMeshBulkData.change_entity_parts(entity, edge_topo_parts, empty_parts);
          entity_topo = stk::topology::LINE_2;
          break;
        case stk::topology::FACE_RANK:
          stkMeshBulkData.change_entity_parts(entity, face_topo_parts, empty_parts);
          entity_topo = stk::topology::QUAD_4;
          break;
        case stk::topology::ELEMENT_RANK:
          stkMeshBulkData.change_entity_parts(entity, element_topo_parts, empty_parts);
          entity_topo = stk::topology::HEX_8;
          break;
        default:
          break;
        }
        if (entity_rank != stk::topology::NODE_RANK)
        {
          for (unsigned ord = 0; ord < entity_topo.num_nodes(); ++ord)
          {
            stkMeshBulkData.declare_relation(entity, node, ord);
          }
        }
      }
    }

    stkMeshBulkData.modification_end();

    EXPECT_EQ(totalCount, requested_entities.size());

    size_t numNodesInMesh = (2+1)*(2+1)*(2+1);
    size_t numEdgesInMesh = 0;
    size_t numElementsInMesh = 2*2*2;

    size_t nodeCounter=0;
    size_t edgeCounter=0;
    size_t elemCounter=0;

    for(size_t i = 0; i < requested_entities.size(); i++)
    {
      if(myProc == 0)
      {
        if (stkMeshBulkData.entity_rank(requested_entities[i]) == stk::topology::NODE_RANK )
        {
          nodeCounter++;
          EXPECT_EQ(numNodesInMesh+nodeCounter, stkMeshBulkData.identifier(requested_entities[i]));
        }
        else if (stkMeshBulkData.entity_rank(requested_entities[i]) == stk::topology::EDGE_RANK )
        {
          edgeCounter++;
          EXPECT_EQ(numEdgesInMesh+edgeCounter, stkMeshBulkData.identifier(requested_entities[i]));
        }
        else if (stkMeshBulkData.entity_rank(requested_entities[i]) == stk::topology::ELEMENT_RANK )
        {
          elemCounter++;
          EXPECT_EQ(numElementsInMesh+elemCounter, stkMeshBulkData.identifier(requested_entities[i]));
        }
      }
    }

    EXPECT_EQ(requests[stk::topology::NODE_RANK], nodeCounter);
    EXPECT_EQ(requests[stk::topology::EDGE_RANK], edgeCounter);
    EXPECT_EQ(requests[stk::topology::ELEMENT_RANK], elemCounter);

    MPI_Barrier(communicator);
  }
}

TEST(UnderstandingDistributedIndex, TestSharedAndGhostedAndOwnedEntitiesWithoutAnyModification)
{
  stk::ParallelMachine communicator = MPI_COMM_WORLD;
  int procCount = stk::parallel_machine_size(communicator);
  int myProc = stk::parallel_machine_rank(communicator);

  if(procCount == 2)
  {
    const std::string generatedMeshSpec = "generated:2x2x2|sideset:xXyYzZ|nodeset:xXyYzZ";
    stk::unit_test_util::StkMeshCreator stkMesh(generatedMeshSpec, communicator);

    stk::unit_test_util::BulkDataTester &stkMeshBulkData = *stkMesh.getBulkData();

    int sharedNodeIds[] = { 10, 11, 12, 13, 14, 15, 16, 17, 18 };

    int otherProcId = 1;
    if ( myProc == 1 )
    {
      otherProcId = 0;
    }

    size_t numSharedNodes=9;
    for (size_t i=0;i<numSharedNodes;i++)
    {
      stk::mesh::Entity sharedNode = stkMeshBulkData.get_entity(stk::topology::NODE_RANK, sharedNodeIds[i]);
      stk::mesh::Bucket& bucket = stkMeshBulkData.bucket(sharedNode);
      EXPECT_TRUE(bucket.shared());
      stk::mesh::PairIterEntityComm commStuff = stkMeshBulkData.my_internal_entity_comm_map_shared(stkMeshBulkData.entity_key(sharedNode));
      EXPECT_TRUE(commStuff.size() == 1);
      EXPECT_TRUE((*commStuff.first).proc == otherProcId);
      size_t sharedButNotGhostedId = 0;
      EXPECT_TRUE((*commStuff.first).ghost_id == sharedButNotGhostedId);
    }

    int ghostedNodeIds_onProc1[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    int ghostedNodeIds_onProc0[] = { 19, 20, 21, 22, 23, 24, 25, 26, 27 };

    int *ghostedNodeIds=ghostedNodeIds_onProc1;
    if ( myProc == 0 )
    {
      ghostedNodeIds = ghostedNodeIds_onProc0;
    }

    size_t numGhostedNodes=9;
    for (size_t i=0;i<numGhostedNodes;i++)
    {
      stk::mesh::Entity ghostedNode = stkMeshBulkData.get_entity(stk::topology::NODE_RANK, ghostedNodeIds[i]);
      stk::mesh::Bucket& bucket = stkMeshBulkData.bucket(ghostedNode);
      bool isGhosted = !bucket.shared() && !bucket.owned();
      EXPECT_TRUE(isGhosted);
      stk::mesh::PairIterEntityComm commStuff = stkMeshBulkData.my_internal_entity_comm_map(stkMeshBulkData.entity_key(ghostedNode));
      size_t numProcsOtherThanMeOnWhichEntityIsAGhost = 1;
      EXPECT_TRUE(commStuff.size() == numProcsOtherThanMeOnWhichEntityIsAGhost);
      EXPECT_TRUE((*commStuff.first).proc == otherProcId);
      size_t auraGhostedId = 1; // aura is not custom ghosted (automatically generated one layer of ghosting)
      EXPECT_TRUE((*commStuff.first).ghost_id == auraGhostedId);
    }

    int elementsOnProc0[] = { 1, 2, 3, 4 }; // ghosted on proc 1
    int elementsOnProc1[] = { 5, 6, 7, 8 }; // ghosted on proc 0
    int *ghostedElements=elementsOnProc0;
    if ( myProc == 0 )
    {
      ghostedElements = elementsOnProc1;
    }

    size_t numGhostedElements=4;
    for (size_t i=0;i<numGhostedElements;i++)
    {
      stk::mesh::Entity ghostedElement = stkMeshBulkData.get_entity(stk::topology::ELEMENT_RANK, ghostedElements[i]);
      stk::mesh::Bucket& bucket = stkMeshBulkData.bucket(ghostedElement);
      bool isGhosted = !bucket.shared() && !bucket.owned();
      EXPECT_TRUE(isGhosted);
      stk::mesh::PairIterEntityComm commStuff = stkMeshBulkData.my_internal_entity_comm_map(stkMeshBulkData.entity_key(ghostedElement));
      size_t numProcsOtherThanMeOnWhichEntityIsAGhost = 1;
      EXPECT_TRUE(commStuff.size() == numProcsOtherThanMeOnWhichEntityIsAGhost);
      EXPECT_TRUE((*commStuff.first).proc == otherProcId);
      size_t auraGhostedId = 1; // aura is not custom ghosted (automatically generated one layer of ghosting)
      EXPECT_TRUE((*commStuff.first).ghost_id == auraGhostedId);
    }

    int *ownedElements = elementsOnProc1;
    if ( myProc == 0 )
    {
      ownedElements = elementsOnProc0;
    }

    size_t numOwnedElements=4;
    for (size_t i=0;i<numOwnedElements;i++)
    {
      stk::mesh::Entity ownedElement = stkMeshBulkData.get_entity(stk::topology::ELEMENT_RANK, ownedElements[i]);
      stk::mesh::Bucket& bucket = stkMeshBulkData.bucket(ownedElement);
      EXPECT_TRUE(bucket.owned());
      stk::mesh::PairIterEntityComm commStuff = stkMeshBulkData.my_internal_entity_comm_map(stkMeshBulkData.entity_key(ownedElement));
      size_t numProcsOtherThanMeOnWhichEntityIsAGhost = 1;
      EXPECT_TRUE(commStuff.size() == numProcsOtherThanMeOnWhichEntityIsAGhost);
      EXPECT_TRUE((*commStuff.first).proc == otherProcId);
      size_t auraGhostedId = 1; // aura is not custom ghosted (automatically generated one layer of ghosting)
      EXPECT_TRUE((*commStuff.first).ghost_id == auraGhostedId);
    }
  }
}

void testSharedNodesFor2x2x4MeshForTwoProcs(const int myProc, const stk::unit_test_util::BulkDataTester &stkMeshBulkData)
{
  int sharedNodeIds[] = { 19, 20, 21, 22, 23, 24, 25, 26, 27 };

  int otherProcId = 1;
  if ( myProc == 1 )
  {
    otherProcId = 0;
  }

  size_t numSharedNodes=9;
  for (size_t i=0;i<numSharedNodes;i++)
  {
    stk::mesh::Entity sharedNode = stkMeshBulkData.get_entity(stk::topology::NODE_RANK, sharedNodeIds[i]);
    stk::mesh::Bucket& bucket = stkMeshBulkData.bucket(sharedNode);
    EXPECT_TRUE(bucket.shared());
    stk::mesh::PairIterEntityComm commStuff = stkMeshBulkData.my_internal_entity_comm_map_shared(stkMeshBulkData.entity_key(sharedNode));
    EXPECT_TRUE(commStuff.size() == 1);
    EXPECT_TRUE((*commStuff.first).proc == otherProcId);
    size_t sharedButNotGhostedId = 0;
    EXPECT_TRUE((*commStuff.first).ghost_id == sharedButNotGhostedId);
  }
}

TEST(UnderstandingDistributedIndex, GhostAnElement)
{
  stk::ParallelMachine communicator = MPI_COMM_WORLD;
  int procCount = stk::parallel_machine_size(communicator);
  int myProc = stk::parallel_machine_rank(communicator);

  if(procCount == 2)
  {
    const std::string generatedMeshSpec = "generated:2x2x4|sideset:xXyYzZ|nodeset:xXyYzZ";
    stk::unit_test_util::StkMeshCreator stkMesh(generatedMeshSpec, communicator);

    stk::unit_test_util::BulkDataTester &stkMeshBulkData = *stkMesh.getBulkData();

    testSharedNodesFor2x2x4MeshForTwoProcs(myProc, stkMeshBulkData);

    int otherProcId = 1;
    if ( myProc == 1 )
    {
      otherProcId = 0;
    }

    // elements 1-8 are proc 0
    // elements 9-16 are proc 1
    // 5-8 are ghosted on proc 1
    // 9-12 are ghosted on proc 0

    // proc 1 needs to ghost element 13 to proc 0

    //BEGIN_DOC_FOR_GHOSTING_AN_ELEMENT
    int ghostedElementId = 13;
    int processorThatOwnsElement13 = 1;
    int processorThatReceivesGhostedElement13 = 0;

    stkMeshBulkData.modification_begin();
    stk::mesh::Ghosting &ghosting = stkMeshBulkData.create_ghosting("Ghost Element 13");
    std::vector< std::pair<stk::mesh::Entity, int> > ghostingStruct;
    if ( myProc == processorThatOwnsElement13 )
    {
      stk::mesh::Entity element = stkMeshBulkData.get_entity(stk::topology::ELEMENT_RANK, ghostedElementId);
      ghostingStruct.push_back(std::make_pair(element, processorThatReceivesGhostedElement13));
    }
    stkMeshBulkData.change_ghosting(ghosting, ghostingStruct);
    stkMeshBulkData.modification_end();

    if ( myProc == processorThatReceivesGhostedElement13 )
    {
      stk::mesh::Entity element = stkMeshBulkData.get_entity(stk::topology::ELEMENT_RANK, ghostedElementId);
      stk::mesh::Bucket& bucket = stkMeshBulkData.bucket(element);
      bool isGhosted = !bucket.shared() && !bucket.owned();
      EXPECT_TRUE(isGhosted);

      int nodeConnectedToGhostedElement = 37;
      stk::mesh::Entity node = stkMeshBulkData.get_entity(stk::topology::NODE_RANK, nodeConnectedToGhostedElement);
      stk::mesh::Bucket& nodeBucket = stkMeshBulkData.bucket(node);
      isGhosted = !nodeBucket.shared() && !nodeBucket.owned();
      EXPECT_TRUE(isGhosted);
    }
    //END_DOC_FOR_GHOSTING_AN_ELEMENT

    if ( myProc == processorThatReceivesGhostedElement13 )
    {
      stk::mesh::Entity element = stkMeshBulkData.get_entity(stk::topology::ELEMENT_RANK, ghostedElementId);
      stk::mesh::PairIterEntityComm commStuff = stkMeshBulkData.my_internal_entity_comm_map(stkMeshBulkData.entity_key(element));
      size_t numProcsToCommunicateWithForEntity = 1;
      EXPECT_TRUE(commStuff.size() == numProcsToCommunicateWithForEntity);
      EXPECT_TRUE((*commStuff.first).proc == otherProcId);

      size_t customGhosting = ghosting.ordinal();
      EXPECT_TRUE((*commStuff.first).ghost_id == customGhosting);
    }
  }
}

TEST(UnderstandingDistributedIndex, KillAGhostedElement)
{
  stk::ParallelMachine communicator = MPI_COMM_WORLD;
  int procCount = stk::parallel_machine_size(communicator);
  int myProc = stk::parallel_machine_rank(communicator);

  if(procCount == 2)
  {
    const std::string generatedMeshSpec = "generated:2x2x4|sideset:xXyYzZ|nodeset:xXyYzZ";
    stk::unit_test_util::StkMeshCreator stkMesh(generatedMeshSpec, communicator);

    stk::unit_test_util::BulkDataTester &stkMeshBulkData = *stkMesh.getBulkData();

    testSharedNodesFor2x2x4MeshForTwoProcs(myProc, stkMeshBulkData);

    int otherProcId = 1;
    if ( myProc == 1 )
    {
      otherProcId = 0;
    }

    // elements 1-8 are proc 0
    // elements 9-16 are proc 1
    // 5-8 are ghosted on proc 1
    // 9-12 are ghosted on proc 0

    // proc 1 needs to ghost element 13 to proc 0
    int elementIdToKill = 8;
    int owningProc = 0;
    int ghostedToProc = 1;

    stk::mesh::Entity elementToKill = stkMeshBulkData.get_entity(stk::topology::ELEMENT_RANK, elementIdToKill);

    stkMeshBulkData.modification_begin();
    if ( myProc == ghostedToProc )
    {
      stkMeshBulkData.destroy_entity(elementToKill);
    }
    stkMeshBulkData.modification_end();

    if ( myProc == ghostedToProc )
    {
      stk::mesh::Entity element = stkMeshBulkData.get_entity(stk::topology::ELEMENT_RANK, elementIdToKill);
      stk::mesh::Bucket& bucket = stkMeshBulkData.bucket(element);
      bool isGhosted = !bucket.shared() && !bucket.owned();
      EXPECT_TRUE(isGhosted);

      stk::mesh::PairIterEntityComm commStuff = stkMeshBulkData.my_internal_entity_comm_map(stkMeshBulkData.entity_key(element));
      size_t numProcsToCommunicateWithForEntity = 1;
      EXPECT_TRUE(commStuff.size() == numProcsToCommunicateWithForEntity);
      EXPECT_TRUE((*commStuff.first).proc == otherProcId);

      size_t auraGhostingId = 1;
      EXPECT_TRUE((*commStuff.first).ghost_id == auraGhostingId);
    }

    stkMeshBulkData.modification_begin();
    if ( myProc == owningProc )
    {
      stkMeshBulkData.destroy_entity(elementToKill);
    }
    stkMeshBulkData.modification_end();

    stk::mesh::Entity element = stkMeshBulkData.get_entity(stk::topology::ELEMENT_RANK, elementIdToKill);
    EXPECT_FALSE(stkMeshBulkData.is_valid(element));
  }
}

TEST(UnderstandingDistributedIndex, CreateDisconnectedElement)
{
  stk::ParallelMachine communicator = MPI_COMM_WORLD;
  int procCount = stk::parallel_machine_size(communicator);
  int myProc = stk::parallel_machine_rank(communicator);

  if(procCount == 2)
  {
    const std::string generatedMeshSpec = "generated:2x2x4|sideset:xXyYzZ|nodeset:xXyYzZ";
    stk::unit_test_util::StkMeshCreator stkMesh(generatedMeshSpec, communicator);

    stk::unit_test_util::BulkDataTester &stkMeshBulkData = *stkMesh.getBulkData();

    testSharedNodesFor2x2x4MeshForTwoProcs(myProc, stkMeshBulkData);

    // elements 1-8 are proc 0
    // elements 9-16 are proc 1
    // 5-8 are ghosted on proc 1
    // 9-12 are ghosted on proc 0

    int owningProc = 0;

    stk::mesh::MetaData &stkMeshMetaData = *stkMesh.getMetaData();

    stkMeshBulkData.modification_begin();
    std::vector<size_t> requestsForNewEntities(4, 0);
    std::vector<stk::mesh::Entity> generatedEntities;
    if(myProc == owningProc)
    {
      requestsForNewEntities[stk::topology::NODE_RANK] = 8;
      requestsForNewEntities[stk::topology::ELEMENT_RANK] = 1;
    }
    stkMeshBulkData.generate_new_entities(requestsForNewEntities, generatedEntities);

    stk::mesh::EntityId elementId = static_cast<stk::mesh::EntityId>(-1);
    if(myProc == owningProc)
    {
      stk::mesh::EntityIdVector nodeIds;
      for(size_t i=0; i < generatedEntities.size(); i++)
      {
        stk::mesh::Entity current = generatedEntities[i];
        stk::mesh::EntityRank rank = stkMeshBulkData.entity_rank(current);
        stk::mesh::EntityId id = stkMeshBulkData.identifier(current);
        if(rank == stk::topology::NODE_RANK)
        {
          nodeIds.push_back(id);
        }
        else if(rank == stk::topology::ELEMENT_RANK)
        {
          elementId = id;
        }
      }
      stk::mesh::Part *block1 = stkMeshMetaData.get_part("block_1");
      stk::mesh::PartVector justBlock1Really(1, block1);

      stkMeshBulkData.change_entity_parts(generatedEntities[0], justBlock1Really);
      stk::mesh::declare_element(stkMeshBulkData, justBlock1Really, elementId, nodeIds);
    }
    stkMeshBulkData.modification_end();

    if ( myProc == owningProc )
    {
      stk::mesh::Entity ownedElement = stkMeshBulkData.get_entity(stk::topology::ELEMENT_RANK, elementId);
      stk::mesh::Bucket& bucket = stkMeshBulkData.bucket(ownedElement);
      EXPECT_TRUE(bucket.owned());
    }
    else
    {
      stk::mesh::Entity ownedElement = stkMeshBulkData.get_entity(stk::topology::ELEMENT_RANK, elementId);
      EXPECT_FALSE(stkMeshBulkData.is_valid(ownedElement));
    }
  }
}

void testElementMove(int fromProc, int toProc, int myProc, int elementToMoveId, stk::mesh::BulkData &stkMeshBulkData)
{
  //BEGIN_DOC_FOR_ELEMENT_MOVE
  stk::mesh::Entity elementToMove = stkMeshBulkData.get_entity(stk::topology::ELEMENT_RANK, elementToMoveId);

  std::vector<std::pair<stk::mesh::Entity, int> > entityProcPairs;
  if(myProc == fromProc)
  {
    entityProcPairs.push_back(std::make_pair(elementToMove, toProc));
  }
  stkMeshBulkData.change_entity_owner(entityProcPairs);
  //END_DOC_FOR_ELEMENT_MOVE
}

TEST(UnderstandingDistributedIndex, MoveAnElement)
{
  stk::ParallelMachine communicator = MPI_COMM_WORLD;
  int procCount = stk::parallel_machine_size(communicator);
  int myProc = stk::parallel_machine_rank(communicator);

  if(procCount == 2)
  {
    const std::string generatedMeshSpec = "generated:2x2x4|sideset:xXyYzZ|nodeset:xXyYzZ";
    stk::unit_test_util::StkMeshCreator stkMesh(generatedMeshSpec, communicator);

    stk::unit_test_util::BulkDataTester &stkMeshBulkData = *stkMesh.getBulkData();

    testSharedNodesFor2x2x4MeshForTwoProcs(myProc, stkMeshBulkData);

    // elements 1-8 are proc 0
    // elements 9-16 are proc 1
    // 5-8 are ghosted on proc 1
    // 9-12 are ghosted on proc 0

    int fromProc = 0;
    int toProc = 1;
    int elementToMoveId = 2;

    testElementMove(fromProc, toProc, myProc, elementToMoveId, stkMeshBulkData);
    stk::mesh::Entity elementToMove = stkMeshBulkData.get_entity(stk::topology::ELEMENT_RANK, elementToMoveId);

    if ( myProc == fromProc )
    {
      stk::mesh::Bucket& bucket = stkMeshBulkData.bucket(elementToMove);
      bool isGhosted = !bucket.owned() && !bucket.shared();
      EXPECT_TRUE(isGhosted);

      stk::mesh::PairIterEntityComm commStuff = stkMeshBulkData.my_internal_entity_comm_map(stkMeshBulkData.entity_key(elementToMove));
      size_t numProcsToCommunicateWithForEntity = 1;
      EXPECT_TRUE(commStuff.size() == numProcsToCommunicateWithForEntity);
      EXPECT_TRUE((*commStuff.first).proc == toProc);

      size_t auraGhostingId = 1;
      EXPECT_TRUE((*commStuff.first).ghost_id == auraGhostingId);
    }
    else
    {
      stk::mesh::Bucket& bucket = stkMeshBulkData.bucket(elementToMove);
      EXPECT_TRUE(bucket.owned());

      stk::mesh::PairIterEntityComm commStuff = stkMeshBulkData.my_internal_entity_comm_map(stkMeshBulkData.entity_key(elementToMove));
      size_t numProcsToCommunicateWithForEntity = 1;
      EXPECT_TRUE(commStuff.size() == numProcsToCommunicateWithForEntity);
      EXPECT_TRUE((*commStuff.first).proc == fromProc);

      size_t auraGhostingId = 1;
      EXPECT_TRUE((*commStuff.first).ghost_id == auraGhostingId);
    }

    elementToMoveId = 9;
    fromProc = 1;
    toProc = 0;

    testElementMove(fromProc, toProc, myProc, elementToMoveId, stkMeshBulkData);
    elementToMove = stkMeshBulkData.get_entity(stk::topology::ELEMENT_RANK, elementToMoveId);

    if ( myProc == fromProc )
    {
      stk::mesh::Bucket& bucket = stkMeshBulkData.bucket(elementToMove);
      bool isGhosted = !bucket.owned() && !bucket.shared();
      EXPECT_TRUE(isGhosted);

      stk::mesh::PairIterEntityComm commStuff = stkMeshBulkData.my_internal_entity_comm_map(stkMeshBulkData.entity_key(elementToMove));
      size_t numProcsToCommunicateWithForEntity = 1;
      EXPECT_TRUE(commStuff.size() == numProcsToCommunicateWithForEntity);
      EXPECT_TRUE((*commStuff.first).proc == toProc);

      size_t auraGhostingId = 1;
      EXPECT_TRUE((*commStuff.first).ghost_id == auraGhostingId);
    }
    else
    {
      stk::mesh::Bucket& bucket = stkMeshBulkData.bucket(elementToMove);
      EXPECT_TRUE(bucket.owned());
      stk::mesh::PairIterEntityComm commStuff = stkMeshBulkData.my_internal_entity_comm_map(stkMeshBulkData.entity_key(elementToMove));
      size_t numProcsToCommunicateWithForEntity = 1;
      EXPECT_TRUE(commStuff.size() == numProcsToCommunicateWithForEntity);
      EXPECT_TRUE((*commStuff.first).proc == fromProc);

      size_t auraGhostingId = 1;
      EXPECT_TRUE((*commStuff.first).ghost_id == auraGhostingId);

    }
  }
}

TEST(UnderstandingDistributedIndex, GhostANode)
{
  stk::ParallelMachine communicator = MPI_COMM_WORLD;
  int procCount = stk::parallel_machine_size(communicator);
  int myProc = stk::parallel_machine_rank(communicator);

  if(procCount == 2)
  {
    const std::string generatedMeshSpec = "generated:2x2x4|sideset:xXyYzZ|nodeset:xXyYzZ";
    stk::unit_test_util::StkMeshCreator stkMesh(generatedMeshSpec, communicator);

    stk::unit_test_util::BulkDataTester &stkMeshBulkData = *stkMesh.getBulkData();

    testSharedNodesFor2x2x4MeshForTwoProcs(myProc, stkMeshBulkData);

    int otherProcId = 1;
    if ( myProc == 1 )
    {
      otherProcId = 0;
    }

    // elements 1-8 are proc 0
    // elements 9-16 are proc 1
    // 5-8 are ghosted on proc 1
    // 9-12 are ghosted on proc 0

    int node45 = 45;
    int procThatOwnsNode45 = 1;
    int procThatReceivesGhostedNode45 = 0;

    if (myProc == procThatOwnsNode45)
    {
      stk::mesh::Entity ghostedNode = stkMeshBulkData.get_entity(stk::topology::NODE_RANK, node45);
      stk::mesh::Bucket& bucket = stkMeshBulkData.bucket(ghostedNode);
      bool isGhosted = !bucket.shared() && !bucket.owned();
      EXPECT_FALSE(isGhosted);

      stk::mesh::PairIterEntityComm commStuff = stkMeshBulkData.my_internal_entity_comm_map(stkMeshBulkData.entity_key(ghostedNode));
      size_t numProcsToCommunicateWithForEntity = 0;
      EXPECT_TRUE(commStuff.size() == numProcsToCommunicateWithForEntity);
    }

    //BEGIN_DOC_FOR_GHOSTED_NODE
    stkMeshBulkData.modification_begin();
    stk::mesh::Ghosting &ghosting = stkMeshBulkData.create_ghosting("Ghost Node 45");
    std::vector< std::pair<stk::mesh::Entity, int> > ghostingStruct;
    if ( myProc == procThatOwnsNode45 )
    {
      stk::mesh::Entity nodeToGhost = stkMeshBulkData.get_entity(stk::topology::NODE_RANK, node45);
      ghostingStruct.push_back(std::make_pair(nodeToGhost, procThatReceivesGhostedNode45));
    }
    stkMeshBulkData.change_ghosting(ghosting, ghostingStruct);
    stkMeshBulkData.modification_end();

    if ( myProc == procThatReceivesGhostedNode45 )
    {
      int elementIdConnectedToGhostedNode45 = 16;
      stk::mesh::Entity ghostedElement = stkMeshBulkData.get_entity(stk::topology::ELEMENT_RANK, elementIdConnectedToGhostedNode45);
      EXPECT_FALSE(stkMeshBulkData.is_valid(ghostedElement));
    }
    //END_DOC_FOR_GHOSTED_NODE

    if ( myProc == procThatReceivesGhostedNode45 )
    {
      stk::mesh::Entity ghostedNode = stkMeshBulkData.get_entity(stk::topology::NODE_RANK, node45);
      stk::mesh::Bucket& bucket = stkMeshBulkData.bucket(ghostedNode);
      bool isGhosted = !bucket.shared() && !bucket.owned();
      EXPECT_TRUE(isGhosted);

      stk::mesh::PairIterEntityComm commStuff = stkMeshBulkData.my_internal_entity_comm_map(stkMeshBulkData.entity_key(ghostedNode));
      size_t numProcsToCommunicateWithForEntity = 1;
      EXPECT_TRUE(commStuff.size() == numProcsToCommunicateWithForEntity);
      EXPECT_TRUE((*commStuff.first).proc == otherProcId);

      size_t customGhosting = ghosting.ordinal();
      EXPECT_TRUE((*commStuff.first).ghost_id == customGhosting);
    }
    else
    {
      stk::mesh::Entity ghostedNode = stkMeshBulkData.get_entity(stk::topology::NODE_RANK, node45);
      stk::mesh::Bucket& bucket = stkMeshBulkData.bucket(ghostedNode);
      bool isThisAGhostOnThisProc = !bucket.shared() && !bucket.owned();
      EXPECT_FALSE(isThisAGhostOnThisProc);

      stk::mesh::PairIterEntityComm commStuff = stkMeshBulkData.my_internal_entity_comm_map(stkMeshBulkData.entity_key(ghostedNode));
      size_t numProcsToCommunicateWithForEntity = 1;
      EXPECT_TRUE(commStuff.size() == numProcsToCommunicateWithForEntity);
    }
  }
}

void testReceivingProcHasOneNodeGhosted(stk::unit_test_util::BulkDataTester &stkMeshBulkData,
                                        stk::mesh::Ghosting &ghosting1,
                                        stk::mesh::Ghosting &ghosting2,
                                        int nodeIdToGhost, int owningProc)
{
  stk::mesh::Entity ghostedNode = stkMeshBulkData.get_entity(stk::topology::NODE_RANK, nodeIdToGhost);
  stk::mesh::Bucket& bucket = stkMeshBulkData.bucket(ghostedNode);
  bool isGhosted = !bucket.shared() && !bucket.owned();
  EXPECT_TRUE(isGhosted);

  stk::mesh::PairIterEntityComm commStuff = stkMeshBulkData.my_internal_entity_comm_map(stkMeshBulkData.entity_key(ghostedNode));
  size_t numGhostingsInvolvingEntity = 2;
  EXPECT_EQ(numGhostingsInvolvingEntity, commStuff.size());
  EXPECT_EQ(owningProc, (*commStuff.first).proc);
  EXPECT_EQ(ghosting1.ordinal(), (*commStuff.first).ghost_id);
  commStuff++;
  EXPECT_EQ(owningProc, (*commStuff.first).proc);
  EXPECT_EQ(ghosting2.ordinal(), (*commStuff.first).ghost_id);
}

void testOwningProcHasOneNodeGhosted(stk::unit_test_util::BulkDataTester &stkMeshBulkData,
                                     stk::mesh::Ghosting &ghosting1,
                                     stk::mesh::Ghosting &ghosting2,
                                     int nodeIdToGhost, int ghostReceivingProc)
{
  stk::mesh::Entity ghostedNode = stkMeshBulkData.get_entity(stk::topology::NODE_RANK, nodeIdToGhost);
  stk::mesh::PairIterEntityComm commStuff = stkMeshBulkData.my_internal_entity_comm_map(stkMeshBulkData.entity_key(ghostedNode));
  size_t numGhostingsInvolvingEntity = 2;
  EXPECT_EQ(numGhostingsInvolvingEntity, commStuff.size());
  EXPECT_EQ(ghostReceivingProc, (*commStuff.first).proc);
  EXPECT_EQ(ghosting1.ordinal(), (*commStuff.first).ghost_id);
  commStuff++;
  EXPECT_EQ(ghostReceivingProc, (*commStuff.first).proc);
  EXPECT_EQ(ghosting2.ordinal(), (*commStuff.first).ghost_id);
}

void testReceivingProcHasNoGhosts(stk::unit_test_util::BulkDataTester &stkMeshBulkData,
                                  stk::mesh::Ghosting &ghosting1,
                                  stk::mesh::Ghosting &ghosting2,
                                  int nodeIdToGhost, int owningProc)
{
  stk::mesh::Entity ghostedNode = stkMeshBulkData.get_entity(stk::topology::NODE_RANK, nodeIdToGhost);
  EXPECT_FALSE(stkMeshBulkData.is_valid(ghostedNode));
}

void testOwningProcHasNoGhosts(stk::unit_test_util::BulkDataTester &stkMeshBulkData,
                               stk::mesh::Ghosting &ghosting1,
                               stk::mesh::Ghosting &ghosting2,
                               int nodeIdToGhost, int ghostReceivingProc)
{
  stk::mesh::Entity ghostedNode = stkMeshBulkData.get_entity(stk::topology::NODE_RANK, nodeIdToGhost);
  stk::mesh::PairIterEntityComm commStuff = stkMeshBulkData.my_internal_entity_comm_map(stkMeshBulkData.entity_key(ghostedNode));
  size_t numGhostingsInvolvingEntity = 0;
  EXPECT_EQ(numGhostingsInvolvingEntity, commStuff.size());
}

void testReceivingProcAfterOneGhostingDestroyed(stk::unit_test_util::BulkDataTester &stkMeshBulkData,
                                                stk::mesh::Ghosting &ghosting1,
                                                stk::mesh::Ghosting &ghosting2,
                                                int nodeIdToGhost, int owningProc)
{
  stk::mesh::Entity ghostedNode = stkMeshBulkData.get_entity(stk::topology::NODE_RANK, nodeIdToGhost);
  stk::mesh::Bucket& bucket = stkMeshBulkData.bucket(ghostedNode);
  bool isGhosted = !bucket.shared() && !bucket.owned();
  EXPECT_TRUE(isGhosted);

  stk::mesh::PairIterEntityComm commStuff = stkMeshBulkData.my_internal_entity_comm_map(stkMeshBulkData.entity_key(ghostedNode));
  size_t numGhostingsInvolvingEntity = 1;
  EXPECT_EQ(numGhostingsInvolvingEntity, commStuff.size());
  EXPECT_EQ(owningProc, (*commStuff.first).proc);
  EXPECT_EQ(ghosting2.ordinal(), (*commStuff.first).ghost_id);
}

void testOwningProcAfterOneGhostingDestroyed(stk::unit_test_util::BulkDataTester &stkMeshBulkData,
                                             stk::mesh::Ghosting &ghosting1,
                                             stk::mesh::Ghosting &ghosting2,
                                             int nodeIdToGhost, int ghostReceivingProc)
{
  stk::mesh::Entity ghostedNode = stkMeshBulkData.get_entity(stk::topology::NODE_RANK, nodeIdToGhost);
  stk::mesh::PairIterEntityComm commStuff = stkMeshBulkData.my_internal_entity_comm_map(stkMeshBulkData.entity_key(ghostedNode));
  size_t numGhostingsInvolvingEntity = 1;
  EXPECT_EQ(numGhostingsInvolvingEntity, commStuff.size());
  EXPECT_EQ(ghostReceivingProc, (*commStuff.first).proc);
  EXPECT_EQ(ghosting2.ordinal(), (*commStuff.first).ghost_id);
}

void testOneNodeGhostedTwice(stk::unit_test_util::BulkDataTester &stkMeshBulkData,
                             stk::mesh::Ghosting &ghosting1,
                             stk::mesh::Ghosting &ghosting2,
                             int nodeIdToGhost,
                             int owningProc,
                             int ghostReceivingProc,
                             int myProc)
{
  if(myProc == ghostReceivingProc)
  {
    testReceivingProcHasOneNodeGhosted(stkMeshBulkData, ghosting1, ghosting2, nodeIdToGhost, owningProc);
  }
  else
  {
    testOwningProcHasOneNodeGhosted(stkMeshBulkData, ghosting1, ghosting2, nodeIdToGhost, ghostReceivingProc);
  }
}

void testNoGhosts(stk::unit_test_util::BulkDataTester &stkMeshBulkData,
                  stk::mesh::Ghosting &ghosting1,
                  stk::mesh::Ghosting &ghosting2,
                  int nodeIdToGhost,
                  int owningProc,
                  int ghostReceivingProc,
                  int myProc)
{
  if(myProc == ghostReceivingProc)
  {
    testReceivingProcHasNoGhosts(stkMeshBulkData, ghosting1, ghosting2, nodeIdToGhost, owningProc);
  }
  else
  {
    testOwningProcHasNoGhosts(stkMeshBulkData, ghosting1, ghosting2, nodeIdToGhost, ghostReceivingProc);
  }
}

void testOneGhostingDestroyed(stk::unit_test_util::BulkDataTester &stkMeshBulkData,
                              stk::mesh::Ghosting &ghosting1,
                              stk::mesh::Ghosting &ghosting2,
                              int nodeIdToGhost,
                              int owningProc,
                              int ghostReceivingProc,
                              int myProc)
{
  if(myProc == ghostReceivingProc)
  {
    testReceivingProcAfterOneGhostingDestroyed(stkMeshBulkData, ghosting1, ghosting2, nodeIdToGhost, owningProc);
  }
  else
  {
    testOwningProcAfterOneGhostingDestroyed(stkMeshBulkData, ghosting1, ghosting2, nodeIdToGhost, ghostReceivingProc);
  }
}

TEST(UnderstandingDistributedIndex, MultipleCustomGhostings)
{
  stk::ParallelMachine communicator = MPI_COMM_WORLD;
  int procCount = stk::parallel_machine_size(communicator);
  int myProc = stk::parallel_machine_rank(communicator);

  if(procCount == 2)
  {
    const std::string generatedMeshSpec = "generated:2x2x4|sideset:xXyYzZ|nodeset:xXyYzZ";
    stk::unit_test_util::StkMeshCreator stkMesh(generatedMeshSpec, communicator);

    int owningProc = 0;
    int ghostReceivingProc = 1;
    int nodeIdToGhost = 2;

    stk::unit_test_util::BulkDataTester &stkMeshBulkData = *stkMesh.getBulkData();

    stkMeshBulkData.modification_begin();
    stk::mesh::Ghosting &ghosting1 = stkMeshBulkData.create_ghosting("Ghosting For Algorithm 1");
    stk::mesh::Ghosting &ghosting2 = stkMeshBulkData.create_ghosting("Ghosting For Algorithm 2");
    std::vector<std::pair<stk::mesh::Entity, int> > ghostsToSend;
    if(myProc == owningProc)
    {
      stk::mesh::Entity nodeToGhost = stkMeshBulkData.get_entity(stk::topology::NODE_RANK, nodeIdToGhost);
      ghostsToSend.push_back(std::make_pair(nodeToGhost, ghostReceivingProc));
    }
    stkMeshBulkData.change_ghosting(ghosting1, ghostsToSend);
    stkMeshBulkData.change_ghosting(ghosting2, ghostsToSend);
    stkMeshBulkData.modification_end();

    testOneNodeGhostedTwice(stkMeshBulkData, ghosting1, ghosting2, nodeIdToGhost, owningProc, ghostReceivingProc, myProc);

    stkMeshBulkData.modification_begin();
    std::vector<std::pair<stk::mesh::Entity, int> > noGhostsToSend;
    std::vector<stk::mesh::EntityKey> receivedGhostsToRemove;
    if(myProc == ghostReceivingProc)
    {
      receivedGhostsToRemove.push_back(stk::mesh::EntityKey(stk::topology::NODE_RANK, nodeIdToGhost));
    }
    stkMeshBulkData.change_ghosting(ghosting1, noGhostsToSend, receivedGhostsToRemove);
    stkMeshBulkData.modification_end();

    testOneGhostingDestroyed(stkMeshBulkData, ghosting1, ghosting2, nodeIdToGhost, owningProc, ghostReceivingProc, myProc);
  }
}

TEST(UnderstandingDistributedIndex, MultipleCustomGhostingsWithDestroy)
{
  stk::ParallelMachine communicator = MPI_COMM_WORLD;
  int procCount = stk::parallel_machine_size(communicator);
  int myProc = stk::parallel_machine_rank(communicator);

  if(procCount == 2)
  {
    const std::string generatedMeshSpec = "generated:2x2x4|sideset:xXyYzZ|nodeset:xXyYzZ";
    stk::unit_test_util::StkMeshCreator stkMesh(generatedMeshSpec, communicator);

    int owningProc = 0;
    int ghostReceivingProc = 1;
    int nodeIdToGhost = 2;

    stk::unit_test_util::BulkDataTester &stkMeshBulkData = *stkMesh.getBulkData();

    stkMeshBulkData.modification_begin();
    stk::mesh::Ghosting &ghosting1 = stkMeshBulkData.create_ghosting("Ghosting For Algorithm 1");
    stk::mesh::Ghosting &ghosting2 = stkMeshBulkData.create_ghosting("Ghosting For Algorithm 2");
    std::vector<std::pair<stk::mesh::Entity, int> > ghostsToSend;
    if(myProc == owningProc)
    {
      stk::mesh::Entity nodeToGhost = stkMeshBulkData.get_entity(stk::topology::NODE_RANK, nodeIdToGhost);
      ghostsToSend.push_back(std::make_pair(nodeToGhost, ghostReceivingProc));
    }

    stkMeshBulkData.change_ghosting(ghosting1, ghostsToSend);
    stkMeshBulkData.change_ghosting(ghosting2, ghostsToSend);
    testOneNodeGhostedTwice(stkMeshBulkData, ghosting1, ghosting2, nodeIdToGhost, owningProc, ghostReceivingProc, myProc);

    stkMeshBulkData.destroy_ghosting(ghosting1);
    stkMeshBulkData.destroy_ghosting(ghosting2);
    testNoGhosts(stkMeshBulkData, ghosting1, ghosting2, nodeIdToGhost, owningProc, ghostReceivingProc, myProc);

    stkMeshBulkData.change_ghosting(ghosting1, ghostsToSend);
    stkMeshBulkData.change_ghosting(ghosting2, ghostsToSend);
    testOneNodeGhostedTwice(stkMeshBulkData, ghosting1, ghosting2, nodeIdToGhost, owningProc, ghostReceivingProc, myProc);

    stkMeshBulkData.destroy_ghosting(ghosting1);
    stkMeshBulkData.modification_end();
    testOneGhostingDestroyed(stkMeshBulkData, ghosting1, ghosting2, nodeIdToGhost, owningProc, ghostReceivingProc, myProc);
  }
}

// Test multiple ghostings (delete a node from one ghosting, but not the other)
// Try node only storage in DI - then use create_edges and see what happens()
//          - add test for legacy capability that tests correct sharing of edges
//          - use that test to try new things?

}
