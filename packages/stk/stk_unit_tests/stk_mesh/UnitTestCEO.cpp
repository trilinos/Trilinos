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

#include <gtest/gtest.h>                // for AssertHelper, EXPECT_TRUE, etc
#include <iostream>                     // for basic_ostream::operator<<
#include <stk_mesh/base/BulkData.hpp>   // for BulkData, etc
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/MeshUtils.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_io/FillMesh.hpp>
#include <stk_util/parallel/Parallel.hpp>  // for parallel_machine_rank, etc
#include <string>                       // for string
#include <stk_unit_test_utils/BulkDataTester.hpp>
#include <utility>                      // for pair, make_pair
#include <vector>                       // for vector

#include "UnitTestCEOCommonUtils.hpp"
#include "UnitTestCEO2Elem.hpp"
#include "UnitTestCEO3Elem.hpp"
#include "UnitTestCEO4ElemEdge.hpp"
#include "UnitTestCEO4ElemRotate.hpp"
#include "UnitTestCEO8Elem.hpp"
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/EntityKey.hpp"  // for EntityKey
#include "stk_mesh/base/MetaData.hpp"   // for MetaData
#include "stk_mesh/base/Types.hpp"      // for EntityProcVec, EntityProc, etc
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_unit_test_utils/GenerateALefRAMesh.hpp"
namespace stk { namespace mesh { class Ghosting; } }
namespace stk { namespace mesh { class Part; } }
namespace stk { namespace mesh { class Selector; } }
namespace stk { namespace mesh { namespace fixtures { class BoxFixture; } } }
namespace stk { namespace mesh { namespace fixtures { class RingFixture; } } }

namespace stk
{
namespace mesh
{
class FieldBase;
}
}

using stk::mesh::Part;
using stk::mesh::MetaData;
using stk::mesh::BulkData;
using stk::mesh::Selector;
using stk::mesh::PartVector;
using stk::mesh::PairIterRelation;
using stk::mesh::EntityProc;
using stk::mesh::Entity;
using stk::mesh::EntityId;
using stk::mesh::EntityKey;
using stk::mesh::EntityVector;
using stk::mesh::EntityRank;
using stk::mesh::fixtures::RingFixture;
using stk::mesh::fixtures::BoxFixture;

namespace
{
//==============================================================================

TEST(CEO, change_entity_owner_2Elem2ProcMove)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int p_rank = stk::parallel_machine_rank(pm);
  const int p_size = stk::parallel_machine_size(pm);

  if(p_size != 2)
  {
    return;
  }

  const int spatial_dimension = 2;
  stk::mesh::MetaData meta(spatial_dimension);
  stk::unit_test_util::BulkDataTester bulk(meta, pm);

  stk::mesh::EntityVector elems;
  CEOUtils::fillMeshfor2Elem2ProcMoveAndTest(bulk, meta, elems);

  stk::mesh::EntityProcVec entity_procs;
  if(p_rank == 0)
  {
    entity_procs.push_back(stk::mesh::EntityProc(elems[1], 1));
    entity_procs.push_back(stk::mesh::EntityProc(bulk.get_entity(stk::topology::NODE_RANK, 4), 1));
    entity_procs.push_back(stk::mesh::EntityProc(bulk.get_entity(stk::topology::NODE_RANK, 5), 1));
    entity_procs.push_back(stk::mesh::EntityProc(bulk.get_entity(stk::topology::NODE_RANK, 6), 1));
  }

  bulk.modification_begin("testing CEO");
  bulk.my_internal_change_entity_owner(entity_procs);

  CEOUtils::checkStatesAfterCEO_2Elem2ProcMove(bulk);
}

void check_sideset_orientation(const stk::mesh::BulkData& bulk,
                               const std::vector<stk::mesh::SideSet*> & sidesets,
                               const stk::mesh::EntityId expectedId,
                               const stk::mesh::ConnectivityOrdinal expectedOrdinal)
{
  stk::mesh::SideSet sideSet = *sidesets[0];
  stk::mesh::SideSetEntry sideSetEntry;
  if (sideSet.size() > 0) {
    sideSetEntry = sideSet[0];
  }
  EXPECT_EQ(expectedId, bulk.identifier(sideSetEntry.element));
  EXPECT_EQ(expectedOrdinal, sideSetEntry.side);
}

TEST(CEO, change_entity_owner_2ElemWithSideset) {
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int pRank = stk::parallel_machine_rank(pm);
  const int pSize = stk::parallel_machine_size(pm);

  if(pSize != 2) return;

  const std::string outputMeshName = "2ElemWithSideset.e";
  const stk::mesh::EntityId expectedId = 2;
  const stk::mesh::ConnectivityOrdinal expectedOrdinal = 5;

  stk::mesh::MeshBuilder builder(pm);
  builder.set_spatial_dimension(3);
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = builder.create();
  stk::mesh::BulkData& bulk = *bulkPtr;

  stk::unit_test_util::create_AB_mesh_with_sideset_and_field(bulk, stk::unit_test_util::LEFT, stk::unit_test_util::DECREASING, "dummyField");

  if (pRank == 0)
  {
    std::vector<stk::mesh::SideSet *> sidesets = bulk.get_sidesets();
    ASSERT_EQ(1u, sidesets.size());
    EXPECT_EQ(1u, sidesets[0]->size());

    check_sideset_orientation(bulk, sidesets, expectedId, expectedOrdinal);
  }

  stk::mesh::EntityProcVec entityProc;
  if (pRank == 0)
  {
    entityProc.push_back(stk::mesh::EntityProc(bulk.get_entity(stk::topology::ELEMENT_RANK, expectedId), 1));
  }
  bulk.change_entity_owner(entityProc);

  if (pRank == 0)
  {
    std::vector<stk::mesh::SideSet *> sidesets = bulk.get_sidesets();
    EXPECT_EQ(1u, sidesets.size());
    EXPECT_EQ(0u, sidesets[0]->size());
  }
  if (pRank == 1)
  {
    std::vector<stk::mesh::SideSet *> sidesets = bulk.get_sidesets();
    check_sideset_orientation(bulk, sidesets, expectedId, expectedOrdinal);
  }
}

void test_change_entity_owner_3Elem3Proc_WithCustomGhosts(stk::mesh::BulkData::AutomaticAuraOption autoAuraOption)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  int psize = stk::parallel_machine_size(communicator);
  int prank = stk::parallel_machine_rank(communicator);
  if(psize == 3)
  { // Skip unless we're on 3 processors

    const int spatialDim = 3;
    std::vector<std::string> rankNames;
    rankNames.push_back("node");
    rankNames.push_back("edge");
    rankNames.push_back("face");
    rankNames.push_back("elem");
    rankNames.push_back("const");

    stk::mesh::MeshBuilder builder(communicator);
    builder.set_spatial_dimension(spatialDim);
    builder.set_entity_rank_names(rankNames);
    std::shared_ptr<stk::mesh::BulkData> bulkPtr = builder.create();
    stk::mesh::BulkData& stkMeshBulkData = *bulkPtr;

    const std::string generatedMeshSpecification = "generated:1x1x6";

    stk::io::fill_mesh(generatedMeshSpecification, stkMeshBulkData);

    /////////////////////////
    stkMeshBulkData.modification_begin();
    std::vector< std::pair<stk::mesh::Entity, int> > ghostingStruct;
    stk::mesh::Ghosting &ghosting = stkMeshBulkData.create_ghosting("Ghost Node 1");
    if(prank == 0)
    {
      int proc2 = 2;
      stk::mesh::Entity node1 = stkMeshBulkData.get_entity(stk::topology::NODE_RANK, 1);
      ghostingStruct.push_back(std::make_pair(node1, proc2));
    }
    stkMeshBulkData.change_ghosting(ghosting, ghostingStruct);
    stkMeshBulkData.modification_end();

    if (prank==2)
    {
      stk::mesh::Entity node1 = stkMeshBulkData.get_entity(stk::topology::NODE_RANK, 1);
      ASSERT_TRUE(stkMeshBulkData.is_valid(node1));
      EXPECT_EQ(0, stkMeshBulkData.parallel_owner_rank(node1));
    }
    /////////////////////////

    stkMeshBulkData.modification_begin();

    if (prank==2)
    {
      stk::mesh::Entity node1 = stkMeshBulkData.get_entity(stk::topology::NODE_RANK, 1);
      stk::mesh::Entity constraint = stkMeshBulkData.declare_constraint(1);
      stkMeshBulkData.declare_relation(constraint,node1,0);
    }
    stk::mesh::fixup_ghosted_to_shared_nodes(stkMeshBulkData);
    stkMeshBulkData.modification_end();

    if(prank==2)
    {
      stk::mesh::Entity node1 = stkMeshBulkData.get_entity(stk::topology::NODE_RANK, 1);
      EXPECT_TRUE(stkMeshBulkData.bucket(node1).shared());
      EXPECT_TRUE(stkMeshBulkData.in_receive_ghost(ghosting, node1));
    }

    if(prank==0 || prank==2)
    {
      stk::mesh::Entity node1 = stkMeshBulkData.get_entity(stk::topology::NODE_RANK, 1);
      stk::mesh::EntityKey key = stkMeshBulkData.entity_key(node1);
      int otherProc = 2-prank;
      EXPECT_TRUE(stkMeshBulkData.in_ghost(ghosting, key, otherProc));
    }

    /////////////////////////

    stk::mesh::EntityProcVec entity_procs;

    if(prank == 0)
    {
      int proc2 = 2;
      stk::mesh::Entity elem1 = stkMeshBulkData.get_entity(stk::topology::ELEM_RANK, 1);
      entity_procs.push_back(stk::mesh::EntityProc(elem1, proc2));
      stk::mesh::Entity node1 = stkMeshBulkData.get_entity(stk::topology::NODE_RANK, 1);
      entity_procs.push_back(stk::mesh::EntityProc(node1, proc2));
    }

    stkMeshBulkData.change_entity_owner(entity_procs);

    if (prank == 2)
    {
      stk::mesh::Entity elem1 = stkMeshBulkData.get_entity(stk::topology::ELEM_RANK, 1);
      EXPECT_TRUE(stkMeshBulkData.parallel_owner_rank(elem1) == 2);
      stk::mesh::Entity node1 = stkMeshBulkData.get_entity(stk::topology::NODE_RANK, 1);
      ASSERT_TRUE(stkMeshBulkData.is_valid(node1));
      EXPECT_FALSE(stkMeshBulkData.in_receive_ghost(ghosting, node1));
      stk::mesh::Part& ghostingPart = stkMeshBulkData.ghosting_part(ghosting);
      EXPECT_FALSE(stkMeshBulkData.bucket(node1).member(ghostingPart));
    }
  }
}

TEST(CEO, change_entity_owner_3Elem3Proc_WithCustomGhosts_WithAura)
{
  test_change_entity_owner_3Elem3Proc_WithCustomGhosts(stk::mesh::BulkData::AUTO_AURA);
}

TEST(CEO, change_entity_owner_3Elem3Proc_WithCustomGhosts_WithoutAura)
{
  test_change_entity_owner_3Elem3Proc_WithCustomGhosts(stk::mesh::BulkData::NO_AUTO_AURA);
}

TEST(CEO,moveElem_fieldDataOfNodes)
{
  stk::ParallelMachine comm{MPI_COMM_WORLD};
  if(stk::parallel_machine_size(comm) == 2)
  {
    stk::mesh::MeshBuilder builder(comm);
    std::shared_ptr<stk::mesh::BulkData> bulkPtr = builder.create();
    stk::mesh::BulkData& bulk = *bulkPtr;
    stk::mesh::MetaData& meta = bulk.mesh_meta_data();
    auto &field1 = meta.declare_field<int>(stk::topology::NODE_RANK, "field1");
    stk::mesh::put_field_on_entire_mesh(field1);
    stk::io::fill_mesh("generated:1x1x2", bulk);

    stk::mesh::EntityVector sharedNodes;
    stk::mesh::get_selected_entities(meta.globally_shared_part(), bulk.buckets(stk::topology::NODE_RANK), sharedNodes);
    for(stk::mesh::Entity node : sharedNodes)
    {
      int *data = stk::mesh::field_data(field1, node);
      *data = stk::parallel_machine_rank(comm);
    }

    stk::mesh::EntityProcVec thingsToChange;
    if(stk::parallel_machine_rank(comm) == 1)
    {
      stk::mesh::EntityVector ownedElems;
      stk::mesh::get_selected_entities(meta.locally_owned_part(), bulk.buckets(stk::topology::ELEM_RANK), ownedElems);

      ASSERT_EQ(1u, ownedElems.size());
      stk::mesh::Entity elem{ownedElems[0]};
      thingsToChange.push_back(stk::mesh::EntityProc(elem, 0));
    }

    if(stk::parallel_machine_rank(comm) == 0)
    {
      for(stk::mesh::Entity node : sharedNodes)
      {
        int *data = stk::mesh::field_data(field1, node);
        EXPECT_EQ(0, *data);
      }
    }

    bulk.change_entity_owner(thingsToChange);

    if(stk::parallel_machine_rank(comm) == 0)
    {
      for(stk::mesh::Entity node : sharedNodes)
      {
        int *data = stk::mesh::field_data(field1, node);
        EXPECT_EQ(0, *data);
      }
    }
  }
}

}
