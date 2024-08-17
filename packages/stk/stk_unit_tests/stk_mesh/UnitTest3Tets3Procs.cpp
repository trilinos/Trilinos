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
#include <stk_util/parallel/Parallel.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>  // for declare_element
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/DestroyElements.hpp>
#include <stk_mesh/baseImpl/CommEntityMods.hpp>
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/MetaData.hpp"   // for MetaData
#include "stk_mesh/base/Types.hpp"      // for PartVector, EntityId, etc
#include "stk_mesh/base/Comm.hpp"
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_io/FillMesh.hpp"
#include "stk_mesh/base/MeshBuilder.hpp"
#include "stk_unit_test_utils/BuildMesh.hpp"
#include "UnitTestTextMeshFixture.hpp"

namespace {

TEST(MakeUnique, entityParallelStates_empty)
{
  std::vector<stk::mesh::EntityParallelState> empty;
  stk::mesh::impl::make_unique(empty);
  EXPECT_TRUE(empty.empty());
}

stk::mesh::EntityKey get_node_key(stk::mesh::EntityId nodeId)
{
  return stk::mesh::EntityKey(stk::topology::NODE_RANK, nodeId);
}

stk::mesh::EntityParallelState get_pll_state(stk::mesh::EntityId id,
                                             int proc,
                                             stk::mesh::EntityState state)
{
  return {proc, state, {get_node_key(id), stk::mesh::Entity(), -1}, false};
}

TEST(MakeUnique, entityParallelStates_alreadyUnique)
{
  std::vector<stk::mesh::EntityParallelState> pllStates;
  pllStates.push_back(get_pll_state(1, 0, stk::mesh::Modified));
  pllStates.push_back(get_pll_state(2, 1, stk::mesh::Modified));
  pllStates.push_back(get_pll_state(3, 2, stk::mesh::Modified));

  stk::mesh::impl::make_unique(pllStates);
  EXPECT_EQ(3u, pllStates.size());
}

TEST(MakeUnique, entityParallelStates_duplicate)
{
  std::vector<stk::mesh::EntityParallelState> pllStates;
  pllStates.push_back(get_pll_state(1, 0, stk::mesh::Modified));
  pllStates.push_back(get_pll_state(2, 1, stk::mesh::Modified));
  pllStates.push_back(get_pll_state(2, 1, stk::mesh::Deleted));
  pllStates.push_back(get_pll_state(3, 2, stk::mesh::Modified));

  stk::mesh::impl::make_unique(pllStates);
  EXPECT_EQ(3u, pllStates.size());
  EXPECT_EQ(stk::mesh::Deleted, pllStates[1].state);
}

TEST(MakeUnique, entityParallelStates_duplicate_sort_then_unique)
{
  std::vector<stk::mesh::EntityParallelState> pllStates;
  pllStates.push_back(get_pll_state(2, 1, stk::mesh::Deleted));
  pllStates.push_back(get_pll_state(1, 0, stk::mesh::Modified));
  pllStates.push_back(get_pll_state(2, 1, stk::mesh::Modified));
  pllStates.push_back(get_pll_state(3, 2, stk::mesh::Modified));
  pllStates.push_back(get_pll_state(1, 0, stk::mesh::Deleted));

  std::sort(pllStates.begin(), pllStates.end());
  stk::mesh::impl::make_unique(pllStates);
  EXPECT_EQ(3u, pllStates.size());
  EXPECT_EQ(stk::mesh::Deleted, pllStates[0].state);
  EXPECT_EQ(stk::mesh::Deleted, pllStates[1].state);
}

TEST(MakeUnique, entityParallelStates_3duplicates_sort_then_unique)
{
  std::vector<stk::mesh::EntityParallelState> pllStates;
  pllStates.push_back(get_pll_state(2, 1, stk::mesh::Deleted));
  pllStates.push_back(get_pll_state(1, 0, stk::mesh::Modified));
  pllStates.push_back(get_pll_state(2, 1, stk::mesh::Modified));
  pllStates.push_back(get_pll_state(3, 2, stk::mesh::Modified));
  pllStates.push_back(get_pll_state(2, 1, stk::mesh::Created));
  pllStates.push_back(get_pll_state(1, 0, stk::mesh::Deleted));

  std::sort(pllStates.begin(), pllStates.end());
  stk::mesh::impl::make_unique(pllStates);
  EXPECT_EQ(3u, pllStates.size());
  EXPECT_EQ(stk::mesh::Deleted, pllStates[0].state);
  EXPECT_EQ(stk::mesh::Deleted, pllStates[1].state);
}

TEST(ThreeTet4sOn3Procs, deleteMiddleTetOnP1AndRecreateOnP2_works)
{
    if (stk::parallel_machine_size(MPI_COMM_WORLD) != 3) { GTEST_SKIP(); }

    std::string meshDesc = "0,1,TET_4,5,7,2,1\n"
                           "1,2,TET_4,6,7,1,4\n"
                           "2,3,TET_4,6,8,4,3\n";

    std::vector<double> coordinates {
        0.082356, -0.411779, 0.082356,
        0.082356, -0.247067, 0.082356,
        0.082356, -0.576490, 0.247067,
       0.082356, -0.411779, 0.247067,
        0.164711, -0.329423, 0.000000,
        0.164711, -0.494134, 0.164711,
        0.164711, -0.329423, 0.164711,
        0.164711, -0.494134, 0.329423,
    };

    stk::mesh::MeshBuilder builder(MPI_COMM_WORLD);
    builder.set_spatial_dimension(3);
    std::unique_ptr<stk::mesh::BulkData> bulkPtr = builder.create();
    stk::mesh::BulkData &bulk = *bulkPtr;

    stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

    const int proc = bulk.parallel_rank();
    stk::mesh::Part *partWithTopology = bulk.mesh_meta_data().get_part("block_TETRAHEDRON_4");

    const stk::mesh::Entity node1 = bulk.get_entity(stk::topology::NODE_RANK, 1);
    const stk::mesh::Entity node7 = bulk.get_entity(stk::topology::NODE_RANK, 7);

    bulk.modification_begin();
    stk::mesh::EntityVector elem2Vec;
    if(proc == 0)
    {
        bulk.add_node_sharing(node1, 2);
        bulk.add_node_sharing(node7, 2);
    }
    else if(proc == 1)
    {
        elem2Vec = {bulk.get_entity(stk::topology::ELEM_RANK, 2)};
    }
    else if(proc == 2)
    {
        stk::mesh::declare_element(bulk, *partWithTopology, 11, {6,7,1,4});
        bulk.add_node_sharing(node1, 0);
        bulk.add_node_sharing(node7, 0);
    }

    stk::mesh::destroy_elements_no_mod_cycle(bulk, elem2Vec, bulk.mesh_meta_data().universal_part());

    EXPECT_NO_THROW(bulk.modification_end());
}

TEST(ThreeTet10sOn3Procs, deleteMiddleTetOnP1AndRecreateOnP2_works)
{
    if (stk::parallel_machine_size(MPI_COMM_WORLD) != 3) { GTEST_SKIP(); }

    std::string meshDesc = "0,1,TET_10,5,7,2,1,9,10,11,12,13,14\n"
                           "1,2,TET_10,6,7,1,4,15,13,16,17,18,19\n"
                           "2,3,TET_10,6,8,4,3,20,21,17,22,23,24\n";

    std::vector<double> coordinates {
        0.082356, -0.411779, 0.082356,
        0.082356, -0.247067, 0.082356,
        0.082356, -0.576490, 0.247067,
        0.082356, -0.411779, 0.247067,
        0.164711, -0.329423, 0.000000,
        0.164711, -0.494134, 0.164711,
        0.164711, -0.329423, 0.164711,
        0.164711, -0.494134, 0.329423,
        0.164711, -0.329423, 0.082356,
        0.123534, -0.288245, 0.123534,
        0.123534, -0.288245, 0.041178,
        0.123534, -0.370601, 0.041178,
        0.123534, -0.370601, 0.123534,
        0.082356, -0.329423, 0.082356,
        0.164711, -0.411779, 0.164711,
        0.123534, -0.452956, 0.123534,
        0.123534, -0.452956, 0.205889,
        0.123534, -0.370601, 0.205889,
        0.082356, -0.411779, 0.164711,
        0.164711, -0.494134, 0.247067,
        0.123534, -0.452956, 0.288245,
        0.123534, -0.535312, 0.205889,
        0.123534, -0.535312, 0.288245,
        0.082356, -0.494134, 0.247067,
    };

    stk::mesh::MeshBuilder builder(MPI_COMM_WORLD);
    builder.set_spatial_dimension(3);
    std::unique_ptr<stk::mesh::BulkData> bulkPtr = builder.create();
    stk::mesh::BulkData &bulk = *bulkPtr;

    stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

    const int proc = bulk.parallel_rank();
    stk::mesh::Part *partWithTopology = bulk.mesh_meta_data().get_part("block_TETRAHEDRON_10");

    const stk::mesh::Entity node1 = bulk.get_entity(stk::topology::NODE_RANK, 1);
    const stk::mesh::Entity node7 = bulk.get_entity(stk::topology::NODE_RANK, 7);
    stk::mesh::Entity node13 = bulk.get_entity(stk::topology::NODE_RANK, 13);

    bulk.modification_begin();
    stk::mesh::EntityVector elem2Vec;
    if(proc == 0)
    {
        bulk.add_node_sharing(node1, 2);
        bulk.add_node_sharing(node7, 2);
        bulk.add_node_sharing(node13, 2);
    }
    else if(proc == 1)
    {
        elem2Vec = {bulk.get_entity(stk::topology::ELEM_RANK, 2)};
    }
    else if(proc == 2)
    {
        stk::mesh::declare_element(bulk, *partWithTopology, 11, {6,7,1,4,15,13,16,17,18,19});
        bulk.add_node_sharing(node1, 0);
        bulk.add_node_sharing(node7, 0);
        bulk.add_node_sharing(node13, 0);
    }

    stk::mesh::destroy_elements_no_mod_cycle(bulk, elem2Vec, bulk.mesh_meta_data().universal_part());

    EXPECT_NO_THROW(bulk.modification_end());

    node13 = bulk.get_entity(stk::topology::NODE_RANK, 13);
    if (bulk.parallel_rank() != 1) {
      EXPECT_TRUE(bulk.is_valid(node13));
      EXPECT_TRUE(bulk.bucket(node13).shared());
    }
    else {
      EXPECT_FALSE(bulk.is_valid(node13));
    }
}

} // empty namespace

