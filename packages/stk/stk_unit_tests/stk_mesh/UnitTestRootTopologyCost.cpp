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
#include <algorithm>

#include <stk_mesh/base/FEMHelpers.hpp>  // for declare_element
#include <stk_unit_test_utils/BuildMesh.hpp>

using stk::unit_test_util::build_mesh;

void make_small_hybrid_mesh(stk::mesh::MetaData &meta, stk::mesh::BulkData &mesh,
                            bool user_attempt_no_induce = false, bool user_parts_force_no_induce = true)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  int p_size = stk::parallel_machine_size(pm);

  if(p_size > 2)
  {
    return;
  }

  const unsigned p_rank = mesh.parallel_rank();

  stk::mesh::Part * hexPart = &meta.get_topology_root_part(stk::topology::HEX_8);
  stk::mesh::Part * pyrPart = &meta.get_topology_root_part(stk::topology::PYRAMID_5);
  stk::mesh::Part * tetPart = &meta.get_topology_root_part(stk::topology::TET_4);

  if (user_attempt_no_induce)
  {
    hexPart = &meta.declare_part_with_topology("my_hex_part",stk::topology::HEX_8, user_parts_force_no_induce);
    pyrPart = &meta.declare_part_with_topology("my_pyr_part",stk::topology::PYRAMID_5, user_parts_force_no_induce);
    tetPart = &meta.declare_part_with_topology("my_tet_part",stk::topology::TET_4, user_parts_force_no_induce);

    EXPECT_EQ(user_parts_force_no_induce, hexPart->force_no_induce());
    EXPECT_EQ(user_parts_force_no_induce, pyrPart->force_no_induce());
    EXPECT_EQ(user_parts_force_no_induce, tetPart->force_no_induce());
  }

  meta.commit();

  const size_t numHex = 1;
  stk::mesh::EntityIdVector hexNodeIDs[] {
    { 1, 2, 3, 4, 5, 6, 7, 8 }
  };
  stk::mesh::EntityId hexElemIDs[] = { 1 };

  const size_t numPyr = 1;
  stk::mesh::EntityIdVector pyrNodeIDs[] {
    { 5, 6, 7, 8, 9 }
  };
  stk::mesh::EntityId pyrElemIDs[] = { 2 };

  const size_t numTet = 4;
  stk::mesh::EntityIdVector tetNodeIDs[] {
    { 7, 8, 9, 12 },
    { 6, 9, 10, 7 },
    { 7, 9, 10, 12 },
    { 7, 12, 10, 11 }
  };
  stk::mesh::EntityId tetElemIDs[] = { 3, 4, 5, 6 };

  // list of triplets: (owner-proc, shared-nodeID, sharing-proc)
  std::vector< std::vector<unsigned> > shared_nodeIDs_and_procs
  {
    { 0, 5, 1 },  // proc 0
    { 0, 6, 1 },
    { 0, 7, 1 },
    { 0, 8, 1 },
    { 1, 5, 0 },  // proc 1
    { 1, 6, 0 },
    { 1, 7, 0 },
    { 1, 8, 0 }
  };

  mesh.modification_begin();

  if (0 == p_rank) {
    for (size_t i = 0; i < numHex; ++i) {
      stk::mesh::declare_element(mesh, *hexPart, hexElemIDs[i], hexNodeIDs[i]);
    }
  }
  if ( (1 == p_rank) || (1 == p_size) )  { // setup the pyramids/tets for either np 2 or serial
    for (size_t i = 0; i < numPyr; ++i) {
      stk::mesh::declare_element(mesh, *pyrPart, pyrElemIDs[i], pyrNodeIDs[i]);
    }
    for (size_t i = 0; i < numTet; ++i) {
      stk::mesh::declare_element(mesh, *tetPart, tetElemIDs[i], tetNodeIDs[i]);
    }
  }

  if (p_size > 1)
  {
    for (size_t nodeIdx = 0, end = shared_nodeIDs_and_procs.size(); nodeIdx < end; ++nodeIdx) {
      if (p_rank == shared_nodeIDs_and_procs[nodeIdx][0]) {
        stk::mesh::EntityId nodeID = shared_nodeIDs_and_procs[nodeIdx][1];
        int sharingProc = shared_nodeIDs_and_procs[nodeIdx][2];
        stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, nodeID);
        mesh.add_node_sharing(node, sharingProc);
      }
    }
  }

  mesh.modification_end();
}


TEST(UnitTestRootTopology, autoInduceFragmentation)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  int p_size = stk::parallel_machine_size(pm);

  if (p_size != 1)
  {
    return;
  }

  const unsigned spatialDim = 3;
  {
    // Root topology parts do not induce.  No user parts.  No bucket fragmentation.
    std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, pm);
    stk::mesh::MetaData& meta_no_induce = bulkPtr->mesh_meta_data();
    stk::mesh::BulkData& mesh = *bulkPtr;
    stk::mesh::PartVector parts = meta_no_induce.get_parts();
    for (auto part : parts)
    {
      meta_no_induce.force_no_induce(*part);
    }
    make_small_hybrid_mesh(meta_no_induce, mesh);

    const stk::mesh::BucketVector &node_buckets = mesh.buckets(stk::topology::NODE_RANK);
    unsigned num_buckets = node_buckets.size();
    EXPECT_EQ(num_buckets, 1u);
  }

  {
    // Root topology parts do induce.  No user parts.  Bucket fragmentation.
    std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, pm);
    stk::mesh::MetaData& meta_induce = bulkPtr->mesh_meta_data();
    stk::mesh::BulkData& mesh = *bulkPtr;
    make_small_hybrid_mesh(meta_induce, mesh);

    const stk::mesh::BucketVector &node_buckets = mesh.buckets(stk::topology::NODE_RANK);
    unsigned num_buckets = node_buckets.size();
    EXPECT_EQ(num_buckets, 5u);
  }

  {
    // Root topology parts do not induce.  User parts don't induce.  No bucket fragmentation.
    std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, pm);
    stk::mesh::MetaData& meta_no_induce = bulkPtr->mesh_meta_data();
    stk::mesh::BulkData& mesh = *bulkPtr;
    stk::mesh::PartVector parts = meta_no_induce.get_parts();
    for (auto part : parts)
    {
      meta_no_induce.force_no_induce(*part);
    }
    make_small_hybrid_mesh(meta_no_induce, mesh, true);

    const stk::mesh::BucketVector &node_buckets = mesh.buckets(stk::topology::NODE_RANK);
    unsigned num_buckets = node_buckets.size();
    EXPECT_EQ(num_buckets, 1u);
  }

  {
    // Root topology parts do induce.  User parts induce.  Bucket fragmentation.
    std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, pm);
    stk::mesh::MetaData& meta_induce = bulkPtr->mesh_meta_data();
    stk::mesh::BulkData& mesh = *bulkPtr;
    make_small_hybrid_mesh(meta_induce, mesh, true);

    const stk::mesh::BucketVector &node_buckets = mesh.buckets(stk::topology::NODE_RANK);
    unsigned num_buckets = node_buckets.size();
    EXPECT_EQ(num_buckets, 5u);
  }

  {
    // Root topology parts do not induce.  User parts induce.  No bucket fragmentation.
    std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, pm);
    stk::mesh::MetaData& meta_no_induce = bulkPtr->mesh_meta_data();
    stk::mesh::BulkData& mesh = *bulkPtr;
    stk::mesh::PartVector parts = meta_no_induce.get_parts();
    for (auto part : parts)
    {
      meta_no_induce.force_no_induce(*part);
    }
    make_small_hybrid_mesh(meta_no_induce, mesh, true, false);

    const stk::mesh::BucketVector &node_buckets = mesh.buckets(stk::topology::NODE_RANK);
    unsigned num_buckets = node_buckets.size();
    EXPECT_EQ(num_buckets, 5u);
  }

  {
    // Root topology parts do induce.  User parts induce.  Bucket fragmentation.
    std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, pm);
    stk::mesh::MetaData& meta_induce = bulkPtr->mesh_meta_data();
    stk::mesh::BulkData& mesh = *bulkPtr;
    make_small_hybrid_mesh(meta_induce, mesh, true, false);

    const stk::mesh::BucketVector &node_buckets = mesh.buckets(stk::topology::NODE_RANK);
    unsigned num_buckets = node_buckets.size();
    EXPECT_EQ(num_buckets, 5u);
  }
}

void check_parts_allowed(const stk::mesh::BucketVector &buckets, const stk::mesh::PartVector &parts_allowed)
{
  for (unsigned i = 0; i < buckets.size(); ++i)
  {
    const stk::mesh::Bucket &bucket = *buckets[i];
    const stk::mesh::PartVector &parts_found = bucket.supersets();

    for (unsigned j = 0; j < parts_found.size(); ++j)
    {
      const stk::mesh::Part *part = parts_found[j];
      bool found = (std::find(parts_allowed.begin(), parts_allowed.end(), part) != parts_allowed.end());
      EXPECT_TRUE(found);
    }
  }
}


TEST(UnitTestRootTopology, unusedPartsNotInBuckets)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  int p_size = stk::parallel_machine_size(pm);

  if (p_size != 1)
  {
    return;
  }

  const unsigned spatialDim = 3;

  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, pm);
  stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();
  stk::mesh::BulkData& mesh = *bulkPtr;
  make_small_hybrid_mesh(meta, mesh);

  stk::mesh::PartVector parts_allowed;
  parts_allowed.push_back(&meta.universal_part());
  parts_allowed.push_back(&meta.locally_owned_part());
  parts_allowed.push_back(&meta.get_topology_root_part(stk::topology::NODE));
  parts_allowed.push_back(&meta.get_topology_root_part(stk::topology::HEXAHEDRON_8));
  parts_allowed.push_back(&meta.get_topology_root_part(stk::topology::TETRAHEDRON_4));
  parts_allowed.push_back(&meta.get_topology_root_part(stk::topology::PYRAMID_5));

  size_t num_entities = 0;
  for (stk::mesh::EntityRank e_rank = stk::topology::NODE_RANK; e_rank <= stk::topology::ELEM_RANK; ++e_rank)
  {
    const stk::mesh::BucketVector &buckets = mesh.buckets(e_rank);
    check_parts_allowed(buckets, parts_allowed);

    for (unsigned bkt_j = 0; bkt_j < buckets.size(); ++bkt_j)
    {
      num_entities += buckets[bkt_j]->size();
    }
  }
  EXPECT_EQ(num_entities, 18u);
}

