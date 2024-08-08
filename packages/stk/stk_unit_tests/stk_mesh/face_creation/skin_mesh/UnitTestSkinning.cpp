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
#include <map>                          // for map, operator==, etc
#include <set>                          // for set, operator==
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/Entity.hpp>     // for Entity
#include <stk_mesh/base/GetEntities.hpp>  // for count_entities, etc
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/Selector.hpp>   // for Selector
#include <stk_mesh/base/SkinMesh.hpp>   // for skin_mesh
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine, etc
#include <vector>                       // for vector, etc
#include "mpi.h"                        // for ompi_communicator_t, etc

#include "stk_mesh/base/Types.hpp"      // for PartVector, BucketVector, etc
#include "stk_topology/topology.hpp"    // for topology, topology::rank_t, etc
#include "stk_unit_test_utils/stk_mesh_fixtures/GridFixture.hpp"  // for GridFixture
namespace stk { namespace mesh { class Part; } }
namespace stk { namespace mesh { class Bucket; } }

static const stk::topology::rank_t NODE_RANK = stk::topology::NODE_RANK;

using stk::mesh::MetaData;

class UnitTestStkMeshSkinning {
public:
  UnitTestStkMeshSkinning(stk::ParallelMachine pm) : m_comm(pm),  m_num_procs(0), m_rank(0)
  {
    m_num_procs = stk::parallel_machine_size( m_comm );
    m_rank = stk::parallel_machine_rank( m_comm );
  }

  void test_skinning();

  stk::ParallelMachine m_comm;
  int m_num_procs;
  int m_rank;
};

void UnitTestStkMeshSkinning::test_skinning()
{
  // This test will only work for np=1
  if (m_num_procs > 1) {
    return;
  }

  stk::mesh::fixtures::GridFixture grid_mesh(MPI_COMM_WORLD);

  stk::mesh::BulkData& bulk_data = grid_mesh.bulk_data();
  stk::mesh::MetaData& fem_meta = grid_mesh.fem_meta();
  const stk::mesh::EntityRank element_rank = stk::topology::ELEMENT_RANK;

  // Create a part for the skin and the shells
  stk::mesh::Part & skin_part = fem_meta.declare_part("skin_part");
  stk::mesh::Part & shell_part = fem_meta.declare_part_with_topology("shell_part", stk::topology::SHELL_LINE_2);
  fem_meta.commit();

  // Begin modification cycle
  grid_mesh.bulk_data().modification_begin();

  // Generate the plain grid
  grid_mesh.generate_grid();

  // Add the shells
  std::vector<size_t> count;
  stk::mesh::Selector locally_owned(fem_meta.locally_owned_part());
  stk::mesh::count_entities(locally_owned, bulk_data, count);
  const unsigned num_shell_1_faces = 4;
  const unsigned num_shell_2_faces = 2;
  const unsigned num_shell_faces = num_shell_1_faces + num_shell_2_faces;
  const unsigned num_entities = count[NODE_RANK] +
      count[element_rank];

  stk::mesh::PartVector shell_parts;
  shell_parts.push_back(&shell_part);

  std::vector<stk::mesh::Entity> shell_faces;
  for (unsigned i = 1; i <= num_shell_faces; ++i) {
    stk::mesh::Entity new_shell = bulk_data.declare_element(num_entities + i,
                                                            shell_parts);
    shell_faces.push_back(new_shell);
  }

  // Set up relationships for shells

  // declare shell relationships for first shell
  unsigned node_list_1[5] = {21, 26, 31, 36, 41};
  for (unsigned i = 0; i < num_shell_1_faces; ++i) {
    stk::mesh::Entity shell = shell_faces[i];
    stk::mesh::Entity node1 = bulk_data.get_entity(NODE_RANK, node_list_1[i]);
    stk::mesh::Entity node2 = bulk_data.get_entity(NODE_RANK, node_list_1[i+1]);
    bulk_data.declare_relation(shell, node1, 0);
    bulk_data.declare_relation(shell, node2, 1);
  }

  // declare shell relationships for second shell
  unsigned node_list_2[3] = {31, 36, 41};
  for (unsigned i = 0; i < num_shell_2_faces; ++i) {
    stk::mesh::Entity shell = shell_faces[i + num_shell_1_faces];
    stk::mesh::Entity node1 = bulk_data.get_entity(NODE_RANK, node_list_2[i]);
    stk::mesh::Entity node2 = bulk_data.get_entity(NODE_RANK, node_list_2[i+1]);
    bulk_data.declare_relation(shell, node1, 0);
    bulk_data.declare_relation(shell, node2, 1);
  }

  grid_mesh.bulk_data().modification_end();

  // skin the boundary
  {
    stk::mesh::PartVector add_parts(1,&skin_part);
    stk::mesh::skin_mesh(bulk_data, add_parts);
  }

  // Grab the skin entities
  stk::mesh::Selector skin_selector(skin_part);
  const stk::mesh::BucketVector& edge_buckets = bulk_data.buckets(stk::topology::EDGE_RANK);
  std::vector<stk::mesh::Entity> skin_entities;
  stk::mesh::get_selected_entities(skin_selector, edge_buckets, skin_entities);

  unsigned num_expected_skin_entites = 16;
  EXPECT_EQ(num_expected_skin_entites, skin_entities.size());

  // Map the element id to the number of skins associated with that element

  std::map<unsigned, unsigned> results;
  std::map<unsigned, unsigned> expected_results;

  expected_results[1] = 2;
  expected_results[2] = 1;
  expected_results[3] = 1;
  expected_results[4] = 1;
  expected_results[5] = 1;
  expected_results[9] = 1;
  expected_results[13] = 2;
  expected_results[14] = 1;
  expected_results[15] = 1;
  expected_results[16] = 1;
  expected_results[42] = 1;
  expected_results[43] = 1;
  expected_results[44] = 1;
  expected_results[45] = 1;
  expected_results[46] = 1;
  expected_results[47] = 1;

  // Vector of of vector of entities (shells) that are expected to share a skin

  std::set<std::set<unsigned> > sharing;
  std::set<std::set<unsigned> > expected_sharing;

  std::set<unsigned> temp;
  temp.insert(44);
  temp.insert(46);
  expected_sharing.insert(temp);

  temp.clear();
  temp.insert(45);
  temp.insert(47);
  expected_sharing.insert(temp);

  // map skin-id to ids of elements it is attached to; we will use this to
  // compute sharing
  for (std::vector<stk::mesh::Entity>::const_iterator
       itr = skin_entities.begin(); itr != skin_entities.end(); ++itr)
  {
    stk::mesh::Entity const *elem_itr = bulk_data.begin_elements(*itr);
    stk::mesh::Entity const *elem_end = bulk_data.end_elements(*itr);
    bool has_multiple = bulk_data.num_elements(*itr) > 1;
    std::set<unsigned> sharing_elements;
    for ( ; elem_itr != elem_end ; ++elem_itr )
    {
      unsigned elem_id = bulk_data.identifier(*elem_itr);

      if (results.find(elem_id) != results.end()) {
        ++results[elem_id];
      }
      else {
        results[elem_id] = 1;
      }

      if (has_multiple) {
        sharing_elements.insert(elem_id);
      }
    }
    if (has_multiple) {
      sharing.insert(sharing_elements);
    }
  }

  EXPECT_TRUE(results == expected_results);
  EXPECT_TRUE(sharing == expected_sharing);
}

