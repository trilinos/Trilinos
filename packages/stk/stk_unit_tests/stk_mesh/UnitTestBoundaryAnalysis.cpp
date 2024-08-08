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

#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stddef.h>                     // for size_t
#include <algorithm>                    // for sort
#include <stk_mesh/base/BoundaryAnalysis.hpp>  // for EntitySideComponent, etc
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/Entity.hpp>     // for Entity
#include <stk_mesh/base/FEMHelpers.hpp>  // for get_entity_subcell_id, etc
#include <stk_mesh/base/GetEntities.hpp>  // for count_entities
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/Selector.hpp>   // for Selector
#include <string>                       // for string
#include <utility>                      // for pair, operator==
#include <vector>                       // for vector, vector<>::iterator
#include "mpi.h"                        // for MPI_COMM_WORLD, etc

#include "stk_io/DatabasePurpose.hpp"   // for DatabasePurpose::READ_MESH
#include "stk_io/StkMeshIoBroker.hpp"   // for StkMeshIoBroker
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/EntityLess.hpp"  // for EntityLess
#include "stk_mesh/base/Types.hpp"      // for EntityId, Ordinal, etc
#include "stk_topology/topology.hpp"    // for topology::num_faces
#include "stk_unit_test_utils/stk_mesh_fixtures/GridFixture.hpp"  // for GridFixture
#include "stk_util/parallel/Parallel.hpp"  // for parallel_machine_size, etc
#include "stk_util/util/NamedPair.hpp"
namespace stk { namespace mesh { class Part; } }

static const stk::mesh::EntityRank NODE_RANK = stk::topology::NODE_RANK;

using stk::mesh::MetaData;

class UnitTestStkMeshBoundaryAnalysis {
public:
  UnitTestStkMeshBoundaryAnalysis(stk::ParallelMachine pm) : m_comm(pm),  m_num_procs(0), m_rank(0)
  {
    m_num_procs = stk::parallel_machine_size( m_comm );
    m_rank = stk::parallel_machine_rank( m_comm );
  }

  void test_boundary_analysis();

  stk::ParallelMachine m_comm;
  int m_num_procs;
  int m_rank;
};

namespace {

TEST( UnitTestStkMeshBoundaryAnalysis , testUnit )
{
  UnitTestStkMeshBoundaryAnalysis unit(MPI_COMM_WORLD);
  unit.test_boundary_analysis();
}

} //end namespace

void UnitTestStkMeshBoundaryAnalysis::test_boundary_analysis()
{
  // Test the boundary_analysis algorithm in stk_mesh/base/BoundaryAnalysis.hpp
  // with a boundary that is both shelled and non-shelled.
  //
  // We will be testing the algorithm on the following 2D mesh:
  //
  //  17---18---19---20---21
  //  |  1 |  2 |  3 || 4 |
  //  22---23---24---25---26
  //  |  5 |  6 |  7 || 8 |
  //  27---28---29---30---31
  //  |  9 | 10 | 11 ||12 |
  //  32---33---34---35---36
  //  | 13 | 14 | 15 ||16 |
  //  37---38---39---40---41
  //
  // Note the shells along nodes 20-40.
  //
  // We will be computing the boundary of the closure of
  // elements: 6, 7, 10, 11, 14, 15

  // This test will only work for np=1
  if (m_num_procs > 1) {
    return;
  }

  // set up grid_mesh
  stk::mesh::fixtures::GridFixture grid_mesh(MPI_COMM_WORLD);

  stk::mesh::MetaData& fem_meta = grid_mesh.fem_meta();
  stk::mesh::BulkData& bulk_data = grid_mesh.bulk_data();

  const stk::mesh::EntityRank element_rank = stk::topology::ELEMENT_RANK;

  // make shell part
  stk::mesh::Part& shell_part = fem_meta.declare_part_with_topology("shell_part", stk::topology::SHELL_LINE_2);

  fem_meta.commit();

  bulk_data.modification_begin();
  grid_mesh.generate_grid();

  // Add some shells
  const unsigned num_shells = 4;

  // get a count of entities that have already been created
  std::vector<size_t> count;
  stk::mesh::Selector locally_owned(fem_meta.locally_owned_part());
  stk::mesh::count_entities(locally_owned, bulk_data, count);
  const unsigned num_entities = count[NODE_RANK] + count[element_rank];

  // Declare the shell entities, placing them in the shell part
  std::vector<stk::mesh::Entity> shells;
  stk::mesh::PartVector shell_parts;
  shell_parts.push_back(&shell_part);
  for (unsigned i = 1; i <= num_shells; ++i) {
    stk::mesh::Entity new_shell = bulk_data.declare_element(num_entities + i, shell_parts);
    shells.push_back(new_shell);
  }

  // declare shell relationships
  unsigned node_list[5] = {20, 25, 30, 35, 40};
  for (unsigned i = 0; i < num_shells; ++i) {
    stk::mesh::Entity shell = shells[i];
    stk::mesh::Entity node1 = bulk_data.get_entity(NODE_RANK, node_list[i]);
    stk::mesh::Entity node2 = bulk_data.get_entity(NODE_RANK, node_list[i+1]);
    bulk_data.declare_relation(shell, node1, 0);
    bulk_data.declare_relation(shell, node2, 1);
  }

  bulk_data.modification_end();

  // create the closure we want to analyze
  std::vector<stk::mesh::Entity> closure;
  unsigned num_elems_in_closure = 6;
  stk::mesh::EntityId ids_of_entities_in_closure[] =
  {6, 7, 10, 11, 14, 15, 23, 24, 25, 28, 29, 30, 33, 34, 35, 38, 39, 40};
  for (unsigned i = 0;
       i < sizeof(ids_of_entities_in_closure)/sizeof(stk::mesh::EntityId);
       ++i) {
    stk::mesh::EntityRank rank_of_entity=stk::topology::NODE_RANK;
    if (i < num_elems_in_closure) {
      rank_of_entity = element_rank;
    }
    else {
      rank_of_entity = NODE_RANK;
    }
    stk::mesh::Entity closure_entity =
        bulk_data.get_entity(rank_of_entity, ids_of_entities_in_closure[i]);
    closure.push_back(closure_entity);
  }
  // sort the closure (boundary analysis expects it this way)
  std::sort(closure.begin(), closure.end(), stk::mesh::EntityLess(bulk_data));

  // Run the boundary analysis!
  stk::mesh::EntitySideVector boundary;
  stk::mesh::boundary_analysis(bulk_data, closure, element_rank, boundary);
  EXPECT_TRUE(!boundary.empty());

  // Prepare the expected-results as a vector of pairs of pairs representing
  // ( inside, outside ) where
  // inside is ( element-id, side-ordinal ) of element inside the closure
  // outside is ( element-id, side-ordinal ) of element outside the closure
  // and inside and outside border each other

  typedef std::pair<stk::mesh::EntityId, stk::mesh::Ordinal> BoundaryItem;
  typedef std::pair<BoundaryItem, BoundaryItem>              BoundaryPair;

  // Note that certain sides of elements 7, 11, 15 have a boundary with
  // a shell AND the adjacent element outside the closure.

  BoundaryPair results[] = {
    BoundaryPair(BoundaryItem(6,  0), BoundaryItem(5,  2)),

    BoundaryPair(BoundaryItem(6,  3), BoundaryItem(2,  1)),

    BoundaryPair(BoundaryItem(7,  2), BoundaryItem(8,  0)),
    BoundaryPair(BoundaryItem(7,  2), BoundaryItem(43, 0)),

    BoundaryPair(BoundaryItem(7,  3), BoundaryItem(3,  1)),

    BoundaryPair(BoundaryItem(10, 0), BoundaryItem(9,  2)),

    BoundaryPair(BoundaryItem(11, 2), BoundaryItem(12, 0)),
    BoundaryPair(BoundaryItem(11, 2), BoundaryItem(44, 0)),

    BoundaryPair(BoundaryItem(14, 0), BoundaryItem(13, 2)),

    BoundaryPair(BoundaryItem(14, 1), BoundaryItem(0,  0)),

    BoundaryPair(BoundaryItem(15, 1), BoundaryItem(0,  0)),

    BoundaryPair(BoundaryItem(15, 2), BoundaryItem(16, 0)),
    BoundaryPair(BoundaryItem(15, 2), BoundaryItem(45, 0))
  };

  // Convert the boundary returned by boundary_analysis into a data-structure
  // comparable to expected_results

  BoundaryPair expected_results[sizeof(results)/sizeof(BoundaryPair)];

  unsigned i = 0;
  stk::mesh::EntitySideVector::iterator itr = boundary.begin();

  for (; itr != boundary.end(); ++itr, ++i)
  {
    stk::mesh::EntitySide& side = *itr;
    stk::mesh::EntitySideComponent& inside_closure = side.inside;
    stk::mesh::EntityId inside_id =
        bulk_data.is_valid(inside_closure.entity) ? bulk_data.identifier(inside_closure.entity) : 0;
    stk::mesh::EntityId inside_side =
        bulk_data.is_valid(inside_closure.entity) ? inside_closure.side_ordinal : 0;
    stk::mesh::EntitySideComponent& outside_closure = side.outside;
    stk::mesh::EntityId outside_id =
        bulk_data.is_valid(outside_closure.entity) ? bulk_data.identifier(outside_closure.entity) : 0;
    stk::mesh::EntityId outside_side =
        bulk_data.is_valid(outside_closure.entity) ? outside_closure.side_ordinal : 0;

    expected_results[i] = BoundaryPair(BoundaryItem(inside_id, inside_side),
                                       BoundaryItem(outside_id, outside_side));
  }

  // Check that results match expected results

  EXPECT_EQ(sizeof(results), sizeof(expected_results));

  for (i = 0; i < sizeof(results)/sizeof(BoundaryPair); ++i) {
    EXPECT_TRUE(results[i] == expected_results[i]);
  }
}

TEST(BoundaryAnalysis, get_adjacent_entities)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int numProcs = stk::parallel_machine_size(comm);
  if (numProcs > 1)
  {
    return;
  }

  stk::io::StkMeshIoBroker meshReader(comm);
  std::string mesh_spec("generated:3x3x3");
  meshReader.add_mesh_database(mesh_spec, stk::io::READ_MESH);
  meshReader.create_input_mesh();
  meshReader.populate_bulk_data();

  stk::mesh::BulkData& stkMeshBulkData = meshReader.bulk_data();

  unsigned numEntitiesToTest = 4;
  stk::mesh::EntityId ids[] = { 14, 1, 2, 5 };
  unsigned goldNumConnectedEntities[] = { 6, 3, 4, 5 };
  stk::mesh::EntityId nodeIdsForFace[][4] = {
    {22,       23,      39,      38},
    {23,       27,      43,      39},
    {27,       26,      42,      43},
    {22,       38,      42,      26},
    {22,       26,      27,      23},
    {38,       39,      43,      42},
    {1 ,       2 ,      18,      17},
    {2 ,       6 ,      22,      18},
    {6 ,       5 ,      21,      22},
    {1 ,       17,      21,      5 },
    {1 ,       5 ,      6 ,      2 },
    {17,       18,      22,      21},
    {2 ,       3 ,      19,      18},
    {3 ,       7 ,      23,      19},
    {7 ,       6 ,      22,      23},
    {2 ,       18,      22,      6 },
    {2 ,       6 ,      7 ,      3 },
    {18,       19,      23,      22},
    {6 ,       7 ,      23,      22},
    {7 ,       11,      27,      23},
    {11,       10,      26,      27},
    {6 ,       22,      26,      10},
    {6 ,       10,      11,      7 },
    {22,       23,      27,      26}
  };

  unsigned testCounter=0;
  for (unsigned int i=0;i<numEntitiesToTest;i++)
  {
    stk::mesh::Entity element = stkMeshBulkData.get_entity(stk::topology::ELEM_RANK, ids[i]);

    std::vector<stk::mesh::EntitySideComponent> adjacent_entities;
    stk::topology elemTopology = stkMeshBulkData.bucket(element).topology();

    unsigned numConnectedEntities = 0;
    for(unsigned faceIndex = 0; faceIndex < elemTopology.num_faces(); ++faceIndex)
    {
      stk::mesh::get_adjacent_entities(stkMeshBulkData,
                                       element,
                                       stk::topology::FACE_RANK,
                                       faceIndex,
                                       adjacent_entities);

      numConnectedEntities += adjacent_entities.size();
      stk::mesh::EntityVector subcell_nodes;
      stk::topology stksubcell_topo = stk::mesh::get_subcell_nodes(stkMeshBulkData,
                                                                   element,
                                                                   stk::topology::FACE_RANK,
                                                                   faceIndex,
                                                                   subcell_nodes);
      EXPECT_FALSE(stksubcell_topo == stk::topology::INVALID_TOPOLOGY);

      for(size_t j = 0; j < subcell_nodes.size(); j++)
      {
        EXPECT_EQ(nodeIdsForFace[testCounter][j], stkMeshBulkData.identifier(subcell_nodes[j]));
      }

      size_t local_subcell_num = stk::mesh::get_entity_subcell_id(stkMeshBulkData,
                                                                  element,
                                                                  stk::topology::FACE_RANK,
                                                                  stksubcell_topo,
                                                                  subcell_nodes);
      EXPECT_EQ(faceIndex, local_subcell_num);
      testCounter++;
    }
    EXPECT_EQ(goldNumConnectedEntities[i], numConnectedEntities);
  }
}
