/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <algorithm>                    // for sort
#include <stk_mesh/base/BoundaryAnalysis.hpp>  // for EntitySideComponent, etc
#include <stk_mesh/base/BulkData.hpp>   // for EntityLess, BulkData
#include <stk_mesh/base/Entity.hpp>     // for Entity
#include <stk_mesh/base/GetEntities.hpp>  // for count_entities
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/Selector.hpp>   // for Selector
#include <stk_mesh/fixtures/GridFixture.hpp>  // for GridFixture
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine, etc
#include <stk_util/unit_test_support/stk_utest_macros.hpp>
#include <utility>                      // for pair, operator==
#include <vector>                       // for vector, vector<>::iterator
#include "mpi.h"                        // for MPI_COMM_WORLD, etc
#include "stk_mesh/base/Types.hpp"      // for EntityId, Ordinal, etc
#include "stk_topology/topology.hpp"    // for topology, etc
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
  void test_boundary_analysis_null_topology();

  stk::ParallelMachine m_comm;
  int m_num_procs;
  int m_rank;
};

namespace {

STKUNIT_UNIT_TEST( UnitTestStkMeshBoundaryAnalysis , testUnit )
{
  UnitTestStkMeshBoundaryAnalysis unit(MPI_COMM_WORLD);
  unit.test_boundary_analysis();
}

STKUNIT_UNIT_TEST( UnitTestStkMeshBoundaryAnalysis , testNullTopology )
{
  UnitTestStkMeshBoundaryAnalysis unit(MPI_COMM_WORLD);
  unit.test_boundary_analysis_null_topology();
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
  std::vector<unsigned> count;
  stk::mesh::Selector locally_owned(fem_meta.locally_owned_part());
  stk::mesh::count_entities(locally_owned, bulk_data, count);
  const unsigned num_entities = count[NODE_RANK] + count[element_rank];

  // Declare the shell entities, placing them in the shell part
  std::vector<stk::mesh::Entity> shells;
  stk::mesh::PartVector shell_parts;
  shell_parts.push_back(&shell_part);
  for (unsigned i = 1; i <= num_shells; ++i) {
    stk::mesh::Entity new_shell = bulk_data.declare_entity(element_rank,
                                                            num_entities + i,
                                                            shell_parts);
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

  // Run the bounary analysis!
  stk::mesh::EntitySideVector boundary;
  stk::mesh::boundary_analysis(bulk_data, closure, element_rank, boundary);
  STKUNIT_EXPECT_TRUE(!boundary.empty());

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

  STKUNIT_EXPECT_EQ(sizeof(results), sizeof(expected_results));

  for (i = 0; i < sizeof(results)/sizeof(BoundaryPair); ++i) {
    STKUNIT_EXPECT_TRUE(results[i] == expected_results[i]);
  }
}

void UnitTestStkMeshBoundaryAnalysis::test_boundary_analysis_null_topology()
{
  //test on boundary_analysis for closure with a NULL topology - coverage of lines 39-40 of BoundaryAnalysis.cpp

  //create new fem_meta, bulk and boundary for this test
  const int spatial_dimension = 3;
  stk::mesh::MetaData fem_meta(spatial_dimension);

  const stk::mesh::EntityRank side_rank = fem_meta.side_rank();

  //declare part with topology = NULL
  stk::mesh::Part & quad_part = fem_meta.declare_part("quad_part", side_rank);
  fem_meta.commit();

  stk::ParallelMachine comm(MPI_COMM_WORLD);
  stk::mesh::BulkData bulk ( fem_meta , comm , 100 );

  stk::mesh::EntitySideVector boundary;
  std::vector<stk::mesh::Entity> newclosure;

  stk::mesh::PartVector face_parts;
  face_parts.push_back(&quad_part);

  bulk.modification_begin();
  if (m_rank == 0) {
    stk::mesh::Entity new_face = bulk.declare_entity(side_rank, 1, face_parts);
    newclosure.push_back(new_face);
  }

  stk::mesh::boundary_analysis(bulk, newclosure, side_rank, boundary);
  /*
  STKUNIT_EXPECT_TRUE(!boundary.empty());
  */

  bulk.modification_end();
}
