/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <unit_tests/stk_utest_macros.hpp>
#include <Shards_BasicTopologies.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/GetBuckets.hpp>

#include <stk_mesh/fem/BoundaryAnalysis.hpp>
#include <stk_mesh/fem/EntityRanks.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>
#include <stk_mesh/fem/EntityRanks.hpp>

#include <unit_tests/UnitTestGridMeshFixture.hpp>

#include <iomanip>
#include <algorithm>

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

  GridMeshFixture grid_mesh(MPI_COMM_WORLD);

  if (m_rank > 0) return;

  stk::mesh::BulkData& bulk_data = grid_mesh.bulk_data();

  const std::vector<stk::mesh::Entity*>& closure = grid_mesh.get_predefined_closure();

  stk::mesh::EntitySideVector boundary;
  stk::mesh::boundary_analysis(bulk_data, closure, stk::mesh::Face, boundary);
  STKUNIT_EXPECT_TRUE(!boundary.empty());

  std::vector<std::pair<std::pair<unsigned, unsigned>, std::pair<unsigned, unsigned> > > results;
  std::vector<std::pair<std::pair<unsigned, unsigned>, std::pair<unsigned, unsigned> > > expected_results;

  {
    std::pair<unsigned, unsigned> inside(6, 0);
    std::pair<unsigned, unsigned> outside(5, 2);
    expected_results.push_back(
      std::pair<std::pair<unsigned, unsigned>, std::pair<unsigned, unsigned> >(inside, outside));
  }

  {
    std::pair<unsigned, unsigned> inside(6, 3);
    std::pair<unsigned, unsigned> outside(2, 1);
    expected_results.push_back(
      std::pair<std::pair<unsigned, unsigned>, std::pair<unsigned, unsigned> >(inside, outside));
  }

  {
    std::pair<unsigned, unsigned> inside(7, 2);
    std::pair<unsigned, unsigned> outside(8, 0);
    expected_results.push_back(
      std::pair<std::pair<unsigned, unsigned>, std::pair<unsigned, unsigned> >(inside, outside));
  }

  {
    std::pair<unsigned, unsigned> inside(7, 2);
    std::pair<unsigned, unsigned> outside(43, 0);
    expected_results.push_back(
      std::pair<std::pair<unsigned, unsigned>, std::pair<unsigned, unsigned> >(inside, outside));
  }

  {
    std::pair<unsigned, unsigned> inside(7, 3);
    std::pair<unsigned, unsigned> outside(3, 1);
    expected_results.push_back(
      std::pair<std::pair<unsigned, unsigned>, std::pair<unsigned, unsigned> >(inside, outside));
  }

  {
    std::pair<unsigned, unsigned> inside(10, 0);
    std::pair<unsigned, unsigned> outside(9, 2);
    expected_results.push_back(
      std::pair<std::pair<unsigned, unsigned>, std::pair<unsigned, unsigned> >(inside, outside));
  }

  {
    std::pair<unsigned, unsigned> inside(11, 2);
    std::pair<unsigned, unsigned> outside(12, 0);
    expected_results.push_back(
      std::pair<std::pair<unsigned, unsigned>, std::pair<unsigned, unsigned> >(inside, outside));
  }

  {
    std::pair<unsigned, unsigned> inside(11, 2);
    std::pair<unsigned, unsigned> outside(44, 0);
    expected_results.push_back(
      std::pair<std::pair<unsigned, unsigned>, std::pair<unsigned, unsigned> >(inside, outside));
  }

  {
    std::pair<unsigned, unsigned> inside(14, 0);
    std::pair<unsigned, unsigned> outside(13, 2);
    expected_results.push_back(
      std::pair<std::pair<unsigned, unsigned>, std::pair<unsigned, unsigned> >(inside, outside));
  }

  {
    std::pair<unsigned, unsigned> inside(14, 1);
    std::pair<unsigned, unsigned> outside(0, 0);
    expected_results.push_back(
      std::pair<std::pair<unsigned, unsigned>, std::pair<unsigned, unsigned> >(inside, outside));
  }

  {
    std::pair<unsigned, unsigned> inside(15, 1);
    std::pair<unsigned, unsigned> outside(0, 0);
    expected_results.push_back(
      std::pair<std::pair<unsigned, unsigned>, std::pair<unsigned, unsigned> >(inside, outside));
  }

  {
    std::pair<unsigned, unsigned> inside(15, 2);
    std::pair<unsigned, unsigned> outside(16, 0);
    expected_results.push_back(
      std::pair<std::pair<unsigned, unsigned>, std::pair<unsigned, unsigned> >(inside, outside));
  }

  {
    std::pair<unsigned, unsigned> inside(15, 2);
    std::pair<unsigned, unsigned> outside(45, 0);
    expected_results.push_back(
      std::pair<std::pair<unsigned, unsigned>, std::pair<unsigned, unsigned> >(inside, outside));
  }


  for (stk::mesh::EntitySideVector::iterator itr = boundary.begin(); itr != boundary.end(); ++itr)
  {
    stk::mesh::EntitySide& side = *itr;
    stk::mesh::EntitySideComponent& inside_closure = side.inside;
    stk::mesh::EntityId inside_id = inside_closure.entity != NULL ? inside_closure.entity->identifier() : 0;
    stk::mesh::EntityId inside_side = inside_closure.entity != NULL ? inside_closure.side_ordinal : 0;
    stk::mesh::EntitySideComponent& outside_closure = side.outside;
    stk::mesh::EntityId outside_id = outside_closure.entity != NULL ? outside_closure.entity->identifier() : 0;
    stk::mesh::EntityId outside_side = outside_closure.entity != NULL ? outside_closure.side_ordinal : 0;

    std::pair<unsigned, unsigned> inside(inside_id, inside_side);
    std::pair<unsigned, unsigned> outside(outside_id, outside_side);
    results.push_back(std::pair<std::pair<unsigned, unsigned>, std::pair<unsigned, unsigned> >(inside, outside));
  }
  STKUNIT_EXPECT_TRUE(results == expected_results);
}

void UnitTestStkMeshBoundaryAnalysis::test_boundary_analysis_null_topology()
{
  //test on boundary_analysis for closure with a NULL topology - coverage of lines 39-40 of BoundaryAnalysis.cpp

  //create new meta, bulk and boundary for this test
  stk::mesh::MetaData meta( stk::mesh::fem_entity_rank_names() );

  //declare part with topology = NULL
  stk::mesh::Part & quad_part = meta.declare_part("quad_part", stk::mesh::Face);
  meta.commit();

  stk::ParallelMachine comm(MPI_COMM_WORLD);
  stk::mesh::BulkData bulk ( meta , comm , 100 );

  stk::mesh::EntitySideVector boundary;
  std::vector<stk::mesh::Entity*> newclosure;

  stk::mesh::PartVector face_parts;
  face_parts.push_back(&quad_part);

  bulk.modification_begin();
  if (m_rank == 0) {
    stk::mesh::Entity & new_face = bulk.declare_entity(stk::mesh::Face, 1, face_parts);
    newclosure.push_back(&new_face);
  }

  stk::mesh::boundary_analysis(bulk, newclosure, stk::mesh::Face, boundary);
  /*
  STKUNIT_EXPECT_TRUE(!boundary.empty());
  */

  bulk.modification_end();
}
