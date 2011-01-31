
/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <Shards_BasicTopologies.hpp>

#include <stk_mesh/fixtures/HexFixture.hpp>

#include <stk_mesh/base/BulkModification.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>

#include <stk_mesh/fem/TopologyHelpers.hpp>

#include <stk_util/parallel/ParallelReduce.hpp>

STKUNIT_UNIT_TEST( UnitTestDeclareElement , inject_shell )
{
  // This tests creates a small HexFixture with two hexes then, in a separate
  // modification cycle, inserts a shell between the two elements.

  stk::ParallelMachine pm = MPI_COMM_WORLD ;

  // Create the fixture, adding a part for the shell

  stk::mesh::fixtures::HexFixture fixture( pm , 2 , 1 , 1 );

  const unsigned p_rank = fixture.m_bulk_data.parallel_rank();

  stk::mesh::Part & shell_part = stk::mesh::declare_part<shards::ShellQuadrilateral<4> >( fixture.m_meta_data, "shell_part");

  fixture.m_meta_data.commit();

  fixture.generate_mesh();

  stk::mesh::Entity * elem = fixture.elem( 0 , 0 , 0 );

  fixture.m_bulk_data.modification_begin();

  bool no_throw = true;

  // Whoever owns the 0,0,0 element create the shell and insert it between
  // the two elements.
  if ( elem != NULL && p_rank == elem->owner_rank() ) {
    stk::mesh::EntityId elem_node[4] ;
    elem_node[0] = fixture.node_id( 1, 0, 0 );
    elem_node[1] = fixture.node_id( 1, 1, 0 );
    elem_node[2] = fixture.node_id( 1, 1, 1 );
    elem_node[3] = fixture.node_id( 1, 0, 1 );

    stk::mesh::EntityId elem_id = 3;

    try {
      stk::mesh::declare_element( fixture.m_bulk_data, shell_part, elem_id, elem_node);
    }
    catch (...) {
      no_throw = false;
    }

  }
  fixture.m_bulk_data.modification_end();

  STKUNIT_EXPECT_TRUE(no_throw);
}
