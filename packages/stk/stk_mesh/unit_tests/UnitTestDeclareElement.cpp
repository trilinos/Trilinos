
/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/Entity.hpp>     // for Entity
#include <stk_mesh/base/FEMHelpers.hpp>  // for declare_element
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/fixtures/HexFixture.hpp>  // for HexFixture
#include <stk_util/unit_test_support/stk_utest_macros.hpp>
#include "mpi.h"                        // for MPI_COMM_WORLD, etc
#include "stk_mesh/base/Types.hpp"      // for EntityId
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/parallel/Parallel.hpp"  // for ParallelMachine
namespace stk { namespace mesh { class Part; } }





STKUNIT_UNIT_TEST( UnitTestDeclareElement , inject_shell )
{
  // This tests creates a small HexFixture with two hexes then, in a separate
  // modification cycle, inserts a shell between the two elements.

  stk::ParallelMachine pm = MPI_COMM_WORLD ;

  // Create the fixture, adding a part for the shell

  stk::mesh::fixtures::HexFixture fixture( pm , 2 , 1 , 1 );

  const int p_rank = fixture.m_bulk_data.parallel_rank();

  stk::mesh::Part & shell_part = fixture.m_meta.declare_part_with_topology("shell_part", stk::topology::SHELL_QUAD_4);

  fixture.m_meta.commit();

  fixture.generate_mesh();

  stk::mesh::Entity elem = fixture.elem( 0 , 0 , 0 );

  fixture.m_bulk_data.modification_begin();

  bool no_throw = true;

  // Whoever owns the 0,0,0 element create the shell and insert it between
  // the two elements.
  if ( fixture.m_bulk_data.is_valid(elem) && p_rank == fixture.m_bulk_data.parallel_owner_rank(elem) ) {
    stk::mesh::EntityId elem_node[4];
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
