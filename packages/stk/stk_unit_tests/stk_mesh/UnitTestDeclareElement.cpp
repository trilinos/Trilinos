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
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/Entity.hpp>     // for Entity
#include <stk_mesh/base/FEMHelpers.hpp>  // for declare_element
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include "mpi.h"                        // for MPI_COMM_WORLD, etc

#include "stk_mesh/base/Types.hpp"      // for EntityIdVector, EntityId
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_unit_test_utils/stk_mesh_fixtures/HexFixture.hpp"  // for HexFixture
#include "stk_util/parallel/Parallel.hpp"  // for ParallelMachine
namespace stk { namespace mesh { class Part; } }

TEST( UnitTestDeclareElement , inject_shell )
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
    stk::mesh::EntityIdVector elem_node(4);
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

  EXPECT_TRUE(no_throw);
}
