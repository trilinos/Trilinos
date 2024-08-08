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

#include <gtest/gtest.h>                // for TEST
#include <iostream>                     // for cout, endl
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/Entity.hpp>     // for Entity
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/Types.hpp>      // for EntityRank, PartVector, etc
#include <stk_util/parallel/Parallel.hpp>  // for parallel_machine_rank, etc
#include <vector>                       // for vector
#include "mpi.h"                        // for MPI_COMM_WORLD

#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_unit_test_utils/stk_mesh_fixtures/GridFixture.hpp"  // for GridFixture
namespace stk { namespace mesh { class Part; } }
namespace stk { namespace mesh { class Selector; } }

using stk::mesh::MetaData;
using stk::mesh::Part;
using stk::mesh::PartVector;
using stk::mesh::Entity;
using stk::mesh::EntityVector;
using stk::mesh::EntityRank;
using stk::mesh::Selector;
using stk::mesh::BulkData;
using stk::ParallelMachine;
using std::cout;
using std::endl;

TEST( UnitTestGridFixture, test_gridfixture )
{
  //Coverage of GridFixture, Hexfixture, BoxFixture,QuadFixture
  //and RingFixture in fixture directory for more than one
  //processor.
  stk::mesh::fixtures::GridFixture grid_mesh(MPI_COMM_WORLD);

  stk::mesh::BulkData& bulk_data = grid_mesh.bulk_data();
  stk::mesh::MetaData& fem_meta = grid_mesh.fem_meta();

  int rank = stk::parallel_machine_rank( MPI_COMM_WORLD );
  int size = stk::parallel_machine_size( MPI_COMM_WORLD );

  // Create a part for the shells
  stk::mesh::Part & shell_part = fem_meta.declare_part_with_topology("shell_part", stk::topology::SHELL_LINE_2);

  fem_meta.commit();

  // Generate the plain grid
  bulk_data.modification_begin();
  grid_mesh.generate_grid();
  bulk_data.modification_end();

  // Add the shells
  bulk_data.modification_begin();

  const unsigned num_shell_1_faces = 4*size + rank;
  const unsigned num_shell_2_faces = 2*size + rank;
  const unsigned num_shell_faces = num_shell_1_faces + num_shell_2_faces;

  stk::mesh::PartVector shell_parts;
  shell_parts.push_back(&shell_part);

  std::vector<stk::mesh::Entity> shell_faces;

  unsigned id_base = 0;
  unsigned id_offset = 500; // a safe offset to avoid id overlap
  //Start at 1 so as not to have same element on different processors
  for (id_base = 1; id_base <= num_shell_faces; ++id_base) {

    int new_id = rank * num_shell_faces + id_base;
    stk::mesh::Entity new_shell = bulk_data.declare_element(id_offset + new_id, shell_parts);
    for(unsigned node_ord = 0 ; node_ord < 2; ++node_ord)
    {
      stk::mesh::Entity new_node = bulk_data.declare_node((node_ord+10)*(id_offset + new_id));
      bulk_data.declare_relation( new_shell , new_node , node_ord);
    }
    shell_faces.push_back(new_shell);
  }

  bulk_data.modification_end();
}


