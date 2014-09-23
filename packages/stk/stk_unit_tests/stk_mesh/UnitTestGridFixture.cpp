/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <iostream>                     // for cout, endl
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/Entity.hpp>     // for Entity
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/Types.hpp>      // for PartVector, EntityRank, etc
#include <stk_mesh/fixtures/GridFixture.hpp>  // for GridFixture
#include <stk_util/parallel/Parallel.hpp>  // for parallel_machine_rank, etc
#include <gtest/gtest.h>
#include <vector>                       // for vector
#include "stk_topology/topology.hpp"    // for topology, etc
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
  const stk::mesh::EntityRank elem_rank = stk::topology::ELEMENT_RANK;
  const stk::mesh::EntityRank node_rank = stk::topology::NODE_RANK;

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
    stk::mesh::Entity new_shell = bulk_data.declare_entity(elem_rank,
                                                            id_offset + new_id,
                                                            shell_parts);
    for(unsigned node_ord = 0 ; node_ord < 2; ++node_ord)
    {
      stk::mesh::Entity new_node = bulk_data.declare_entity(node_rank, (node_ord+10)*(id_offset + new_id));
      bulk_data.declare_relation( new_shell , new_node , node_ord);
    }
    shell_faces.push_back(new_shell);
  }

   bulk_data.modification_end();
}


