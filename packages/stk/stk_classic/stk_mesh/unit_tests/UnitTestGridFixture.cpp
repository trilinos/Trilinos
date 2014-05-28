/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stk_util/unit_test_support/stk_utest_macros.hpp>
#include <Shards_BasicTopologies.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/BulkModification.hpp>

#include <stk_mesh/fem/BoundaryAnalysis.hpp>
#include <stk_mesh/fem/SkinMesh.hpp>

#include <stk_mesh/fixtures/GridFixture.hpp>

#include <stk_util/parallel/ParallelReduce.hpp>

#include <iomanip>
#include <algorithm>

using stk::mesh::MetaData;
using stk::mesh::fem::FEMMetaData;
using stk::mesh::Part;
using stk::mesh::PartVector;
using stk::mesh::PartRelation;
using stk::mesh::Entity;
using stk::mesh::EntityVector;
using stk::mesh::EntityRank;
using stk::mesh::Selector;
using stk::mesh::BulkData;
using stk::ParallelMachine;
using std::cout;
using std::endl;

STKUNIT_UNIT_TEST( UnitTestGridFixture, test_gridfixture )
{
  //Coverage of GridFixture, Hexfixture, BoxFixture,QuadFixture
  //and RingFixture in fixture directory for more than one
  //processor.
  stk::mesh::fixtures::GridFixture grid_mesh(MPI_COMM_WORLD);

  stk::mesh::BulkData& bulk_data = grid_mesh.bulk_data();
  stk::mesh::fem::FEMMetaData& fem_meta = grid_mesh.fem_meta();
  const stk::mesh::EntityRank elem_rank = fem_meta.element_rank();
  
  int  size , rank;
  rank = stk::parallel_machine_rank( MPI_COMM_WORLD );
  size = stk::parallel_machine_size( MPI_COMM_WORLD );

  // Create a part for the shells
  stk::mesh::fem::CellTopology line_top(shards::getCellTopologyData<shards::ShellLine<2> >());
  stk::mesh::Part & shell_part = fem_meta.declare_part("shell_part", line_top);

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

  std::vector<stk::mesh::Entity*> shell_faces;

  unsigned id_base = 0;
  unsigned id_offset = 500; // a safe offset to avoid id overlap
  //Start at 1 so as not to have same element on different processors
  for (id_base = 1; id_base <= num_shell_faces; ++id_base) {

    int new_id = rank * num_shell_faces + id_base;
    stk::mesh::Entity& new_shell = bulk_data.declare_entity(elem_rank,
                                                            id_offset + new_id,
                                                            shell_parts);
    shell_faces.push_back(&new_shell);
  }

   bulk_data.modification_end();
}


