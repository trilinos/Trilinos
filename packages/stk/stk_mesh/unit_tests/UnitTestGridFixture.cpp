/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
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
#include <stk_mesh/fem/EntityRanks.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>

#include <stk_mesh/fixtures/GridFixture.hpp>

#include <stk_util/parallel/ParallelReduce.hpp>

#include <iomanip>
#include <algorithm>

using stk::mesh::MetaData;
using stk::mesh::TopologicalMetaData;
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
  stk::mesh::MetaData& meta_data = grid_mesh.meta_data();
  stk::mesh::TopologicalMetaData& top_data = grid_mesh.top_data();

  int  size , rank;
  rank = stk::parallel_machine_rank( MPI_COMM_WORLD );
  size = stk::parallel_machine_size( MPI_COMM_WORLD );

  // Create a part for the shells
  stk::mesh::Part & shell_part = top_data.declare_part<shards::ShellLine<2> >("shell_part");

  meta_data.commit();

  // Begin modification cycle
  bulk_data.modification_begin();

  // Generate the plain grid

  grid_mesh.generate_grid();

  // Add the shells

  std::vector<unsigned> count;
  stk::mesh::Selector locally_owned(meta_data.locally_owned_part());
  stk::mesh::count_entities(locally_owned, bulk_data, count);
  const unsigned num_shell_1_faces = 4*size + rank;
  const unsigned num_shell_2_faces = 2*size + rank;
  const unsigned num_shell_faces = num_shell_1_faces + num_shell_2_faces;

  const unsigned num_entities = count[top_data.node_rank] +
                                count[top_data.element_rank];

  stk::mesh::PartVector shell_parts;
  shell_parts.push_back(&shell_part);

  std::vector<stk::mesh::Entity*> shell_faces;

  unsigned id_base = 0;
  //Start at 1 so as not to have same element on different processors
  for (id_base = 1; id_base <= num_shell_faces; ++id_base) {

    int new_id = size * id_base + rank;
    stk::mesh::Entity& new_shell = bulk_data.declare_entity(top_data.element_rank,
                                                            num_entities + new_id,
                                                            shell_parts);
    shell_faces.push_back(&new_shell);
  }

   bulk_data.modification_end();
}


