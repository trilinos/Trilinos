/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <Shards_BasicTopologies.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>

#include <stk_mesh/fem/TopologyHelpers.hpp>
#include <stk_mesh/fem/DefaultFEM.hpp>

using stk::mesh::fem::NODE_RANK;

STKUNIT_UNIT_TEST( UnitTestStkMeshGenerateNewEntities , testUnit )
{
  // Test BulkData's generate_new_entities method.

  stk::ParallelMachine pm(MPI_COMM_WORLD);

  const int spatial_dimension = 3;
  stk::mesh::MetaData meta_data( stk::mesh::fem::entity_rank_names(spatial_dimension) );
  stk::mesh::BulkData bulk_data( meta_data , pm );
  stk::mesh::DefaultFEM top_data( meta_data, spatial_dimension );

  meta_data.commit();

  const stk::mesh::PartVector no_parts;

  bulk_data.modification_begin();

  bulk_data.declare_entity(NODE_RANK, bulk_data.parallel_rank() + 1, no_parts);

  bulk_data.modification_end();

  // Create a request vector for 2 new nodes on each processor
  size_t num_nodes_requested = 2;
  std::vector<size_t> requests(meta_data.entity_rank_count(), 0);
  requests[0] = num_nodes_requested;

  bulk_data.modification_begin();

  // generate_new_entities creates new blank entities of the requested ranks
  stk::mesh::EntityVector new_nodes;
  bulk_data.generate_new_entities(requests, new_nodes);
  STKUNIT_ASSERT_EQ(new_nodes.size(), num_nodes_requested);

  // confirm that the nodes we created earlier are not in the new entities
  for (stk::mesh::EntityVector::const_iterator itr = new_nodes.begin();
       itr != new_nodes.end(); ++itr) {
    STKUNIT_ASSERT_GT((*itr)->identifier(), bulk_data.parallel_size());
  }

  bulk_data.modification_end();
}
