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

#include <stk_mesh/fem/FEMMetaData.hpp>

namespace {

const stk_classic::mesh::EntityRank NODE_RANK = stk_classic::mesh::fem::FEMMetaData::NODE_RANK;

STKUNIT_UNIT_TEST( UnitTestStkMeshGenerateNewEntities , testUnit )
{
  // Test BulkData's generate_new_entities method.

  stk_classic::ParallelMachine pm(MPI_COMM_WORLD);

  const int spatial_dimension = 3;
  stk_classic::mesh::fem::FEMMetaData meta_data( spatial_dimension );
  stk_classic::mesh::BulkData bulk_data( stk_classic::mesh::fem::FEMMetaData::get_meta_data(meta_data) , pm );

  meta_data.commit();

  const stk_classic::mesh::PartVector no_parts;

  bulk_data.modification_begin();

  bulk_data.declare_entity(NODE_RANK, bulk_data.parallel_rank() + 1, no_parts);

  bulk_data.modification_end();

  // Create a request vector for 2 new nodes on each processor
  size_t num_nodes_requested = 2;
  std::vector<size_t> requests(meta_data.entity_rank_count(), 0);
  requests[0] = num_nodes_requested;

  bulk_data.modification_begin();

  // generate_new_entities creates new blank entities of the requested ranks
  stk_classic::mesh::EntityVector new_nodes;
  bulk_data.generate_new_entities(requests, new_nodes);
  STKUNIT_ASSERT_EQ(new_nodes.size(), num_nodes_requested);

  // confirm that the nodes we created earlier are not in the new entities
  for (stk_classic::mesh::EntityVector::const_iterator itr = new_nodes.begin();
       itr != new_nodes.end(); ++itr) {
    STKUNIT_ASSERT_GT((*itr)->identifier(), bulk_data.parallel_size());
  }

  bulk_data.modification_end();
}

}
