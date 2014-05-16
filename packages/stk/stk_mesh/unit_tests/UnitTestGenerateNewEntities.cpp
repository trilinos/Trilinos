/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stddef.h>                     // for size_t
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine
#include <gtest/gtest.h>
#include <vector>                       // for vector, etc
#include "mpi.h"                        // for MPI_COMM_WORLD, etc
#include "stk_mesh/base/Types.hpp"      // for EntityVector, EntityRank, etc
#include "stk_topology/topology.hpp"    // for topology, etc




namespace {

const stk::mesh::EntityRank NODE_RANK = stk::topology::NODE_RANK;

TEST( UnitTestStkMeshGenerateNewEntities , testUnit )
{
  // Test BulkData's generate_new_entities method.

  stk::ParallelMachine pm(MPI_COMM_WORLD);

  const int spatial_dimension = 3;
  stk::mesh::MetaData meta_data( spatial_dimension );
  stk::mesh::BulkData bulk_data( meta_data , pm );

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
  ASSERT_EQ(new_nodes.size(), num_nodes_requested);

  // confirm that the nodes we created earlier are not in the new entities
  for (stk::mesh::EntityVector::const_iterator itr = new_nodes.begin();
       itr != new_nodes.end(); ++itr) {
    ASSERT_GT(static_cast<int>(bulk_data.identifier(*itr)), bulk_data.parallel_size());
  }

  bulk_data.modification_end();
}

}
