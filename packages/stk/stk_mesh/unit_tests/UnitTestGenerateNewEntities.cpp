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

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>

#include <stk_mesh/fem/EntityRanks.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>
#include <stk_mesh/fem/TopologicalMetaData.hpp>

STKUNIT_UNIT_TEST( UnitTestStkMeshGenerateNewEntities , testUnit ) {

  stk::ParallelMachine pm(MPI_COMM_WORLD);

  const int spatial_dimension = 3;
  stk::mesh::MetaData meta_data( stk::mesh::TopologicalMetaData::entity_rank_names(spatial_dimension) );
  stk::mesh::BulkData bulk_data( meta_data , pm );
  stk::mesh::TopologicalMetaData top_data( meta_data, spatial_dimension );


  meta_data.commit();


  const stk::mesh::PartVector no_parts;

  bulk_data.modification_begin();

  bulk_data.declare_entity(top_data.node_rank, bulk_data.parallel_rank() + 1, no_parts);

  bulk_data.modification_end();

  std::vector<size_t> requests(meta_data.entity_rank_count(), 0);
  requests[0] = 2;

  bulk_data.modification_begin();

  // generate_new_entities creates new blank entities of the requested ranks
  stk::mesh::EntityVector new_nodes;
  //STKUNIT_ASSERT_NO_THROW(
      bulk_data.generate_new_entities(requests, new_nodes);
        //);

  bulk_data.modification_end();
}
