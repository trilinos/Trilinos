/*------------------------------------------------------------------------*/
/*         _        Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <stddef.h>                     // for size_t
#include <sstream>                      // for ostringstream, etc
#include <stk_mesh/base/Bucket.hpp>     // for has_superset, Bucket, etc
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/Field.hpp>      // for Field
#include <stk_mesh/base/GetBuckets.hpp>  // for get_involved_parts
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, put_field
#include <stk_mesh/fixtures/BoxFixture.hpp>  // for BoxFixture
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine, etc
#include <gtest/gtest.h>
#include <string>                       // for string, basic_string, etc
#include <vector>                       // for vector, etc
#include "stk_mesh/base/Types.hpp"      // for PartVector, BucketVector, etc
#include "stk_topology/topology.hpp"    // for topology, etc
namespace stk { namespace mesh { class FieldBase; } }
namespace stk { namespace mesh { class Part; } }
namespace stk { namespace mesh { class Selector; } }
namespace stk { namespace mesh { struct Entity; } }
namespace stk { namespace mesh { class Part; } }

using stk::mesh::MetaData;
using stk::mesh::BulkData;
using stk::mesh::Part;
using stk::mesh::PartVector;
using stk::mesh::EntityRank;
using stk::mesh::EntityId;
using stk::mesh::PairIterEntityComm;
using stk::mesh::Entity;
using stk::mesh::Bucket;
using stk::mesh::BucketIterator;
using stk::mesh::Selector;
using stk::mesh::Field;
using stk::mesh::FieldBase;
using stk::mesh::put_field;
using stk::mesh::BucketVector;

typedef Field<double> ScalarFieldType;

namespace {

TEST(UnitTestingOfBucket, testBucket)
{
  // Unit test the Part functionality in isolation:

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  MPI_Barrier( pm );

  // Create a mesh for testing buckets...

  // Create dummy names for entity ranks to be given to MetaData
  std::vector<std::string> entity_names(5);
  for ( size_t i = 0 ; i < 5 ; ++i ) {
    std::ostringstream name ;
    name << "EntityRank" << i ;
    entity_names[i] = name.str();
  }

  // Create MetaData, BulkData
  unsigned max_bucket_size = 4;
  stk::mesh::fixtures::BoxFixture fixture(pm, max_bucket_size, entity_names);
  MetaData& meta = fixture.fem_meta();
  BulkData& bulk = fixture.bulk_data();
  // Create two scalar fields, temperature and volume. Put temperature
  // on all the nodes and put volume on all the elements.
  unsigned number_of_states = 4;

  ScalarFieldType & temperature =
    meta.declare_field < ScalarFieldType > (stk::topology::NODE_RANK, "temperature" , number_of_states );
  ScalarFieldType & volume =

    meta.declare_field < ScalarFieldType > ( stk::topology::ELEMENT_RANK, "volume" , number_of_states );
  Part & universal     = meta.universal_part ();
  put_field ( temperature , universal );
  put_field ( volume , universal );
  meta.commit();

  // Generate the mesh
  int root_box[3][2] = { { 0,4 } , { 0,5 } , { 0,6 } };
  int local_box[3][2] = { { 0,0 } , { 0,0 } , { 0,0 } };
  bulk.modification_begin();
  fixture.generate_boxes( root_box, local_box );
  ASSERT_TRUE(bulk.modification_end());

  //  First, test for streaming IO;
  {
    std::string gold1;
    gold1 = "Bucket( EntityRank0 : {UNIVERSAL} {OWNS} {FEM_ROOT_CELL_TOPOLOGY_PART_Hexahedron_8} elem_part )";
    Bucket *b1 = bulk.buckets(stk::topology::NODE_RANK)[0];
    std::stringstream  out1_str;
    out1_str << (*b1);
    bool equal = (gold1 == out1_str.str());
    ASSERT_EQ ( equal, true );
  }

  // Second, update state of bucket until circular cue is filled
  {
    /* Need to set some data in state, rotate look for it, rotate 3 more times
       and look for it again */
    for ( size_t i = 0 ; i != 10 ; ++i )
      bulk.update_field_data_states ();
  }

  // next, check has_superset (...) and membership functions
  {
    PartVector tmp(2) ;
    tmp[0] = & meta.universal_part();
    tmp[1] = & meta.locally_owned_part();
    ASSERT_TRUE ( has_superset ( *bulk.buckets(stk::topology::NODE_RANK)[0] , tmp ) );
    ASSERT_TRUE ( bulk.buckets(stk::topology::NODE_RANK)[0]->member_any ( tmp ) );
    ASSERT_TRUE ( bulk.buckets(stk::topology::NODE_RANK)[0]->member_all ( tmp ) );
    ASSERT_TRUE ( bulk.buckets(stk::topology::NODE_RANK)[0]->member ( **meta.get_parts().begin() ) );
  }
}

TEST(UnitTestingOfBucket, bucketSortChangeEntityId)
{
  const unsigned spatialDim=3;
  stk::mesh::MetaData meta(spatialDim);
  stk::mesh::Part& part = meta.declare_part_with_topology("node_part", stk::topology::NODE);
  meta.commit();
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  if (bulk.parallel_size() > 1) {
    return;
  }
  stk::mesh::EntityId nodeID=1;
  bulk.modification_begin();
  bulk.declare_entity(stk::topology::NODE_RANK, nodeID, part);
  nodeID=3;
  bulk.declare_entity(stk::topology::NODE_RANK, nodeID, part);
  nodeID=5;
  bulk.declare_entity(stk::topology::NODE_RANK, nodeID, part);
  bulk.modification_end();

  const stk::mesh::BucketVector& node_buckets_1 = bulk.get_buckets(stk::topology::NODE_RANK, meta.universal_part());
  size_t expected_num_buckets = 1;
  EXPECT_EQ(expected_num_buckets, node_buckets_1.size());
  size_t expected_bucket_size = 3;
  EXPECT_EQ(expected_bucket_size, node_buckets_1[0]->size());

  stk::mesh::Entity node3 = (*node_buckets_1[0])[1];
  stk::mesh::EntityId node3ID = 3;
  EXPECT_EQ(node3ID, bulk.identifier(node3));

  stk::mesh::Entity node5 = (*node_buckets_1[0])[2];

  bulk.modification_begin();
  stk::mesh::EntityId node2ID = 2;
  bulk.change_entity_id(node2ID, node5);
  bulk.modification_end();

  const stk::mesh::BucketVector& node_buckets_2 = bulk.get_buckets(stk::topology::NODE_RANK, meta.universal_part());

  stk::mesh::Entity node2 = (*node_buckets_2[0])[1];
  EXPECT_EQ(node2ID, bulk.identifier(node2));
}

}
