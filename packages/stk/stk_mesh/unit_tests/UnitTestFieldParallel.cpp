/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <stdexcept>
#include <sstream>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/GetEntities.hpp>

#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>

using stk::mesh::MetaData;

namespace {

const stk::mesh::EntityRank NODE_RANK = stk::topology::NODE_RANK;

STKUNIT_UNIT_TEST(FieldParallel, parallel_sum)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD ;

  const int spatial_dimension = 3;
  stk::mesh::MetaData meta_data( spatial_dimension );
  stk::mesh::BulkData bulk_data( meta_data , pm );

  stk::mesh::Field<double>  & field_1 = meta_data.declare_field< stk::mesh::Field<double> >( "field_1" );
  stk::mesh::Field<double>  & field_2 = meta_data.declare_field< stk::mesh::Field<double> >( "field_2" );
  stk::mesh::Field<double>  & field_3 = meta_data.declare_field< stk::mesh::Field<double> >( "field_3" );

  stk::mesh::put_field( field_1 , NODE_RANK , meta_data.universal_part() );
  stk::mesh::put_field( field_2 , NODE_RANK , meta_data.universal_part() );
  stk::mesh::put_field( field_3 , NODE_RANK , meta_data.universal_part() );

  stk::mesh::Part& block_1 = meta_data.declare_part<shards::Hexahedron<8> >("block_1");

  meta_data.commit();

  bulk_data.modification_begin();

  // Declare 8 nodes on each proc
  const int num_nodes_per_proc = 8;
  int this_proc = bulk_data.parallel_rank();
  stk::mesh::EntityId my_first_id = this_proc*4 + 1;
  std::vector<stk::mesh::Entity> nodes;

  int end_id = my_first_id + num_nodes_per_proc;
  for ( int id = my_first_id ; id < end_id ; ++id ) {
    nodes.push_back(bulk_data.declare_entity( NODE_RANK , id));
  }

  stk::mesh::EntityId elem_id = this_proc+1;
  stk::mesh::Entity elem = bulk_data.declare_entity(stk::topology::ELEMENT_RANK, elem_id, block_1);
  const int num_nodes_per_elem = 8;
  for(int i=0; i<num_nodes_per_elem; ++i) {
    bulk_data.declare_relation(elem, nodes[i], i);
  }

  bulk_data.modification_end();

  //first, insist that we have 4 shared nodes.
//  stk::mesh::Selector shared = meta_data.globally_shared_part();
//  unsigned num_shared_nodes = stk::mesh::count_selected_entities(shared, bulk_data.buckets(NODE_RANK));
//  STKUNIT_ASSERT_EQUAL(num_shared_nodes, 4u);
//
//  // Go through node_buckets and initialize the fields
//  const std::vector< stk::mesh::Bucket *> & node_buckets =
//    bulk_data.buckets( NODE_RANK );
//
//  const double value = 1.0;
//
//  for ( size_t i=0; i<node_buckets.size(); ++i) {
//    stk::mesh::Bucket & b = *node_buckets[i];
//
//    double* field_1_ptr = reinterpret_cast<double*>(b.field_data_location(field_1));
//    double* field_2_ptr = reinterpret_cast<double*>(b.field_data_location(field_2));
//    double* field_3_ptr = reinterpret_cast<double*>(b.field_data_location(field_3));
//
//    for(size_t j=0; j<b.size(); ++j) {
//      field_1_ptr[0] = value;
//      field_2_ptr[0] = value;
//      field_3_ptr[0] = value;
//    }
//  }
//
//  std::vector<stk::mesh::FieldBase*> field_vector;
//  field_vector.push_back(&field_1);
//  field_vector.push_back(&field_2);
//  field_vector.push_back(&field_3);
//
//  stk::mesh::parallel_sum(bulk_data, field_vector);

  //now go through the node buckets again and confirm that the field-data values
  //are equal to the number of procs that share the node.
}

} //namespace <anonymous>

