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

void setup_simple_mesh(stk::mesh::BulkData& bulk_data)
{
    stk::mesh::MetaData& meta_data = bulk_data.mesh_meta_data();

    stk::mesh::Field<double>  & field_1 = meta_data.declare_field< stk::mesh::Field<double> >( stk::topology::NODE_RANK, "field_1" );
    stk::mesh::Field<double>  & field_2 = meta_data.declare_field< stk::mesh::Field<double> >( stk::topology::NODE_RANK, "field_2" );
    stk::mesh::Field<double>  & field_3 = meta_data.declare_field< stk::mesh::Field<double> >( stk::topology::NODE_RANK, "field_3" );

    stk::mesh::put_field( field_1 , meta_data.universal_part() );
    stk::mesh::put_field( field_2 , meta_data.universal_part() );
    stk::mesh::put_field( field_3 , meta_data.universal_part() );

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
    STKUNIT_ASSERT_EQUAL(num_nodes_per_elem, static_cast<int>(nodes.size()));

    for(int i=0; i<num_nodes_per_elem; ++i) {
      bulk_data.declare_relation(elem, nodes[i], i);
    }

    bulk_data.modification_end();
}

STKUNIT_UNIT_TEST(FieldParallel, parallel_sum)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD ;

  const int spatial_dimension = 3;
  stk::mesh::MetaData meta_data( spatial_dimension );
  stk::mesh::BulkData bulk_data( meta_data , pm );

  if (bulk_data.parallel_size() != 2) {
      return;
  }

  setup_simple_mesh(bulk_data);

  stk::mesh::Field<double>  & field_1 = *meta_data.get_field< stk::mesh::Field<double> >( "field_1" );
  stk::mesh::Field<double>  & field_2 = *meta_data.get_field< stk::mesh::Field<double> >( "field_2" );
  stk::mesh::Field<double>  & field_3 = *meta_data.get_field< stk::mesh::Field<double> >( "field_3" );

  //first, insist that we have 4 shared nodes.
  stk::mesh::Selector shared = meta_data.globally_shared_part();
  unsigned num_shared_nodes = stk::mesh::count_selected_entities(shared, bulk_data.buckets(NODE_RANK));
  STKUNIT_ASSERT_EQUAL(num_shared_nodes, 4u);

  // Go through node_buckets and initialize the fields
  const std::vector< stk::mesh::Bucket *> & node_buckets =
    bulk_data.buckets( NODE_RANK );

  const double one = 1.0;
  const int this_proc = bulk_data.parallel_rank();

  for ( size_t i=0; i<node_buckets.size(); ++i) {
    stk::mesh::Bucket & b = *node_buckets[i];

    double* field_1_ptr = reinterpret_cast<double*>(b.field_data_location(field_1));
    double* field_2_ptr = reinterpret_cast<double*>(b.field_data_location(field_2));
    double* field_3_ptr = reinterpret_cast<double*>(b.field_data_location(field_3));

    for(size_t j=0; j<b.size(); ++j) {
      field_1_ptr[j] = one;
      field_2_ptr[j] = this_proc;
      field_3_ptr[j] = this_proc;
    }
  }

  std::vector<stk::mesh::FieldBase*> field_vector;
  field_vector.push_back(&field_1);
  field_vector.push_back(&field_2);
  field_vector.push_back(&field_3);

  stk::mesh::parallel_sum(bulk_data, field_vector);

  //now go through the comm nodes and confirm that the field-data values
  //for field_1 are equal to the number of procs that know about the node.

  const stk::mesh::EntityCommListInfoVector& entity_comm_list = bulk_data.comm_list();
  for(size_t i=0; i<entity_comm_list.size(); ++i) {
    stk::mesh::Entity node = entity_comm_list[i].entity;

    if(!bulk_data.is_matching_rank(field_1, node)) continue;

    const double* field_1_ptr = bulk_data.field_data(field_1, node);
    stk::mesh::PairIterEntityComm entity_comm = bulk_data.entity_comm_sharing(bulk_data.entity_key(node));
    if (entity_comm.size() > 0) {
      const double num_sharing_procs = entity_comm.size() + 1;
      std::cout<<"sum: proc "<<this_proc<<", node "<<bulk_data.identifier(node)<<", num-sharing-procs: "<<num_sharing_procs<<std::endl;
      STKUNIT_ASSERT_EQUAL(field_1_ptr[0], num_sharing_procs);
    }
    else {
        std::cout<<"sum: proc "<<this_proc<<" skipping non-shared node "<<bulk_data.identifier(node)<<std::endl;
    }
  }
}

STKUNIT_UNIT_TEST(FieldParallel, parallel_max)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD ;

  const int spatial_dimension = 3;
  stk::mesh::MetaData meta_data( spatial_dimension );
  stk::mesh::BulkData bulk_data( meta_data , pm );

  if (bulk_data.parallel_size() != 2) {
      return;
  }

  setup_simple_mesh(bulk_data);

  stk::mesh::Field<double>  & field_1 = *meta_data.get_field< stk::mesh::Field<double> >( "field_1" );
  stk::mesh::Field<double>  & field_2 = *meta_data.get_field< stk::mesh::Field<double> >( "field_2" );
  stk::mesh::Field<double>  & field_3 = *meta_data.get_field< stk::mesh::Field<double> >( "field_3" );

  //first, insist that we have 4 shared nodes.
  stk::mesh::Selector shared = meta_data.globally_shared_part();
  unsigned num_shared_nodes = stk::mesh::count_selected_entities(shared, bulk_data.buckets(NODE_RANK));
  STKUNIT_ASSERT_EQUAL(num_shared_nodes, 4u);

  // Go through node_buckets and initialize the fields
  const std::vector< stk::mesh::Bucket *> & node_buckets =
    bulk_data.buckets( NODE_RANK );

  const double one = 1.0;
  const int this_proc = bulk_data.parallel_rank();

  for ( size_t i=0; i<node_buckets.size(); ++i) {
    stk::mesh::Bucket & b = *node_buckets[i];



    double* field_1_ptr = reinterpret_cast<double*>(b.field_data_location(field_1));
    double* field_2_ptr = reinterpret_cast<double*>(b.field_data_location(field_2));
    double* field_3_ptr = reinterpret_cast<double*>(b.field_data_location(field_3));

    for(size_t j=0; j<b.size(); ++j) {
      field_1_ptr[j] = one;
      field_2_ptr[j] = this_proc;
      field_3_ptr[j] = this_proc;
    }
  }

  std::vector<stk::mesh::FieldBase*> field_vector;
  field_vector.push_back(&field_1);
  field_vector.push_back(&field_2);
  field_vector.push_back(&field_3);

  stk::mesh::parallel_max(bulk_data, field_vector);

  //now go through the comm nodes and confirm that the field-data values
  //for field_1 are equal to 1.0, and values for field_2 and field_3 are equal to max-proc.

  const double max_proc = bulk_data.parallel_size()-1;

  const stk::mesh::EntityCommListInfoVector& entity_comm_list = bulk_data.comm_list();
  for(size_t i=0; i<entity_comm_list.size(); ++i) {
    stk::mesh::Entity node = entity_comm_list[i].entity;

    if(!bulk_data.is_matching_rank(field_1, node)) continue;

    const double* field_1_ptr = bulk_data.field_data(field_1, node);
    const double* field_2_ptr = bulk_data.field_data(field_2, node);
    const double* field_3_ptr = bulk_data.field_data(field_3, node);



    stk::mesh::PairIterEntityComm entity_comm = bulk_data.entity_comm_sharing(bulk_data.entity_key(node));
    if (entity_comm.size() > 0) {
      const double num_sharing_procs = entity_comm.size() + 1;
      std::cout<<"max: proc "<<this_proc<<", node "<<bulk_data.identifier(node)<<", num-sharing-procs: "<<num_sharing_procs<<std::endl;
      STKUNIT_ASSERT_EQUAL(field_1_ptr[0], one);
      STKUNIT_ASSERT_EQUAL(field_2_ptr[0], max_proc);
      STKUNIT_ASSERT_EQUAL(field_3_ptr[0], max_proc);
    }
    else {
      std::cout<<"max: proc "<<this_proc<<" skipping non-shared node "<<bulk_data.identifier(node)<<std::endl;
    }
  }
}

STKUNIT_UNIT_TEST(FieldParallel, parallel_min)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD ;

  const int spatial_dimension = 3;
  stk::mesh::MetaData meta_data( spatial_dimension );
  stk::mesh::BulkData bulk_data( meta_data , pm );

  if (bulk_data.parallel_size() != 2) {
      return;
  }

  setup_simple_mesh(bulk_data);

  stk::mesh::Field<double>  & field_1 = *meta_data.get_field< stk::mesh::Field<double> >( "field_1" );
  stk::mesh::Field<double>  & field_2 = *meta_data.get_field< stk::mesh::Field<double> >( "field_2" );
  stk::mesh::Field<double>  & field_3 = *meta_data.get_field< stk::mesh::Field<double> >( "field_3" );

  //first, insist that we have 4 shared nodes.
  stk::mesh::Selector shared = meta_data.globally_shared_part();
  unsigned num_shared_nodes = stk::mesh::count_selected_entities(shared, bulk_data.buckets(NODE_RANK));
  STKUNIT_ASSERT_EQUAL(num_shared_nodes, 4u);

  // Go through node_buckets and initialize the fields
  const std::vector< stk::mesh::Bucket *> & node_buckets =
    bulk_data.buckets( NODE_RANK );

  const double one = 1.0;
  const int this_proc = bulk_data.parallel_rank();

  for ( size_t i=0; i<node_buckets.size(); ++i) {
    stk::mesh::Bucket & b = *node_buckets[i];

    double* field_1_ptr = reinterpret_cast<double*>(b.field_data_location(field_1));
    double* field_2_ptr = reinterpret_cast<double*>(b.field_data_location(field_2));
    double* field_3_ptr = reinterpret_cast<double*>(b.field_data_location(field_3));

    for(size_t j=0; j<b.size(); ++j) {
      field_1_ptr[j] = one;
      field_2_ptr[j] = this_proc;
      field_3_ptr[j] = this_proc;
    }
  }

  std::vector<stk::mesh::FieldBase*> field_vector;
  field_vector.push_back(&field_1);
  field_vector.push_back(&field_2);
  field_vector.push_back(&field_3);

  stk::mesh::parallel_min(bulk_data, field_vector);

  //now go through the comm nodes and confirm that the field-data values
  //for field_1 are equal to 1.0, and values for field_2 and field_3 are equal to min-proc == 0.

  const double min_proc = 0;

  const stk::mesh::EntityCommListInfoVector& entity_comm_list = bulk_data.comm_list();
  for(size_t i=0; i<entity_comm_list.size(); ++i) {
    stk::mesh::Entity node = entity_comm_list[i].entity;

    if(!bulk_data.is_matching_rank(field_1, node)) continue;

    const double* field_1_ptr = bulk_data.field_data(field_1, node);
    const double* field_2_ptr = bulk_data.field_data(field_2, node);
    const double* field_3_ptr = bulk_data.field_data(field_3, node);
    stk::mesh::PairIterEntityComm entity_comm = bulk_data.entity_comm_sharing(bulk_data.entity_key(node));
    if (entity_comm.size() > 0) {
      const double num_sharing_procs = entity_comm.size() + 1;
      std::cout<<"min: proc "<<this_proc<<", node "<<bulk_data.identifier(node)<<", num-sharing-procs: "<<num_sharing_procs<<std::endl;
      STKUNIT_ASSERT_EQUAL(field_1_ptr[0], one);
      STKUNIT_ASSERT_EQUAL(field_2_ptr[0], min_proc);
      STKUNIT_ASSERT_EQUAL(field_3_ptr[0], min_proc);
    }
    else {
      stk::mesh::PairIterEntityComm ent_comm = bulk_data.entity_comm(bulk_data.entity_key(node));
      const int num_procs = ent_comm.size()+1;
      std::cout<<"min: proc "<<this_proc<<" skipping non-shared ghost node "<<bulk_data.identifier(node)<<", num-procs: "<<num_procs<<std::endl;
    }
  }
}
} //namespace <anonymous>

