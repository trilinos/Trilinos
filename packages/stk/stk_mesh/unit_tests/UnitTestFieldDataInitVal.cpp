/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <sstream>
#include <stdexcept>
#include <iostream>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/EntityCommDatabase.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>

using stk::mesh::Entity;
using stk::mesh::EntityRank;
using stk::mesh::Part;
using stk::mesh::Field;
using stk::mesh::BulkData;
using stk::mesh::EntityId;
using stk::mesh::MetaData;

namespace {

const EntityRank NODE_RANK = MetaData::NODE_RANK;

STKUNIT_UNIT_TEST(UnitTestFieldDataInitVal, test_scalar_field)
{
  // Test that if an initial-value is set on a scalar field, that value is
  // present the first time field-data is referenced for that field.
  //

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  MPI_Barrier( MPI_COMM_WORLD );

  // Set up meta and bulk data
  const unsigned spatial_dim = 2;
  MetaData meta_data(spatial_dim);

  const unsigned num_states = 1;
  Field<double>& dfield = meta_data.declare_field<Field<double> >(stk::topology::NODE_RANK, "double_scalar", num_states);

  const double initial_value = 99.9;

  stk::mesh::put_field(dfield, meta_data.universal_part(), &initial_value);

  meta_data.commit();

  BulkData mesh(meta_data, pm);
  unsigned p_rank = mesh.parallel_rank();

  // Begin modification cycle so we can create stuff
  mesh.modification_begin();

  // Node will be automatically added to the universal part
  stk::mesh::PartVector empty_parts;

  EntityId node_id = p_rank+1;
  // Create node
  Entity node = mesh.declare_entity(NODE_RANK, node_id, empty_parts);

  mesh.modification_end();

  //now insist that data for dfield on node is equal to the initial-value specified above:

  double* data_ptr = mesh.field_data( dfield, node);

  STKUNIT_ASSERT_EQUAL( *data_ptr, initial_value );
}

STKUNIT_UNIT_TEST(UnitTestFieldDataInitVal, test_vector_field)
{
  // Test that if an initial-value is set on a vector field, that value is
  // present the first time field-data is referenced for that field.
  //

  typedef stk::mesh::Field<double,stk::mesh::Cartesian2d> VectorField;

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  MPI_Barrier( MPI_COMM_WORLD );

  // Set up meta and bulk data
  const unsigned spatial_dim = 2;
  MetaData meta_data(spatial_dim);

  const unsigned num_states = 1;
  VectorField& vfield = meta_data.declare_field<VectorField>(stk::topology::NODE_RANK, "double_vector", num_states);

  const double initial_value[stk::mesh::Cartesian2d::Size] = { 50.0, 99.0 };

  stk::mesh::put_field(vfield, NODE_RANK, meta_data.universal_part(), stk::mesh::Cartesian2d::Size, initial_value);

  meta_data.commit();

  BulkData mesh(meta_data, pm);
  unsigned p_rank = mesh.parallel_rank();

  // Begin modification cycle so we can create stuff
  mesh.modification_begin();

  // Node will be automatically added to the universal part
  stk::mesh::PartVector empty_parts;

  EntityId node_id = p_rank+1;
  // Create node
  Entity node = mesh.declare_entity(NODE_RANK, node_id, empty_parts);

  mesh.modification_end();

  //now insist that data for vfield on node is equal to the initial-value specified above:

  double* data_ptr = mesh.field_data( vfield, node);

  STKUNIT_ASSERT_EQUAL( data_ptr[0], initial_value[0] );
  STKUNIT_ASSERT_EQUAL( data_ptr[1], initial_value[1] );
}

STKUNIT_UNIT_TEST(UnitTestFieldDataInitVal, test_vector_field_move_bucket)
{
  // Test that if an initial-value is set on a vector field, that value is
  // present the first time field-data is referenced for that field, and
  // that the value is present for a node that is moved from a bucket without
  // the field to a new bucket that does have the field.
  //

  typedef stk::mesh::Field<double,stk::mesh::Cartesian2d> VectorField;

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  MPI_Barrier( MPI_COMM_WORLD );

  // Set up meta and bulk data
  const unsigned spatial_dim = 2;
  MetaData meta_data(spatial_dim);

  const unsigned num_states = 1;
  VectorField& vfield = meta_data.declare_field<VectorField>(stk::topology::NODE_RANK, "double_vector", num_states);

  const double initial_value[stk::mesh::Cartesian2d::Size] = { 50.0, 99.0 };

  Part& node_part = meta_data.declare_part<shards::Node>("node_part");

  stk::mesh::put_field(vfield, NODE_RANK, node_part, stk::mesh::Cartesian2d::Size, initial_value);

  meta_data.commit();

  BulkData mesh(meta_data, pm);
  unsigned p_rank = mesh.parallel_rank();

  // Begin modification cycle so we can create stuff
  mesh.modification_begin();

  // Node will be automatically added to the universal part
  stk::mesh::PartVector empty_parts;

  EntityId node_id = p_rank+1;
  // Create node
  Entity node = mesh.declare_entity(NODE_RANK, node_id, empty_parts);

  // need to copy since bucket is going to be deleted during mesh modification
  const stk::mesh::PartVector old_parts = mesh.bucket(node).supersets();

  //Now move the node to the "node_part":
  stk::mesh::PartVector node_part_vec;
  node_part_vec.push_back(&node_part);
  mesh.change_entity_parts(node, node_part_vec);

  mesh.modification_end();

  const stk::mesh::PartVector& new_parts = mesh.bucket(node).supersets();

  //Insist that the node is now in a different bucket:
  bool in_different_bucket = old_parts != new_parts;
  STKUNIT_ASSERT_TRUE(in_different_bucket);

  //now insist that data for vfield on node is equal to the initial-value specified above:

  double* data_ptr = mesh.field_data( vfield, node);

  STKUNIT_ASSERT_EQUAL( data_ptr[0], initial_value[0] );
  STKUNIT_ASSERT_EQUAL( data_ptr[1], initial_value[1] );
}

STKUNIT_UNIT_TEST(UnitTestFieldDataInitVal, test_multi_state_vector_field)
{
  // Test that if an initial-value is set on a multi-state vector field, that value is
  // present the first time field-data is referenced for that field.
  //

  typedef stk::mesh::Field<double,stk::mesh::Cartesian2d> VectorField;

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  MPI_Barrier( MPI_COMM_WORLD );

  // Set up meta and bulk data
  const unsigned spatial_dim = 2;
  MetaData meta_data(spatial_dim);

  const unsigned num_states = 2;
  VectorField& vfield = meta_data.declare_field<VectorField>(stk::topology::NODE_RANK, "double_vector", num_states);

  const double initial_value[stk::mesh::Cartesian2d::Size] = { 50.0, 99.0 };

  stk::mesh::put_field(vfield, NODE_RANK, meta_data.universal_part(), stk::mesh::Cartesian2d::Size, initial_value);

  meta_data.commit();

  BulkData mesh(meta_data, pm);
  unsigned p_rank = mesh.parallel_rank();

  // Begin modification cycle so we can create stuff
  mesh.modification_begin();

  // Node will be automatically added to the universal part
  stk::mesh::PartVector empty_parts;

  EntityId node_id = p_rank+1;
  // Create node
  Entity node = mesh.declare_entity(NODE_RANK, node_id, empty_parts);

  mesh.modification_end();

  //now insist that data for vfield on node is equal to the initial-value specified above:

  STKUNIT_ASSERT_EQUAL( vfield.number_of_states(), num_states);

  VectorField& vfield_new = vfield.field_of_state(stk::mesh::StateNew);
  VectorField& vfield_old = vfield.field_of_state(stk::mesh::StateOld);

  {
    double* data_ptr_new = mesh.field_data( vfield_new, node);
    double* data_ptr_old = mesh.field_data( vfield_old, node);

    STKUNIT_ASSERT_EQUAL( data_ptr_new[0], initial_value[0] );
    STKUNIT_ASSERT_EQUAL( data_ptr_new[1], initial_value[1] );

    STKUNIT_ASSERT_EQUAL( data_ptr_old[0], initial_value[0] );
    STKUNIT_ASSERT_EQUAL( data_ptr_old[1], initial_value[1] );
  }
}

}
