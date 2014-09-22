/*------------------------------------------------------------------------*/
/*                 Copyright (c) 2013, Sandia Corporation.
/*                 Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
/*                 the U.S. Governement retains certain rights in this software.
/*                 
/*                 Redistribution and use in source and binary forms, with or without
/*                 modification, are permitted provided that the following conditions are
/*                 met:
/*                 
/*                     * Redistributions of source code must retain the above copyright
/*                       notice, this list of conditions and the following disclaimer.
/*                 
/*                     * Redistributions in binary form must reproduce the above
/*                       copyright notice, this list of conditions and the following
/*                       disclaimer in the documentation and/or other materials provided
/*                       with the distribution.
/*                 
/*                     * Neither the name of Sandia Corporation nor the names of its
/*                       contributors may be used to endorse or promote products derived
/*                       from this software without specific prior written permission.
/*                 
/*                 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
/*                 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
/*                 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
/*                 A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
/*                 OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
/*                 SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
/*                 LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
/*                 DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
/*                 THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
/*                 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
/*                 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
/*                 
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/CoordinateSystems.hpp>  // for Cartesian2d, etc
#include <stk_mesh/base/Field.hpp>      // for Field
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, put_field
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine
#include <gtest/gtest.h>
#include <vector>                       // for operator!=
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/FieldBase.hpp"  // for field_data
#include "stk_mesh/base/FieldState.hpp"  // for FieldState::StateNew, etc
#include "stk_mesh/base/Types.hpp"      // for PartVector, EntityId, etc
#include "stk_topology/topology.hpp"    // for topology, etc
namespace stk { namespace mesh { class Part; } }


using stk::mesh::Entity;
using stk::mesh::EntityRank;
using stk::mesh::Part;
using stk::mesh::Field;
using stk::mesh::BulkData;
using stk::mesh::EntityId;
using stk::mesh::MetaData;

namespace {

const EntityRank NODE_RANK = stk::topology::NODE_RANK;

TEST(UnitTestFieldDataInitVal, test_scalar_field)
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

  double* data_ptr = stk::mesh::field_data( dfield, node);

  ASSERT_EQ( *data_ptr, initial_value );
}

TEST(UnitTestFieldDataInitVal, test_vector_field)
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

  stk::mesh::put_field(vfield, meta_data.universal_part(), stk::mesh::Cartesian2d::Size, initial_value);

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

  double* data_ptr = stk::mesh::field_data( vfield, node);

  ASSERT_EQ( data_ptr[0], initial_value[0] );
  ASSERT_EQ( data_ptr[1], initial_value[1] );
}

TEST(UnitTestFieldDataInitVal, test_vector_field_move_bucket)
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

  Part& node_part = meta_data.declare_part_with_topology("node_part", stk::topology::NODE);

  stk::mesh::put_field(vfield, node_part, stk::mesh::Cartesian2d::Size, initial_value);

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
  ASSERT_TRUE(in_different_bucket);

  //now insist that data for vfield on node is equal to the initial-value specified above:

  double* data_ptr = stk::mesh::field_data( vfield, node);

  ASSERT_EQ( data_ptr[0], initial_value[0] );
  ASSERT_EQ( data_ptr[1], initial_value[1] );
}

TEST(UnitTestFieldDataInitVal, test_multi_state_vector_field)
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

  stk::mesh::put_field(vfield, meta_data.universal_part(), stk::mesh::Cartesian2d::Size, initial_value);

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

  ASSERT_EQ( vfield.number_of_states(), num_states);

  VectorField& vfield_new = vfield.field_of_state(stk::mesh::StateNew);
  VectorField& vfield_old = vfield.field_of_state(stk::mesh::StateOld);

  {
    double* data_ptr_new = stk::mesh::field_data( vfield_new, node);
    double* data_ptr_old = stk::mesh::field_data( vfield_old, node);

    ASSERT_EQ( data_ptr_new[0], initial_value[0] );
    ASSERT_EQ( data_ptr_new[1], initial_value[1] );

    ASSERT_EQ( data_ptr_old[0], initial_value[0] );
    ASSERT_EQ( data_ptr_old[1], initial_value[1] );
  }
}

}
