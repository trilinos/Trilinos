/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <boost/foreach.hpp>            // for auto_any_base, etc
#include <iostream>                     // for ostream, operator<<, etc
#include <stdexcept>                    // for runtime_error
#include <stk_mesh/base/CoordinateSystems.hpp>  // for Cartesian, etc
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, put_field
#include <gtest/gtest.h>
#include <string>                       // for operator==, string, etc
#include <vector>                       // for vector
#include "stk_mesh/base/Field.hpp"      // for Field
#include "stk_mesh/base/FieldBase.hpp"  // for field_bytes_per_entity, etc
#include "stk_mesh/base/Types.hpp"      // for PartVector, BucketVector, etc
#include "stk_topology/topology.hpp"    // for topology, etc

namespace {

TEST(UnitTestGetFieldByName, test1)
{
  size_t spatialDim = 3;
  stk::mesh::MetaData meta(spatialDim);

  //declare fields on different ranks with names that are unique within a rank but not unique overall:
  //
  stk::mesh::Field<double>& nodeField1 = meta.declare_field<stk::mesh::Field<double> >(stk::topology::NODE_RANK, "field1");
  stk::mesh::Field<double,stk::mesh::Cartesian>& nodeField2 = meta.declare_field<stk::mesh::Field<double,stk::mesh::Cartesian> >(stk::topology::NODE_RANK, "field2");

  stk::mesh::Field<double>& faceField1 = meta.declare_field<stk::mesh::Field<double> >(stk::topology::FACE_RANK, "field1");
  stk::mesh::Field<double,stk::mesh::Cartesian>& faceField2 = meta.declare_field<stk::mesh::Field<double,stk::mesh::Cartesian> >(stk::topology::FACE_RANK, "field2");

  stk::mesh::Field<double>& elemField1 = meta.declare_field<stk::mesh::Field<double> >(stk::topology::ELEM_RANK, "field1");
  stk::mesh::Field<double,stk::mesh::Cartesian>& elemField2 = meta.declare_field<stk::mesh::Field<double,stk::mesh::Cartesian> >(stk::topology::ELEM_RANK, "field2");

  stk::mesh::EntityRank side_rank = meta.side_rank();
  EXPECT_EQ(stk::topology::FACE_RANK, static_cast<stk::topology::rank_t>(side_rank));

  //test the get_field method:

  //node fields:
  stk::mesh::Field<double>* get_field_nodeField1 = meta.get_field<stk::mesh::Field<double> >(stk::topology::NODE_RANK, "field1");
  EXPECT_EQ(nodeField1.mesh_meta_data_ordinal(), get_field_nodeField1->mesh_meta_data_ordinal());

  stk::mesh::FieldBase* get_field_nodeField2 = meta.get_field(stk::topology::NODE_RANK, "field2");
  EXPECT_EQ(nodeField2.mesh_meta_data_ordinal(), get_field_nodeField2->mesh_meta_data_ordinal());

  //side/face fields:
  stk::mesh::Field<double>* get_field_sideField1 = meta.get_field<stk::mesh::Field<double> >(side_rank, "field1");
  EXPECT_EQ(faceField1.mesh_meta_data_ordinal(), get_field_sideField1->mesh_meta_data_ordinal());

  stk::mesh::FieldBase* get_field_sideField2 = meta.get_field(stk::topology::FACE_RANK, "field2");
  EXPECT_EQ(faceField2.mesh_meta_data_ordinal(), get_field_sideField2->mesh_meta_data_ordinal());

  //elem fields:
  stk::mesh::Field<double>* get_field_elemField1 = meta.get_field<stk::mesh::Field<double> >(stk::topology::ELEM_RANK, "field1");
  EXPECT_EQ(elemField1.mesh_meta_data_ordinal(), get_field_elemField1->mesh_meta_data_ordinal());

  stk::mesh::FieldBase* get_field_elemField2 = meta.get_field(stk::topology::ELEM_RANK, "field2");
  EXPECT_EQ(elemField2.mesh_meta_data_ordinal(), get_field_elemField2->mesh_meta_data_ordinal());
}
} //namespace <anonymous>

