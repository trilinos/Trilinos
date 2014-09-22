// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Governement retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

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

