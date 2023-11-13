// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
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
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
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


#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stddef.h>                     // for size_t
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>  // for declare_element
#include <stk_mesh/base/Field.hpp>      // for Field
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, entity_rank_names, etc
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/FieldBase.hpp"  // for field_data, etc
#include "stk_mesh/base/Types.hpp"      // for BucketVector, EntityId
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_io/IossBridge.hpp"
namespace stk { namespace mesh { class Part; } }

namespace {

namespace SpatialDimension { const unsigned three = 3; }

void create_two_tet_element_mesh(stk::mesh::BulkData &bulk)
{
  stk::mesh::MetaData &meta = bulk.mesh_meta_data();
  meta.use_simple_fields();
  stk::mesh::Part &tetPart = meta.declare_part_with_topology("tetElementPart", stk::topology::TET_4);
  meta.commit();

  bulk.modification_begin();
  stk::mesh::EntityId elem1Id = 1;
  stk::mesh::EntityIdVector elem1Nodes {1, 2, 3, 4};
  stk::mesh::declare_element(bulk, tetPart, elem1Id, elem1Nodes);
  stk::mesh::EntityId elem2Id = 2;
  stk::mesh::EntityIdVector elem2Nodes {2, 3, 4, 5};
  stk::mesh::declare_element(bulk, tetPart, elem2Id, elem2Nodes);
  bulk.modification_end();
}

//BEGINUseSimpleFields
TEST(stkMeshHowTo, useSimpleFields)
{
  stk::mesh::MeshBuilder builder(MPI_COMM_WORLD);
  builder.set_spatial_dimension(SpatialDimension::three);
  builder.set_entity_rank_names(stk::mesh::entity_rank_names());
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = builder.create();
  bulkPtr->mesh_meta_data().use_simple_fields();
  stk::mesh::MetaData& metaData = bulkPtr->mesh_meta_data();

  typedef stk::mesh::Field<double> DoubleField;
  DoubleField& pressureField = metaData.declare_field<double>(stk::topology::ELEM_RANK, "pressure");
  DoubleField& displacementsField = metaData.declare_field<double>(stk::topology::NODE_RANK, "displacements");

  double initialPressureValue = 4.4;
  stk::mesh::put_field_on_entire_mesh_with_initial_value(pressureField, &initialPressureValue);
  stk::mesh::put_field_on_mesh(displacementsField, metaData.universal_part(), 3, nullptr);
  stk::io::set_field_output_type(displacementsField, stk::io::FieldOutputType::VECTOR_3D);

  stk::mesh::BulkData& mesh = *bulkPtr;
  create_two_tet_element_mesh(mesh);

  const stk::mesh::BucketVector& nodeBuckets = mesh.buckets(stk::topology::NODE_RANK);
  EXPECT_TRUE(!nodeBuckets.empty());
  for(size_t bucketIndex=0; bucketIndex<nodeBuckets.size(); bucketIndex++)
  {
    const stk::mesh::Bucket& bucket = *nodeBuckets[bucketIndex];
    double* displacementDataForBucket = stk::mesh::field_data(displacementsField, bucket);
    EXPECT_GT(bucket.size(), 0u);
    for(size_t nodeIndex=0; nodeIndex<bucket.size(); nodeIndex++)
    {
      unsigned numValuesPerNode = stk::mesh::field_scalars_per_entity(displacementsField, bucket);
      EXPECT_EQ(SpatialDimension::three, numValuesPerNode);
      for(unsigned i=0; i<numValuesPerNode; i++)
      {
        EXPECT_EQ(0.0, displacementDataForBucket[nodeIndex*numValuesPerNode + i]);
        displacementDataForBucket[nodeIndex*numValuesPerNode + i] = 99.9;
      }
    }
  }

  stk::mesh::Entity elem1 = mesh.get_entity(stk::topology::ELEM_RANK, 1);
  double* pressureFieldDataForElem1 = stk::mesh::field_data(pressureField, elem1);
  EXPECT_EQ(initialPressureValue, *pressureFieldDataForElem1);

  stk::mesh::Entity elem2 = mesh.get_entity(stk::topology::ELEM_RANK, 2);
  double* pressureFieldDataForElem2 = stk::mesh::field_data(pressureField, elem2);
  EXPECT_EQ(initialPressureValue, *pressureFieldDataForElem2);
}
//ENDUseSimpleFields

void create_single_tet_element(stk::mesh::BulkData &bulk)
{
  stk::mesh::MetaData &meta = bulk.mesh_meta_data();
  meta.use_simple_fields();
  stk::mesh::Part &tetPart = meta.declare_part_with_topology("tetElementPart", stk::topology::TET_4);
  meta.commit();
  bulk.modification_begin();
  stk::mesh::EntityId elem1Id = 1;
  stk::mesh::EntityIdVector elem1Nodes {1, 2, 3, 4};
  stk::mesh::declare_element(bulk, tetPart, elem1Id, elem1Nodes);
  bulk.modification_end();
}

//BEGINBADFIELD
TEST(stkMeshHowTo, declareVectorFields_omitOutputType_noSubscriptNaming)
{
  stk::mesh::MeshBuilder builder(MPI_COMM_WORLD);
  builder.set_spatial_dimension(SpatialDimension::three);
  builder.set_entity_rank_names(stk::mesh::entity_rank_names());
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = builder.create();
  stk::mesh::MetaData& metaData = bulkPtr->mesh_meta_data();

  typedef stk::mesh::Field<double> DoubleField;
  DoubleField& velocities = metaData.declare_field<double>(stk::topology::NODE_RANK, "velocities");
  DoubleField& displacements = metaData.declare_field<double>(stk::topology::NODE_RANK, "displacements");

  unsigned fieldLength = 3;
  stk::mesh::put_field_on_mesh(velocities, metaData.universal_part(), fieldLength, nullptr);
  stk::mesh::put_field_on_mesh(displacements, metaData.universal_part(), fieldLength, nullptr);

  stk::io::set_field_output_type(velocities, stk::io::FieldOutputType::VECTOR_3D);

  stk::mesh::BulkData& mesh = *bulkPtr;
  create_single_tet_element(mesh);

  stk::mesh::Entity node1 = mesh.get_entity(stk::topology::NODE_RANK, 1);
  EXPECT_EQ(stk::mesh::field_scalars_per_entity(velocities, node1),
            stk::mesh::field_scalars_per_entity(displacements, node1));
}
//ENDBADFIELD

}
