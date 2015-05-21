// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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


#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stddef.h>                     // for size_t
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/CoordinateSystems.hpp>  // for Cartesian3d, etc
#include <stk_mesh/base/FEMHelpers.hpp>  // for declare_element
#include <stk_mesh/base/Field.hpp>      // for Field
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, entity_rank_names, etc
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/FieldBase.hpp"  // for field_data, etc
#include "stk_mesh/base/Types.hpp"      // for BucketVector, EntityId
#include "stk_topology/topology.hpp"    // for topology, etc
namespace stk { namespace mesh { class Part; } }

namespace {

//BEGIN
TEST(stkMeshHowTo, useSimpleFields)
{
    const unsigned spatialDimension = 3;
    stk::mesh::MetaData metaData(spatialDimension, stk::mesh::entity_rank_names());

    typedef stk::mesh::Field<double> ScalarField;
    typedef stk::mesh::Field<double, stk::mesh::Cartesian3d> VectorField;
    ScalarField& pressureField = metaData.declare_field<ScalarField>(stk::topology::ELEM_RANK, "pressure");
    VectorField& displacementsField = metaData.declare_field<VectorField>(stk::topology::NODE_RANK, "displacements");

    double initialPressureValue = 4.4;
    stk::mesh::put_field_on_entire_mesh_with_initial_value(pressureField, &initialPressureValue);
    stk::mesh::put_field_on_entire_mesh(displacementsField);

    stk::mesh::Part &tetPart = metaData.declare_part_with_topology("tetElementPart", stk::topology::TET_4);

    metaData.commit();
    stk::mesh::BulkData mesh(metaData, MPI_COMM_WORLD);
    mesh.modification_begin();
    stk::mesh::EntityId elem1Id = 1;
    stk::mesh::EntityIdVector elem1Nodes {1, 2, 3, 4};
    stk::mesh::Entity elem1=stk::mesh::declare_element(mesh, tetPart, elem1Id, elem1Nodes);
    stk::mesh::EntityId elem2Id = 2;
    stk::mesh::EntityIdVector elem2Nodes {2, 3, 4, 5};
    stk::mesh::Entity elem2=stk::mesh::declare_element(mesh, tetPart, elem2Id, elem2Nodes);
    mesh.modification_end();

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
            const unsigned cartesian_enum_size = stk::mesh::Cartesian3d::Size;
            EXPECT_EQ(cartesian_enum_size, numValuesPerNode);
            for(unsigned i=0; i<numValuesPerNode; i++)
            {
                EXPECT_EQ(0.0, displacementDataForBucket[nodeIndex*numValuesPerNode + i]);
                displacementDataForBucket[nodeIndex*numValuesPerNode + i] = 99.9;
            }
        }
    }

    double* pressureFieldDataForElem1 = stk::mesh::field_data(pressureField, elem1);
    EXPECT_EQ(initialPressureValue, *pressureFieldDataForElem1);

    double* pressureFieldDataForElem2 = stk::mesh::field_data(pressureField, elem2);
    EXPECT_EQ(initialPressureValue, *pressureFieldDataForElem2);
}
//END

}
