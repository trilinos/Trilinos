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

#include <gtest/gtest.h>                // for TEST
#include <stddef.h>                     // for size_t
#include <set>                          // for set
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, etc
#include <stk_mesh/base/Selector.hpp>   // for Selector
#include <stk_topology/topology.hpp>    // for topology, etc
#include <string>                       // for string
#include "stkMeshTestUtils.hpp"
#include "stk_io/DatabasePurpose.hpp"   // for DatabasePurpose::READ_MESH
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/Field.hpp"      // for Field
#include "stk_mesh/base/FieldBase.hpp"  // for field_data
#include "stk_mesh/base/Types.hpp"      // for EntityId, EntityVector, etc
namespace stk { namespace mesh { class Part; } }

namespace
{
//-BEGIN
TEST(StkMeshHowTo, iterateSidesetNodesMostEfficientlyForFieldDataAccess)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    stk::io::StkMeshIoBroker stkMeshIoBroker(communicator);
    // syntax creates faces for the surface on the positive 'x-side' of the 2x2x2 cube,
    // this part is given the name 'surface_1' when it is created [create_input_mesh()]
    const std::string generatedMeshSpecification = "generated:2x2x2|sideset:X";
    stkMeshIoBroker.add_mesh_database(generatedMeshSpecification, stk::io::READ_MESH);
    stkMeshIoBroker.create_input_mesh();

    stk::mesh::MetaData &stkMeshMetaData = stkMeshIoBroker.meta_data();
    stk::mesh::Field<double> &temperatureField = stkMeshMetaData.declare_field<stk::mesh::Field<double> >(stk::topology::NODE_RANK, "temperature");
    stk::mesh::put_field_on_entire_mesh(temperatureField);
    stkMeshIoBroker.populate_bulk_data();

    stk::mesh::Part &boundaryConditionPart = *stkMeshMetaData.get_part("surface_1");
    stk::mesh::Selector boundaryNodesSelector(boundaryConditionPart);

    stk::mesh::BulkData &stkMeshBulkData = stkMeshIoBroker.bulk_data();
    const stk::mesh::BucketVector &boundaryNodeBuckets = stkMeshBulkData.get_buckets(stk::topology::NODE_RANK, boundaryNodesSelector);

    double prescribedTemperatureValue = 2.0;
    std::set<stk::mesh::EntityId> boundaryNodeIds;
    for (size_t bucketIndex = 0; bucketIndex < boundaryNodeBuckets.size(); ++bucketIndex)
    {
        stk::mesh::Bucket &nodeBucket = *boundaryNodeBuckets[bucketIndex];
        double *temperatureValues = stk::mesh::field_data(temperatureField, nodeBucket);
        for (size_t nodeIndex = 0; nodeIndex < nodeBucket.size(); ++nodeIndex)
        {
            stk::mesh::Entity node = nodeBucket[nodeIndex];
            boundaryNodeIds.insert(stkMeshBulkData.identifier(node));
            temperatureValues[nodeIndex] = prescribedTemperatureValue;
        }
    }

    testTemperatureFieldSetCorrectly(temperatureField, prescribedTemperatureValue, boundaryNodeIds);
}

TEST(StkMeshHowTo, iterateSidesetNodesWithFieldDataAccess)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    stk::io::StkMeshIoBroker stkMeshIoBroker(communicator);
    // syntax creates faces for the surface on the positive 'x-side' of the 2x2x2 cube,
    // this part is given the name 'surface_1' when it is created [create_input_mesh()]
    const std::string generatedMeshSpecification = "generated:2x2x2|sideset:X";
    stkMeshIoBroker.add_mesh_database(generatedMeshSpecification, stk::io::READ_MESH);
    stkMeshIoBroker.create_input_mesh();

    stk::mesh::MetaData &stkMeshMetaData = stkMeshIoBroker.meta_data();
    stk::mesh::Field<double> &temperatureField = stkMeshMetaData.declare_field<stk::mesh::Field<double> >(stk::topology::NODE_RANK, "temperature");
    stk::mesh::put_field_on_entire_mesh(temperatureField);
    stkMeshIoBroker.populate_bulk_data();

    stk::mesh::Part &boundaryConditionPart = *stkMeshMetaData.get_part("surface_1");
    stk::mesh::Selector boundaryNodesSelector(boundaryConditionPart);

    stk::mesh::BulkData &stkMeshBulkData = stkMeshIoBroker.bulk_data();

    stk::mesh::EntityVector nodes;
    stkMeshBulkData.get_selected_nodes(boundaryNodesSelector, nodes);

    double prescribedTemperatureValue = 2.0;
    std::set<stk::mesh::EntityId> boundaryNodeIds;

    for (size_t nodeIndex = 0; nodeIndex < nodes.size(); ++nodeIndex)
    {
        boundaryNodeIds.insert(stkMeshBulkData.identifier(nodes[nodeIndex]));
        double *temperatureValues = stk::mesh::field_data(temperatureField, nodes[nodeIndex]);
        *temperatureValues = prescribedTemperatureValue;
    }

    testTemperatureFieldSetCorrectly(temperatureField, prescribedTemperatureValue, boundaryNodeIds);
}
//-END
}
