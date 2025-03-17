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

#include <gtest/gtest.h>                // for TEST
#include <stddef.h>                     // for size_t
#include <set>                          // for set
#include <stk_io/FillMesh.hpp>
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, etc
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/Selector.hpp>   // for Selector
#include <stk_mesh/base/ForEachEntity.hpp>
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
TEST(StkMeshHowTo, iterateSidesetNodes_BucketLoop_ContiguousFieldDataWithinBucket)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  if (stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::unique_ptr<stk::mesh::BulkData> stkMesh = stk::mesh::MeshBuilder(comm)
                                                   .set_spatial_dimension(3)
                                                   .create();
  stk::mesh::MetaData &stkMeshMeta = stkMesh->mesh_meta_data();
  stk::mesh::Field<double> &temperatureField = stkMeshMeta.declare_field<double>(stk::topology::NODE_RANK, "temperature");
  stk::mesh::put_field_on_entire_mesh(temperatureField);

  // syntax creates faces for the surface on the positive 'x-side' of the 2x2x2 cube,
  // this part is given the name 'surface_1' when it is created.
  const std::string generatedMeshSpecification = "generated:2x2x2|sideset:X";
  stk::io::fill_mesh(generatedMeshSpecification, *stkMesh);

  stk::mesh::Part &boundaryConditionPart = *stkMeshMeta.get_part("surface_1");
  stk::mesh::Selector boundaryNodesSelector(boundaryConditionPart);

  const stk::mesh::BucketVector &boundaryNodeBuckets = stkMesh->get_buckets(stk::topology::NODE_RANK, boundaryNodesSelector);

  constexpr double prescribedTemperatureValue = 2.0;
  for (size_t bucketIndex = 0; bucketIndex < boundaryNodeBuckets.size(); ++bucketIndex) {
    const stk::mesh::Bucket &nodeBucket = *boundaryNodeBuckets[bucketIndex];
    double *temperatureValues = stk::mesh::field_data(temperatureField, nodeBucket);

    for (size_t nodeIndex = 0; nodeIndex < nodeBucket.size(); ++nodeIndex) {
      temperatureValues[nodeIndex] = prescribedTemperatureValue;
    }
  }

  testUtils::testTemperatureFieldSetCorrectly(temperatureField, boundaryNodesSelector, prescribedTemperatureValue);
}

TEST(StkMeshHowTo, iterateSidesetNodes_ForEachEntity_FieldDataAccess)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  if (stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::unique_ptr<stk::mesh::BulkData> stkMesh = stk::mesh::MeshBuilder(comm)
                                                   .set_spatial_dimension(3)
                                                   .create();
  stk::mesh::MetaData &stkMeshMeta = stkMesh->mesh_meta_data();
  stk::mesh::Field<double> &temperatureField = stkMeshMeta.declare_field<double>(stk::topology::NODE_RANK, "temperature");
  stk::mesh::put_field_on_entire_mesh(temperatureField);

  // syntax creates faces for the surface on the positive 'x-side' of the 2x2x2 cube,
  // this part is given the name 'surface_1' when it is created.
  const std::string generatedMeshSpecification = "generated:2x2x2|sideset:X";
  stk::io::fill_mesh(generatedMeshSpecification, *stkMesh);

  stk::mesh::Part &boundaryConditionPart = *stkMeshMeta.get_part("surface_1");
  stk::mesh::Selector boundaryNodesSelector(boundaryConditionPart);

  constexpr double prescribedTemperatureValue = 2.0;

  stk::mesh::for_each_entity_run(*stkMesh, stk::topology::NODE_RANK, boundaryNodesSelector,
    [&](const stk::mesh::BulkData& /*bulk*/, stk::mesh::Entity node) {
      double *temperatureValues = stk::mesh::field_data(temperatureField, node);
      *temperatureValues = prescribedTemperatureValue;
    });

  testUtils::testTemperatureFieldSetCorrectly(temperatureField, boundaryNodesSelector, prescribedTemperatureValue);
}
//-END
}
