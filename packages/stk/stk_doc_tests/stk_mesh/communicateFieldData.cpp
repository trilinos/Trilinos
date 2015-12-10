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
#include <stk_mesh/base/FieldParallel.hpp>  // for communicate_field_data
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine
#include <string>                       // for string
#include <vector>                       // for vector
#include "mpi.h"                        // for MPI_COMM_WORLD, etc
#include "stk_io/DatabasePurpose.hpp"   // for DatabasePurpose::READ_MESH
#include "stk_io/StkMeshIoBroker.hpp"   // for StkMeshIoBroker
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/BulkDataInlinedMethods.hpp"
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/Field.hpp"      // for Field
#include "stk_mesh/base/FieldBase.hpp"  // for field_data, etc
#include "stk_mesh/base/MetaData.hpp"   // for MetaData, etc
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_mesh/base/Selector.hpp"   // for operator!, Selector
#include "stk_mesh/base/Types.hpp"      // for BucketVector
#include "stk_topology/topology.hpp"    // for topology, etc

//BEGIN
TEST(CommunicateFieldData, CommunicateMultipleGhostings)
{
  stk::ParallelMachine communicator = MPI_COMM_WORLD;

  const std::string fileName = "generated:8x8x8";
  stk::io::StkMeshIoBroker meshReader(communicator);
  meshReader.add_mesh_database(fileName, stk::io::READ_MESH);
  meshReader.create_input_mesh();

  stk::mesh::MetaData& meta = meshReader.meta_data();
  typedef stk::mesh::Field<double> ScalarField;
  ScalarField& temperatureField = meta.declare_field<ScalarField>(stk::topology::NODE_RANK, "temperature");

  double initialTemperatureValue = 25.0;
  stk::mesh::put_field_on_entire_mesh_with_initial_value(temperatureField, &initialTemperatureValue);

  meshReader.populate_bulk_data();

  stk::mesh::BulkData& bulk = meshReader.bulk_data();

  stk::mesh::Selector select_not_owned = !meta.locally_owned_part();
  const stk::mesh::BucketVector& buckets_not_owned = bulk.get_buckets(stk::topology::NODE_RANK,select_not_owned);

  for(size_t i=0; i<buckets_not_owned.size(); ++i) {
    stk::mesh::Bucket& bucket = *buckets_not_owned[i];
    for(size_t j=0; j<bucket.size(); ++j) {
        stk::mesh::Entity node = bucket[j];
        double* data = stk::mesh::field_data(temperatureField, node);
        double garbage = -1.2345;
        *data = garbage;
    }
  }

  std::vector<const stk::mesh::FieldBase*> fields(1, &temperatureField);
  stk::mesh::communicate_field_data(bulk, fields);

  for(size_t i=0; i<buckets_not_owned.size(); ++i) {
      stk::mesh::Bucket& bucket = *buckets_not_owned[i];
      for(size_t j=0; j<bucket.size(); ++j) {
          stk::mesh::Entity node = bucket[j];
          double* data = stk::mesh::field_data(temperatureField, node);
          EXPECT_EQ(initialTemperatureValue, *data);
      }
  }
}
//END
