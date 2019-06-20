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
#include <stk_unit_test_utils/MeshFixture.hpp>
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

class ParallelHowTo : public stk::unit_test_util::MeshFixture {};

//BEGINCommuniateFieldData
TEST_F(ParallelHowTo, communicateFieldDataForSharedAndAura)
{
    auto& field = get_meta().declare_field<stk::mesh::Field<double>>(stk::topology::NODE_RANK, "temperature");

    double initialValue = 25.0;
    stk::mesh::put_field_on_entire_mesh_with_initial_value(field, &initialValue);

    setup_mesh("generated:8x8x8", stk::mesh::BulkData::AUTO_AURA);

    const stk::mesh::BucketVector& notOwnedBuckets = get_bulk().get_buckets(stk::topology::NODE_RANK,
                                                                            !get_meta().locally_owned_part());

    for(const stk::mesh::Bucket *bucket : notOwnedBuckets)
        for(stk::mesh::Entity node : *bucket)
            *stk::mesh::field_data(field, node) = -1.2345;

    stk::mesh::communicate_field_data(get_bulk(), {&field});

    for(const stk::mesh::Bucket *bucket : notOwnedBuckets)
        for(stk::mesh::Entity node : *bucket)
            EXPECT_EQ(initialValue, *stk::mesh::field_data(field, node));
}
//ENDCommuniateFieldData

//BEGINSum
void expect_field_has_value(const stk::mesh::BucketVector& buckets,
                            const stk::mesh::Field<double> &field,
                            double value)
{
    for(const stk::mesh::Bucket *bucket : buckets)
        for(stk::mesh::Entity node : *bucket)
            EXPECT_EQ(value, *stk::mesh::field_data(field, node));
}

TEST_F(ParallelHowTo, computeParallelSum)
{
    auto& field = get_meta().declare_field<stk::mesh::Field<double>>(stk::topology::NODE_RANK, "temperature");

    double initialValue = 25.0;
    stk::mesh::put_field_on_entire_mesh_with_initial_value(field, &initialValue);

    setup_mesh("generated:8x8x8", stk::mesh::BulkData::AUTO_AURA);

    const stk::mesh::BucketVector& shared = get_bulk().get_buckets(stk::topology::NODE_RANK,
                                                                   get_meta().globally_shared_part());
    const stk::mesh::BucketVector& notShared = get_bulk().get_buckets(stk::topology::NODE_RANK,
                                                                      !get_meta().globally_shared_part());
    expect_field_has_value(shared, field, initialValue);
    expect_field_has_value(notShared, field, initialValue);

    stk::mesh::parallel_sum(get_bulk(), {&field});

    expect_field_has_value(shared, field, 2*initialValue);
    expect_field_has_value(notShared, field, initialValue);
}
//ENDSum
