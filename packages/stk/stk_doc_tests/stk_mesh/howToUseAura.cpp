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
#include <sstream>                      // for ostringstream, etc
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/GetEntities.hpp>
#include "stk_mesh/base/Types.hpp"      // for BucketVector, PartVector
#include "stk_unit_test_utils/ioUtils.hpp"
namespace stk { namespace mesh { class Part; } }

namespace
{
//BEGIN_AURA_EXAMPLES
void expectNumElementsInAura(stk::mesh::BulkData::AutomaticAuraOption autoAuraOption,
                             unsigned numExpectedElementsInAura)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    if (stk::parallel_machine_size(communicator) == 2)
    {
        stk::mesh::MetaData meta(3);
        stk::mesh::BulkData bulk(meta, communicator, autoAuraOption);
        stk::unit_test_util::fill_mesh_using_stk_io("generated:1x1x2", bulk, communicator);

        EXPECT_EQ(numExpectedElementsInAura,
                  stk::mesh::count_selected_entities(meta.aura_part(), bulk.buckets(stk::topology::ELEMENT_RANK)));
    }
}
TEST(StkMeshHowTo, useNoAura)
{
    expectNumElementsInAura(stk::mesh::BulkData::NO_AUTO_AURA, 0);
}
TEST(StkMeshHowTo, useAutomaticGeneratedAura)
{
    expectNumElementsInAura(stk::mesh::BulkData::AUTO_AURA, 1);
}
TEST(StkMeshHowTo, useAuraDefaultBehavior)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    if (stk::parallel_machine_size(communicator) == 2)
    {
        stk::mesh::MetaData meta(3);
        stk::mesh::BulkData bulk(meta, communicator);
        stk::unit_test_util::fill_mesh_using_stk_io("generated:1x1x2", bulk, communicator);

        EXPECT_EQ(1u, stk::mesh::count_selected_entities(meta.aura_part(), bulk.buckets(stk::topology::ELEMENT_RANK)));
    }
}
//END_AURA_EXAMPLES
}
