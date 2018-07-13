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
#include <stk_topology/topology.hpp>    // for topology, etc
#include <string>                       // for string
#include "stk_io/DatabasePurpose.hpp"   // for DatabasePurpose::READ_MESH
#include "stk_mesh/base/Types.hpp"      // for BucketVector, PartVector
#include <stk_unit_test_utils/ioUtils.hpp>  // for fill_mesh_using_stk_io
#include "stkMeshTestUtils.hpp"

namespace stk { namespace mesh { class Part; } }

namespace
{

void verify_elem_is_owned_on_p0_and_valid_as_aura_on_p1(const stk::mesh::BulkData &bulkData, stk::mesh::Entity elem)
{
    EXPECT_TRUE(bulkData.is_valid(elem));
    EXPECT_EQ(0, bulkData.parallel_owner_rank(elem));
}

void verify_elem_is_now_owned_on_p1(const stk::mesh::BulkData &bulkData, stk::mesh::EntityId elemId)
{
    stk::mesh::Entity elem = bulkData.get_entity(stk::topology::ELEM_RANK, elemId);
    EXPECT_TRUE(bulkData.is_valid(elem));
    EXPECT_EQ(1, bulkData.parallel_owner_rank(elem));
}

//BEGIN
TEST(StkMeshHowTo, changeEntityOwner)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    if (stk::parallel_machine_size(communicator) == 2)
    {
        stk::mesh::MetaData metaData;
        stk::mesh::BulkData bulkData(metaData, communicator);
        stk::io::fill_mesh("generated:1x1x4", bulkData);

        stk::mesh::EntityId elem2Id = 2;
        stk::mesh::Entity elem2 = bulkData.get_entity(stk::topology::ELEM_RANK, elem2Id);
        verify_elem_is_owned_on_p0_and_valid_as_aura_on_p1(bulkData, elem2);

        std::vector<std::pair<stk::mesh::Entity, int> > elemProcPairs;
        if (bulkData.parallel_rank() == 0)
          elemProcPairs.push_back(std::make_pair(elem2, testUtils::get_other_proc(bulkData.parallel_rank())));

        bulkData.change_entity_owner(elemProcPairs);

        verify_elem_is_now_owned_on_p1(bulkData, elem2Id);
    }
}
//END
}
