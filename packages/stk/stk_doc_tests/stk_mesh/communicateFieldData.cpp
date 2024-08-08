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
#include <stk_mesh/base/ForEachEntity.hpp>
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine
#include <string>                       // for string
#include <vector>                       // for vector
#include "stk_io/FillMesh.hpp"
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/Field.hpp"      // for Field
#include "stk_mesh/base/MetaData.hpp"   // for MetaData, etc
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_mesh/base/Selector.hpp"   // for operator!, Selector
#include "stk_mesh/base/Types.hpp"      // for BucketVector
#include "stk_topology/topology.hpp"    // for topology, etc

class ParallelHowTo : public stk::unit_test_util::MeshFixture {};

//BEGINCommuniateFieldData
TEST_F(ParallelHowTo, communicateFieldDataForSharedAndAura)
{
  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  auto& field = get_meta().declare_field<double>(stk::topology::NODE_RANK, "temperature");

  double initialValue = 25.0;
  stk::mesh::put_field_on_entire_mesh_with_initial_value(field, &initialValue);

  stk::io::fill_mesh("generated:8x8x8", get_bulk());

  stk::mesh::Selector notOwned = !get_meta().locally_owned_part();
  stk::mesh::for_each_entity_run(get_bulk(), stk::topology::NODE_RANK, notOwned,
    [&](const stk::mesh::BulkData& bulk, stk::mesh::Entity node) {
      *stk::mesh::field_data(field, node) = -1.2345;
    });

  stk::mesh::communicate_field_data(get_bulk(), {&field});

  stk::mesh::for_each_entity_run(get_bulk(), stk::topology::NODE_RANK, notOwned,
    [&](const stk::mesh::BulkData& bulk, stk::mesh::Entity node) {
      EXPECT_EQ(initialValue, *stk::mesh::field_data(field, node));
    });
}
//ENDCommuniateFieldData

//BEGINSum
void expect_nodal_field_has_value(const stk::mesh::Selector& selector,
                                  const stk::mesh::Field<double> &field,
                                  const double value)
{
  stk::mesh::for_each_entity_run(field.get_mesh(), stk::topology::NODE_RANK, selector,
    [&](const stk::mesh::BulkData& bulk, stk::mesh::Entity node) {
      EXPECT_EQ(value, *stk::mesh::field_data(field, node));
    });
}

TEST_F(ParallelHowTo, computeParallelSum)
{
  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  auto& field = get_meta().declare_field<double>(stk::topology::NODE_RANK, "temperature");

  double initialValue = 25.0;
  stk::mesh::put_field_on_entire_mesh_with_initial_value(field, &initialValue);

  stk::io::fill_mesh("generated:8x8x8", get_bulk());

  expect_nodal_field_has_value(get_meta().globally_shared_part(), field, initialValue);
  expect_nodal_field_has_value(!get_meta().globally_shared_part(), field, initialValue);

  stk::mesh::parallel_sum(get_bulk(), {&field});

  expect_nodal_field_has_value(get_meta().globally_shared_part(), field, 2*initialValue);
  expect_nodal_field_has_value(!get_meta().globally_shared_part(), field, initialValue);
}
//ENDSum
