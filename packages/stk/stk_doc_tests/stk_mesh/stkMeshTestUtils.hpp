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

#ifndef stkMeshTestUtilsHpp
#define stkMeshTestUtilsHpp

#include <gtest/gtest.h>
#include <stk_mesh/base/ForEachEntity.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_topology/topology.hpp>

namespace stk { namespace mesh { class BulkData; } }

namespace testUtils
{

inline
int get_other_proc(int myproc)
{
  return myproc == 0 ? 1 : 0;
}

inline
void testTemperatureFieldSetCorrectly(const stk::mesh::Field<double> &temperatureField,
                                      const stk::mesh::Selector& boundaryNodesSelector,
                                      double prescribedTemperatureValue)
{
  const stk::mesh::BulkData &stkMeshBulkData = temperatureField.get_mesh();

  stk::mesh::for_each_entity_run(stkMeshBulkData, stk::topology::NODE_RANK, boundaryNodesSelector,
    [&](const stk::mesh::BulkData& bulk, stk::mesh::Entity node) {
      const double *temperature = stk::mesh::field_data(temperatureField, node);
      EXPECT_EQ(prescribedTemperatureValue, *temperature);
    });

  stk::mesh::Selector nonBoundaryNodes = !boundaryNodesSelector;

  stk::mesh::for_each_entity_run(stkMeshBulkData, stk::topology::NODE_RANK, nonBoundaryNodes,
    [&](const stk::mesh::BulkData& bulk, stk::mesh::Entity node) {
      const double *temperature = stk::mesh::field_data(temperatureField, node);
      EXPECT_EQ(0.0, *temperature);
    });
}
}

#endif
