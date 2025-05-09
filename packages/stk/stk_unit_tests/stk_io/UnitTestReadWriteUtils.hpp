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

#ifndef UNIT_TEST_READ_WRITE_UTILS_HPP
#define UNIT_TEST_READ_WRITE_UTILS_HPP

#include <stk_topology/topology.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/baseImpl/MeshImplUtils.hpp>
#include <stk_mesh/baseImpl/ConnectEdgesImpl.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_unit_test_utils/GetMeshSpec.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_io/FillMesh.hpp>
#include <stk_io/IossBridge.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_util/parallel/Parallel.hpp>

namespace io_test_utils {

struct ExpectedValues
{
  ExpectedValues()
    : numConnectedEdges(0),
      globalEdgeCount(0),
      globalElemCount(0)
  {}

  ExpectedValues(std::vector<unsigned> numLocalEdgesPerProc_, unsigned /*numFaces_*/,
                 unsigned numConnectedEdges_, unsigned globalEdgeCount_, unsigned globalElemCount_)
    : numLocalEdgesPerProc(numLocalEdgesPerProc_),
      numConnectedEdges(numConnectedEdges_),
      globalEdgeCount(globalEdgeCount_),
      globalElemCount(globalElemCount_)
  {}

  std::vector<unsigned> numEdgesPerProc;
  std::vector<unsigned> numLocalEdgesPerProc;
  std::vector<unsigned> numFacesPerProc;
  std::vector<unsigned> numLocalFacesPerProc;
  unsigned numConnectedEdges;
  unsigned globalEdgeCount;
  unsigned globalElemCount;
};

inline bool is_entity1_connected_to_entity2(const stk::mesh::BulkData& bulk, const stk::mesh::Entity entity1, const stk::mesh::Entity entity2)
{
  stk::mesh::EntityRank entityRank = bulk.entity_rank(entity2);

  unsigned numConnection = bulk.num_connectivity(entity1, entityRank);

  const stk::mesh::Entity* connectedEntities = bulk.begin(entity1, entityRank);
  for(unsigned i = 0; i < numConnection; i++) {
    if(connectedEntities[i] == entity2) {
      return true;
    }
  }
  return false;
}

inline bool is_fully_connected(const stk::mesh::BulkData& bulk, const stk::mesh::Entity entity1, const stk::mesh::Entity entity2)
{
  bool entity1IsConnectedToEntity2 = is_entity1_connected_to_entity2(bulk, entity1, entity2);
  bool entity2IsConnectedToEntity1 = is_entity1_connected_to_entity2(bulk, entity2, entity1);

  return entity1IsConnectedToEntity2 && entity2IsConnectedToEntity1;
}

}//io_test_utils

#endif
