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

#include <gtest/gtest.h>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/ForEachEntity.hpp>
#include <stk_topology/topology.hpp>
#include <string>
#include <vector>
#include "stk_io/FillMesh.hpp"

namespace {

//BEGIN_VISIT_EDGE_NODES
TEST(StkMeshHowTo, VisitEdgeNodes)
{
  // ============================================================
  // INITIALIZATION
  MPI_Comm comm = MPI_COMM_WORLD;
  if (stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::shared_ptr<stk::mesh::BulkData> bulk = stk::mesh::MeshBuilder(comm).create();
  stk::mesh::MetaData& meta = bulk->mesh_meta_data();

  const std::string generatedFileName = "generated:1x1x1";
  stk::io::fill_mesh(generatedFileName, *bulk);

  stk::mesh::EntityVector edgeNodes(2);
  unsigned edgeCount = 0;

  stk::mesh::for_each_entity_run(*bulk, stk::topology::ELEM_RANK, meta.locally_owned_part(),
  [&](const stk::mesh::BulkData& mesh, stk::mesh::Entity elem) {
    stk::topology elemTopo = mesh.bucket(elem).topology();
    const unsigned numEdgesPerElem = elemTopo.num_edges();
    const stk::mesh::Entity* elemNodes = mesh.begin(elem, stk::topology::NODE_RANK);
    for(unsigned i=0; i<numEdgesPerElem; ++i) {
      elemTopo.edge_nodes(&elemNodes[0], i, edgeNodes.data());
      std::cout<<"elem "<<mesh.identifier(elem)<<" edge "<<i<<", nodes: ";
      for(stk::mesh::Entity node : edgeNodes) {
        std::cout<<mesh.identifier(node)<<" ";
      }
      std::cout<<std::endl;
    }
    edgeCount += numEdgesPerElem;
  });

  EXPECT_EQ(12u, edgeCount);
}
//END_VISIT_EDGE_NODES
//
}
