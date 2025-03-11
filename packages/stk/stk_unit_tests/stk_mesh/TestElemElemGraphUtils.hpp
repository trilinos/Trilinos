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

#ifndef stk_unit_tests_elem_elem_graph_utils_hpp
#define stk_unit_tests_elem_elem_graph_utils_hpp

#include "stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp"
#include "stk_mesh/base/ForEachEntity.hpp"

namespace stk {
namespace unit_test {

inline
void verify_no_graph_edges(const stk::mesh::BulkData& bulk)
{
  const stk::mesh::ElemElemGraph& eeGraph = bulk.get_face_adjacent_element_graph();
  stk::mesh::for_each_entity_run(bulk, stk::topology::ELEM_RANK, bulk.mesh_meta_data().locally_owned_part(),
  [&eeGraph](const stk::mesh::BulkData& /*mesh*/, stk::mesh::Entity elem) {
    stk::mesh::impl::LocalId elemLocalId = eeGraph.get_local_element_id(elem);
    stk::mesh::GraphEdgesForElement graphEdges = eeGraph.get_edges_for_element(elemLocalId);
    EXPECT_EQ(0u, graphEdges.size());
  });
}

inline
void verify_graph_edge_exists(const stk::mesh::ElemElemGraph& eeGraph,
                              stk::mesh::impl::LocalId elem1LocalId,
                              stk::mesh::impl::LocalId elem2LocalId)
{
  stk::mesh::GraphEdgesForElement graphEdges = eeGraph.get_edges_for_element(elem1LocalId);
  bool foundIt = false;
  for(const stk::mesh::GraphEdge& graphEdge : graphEdges) {
    if (graphEdge.elem2() == elem2LocalId) {
      foundIt = true;
      break;
    }
  }
  EXPECT_TRUE(foundIt);
}

inline
void verify_graph_edge_between_elems(const stk::mesh::BulkData& bulk,
                                     stk::mesh::Entity elem1, stk::mesh::Entity elem2)
{
  const stk::mesh::ElemElemGraph& eeGraph = bulk.get_face_adjacent_element_graph();
  stk::mesh::impl::LocalId elem1LocalId = eeGraph.get_local_element_id(elem1);
  stk::mesh::impl::LocalId elem2LocalId = eeGraph.get_local_element_id(elem2);

  verify_graph_edge_exists(eeGraph, elem1LocalId, elem2LocalId);
  verify_graph_edge_exists(eeGraph, elem2LocalId, elem1LocalId);
}

} // namespace unit_test
} // namespace stk

#endif
