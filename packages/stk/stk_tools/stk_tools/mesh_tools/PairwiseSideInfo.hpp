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

#ifndef _stk_tools_PairwiseSideInfo_hpp_
#define _stk_tools_PairwiseSideInfo_hpp_

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <vector>
#include <utility>

namespace stk {
namespace mesh { class BulkData; }
namespace tools {
namespace impl {

std::pair<stk::mesh::EntityVector,bool> get_pairwise_common_nodes(const stk::mesh::BulkData& bulk, stk::mesh::Entity elem1, stk::mesh::Entity elem2);

class PairwiseSideInfo {

public:
  PairwiseSideInfo(const stk::mesh::BulkData& bulk_, stk::mesh::Entity elem1_, stk::mesh::Entity elem2_)
    : bulk(&bulk_), elem1(elem1_), elem2(elem2_)
  {
    std::tie(commonNodes, hasAdjacentFace) = get_pairwise_common_nodes(*bulk, elem1, elem2);
  }

  stk::mesh::Entity get_element1() const { return elem1; }
  stk::mesh::Entity get_element2() const { return elem2; }
  const stk::mesh::EntityVector& get_common_nodes() { return commonNodes; }
  bool is_adjacent() const { return hasAdjacentFace; }
  void set_adjacency(bool adjacent) { hasAdjacentFace = adjacent; }

  const stk::mesh::BulkData& get_bulk() const {
    STK_ThrowRequire(nullptr != bulk);
    return *bulk;
  }

private:
  const stk::mesh::BulkData* bulk = nullptr;
  stk::mesh::Entity elem1;
  stk::mesh::Entity elem2;
  stk::mesh::EntityVector commonNodes;
  bool hasAdjacentFace = false;
};

typedef std::vector<PairwiseSideInfo> PairwiseSideInfoVector;

} // namespace impl

}} // namespace stk::tools

#endif
