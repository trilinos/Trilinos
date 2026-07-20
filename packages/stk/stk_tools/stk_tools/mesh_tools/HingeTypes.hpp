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

#ifndef _stk_tools_HingeTypes_hpp_
#define _stk_tools_HingeTypes_hpp_

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_tools/mesh_tools/PairwiseSideInfo.hpp>
#include <vector>
#include <utility>

namespace stk {
namespace tools {

class HingeNode {
public:
  HingeNode() { }
  HingeNode(const HingeNode& h) : node(h.node), info(h.info), isAHinge(h.isAHinge), isOwned(h.isOwned) {}
  HingeNode(stk::mesh::Entity node_, const impl::PairwiseSideInfoVector& info_) : node(node_), info(info_), isAHinge(true), isOwned(false) { }

  stk::mesh::Entity get_node() const { return node; }
  bool is_a_hinge() const { return isAHinge; }
  const impl::PairwiseSideInfoVector& get_info() const { return info; }
  bool is_owned() const { return isOwned; }

  void set_is_owned(bool owned) { isOwned = owned; }

  bool operator() (const HingeNode& h1, const HingeNode& h2) const { return h1.node < h2.node; }
  bool operator== (const stk::mesh::Entity entity) const { return node == entity; }
  bool operator== (const HingeNode& h) const { return node == h.node; }
  bool operator<  (const HingeNode& h) const { return node < h.node; }
  HingeNode& operator=  (const HingeNode& h) { node = h.node; info = h.info; isAHinge = h.isAHinge; isOwned = h.isOwned; return *this; }

private:
  stk::mesh::Entity node;
  impl::PairwiseSideInfoVector info;
  bool isAHinge = false;
  bool isOwned = false;
};

typedef std::vector<HingeNode> HingeNodeVector;
typedef std::pair<HingeNode, HingeNode> HingeEdge;
typedef std::vector<HingeEdge> HingeEdgeVector;
typedef stk::mesh::EntityVector HingeGroup;
typedef std::vector<HingeGroup> HingeGroupVector;

}} // namespace stk::tools

#endif
