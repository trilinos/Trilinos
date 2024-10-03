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
#ifndef _DISCONNECT_TYPES_HPP_
#define _DISCONNECT_TYPES_HPP_

#include <utility>
#include <vector>
#include "stk_util/util/ReportHandler.hpp"

namespace stk { namespace mesh { class Part; } }

namespace stk {
namespace tools {

struct BlockPairIdGetter {

  unsigned operator()(const stk::mesh::Part* part) {
    STK_ThrowRequire(part != nullptr);
    return part->mesh_meta_data_ordinal();
  }
};

struct BlockPair {
  typedef stk::mesh::Part* UNIT;

  BlockPair() : first(nullptr), second(nullptr) { }

  BlockPair(stk::mesh::Part* first_, stk::mesh::Part* second_)
    : first(first_), second(second_)
  { }

  BlockPair(const BlockPair& rhs)
  {
    first = rhs.first;
    second = rhs.second;
  }

  bool operator!=(const BlockPair& rhs) const
  {
    return first->mesh_meta_data_ordinal() != rhs.first->mesh_meta_data_ordinal() ||
        second->mesh_meta_data_ordinal() != rhs.second->mesh_meta_data_ordinal();
  }

  bool operator==(const BlockPair& rhs) const
  {
    return first->mesh_meta_data_ordinal() == rhs.first->mesh_meta_data_ordinal() &&
        second->mesh_meta_data_ordinal() == rhs.second->mesh_meta_data_ordinal();
  }

  bool operator<(const BlockPair& rhs) const
  {
    if(first->mesh_meta_data_ordinal() < rhs.first->mesh_meta_data_ordinal()) {
      return true;
    }
    else if(first->mesh_meta_data_ordinal() == rhs.first->mesh_meta_data_ordinal()) {
      return second->mesh_meta_data_ordinal() < rhs.second->mesh_meta_data_ordinal();
    }
    return false;
  }

  const stk::mesh::Part* get_first() const { return first; }
  const stk::mesh::Part* get_second() const { return second; }
  bool is_adjacent() const { return true; }
  bool is_valid() const {return (first != nullptr) && (second != nullptr) && (first != second); }

  stk::mesh::Part* first;
  stk::mesh::Part* second;
};

using BlockPairVector = std::vector<BlockPair>;
using BlockNamePair = std::pair<std::string, std::string>;
using BlockNamePairVector = std::vector<BlockNamePair>;

}}

#endif
