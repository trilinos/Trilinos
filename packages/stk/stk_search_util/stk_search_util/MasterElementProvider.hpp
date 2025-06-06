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

#ifndef STK_STK_SEARCH_UTIL_STK_SEARCH_UTIL_MASTERELEMENTPROVIDER_HPP_
#define STK_STK_SEARCH_UTIL_STK_SEARCH_UTIL_MASTERELEMENTPROVIDER_HPP_

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include "stk_search/SearchInterface.hpp"  // for ProvideMasterElementInterface
#include "stk_search_util/SearchField.hpp"
#include "stk_search_util/spmd/EntityKeyPair.hpp"
#include <stk_mesh/base/EntityKey.hpp>     // for EntityKey
#include <stk_mesh/base/FieldBase.hpp>     // for FieldBase
#include "stk_topology/topology.hpp"       // for topology
#include <map>                             // for map, map<>::value_compare
#include <vector>                          // for vector

namespace stk { namespace mesh { class Bucket; } }
namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class MetaData; } }
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk {
namespace search {

class SearchTopology
{
public:
  SearchTopology(stk::topology topology,
                 spmd::EntityKeyPair key = spmd::EntityKeyPair(),
                 const stk::mesh::Bucket* bucket = nullptr)
  : m_topology(topology),
    m_key(key),
    m_bucket(bucket) {}

  stk::topology get_topology() const {return m_topology;}
  spmd::EntityKeyPair get_key() const {return m_key;}
  const stk::mesh::Bucket* get_bucket() const {return m_bucket;}

  operator stk::topology() const { return m_topology; }

private:
  stk::topology m_topology;
  spmd::EntityKeyPair m_key;
  const stk::mesh::Bucket* m_bucket{nullptr};

  SearchTopology() = delete;
};

using MasterElementProviderInterface =
    ProvideMasterElementInterface<SearchTopology, spmd::EntityKeyPair, SearchField>;


} // namespace stk
} // namespace search

#endif /* STK_STK_SEARCH_UTIL_STK_SEARCH_UTIL_MASTERELEMENTPROVIDER_HPP_ */
