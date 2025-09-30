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

#ifndef STK_STK_SEARCH_UTIL_STK_SEARCH_UTIL_CACHEDENTITY_HPP_
#define STK_STK_SEARCH_UTIL_STK_SEARCH_UTIL_CACHEDENTITY_HPP_

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include "stk_mesh/base/EntityKey.hpp"                // for EntityKey
#include "stk_mesh/base/Entity.hpp"                   // for Entity
#include "stk_mesh/base/BulkData.hpp"                 // for BulkData
#include "stk_util/util/ReportHandler.hpp"
#include "stk_search_util/spmd/EntityKeyPair.hpp"
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk {
namespace search {

class CachedEntity {
 public:
  CachedEntity(const stk::mesh::BulkData* bulk)
  : m_bulk(bulk)
  { }

  ~CachedEntity() = default;

  stk::mesh::Entity get_entity(const stk::mesh::EntityKey& k) const
  {
    STK_ThrowAssertMsg(nullptr != m_bulk, "BulkData is null");
    stk::mesh::Entity e;

    if (k == m_cachedKey) {
      e = m_cachedEntity;
    }
    else {
      e = m_bulk->get_entity(k);
      m_cachedEntity = e;
      m_cachedKey = k;
    }

    return e;
  }

  stk::mesh::Entity get_entity(const stk::mesh::Entity& e) const
  {
    STK_ThrowAssertMsg(nullptr != m_bulk, "BulkData is null");

      m_cachedEntity = e;
      m_cachedKey = m_bulk->entity_key(e);

    return e;
  }

  stk::mesh::Entity get_entity(const stk::search::spmd::EntityKeyPair& key) const
  {
    STK_ThrowAssertMsg(nullptr != m_bulk, "BulkData is null");

      m_cachedEntity = key;
      m_cachedKey = key;

    return m_cachedEntity;
  }

  void set_bulk_data(const stk::mesh::BulkData* bulk)
  {
    STK_ThrowAssertMsg(nullptr != bulk, "Input BulkData is null");
    m_bulk = bulk;
  }

 private:
  const stk::mesh::BulkData* m_bulk{nullptr};

  mutable stk::mesh::Entity m_cachedEntity;
  mutable stk::mesh::EntityKey m_cachedKey;

  CachedEntity(const CachedEntity&) = delete;
  const CachedEntity& operator()(const CachedEntity&) = delete;
};

} // namespace search
} // namespace stk


#endif /* STK_STK_SEARCH_UTIL_STK_SEARCH_UTIL_CACHEDENTITY_HPP_ */
