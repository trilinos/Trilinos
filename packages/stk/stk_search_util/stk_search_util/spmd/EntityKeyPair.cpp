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

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include "stk_search_util/spmd/EntityKeyPair.hpp"
#include "stk_mesh/base/BulkData.hpp"       // for BulkData
#include "stk_util/util/ReportHandler.hpp"
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk {
namespace search {
namespace spmd {

std::ostream & operator << ( std::ostream &os , const EntityKeyPair &key )
{
  os << key.m_key;
  return os;
}

EntityKeyPair make_entity_key_pair(const stk::mesh::BulkData& bulk, stk::mesh::EntityKey key)
{
  return EntityKeyPair(bulk.get_entity(key), key);
}

EntityKeyPair make_entity_key_pair(const stk::mesh::BulkData& bulk, stk::mesh::Entity entity)
{
  return EntityKeyPair(entity, bulk.entity_key(entity));
}

EntityKeyPair make_entity_key_pair(const std::shared_ptr<stk::mesh::BulkData>& bulk, stk::mesh::EntityKey key)
{
  STK_ThrowRequire(bulk);
  return EntityKeyPair(bulk->get_entity(key), key);
}

EntityKeyPair make_entity_key_pair(const std::shared_ptr<stk::mesh::BulkData>& bulk, stk::mesh::Entity entity)
{
  STK_ThrowRequire(bulk);
  return EntityKeyPair(entity, bulk->entity_key(entity));
}

EntityKeyPair make_entity_key_pair(const stk::mesh::BulkData* bulk, stk::mesh::EntityKey key)
{
  STK_ThrowRequire(bulk != nullptr);
  return EntityKeyPair(bulk->get_entity(key), key);
}

EntityKeyPair make_entity_key_pair(const stk::mesh::BulkData* bulk, stk::mesh::Entity entity)
{
  STK_ThrowRequire(bulk != nullptr);
  return EntityKeyPair(entity, bulk->entity_key(entity));
}

stk::mesh::EntityVector convert_keys_to_entities(const stk::mesh::BulkData &bulk, const std::vector<EntityKeyPair>& keyPairs)
{
  stk::mesh::EntityVector entities(keyPairs.size());
  for (size_t i=0;i<entities.size();++i) {
    entities[i] = bulk.get_entity(keyPairs[i]);
  }
  return entities;
}

} // namespace spmd
} // namespace search
} // namespace stk
