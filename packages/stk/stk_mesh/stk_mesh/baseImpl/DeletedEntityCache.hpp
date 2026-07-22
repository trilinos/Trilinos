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

#ifndef stk_mesh_impl_DeletedEntityCache_hpp
#define stk_mesh_impl_DeletedEntityCache_hpp

#include <stk_mesh/base/Types.hpp>      // for MeshIndex, EntityRank, etc
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/EntityKey.hpp>
#include <stk_mesh/base/HashEntityAndEntityKey.hpp>
#include <unordered_map>

namespace stk {
namespace mesh {

class BulkData;
typedef std::unordered_map<EntityKey, Entity::entity_value_type, std::hash<EntityKey>> GhostReuseMap;


namespace impl {

class DeletedEntityCache
{
  public:
    explicit DeletedEntityCache(BulkData& bulkData) : 
      m_bulkData(bulkData)
    {}

    void mark_entity_as_deleted(Entity entity, bool is_ghost);

    const std::vector<Entity::entity_value_type>& get_deleted_entities_current_mod_cycle() const { return m_deleted_entities_current_modification_cycle; }

    GhostReuseMap& get_ghost_reuse_map() { return m_ghost_reuse_map; }

    const GhostReuseMap& get_ghost_reuse_map() const { return m_ghost_reuse_map; }

    Entity::entity_value_type get_entity_for_reuse();

    void update_deleted_entities_container();

  private:
    BulkData& m_bulkData;
    std::vector<Entity::entity_value_type> m_deleted_entities_current_modification_cycle;
    std::vector<Entity::entity_value_type> m_deleted_entities;
    GhostReuseMap m_ghost_reuse_map;
};

}
}
}

#endif
