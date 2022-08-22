
#include "MeshModification.hpp"
#include <stk_mesh/base/BulkData.hpp>


namespace stk {
namespace mesh {
namespace impl {

void DeletedEntityCache::mark_entity_as_deleted(Entity entity, bool is_ghost)
{
  if (is_ghost) 
  {
    m_ghost_reuse_map[m_bulkData.entity_key(entity)] = entity.local_offset();
  } else 
  {
    m_deleted_entities_current_modification_cycle.push_back(entity.local_offset());
  }
}

Entity::entity_value_type DeletedEntityCache::get_entity_for_reuse()
{
  if (!m_deleted_entities.empty())
  {
    size_t new_local_offset = m_deleted_entities.back();
    m_deleted_entities.pop_back();
    return new_local_offset;
  } else
  {
    return Entity::InvalidEntity;
  }
}

void DeletedEntityCache::update_deleted_entities_container()
{
  m_deleted_entities.insert(m_deleted_entities.end(), m_deleted_entities_current_modification_cycle.begin(), 
                                                      m_deleted_entities_current_modification_cycle.end());
  m_deleted_entities_current_modification_cycle.clear();

  for (auto keyAndOffset : m_ghost_reuse_map) {
    m_deleted_entities.push_back(keyAndOffset.second);
  }
  m_ghost_reuse_map.clear();
}


}
}
}