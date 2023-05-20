#include "entity_sorted_by_owner.hpp"

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

EntitySortedByOwner::EntitySortedByOwner(int numOwnerRanks) :
  m_ownerLocalIds(numOwnerRanks),
  m_localEntities(numOwnerRanks)
{}

void EntitySortedByOwner::insert(const RemoteSharedEntity& owner, MeshEntityPtr localEntity)
{
  auto& ownerLocalIds = m_ownerLocalIds[owner.remoteRank];
  auto& localEntities = m_localEntities[owner.remoteRank];

  if (ownerLocalIds.size() > 0)
  {
    assert(ownerLocalIds.size() > 0 && owner.remoteId > ownerLocalIds.back());
  }
  
  ownerLocalIds.push_back(owner.remoteId);
  localEntities.push_back(localEntity);
}

// returns the entity associated with owner, or nullptr if no such entity;
MeshEntityPtr EntitySortedByOwner::get_entity(const RemoteSharedEntity& owner)
{
  auto& ownerLocalIds = m_ownerLocalIds[owner.remoteRank];
  auto& localEntities = m_localEntities[owner.remoteRank];

  auto it = std::lower_bound(ownerLocalIds.begin(), ownerLocalIds.end(), owner.remoteId);
  MeshEntityPtr entity = nullptr;
  if (it != ownerLocalIds.end() && *it == owner.remoteId)
    entity = localEntities[std::distance(ownerLocalIds.begin(), it)];

  return entity;
}

}
}
}
}
