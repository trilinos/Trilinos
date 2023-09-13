#ifndef STK_MIDDLE_MESH_ENTITY_SORTED_BY_OWNER
#define STK_MIDDLE_MESH_ENTITY_SORTED_BY_OWNER

#include "mesh.hpp"

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

template <typename T>
class ValuesSortedByOwner
{
  public:
    explicit ValuesSortedByOwner(int numOwnerRanks);

    void insert(const RemoteSharedEntity& owner, const T& val);

    // returns the value associated with owner, or defaultVal if no such owner
    const T& get_value(const RemoteSharedEntity& owner, const T& defaultVal = T{});

		const std::vector<T>& get_values(int rank) const { return m_values[rank]; }

  private:
    std::vector<std::vector<int>> m_ownerLocalIds;
    std::vector<std::vector<T>> m_values;
};


template <typename T>
ValuesSortedByOwner<T>::ValuesSortedByOwner(int numOwnerRanks) :
  m_ownerLocalIds(numOwnerRanks),
  m_values(numOwnerRanks)
{}

template <typename T>
void ValuesSortedByOwner<T>::insert(const RemoteSharedEntity& owner, const T& val)
{
  auto& ownerLocalIds = m_ownerLocalIds[owner.remoteRank];
  auto& localEntities = m_values[owner.remoteRank];

  if (ownerLocalIds.size() > 0)
  {
    assert(ownerLocalIds.size() > 0 && owner.remoteId > ownerLocalIds.back());
  }
  
  ownerLocalIds.push_back(owner.remoteId);
  localEntities.push_back(val);
}

template <typename T>
const T& ValuesSortedByOwner<T>::get_value(const RemoteSharedEntity& owner, const T& defaultVal)
{
  auto& ownerLocalIds = m_ownerLocalIds[owner.remoteRank];
  auto& localEntities = m_values[owner.remoteRank];

  auto it = std::lower_bound(ownerLocalIds.begin(), ownerLocalIds.end(), owner.remoteId);
  if (it != ownerLocalIds.end() && *it == owner.remoteId)
    return localEntities[std::distance(ownerLocalIds.begin(), it)];

  return defaultVal;
}

using EntitySortedByOwner = ValuesSortedByOwner<MeshEntityPtr>;

}
}
}
}


#endif
