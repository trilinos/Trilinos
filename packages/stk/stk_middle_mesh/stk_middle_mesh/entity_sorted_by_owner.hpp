#ifndef STK_MIDDLE_MESH_ENTITY_SORTED_BY_OWNER
#define STK_MIDDLE_MESH_ENTITY_SORTED_BY_OWNER

#include "mesh.hpp"

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

class EntitySortedByOwner
{
  public:
    explicit EntitySortedByOwner(int numOwnerRanks);

    void insert(const RemoteSharedEntity& owner, MeshEntityPtr localEntity);

    // returns the entity associated with owner, or nullptr if no such entity;
    MeshEntityPtr get_entity(const RemoteSharedEntity& owner);

		const std::vector<MeshEntityPtr>& get_entities(int rank) const { return m_localEntities[rank]; }

  private:
    std::vector<std::vector<int>> m_ownerLocalIds;
    std::vector<std::vector<MeshEntityPtr>> m_localEntities;
};

}
}
}
}


#endif
