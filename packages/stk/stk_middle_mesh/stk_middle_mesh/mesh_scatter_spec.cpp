#include "mesh_scatter_spec.hpp"

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {  
  
  
void MeshScatterSpec::add_destination(mesh::MeshEntityPtr entity, int destRank)
{
  if(nullptr != entity) {
    MapKey key(entity->get_id(), entity->get_type());

    auto iter = m_entityProcMap.find(key);
    if(iter != m_entityProcMap.end()) {
      stk::util::insert_keep_sorted_and_unique(destRank, iter->second);
    } else {
      m_entityProcMap[key].push_back(destRank);
    }
  }
}

void MeshScatterSpec::get_destinations(mesh::MeshEntityPtr entity, std::vector<int>& destRanks)
{
  if(nullptr != entity) {
    MapKey key(entity->get_id(), entity->get_type());

    auto iter = m_entityProcMap.find(key);
    if(iter != m_entityProcMap.end()) {
      destRanks.reserve(iter->second.size());
      for(int dest : iter->second) {
        destRanks.push_back(dest);
      }
    }
  }
}

}
}
}
}
