#include "mesh_scatter_spec.hpp"

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {  
  
  
void MeshScatterSpec::add_destination(mesh::MeshEntityPtr entity, int destRank)
{
  auto& destRanksField = *m_destRanks;
  auto it = std::find(destRanksField(entity, 0).begin(), destRanksField(entity, 0).end(), destRank);
  if (it == destRanksField(entity, 0).end())
  {
    m_destRanks->insert(entity, 0, destRank);
  }

}

void MeshScatterSpec::get_destinations(mesh::MeshEntityPtr entity, std::vector<int>& destRanks)
{
  auto& destRanksField = *m_destRanks;
  destRanks.assign(destRanksField(entity, 0).begin(), destRanksField(entity, 0).end());
}


}
}
}
}
