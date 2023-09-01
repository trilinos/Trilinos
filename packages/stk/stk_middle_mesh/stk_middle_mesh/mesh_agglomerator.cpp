#include "mesh_agglomerator.hpp"

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

std::vector<int> MeshAgglomerator::get_group_idxs(SetType<MeshEntityPtr>& vertsIn, const int nthres)
{
  std::vector<int> idxs;
  for (int i = 0; i < get_num_groups(); ++i)
  {
    int nfound   = 0;
    auto& vertsI = m_verts[i];
    for (auto& v : vertsIn)
    {
      if (vertsI.count(v) > 0)
      {
        nfound += 1;
      }

      if (nfound == nthres)
      {
        idxs.push_back(i);
        break;
      }
    }
  }

  return idxs;
}

std::ostream& operator<<(std::ostream& os, const MeshAgglomerator& agg)
{
  os << "MeshAgglomerator with " << agg.get_num_groups() << " groups:" << std::endl;
  for (int i = 0; i < agg.get_num_groups(); ++i)
  {
    os << "  group " << i << " verts: ";
    const auto& verts = agg.get_group_verts(i);
    for (auto& v : verts)
      os << v->get_id() << " ";
    os << std::endl;
  }

  return os;
}

} // namespace impl
} // namespace mesh
} // namespace middle_mesh
} // namespace stk
