#include "mesh_layers.hpp"

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

void MeshLayers::initialize_que(std::vector<MeshEntityPtr>& roots, std::queue<MeshEntityPtr>& que)
{
  for (auto& v : roots)
  {
    que.push(v);
    mark_entity_seen(v);
  }
}

void MeshLayers::que_adjacent_verts(MeshEntityPtr v, std::queue<MeshEntityPtr>* que)
{
  for (int j = 0; j < v->count_up(); ++j)
  {
    auto edge            = v->get_up(j);
    MeshEntityPtr vOther = edge->get_down(0) == v ? edge->get_down(1) : edge->get_down(0);
    if (!is_entity_seen(vOther))
    {
      que->push(vOther);
      mark_entity_seen(vOther);
    }
  }
}

void MeshLayers::mark_entity_seen(MeshEntityPtr e)
{
  (*m_seenEntities)(e, 0, 0) = true;
}

bool MeshLayers::is_entity_seen(MeshEntityPtr e)
{
  return (*m_seenEntities)(e, 0, 0);
}

} // namespace impl
} // namespace mesh
} // namespace middle_mesh
} // namespace stk
