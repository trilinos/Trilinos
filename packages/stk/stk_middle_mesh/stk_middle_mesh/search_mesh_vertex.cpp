#include "search_mesh_vertex.hpp"

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

void SearchMeshVertex::fill_bounding_boxes(std::vector<BoundingBox>& boundingBoxes) const
{
  int proc = utils::impl::comm_rank(m_unionComm); // Local search only
  boundingBoxes.reserve(m_mesh->get_vertices().size());
  for (auto& vert : m_mesh->get_vertices())
    if (vert)
    {
      Box box = create_bounding_box(vert);
      EntityProc entityProc(vert->get_id(), proc);
      boundingBoxes.push_back(BoundingBox(box, entityProc));
    }
}

//TODO: what we really want is a line segment class
SearchMeshVertex::Box SearchMeshVertex::create_bounding_box(MeshEntityPtr vert) const 
{
  utils::Point normal = (*m_normalField)(vert, 0, 0);
  utils::Point pt1 =  normal * m_boxNormalFac + vert->get_point_orig(0);
  utils::Point pt2 = -normal * m_boxNormalFac + vert->get_point_orig(0);

  Point minCorner, maxCorner;
  for (int d=0; d < 3; ++d)
  {
    minCorner[d] = std::min(pt1[d], pt2[d]);
    maxCorner[d] = std::max(pt1[d], pt2[d]);
  }

  return Box(minCorner, maxCorner);
}    

}
}
}
}