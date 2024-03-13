#include "search_mesh_element_bounding_box_normal.hpp"

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

void SearchMeshElementBoundingBoxNormal::fill_bounding_boxes(std::vector<BoundingBox>& boundingBoxes) const
{
  int proc;
  MPI_Comm_rank(m_unionComm, &proc);
  auto entities = m_mesh->get_elements();

  for(auto entity : entities)
  {
    if (entity)
    {
      Box box = create_bounding_box(entity);

      EntityProc entityProc(entity->get_id(), proc);

      BoundingBox boundingBox(box, entityProc);
      boundingBoxes.push_back(boundingBox);
    }
  }

  std::sort(boundingBoxes.begin(), boundingBoxes.end(),
            [](const BoundingBox& a, const BoundingBox& b) { return a.second.id() < b.second.id(); });
}

SearchMeshElementBoundingBoxNormal::Box 
SearchMeshElementBoundingBoxNormal::create_bounding_box(mesh::MeshEntityPtr element) const
{
  stk::search::Point<double> minCorner, maxCorner;

  for(unsigned j = 0; j < 3u; ++j)
  {
    minCorner[j] = std::numeric_limits<double>::max();
    maxCorner[j] = -std::numeric_limits<double>::max();
  }    

  std::array<mesh::MeshEntityPtr, mesh::MAX_DOWN> verts;
  int nverts = mesh::get_downward(element, 0, verts.data());
  auto& normalField = *m_averagedNormalField;
  for (int i=0; i < nverts; ++i)
  {
    utils::Point vertCoords = verts[i]->get_point_orig(0);
    utils::Point normal = normalField(verts[i], 0, 0);

    utils::Point vertCoordsPlus = vertCoords + m_normalFac * normal;
    utils::Point vertCoordsMinus = vertCoords - m_normalFac * normal;

    for (int d=0; d < 3; ++d)
    {
      maxCorner[d] = std::max(maxCorner[d], vertCoordsPlus[d]);
      maxCorner[d] = std::max(maxCorner[d], vertCoordsMinus[d]);

      minCorner[d] = std::min(minCorner[d], vertCoordsPlus[d]);
      minCorner[d] = std::min(minCorner[d], vertCoordsMinus[d]); 
    }
  }

  return Box(minCorner, maxCorner);
}

}
}
}
}