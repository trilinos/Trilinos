#include "search_mesh_element_bounding_box.hpp"


namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {



void SearchMeshElementBoundingBox::fill_bounding_boxes(std::vector<BoundingBox>& boundingBoxes) const
{
  int proc;
  MPI_Comm_rank(m_unionComm, &proc);
  auto entities = m_mesh->get_elements();

  stk::search::Point<double> minCorner, maxCorner;

  for(auto entity : entities)
  {
    if (entity)
    {
      fill_bounding_box(entity, minCorner, maxCorner);

      EntityProc entityProc(entity->get_id(), proc);

      BoundingBox boundingBox(Box(minCorner, maxCorner), entityProc);
      boundingBoxes.push_back(boundingBox);
    }
  }

    std::sort(boundingBoxes.begin(), boundingBoxes.end(),
              [](const BoundingBox& a, const BoundingBox& b) { return a.second.id() < b.second.id(); });
  }


void SearchMeshElementBoundingBox::fill_bounding_box(mesh::MeshEntityPtr element,
                                                      stk::search::Point<double>& minCorner,
                                                      stk::search::Point<double>& maxCorner) const
{
  for(unsigned j = 0; j < 3u; ++j) {
    minCorner[j] = std::numeric_limits<double>::max();
    maxCorner[j] = -std::numeric_limits<double>::max();
  }

  int numEdges = element->count_down();
  for(auto i = 0; i < numEdges; ++i) {
    auto edge = element->get_down(i);

    int numNodes = edge->count_down();
    for(auto j = 0; j < numNodes; ++j) {
      auto node = edge->get_down(j);
      auto nodeCoords = node->get_point_orig(0);

      for(unsigned k = 0; k < 3u; ++k) {
        minCorner[k] = std::min(minCorner[k], nodeCoords[k]);
        maxCorner[k] = std::max(maxCorner[k], nodeCoords[k]);
      }
    }
  }
}
}
}
}
}