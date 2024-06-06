#include "CDTInterface.hpp"
#include "CDT.h"

namespace stk {
namespace middle_mesh {
namespace nonconformal4 {
namespace impl {

void CDTInterface::triangulate(const utils::impl::Projection& proj)
{
  assert(m_mesh->get_elements().size() == 0);

  CDT::Triangulation<double> tri(CDT::VertexInsertionOrder::Enum::AsProvided);

  // make vertices
  tri.insertVertices(
      m_mesh->get_vertices().begin(), m_mesh->get_vertices().end(),
      [&](mesh::MeshEntityPtr e) { return proj.project_plane_coords(e->get_point_orig(0)).get_x(); },
      [&](mesh::MeshEntityPtr e) { return proj.project_plane_coords(e->get_point_orig(0)).get_y(); });

  // the vertex indicies in CDT::Triangulation are the same as the indices
  // in mesh_in
  auto getStart = [](mesh::MeshEntityPtr e) -> CDT::VertInd {
    mesh::MeshEntityPtr v1 = e->get_down(0);
    mesh::MeshEntityPtr v2 = e->get_down(1);
    // CDT::Triangulator requires the first vertex to have the lower index
    if (v1->get_id() < v2->get_id())
      return v1->get_id();
    else
      return v2->get_id();
  };

  auto getEnd = [](mesh::MeshEntityPtr e) -> CDT::VertInd {
    mesh::MeshEntityPtr v1 = e->get_down(0);
    mesh::MeshEntityPtr v2 = e->get_down(1);
    // CDT::Triangulator requires the first vertex to have the lower index
    if (v1->get_id() < v2->get_id())
      return v2->get_id();
    else
      return v1->get_id();
  };

  tri.insertEdges(m_mesh->get_edges().begin(), m_mesh->get_edges().end(), getStart, getEnd);

  tri.eraseSuperTriangle();

  if (m_output)
  {
    mesh::impl::print_vert_edges(std::string("mesh_constraints"), m_mesh);
  }

  // write elements from tri back into mesh_in
  auto& verts = m_mesh->get_vertices();
  for (auto& triI : tri.triangles)
  {
    m_mesh->create_triangle_from_verts(verts[triI.vertices[0]], verts[triI.vertices[1]], verts[triI.vertices[2]]);
  }
}
} // namespace impl

} // namespace nonconformal4
} // namespace middle_mesh
} // namespace stk
