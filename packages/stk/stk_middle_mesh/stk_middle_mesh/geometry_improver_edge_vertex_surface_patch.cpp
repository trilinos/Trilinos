#include "geometry_improver_edge_vertex_surface_patch.hpp"

namespace stk {
namespace middle_mesh {
namespace nonconformal4 {
namespace impl {

void GeometryImproverEdgeVertexSurfacePatch::run()
{
  reconstruct_edge_vertex_positions();
}

void GeometryImproverEdgeVertexSurfacePatch::setup_vert_status_field()
{
  auto& meshInVertsFromInputVert = *m_meshInVertsFromInputVert;

  auto& fakeVertsToVertsIn = m_relationalData->fakeVertsToVertsIn;
  auto& verts1ToFakeVerts  = *(m_relationalData->verts1ToFakeVerts);
  for (auto& vert1 : m_mesh1->get_vertices())
    if (vert1)
    {
      auto fv                                = verts1ToFakeVerts(vert1, 0, 0);
      mesh::MeshEntityPtr vertIn             = fakeVertsToVertsIn[fv.id];
      meshInVertsFromInputVert(vertIn, 0, 0) = true;
    }

  auto& verts2ToFakeVerts = *(m_relationalData->verts2ToFakeVerts);
  for (auto& vert2 : m_mesh2->get_vertices())
    if (vert2)
    {
      auto fv                                = verts2ToFakeVerts(vert2, 0, 0);
      mesh::MeshEntityPtr vertIn             = fakeVertsToVertsIn[fv.id];
      meshInVertsFromInputVert(vertIn, 0, 0) = true;
    }
}


utils::Point GeometryImproverEdgeVertexSurfacePatch::get_normal_vector(mesh::MeshEntityPtr edge2, double xi)
{
  auto& normalField    = *m_normalField;
  utils::Point normal1 = normalField(edge2->get_down(0), 0, 0);
  utils::Point normal2 = normalField(edge2->get_down(1), 0, 0);

  normal1 = normal1 / std::sqrt(dot(normal1, normal1));
  normal2 = normal2 / std::sqrt(dot(normal2, normal2));

  return normal1 + xi * (normal2 - normal1);
}

void GeometryImproverEdgeVertexSurfacePatch::reconstruct_edge_vertex_positions()
{
  auto& edges2ToFakeVertsIn = *(m_relationalData->edges2ToFakeVertsIn);
  auto& fakeVertsToVertsIn  = m_relationalData->fakeVertsToVertsIn;
  std::vector<mesh::MeshEntityPtr> closestVerts;
  std::vector<utils::Point> closestPts;
  for (auto& edge2 : m_mesh2->get_edges())
    if (edge2)
    {
      int numVertsIn = edges2ToFakeVertsIn.get_num_comp(edge2, 0);
      if (numVertsIn > 2)
      {
        for (int i = 1; i < numVertsIn - 1; ++i)
        {
          VertOnEdge& vertOnEdge     = edges2ToFakeVertsIn(edge2, 0, i);
          mesh::MeshEntityPtr vertIn = fakeVertsToVertsIn[vertOnEdge.vert.id];

          m_closestPointSearch.search(vertIn, m_numSurfacePatchPoints, closestVerts);

          utils::Point normal(0, 0, 1);
          if (m_normalField)
            normal = get_normal_vector(edge2, vertOnEdge.xi);
          utils::impl::ChangeOfBasis changeOfBasis(utils::impl::compute_basis(normal));

          // apply change of basis to pts
          closestPts.clear();
          for (auto& vert : closestVerts)
            closestPts.push_back(changeOfBasis.project_forward(vert->get_point_orig(0)));
          utils::Point pt0 = changeOfBasis.project_forward(vertIn->get_point_orig(0));

          // compute taylor patch
          m_surfacePatch->construct_patch(closestPts, pt0);

          utils::Point newPt = m_surfacePatch->eval_point(pt0.x, pt0.y);
          auto newPt2        = changeOfBasis.project_back(newPt);
          vertIn->set_point_orig(0, newPt2);
        }
      }
    }
}
} // namespace impl
} // namespace nonconformal4
} // namespace middle_mesh
} // namespace stk
