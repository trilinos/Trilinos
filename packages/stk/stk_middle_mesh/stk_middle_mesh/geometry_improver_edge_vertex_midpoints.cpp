#include "geometry_improver_edge_vertex_midpoints.hpp"
#include "stk_util/util/ReportHandler.hpp"
#include "utils.hpp"

namespace stk {
namespace middle_mesh {
namespace nonconformal4 {
namespace impl {

GeometryImproverEdgeVertexMidPoints::GeometryImproverEdgeVertexMidPoints(
    std::shared_ptr<mesh::Mesh> mesh1, std::shared_ptr<mesh::Mesh> mesh2, std::shared_ptr<mesh::Mesh> meshIn,
    std::shared_ptr<MeshRelationalData> relationalData)
  : GeometryImprover(mesh1, mesh2, meshIn, relationalData)
{
  if (mesh1)
  {
    STK_ThrowRequireMsg(utils::impl::comm_size(mesh1->get_comm()) == 1,
                        "GeometryImproverEdgeVertexMidPoints is not supported in parallel");
  }

  if (mesh2)
  {
    STK_ThrowRequireMsg(utils::impl::comm_size(mesh2->get_comm()) == 1,
                        "GeometryImproverEdgeVertexMidPoints is not supported in parallel");
  }  
}

void GeometryImproverEdgeVertexMidPoints::run()
{
  move_edge_vertices_to_mid_point();
  restore_mesh1_vertices();
}

void GeometryImproverEdgeVertexMidPoints::move_edge_vertices_to_mid_point()
{
  auto& edges2ToFakeVertsIn = *(m_relationalData->edges2ToFakeVertsIn);
  auto& fakeVertsToVertsIn  = m_relationalData->fakeVertsToVertsIn;
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

          utils::Point edge2Pt = mesh::compute_edge_coords_orig(edge2, vertOnEdge.xi);
          utils::Point edge1Pt = vertIn->get_point_orig(0);
          vertIn->set_point_orig(0, (edge1Pt + edge2Pt) / 2);
        }
      }
    }
}

// the first pass can move some mesh1 vertices that are on
// mesh2 edges.  Move them back to their original positions
// (which presumably are on the true geometry)
void GeometryImproverEdgeVertexMidPoints::restore_mesh1_vertices()
{
  auto& verts1ToFakeVerts  = *(m_relationalData->verts1ToFakeVerts);
  auto& fakeVertsToVertsIn = m_relationalData->fakeVertsToVertsIn;

  for (auto& vert1 : m_mesh1->get_vertices())
    if (vert1)
    {
      FakeVert fv                = verts1ToFakeVerts(vert1, 0, 0);
      mesh::MeshEntityPtr vertIn = fakeVertsToVertsIn[fv.id];
      vertIn->set_point_orig(0, vert1->get_point_orig(0));
    }
}

} // namespace impl
} // namespace nonconformal4
} // namespace middle_mesh
} // namespace stk
