#include "geometry_improver_restore_mesh2_verts.hpp"

namespace stk {
namespace middle_mesh {
namespace nonconformal4 {
namespace impl {

GeometryImproverRestoreMesh2Verts::GeometryImproverRestoreMesh2Verts(std::shared_ptr<mesh::Mesh> mesh1,
                                                                     std::shared_ptr<mesh::Mesh> mesh2,
                                                                     std::shared_ptr<mesh::Mesh> meshIn,
                                                                     std::shared_ptr<MeshRelationalData> relationalData)
  : GeometryImprover(mesh1, mesh2, meshIn, relationalData)
{}

void GeometryImproverRestoreMesh2Verts::run()
{
  auto& verts2ToFakeVerts  = *(m_relationalData->verts2ToFakeVerts);
  auto& fakeVertsToVertsIn = m_relationalData->fakeVertsToVertsIn;
  for (auto& vert2 : m_mesh2->get_vertices())
    if (vert2)
    {
      auto fv                    = verts2ToFakeVerts(vert2, 0, 0);
      mesh::MeshEntityPtr vertIn = fakeVertsToVertsIn[fv.id];
      vertIn->set_point_orig(0, vert2->get_point_orig(0));
    }
}

} // namespace impl
} // namespace nonconformal4
} // namespace middle_mesh
} // namespace stk
