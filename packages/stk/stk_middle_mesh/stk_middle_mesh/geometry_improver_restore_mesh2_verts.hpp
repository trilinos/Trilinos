#ifndef GEOMETRY_IMPROVER_RESTORE_MESH2_VERTS_H
#define GEOMETRY_IMPROVER_RESTORE_MESH2_VERTS_H

#include "geometry_improver.hpp"

namespace stk {
namespace middle_mesh {
namespace nonconformal4 {
namespace impl {

// moves the vertices on mesh_in that came from mesh2 back to their original
// locations on mesh2
class GeometryImproverRestoreMesh2Verts : public GeometryImprover
{
  public:
    GeometryImproverRestoreMesh2Verts(std::shared_ptr<mesh::Mesh> mesh1, std::shared_ptr<mesh::Mesh> mesh2,
                                      std::shared_ptr<mesh::Mesh> meshIn,
                                      std::shared_ptr<MeshRelationalData> relationalData);

    void run() override;
};

} // namespace impl

} // namespace nonconformal4
} // namespace middle_mesh
} // namespace stk
#endif