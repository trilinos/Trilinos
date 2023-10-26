#ifndef GEOMETRY_IMPROVER_FACTORY_H
#define GEOMETRY_IMPROVER_FACTORY_H

#include "geometry_improver.hpp"

namespace stk {
namespace middle_mesh {
namespace nonconformal4 {
namespace impl {

enum class GeometryImprovers
{
  RestoreMesh2Verts,
  EdgeVertexMidpoints,
  EdgeVertexTaylorPatchLinear,
  EdgeVertexTaylorPatchQuadratic,
  EdgeVertexTaylorPatchQuadraticZonly, // FOR TESTING ONLY
  EdgeVertexCubicBSplinePatch25Pts
};

std::shared_ptr<GeometryImprover>
geometry_improver_factory(GeometryImprovers improver, std::shared_ptr<mesh::Mesh> mesh1,
                          std::shared_ptr<mesh::Mesh> mesh2, std::shared_ptr<mesh::Mesh> meshIn,
                          std::shared_ptr<MeshRelationalData> relationalData, mesh::FieldPtr<utils::Point> normalField);

} // namespace impl
} // namespace nonconformal4
} // namespace middle_mesh
} // namespace stk

#endif