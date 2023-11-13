#include "geometry_improver_factory.hpp"

#include "bspline_patch.hpp"
#include "geometry_improver_edge_vertex_midpoints.hpp"
#include "geometry_improver_edge_vertex_surface_patch.hpp"
#include "geometry_improver_restore_mesh2_verts.hpp"
#include "taylor_patch.hpp"

namespace stk {
namespace middle_mesh {
namespace nonconformal4 {
namespace impl {

std::shared_ptr<GeometryImprover>
geometry_improver_factory(GeometryImprovers improver, std::shared_ptr<mesh::Mesh> mesh1,
                          std::shared_ptr<mesh::Mesh> mesh2, std::shared_ptr<mesh::Mesh> meshIn,
                          std::shared_ptr<MeshRelationalData> relationalData, mesh::FieldPtr<utils::Point> normalField)
{
  switch (improver)
  {
    case GeometryImprovers::RestoreMesh2Verts: {
      return std::make_shared<GeometryImproverRestoreMesh2Verts>(mesh1, mesh2, meshIn, relationalData);
    }

    case GeometryImprovers::EdgeVertexMidpoints: {
      return std::make_shared<GeometryImproverEdgeVertexMidPoints>(mesh1, mesh2, meshIn, relationalData);
    }

    case GeometryImprovers::EdgeVertexTaylorPatchLinear: {
      auto taylorPatch = std::make_shared<utils::impl::TaylorPatch>();
      return std::make_shared<GeometryImproverEdgeVertexSurfacePatch>(mesh1, mesh2, meshIn, relationalData, normalField,
                                                                      taylorPatch, 3);
    }

    case GeometryImprovers::EdgeVertexTaylorPatchQuadratic: {
      auto taylorPatch = std::make_shared<utils::impl::TaylorPatch>();
      return std::make_shared<GeometryImproverEdgeVertexSurfacePatch>(mesh1, mesh2, meshIn, relationalData, normalField,
                                                                      taylorPatch, 6);
    }

    case GeometryImprovers::EdgeVertexTaylorPatchQuadraticZonly: {
      auto taylorPatch = std::make_shared<utils::impl::TaylorPatch>();
      return std::make_shared<GeometryImproverEdgeVertexSurfacePatch>(mesh1, mesh2, meshIn, relationalData, nullptr,
                                                                      taylorPatch, 6);
    }

    case GeometryImprovers::EdgeVertexCubicBSplinePatch25Pts: {
      auto bsplinePatch = std::make_shared<utils::impl::BSplinePatch>(2);
      return std::make_shared<GeometryImproverEdgeVertexSurfacePatch>(mesh1, mesh2, meshIn, relationalData, normalField,
                                                                      bsplinePatch, 25);
    }
    default:
      throw std::runtime_error("unhandled enum value");
  }
}

} // namespace impl

} // namespace nonconformal4
} // namespace middle_mesh
} // namespace stk
