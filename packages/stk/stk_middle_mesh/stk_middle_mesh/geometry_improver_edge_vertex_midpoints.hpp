#ifndef GEOMETRY_IMPROVER_EDGE_VERTEX_MIDPOINTS
#define GEOMETRY_IMPROVER_EDGE_VERTEX_MIDPOINTS

#include "geometry_improver.hpp"

namespace stk {
namespace middle_mesh {
namespace nonconformal4 {
namespace impl {

// This updates the coordinates of so-called "phantom vertices" (vertices that
// come from edge-edge intersections) to be at the midpoint of the line segment
// defined by the intersection point on the mesh2 edge and the intersection point
// on the mesh1 edge
class GeometryImproverEdgeVertexMidPoints : public GeometryImprover
{
  public:
    GeometryImproverEdgeVertexMidPoints(std::shared_ptr<mesh::Mesh> mesh1, std::shared_ptr<mesh::Mesh> mesh2,
                                        std::shared_ptr<mesh::Mesh> meshIn,
                                        std::shared_ptr<MeshRelationalData> relationalData);

    void run() override;

  private:
    void move_edge_vertices_to_mid_point();

    // the first pass can move some mesh1 vertices that are on
    // mesh2 edges.  Move them back to their original positions
    // (which presumably are on the true geometry)
    void restore_mesh1_vertices();
};

} // namespace impl

} // namespace nonconformal4
} // namespace middle_mesh
} // namespace stk
#endif