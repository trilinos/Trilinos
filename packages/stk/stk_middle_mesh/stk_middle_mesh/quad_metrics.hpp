#ifndef STK_MIDDLE_MESH_QUAD_METRICS
#define STK_MIDDLE_MESH_QUAD_METRICS

#include "projection.hpp"
#include "mesh.hpp"

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

class QuadMetrics
{
  public:
    void compute_det_jacobian(mesh::MeshEntityPtr quad, const std::vector<utils::Point>& xiPts, std::vector<double>& detJ);

  private:

    void compute_det_jacobian(std::array<utils::Point, mesh::MAX_DOWN>& vertCoords, const std::vector<utils::Point>& xiPts,
                              std::vector<double>& detJ);
};

}
}
}
}

#endif