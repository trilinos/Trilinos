#include "util/meshes.hpp"

namespace stk {
namespace middle_mesh {
namespace impl {

const double PI = std::atan(1) * 4;

std::shared_ptr<mesh::Mesh> make_annulus_mesh(const int nelemR, const int nelemTheta, const double rIn,
                                              const double rOut, const double dtheta, MPI_Comm comm, bool createTriangles)
{
  auto f2 = [](const double /*x*/, const double /*y*/) { return 0.0; };
  return make_annulus_mesh(nelemR, nelemTheta, rIn, rOut, dtheta, comm, createTriangles, f2);
}

} // namespace impl
} // namespace middle_mesh
} // namespace stk
