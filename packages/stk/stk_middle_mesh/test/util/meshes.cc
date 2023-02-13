#include "meshes.h"

namespace stk {
namespace middle_mesh {
namespace impl {

const double PI = std::atan(1) * 4;

std::shared_ptr<mesh::Mesh> make_annulus_mesh(const int nelemR, const int nelemTheta, const double rIn,
                                              const double rOut, const double dtheta)
{
  auto f2 = [](const double x, const double y) { return 0.0; };
  return make_annulus_mesh(nelemR, nelemTheta, rIn, rOut, dtheta, f2);
}

} // namespace impl
} // namespace middle_mesh
} // namespace stk
