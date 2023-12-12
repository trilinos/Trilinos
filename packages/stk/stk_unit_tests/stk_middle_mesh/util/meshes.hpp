#ifndef TEST_MESHES_H
#define TEST_MESHES_H

#include "stk_middle_mesh/create_mesh.hpp"

namespace stk {
namespace middle_mesh {
namespace impl {

extern const double PI;

// makes a mesh where the x,y coordinates form an annulus
// dtheta: theta offset applied to all theta values except those
// at theta = 0 (or 2pi)
// func2: is a function z = f(x, y) that determines the z coordinate
template <typename T>
std::shared_ptr<mesh::Mesh> make_annulus_mesh(const int nelemR, const int nelemTheta, const double rIn,
                                              const double rOut, double dtheta, MPI_Comm comm, bool createTriangles,
                                              T func2)
{
  assert(rIn > 0);
  assert(rOut > rIn);

  mesh::impl::MeshSpec spec;
  spec.numelX    = nelemR;
  spec.numelY    = nelemTheta;
  spec.xmin      = rIn;
  spec.xmax      = rOut;
  spec.ymin      = 0;
  spec.ymax      = 2 * PI;
  spec.yPeriodic = true;

  auto func = [&](const utils::Point& pt) {
    // interpret x and y as r and theta
    double r     = pt.x;
    double theta = pt.y;
    if (std::fmod(theta, 2 * PI) > 1e-12)
      theta += dtheta;

    double x = r * std::cos(theta);
    double y = r * std::sin(theta);
    double z = func2(x, y);
    utils::Point pt2(x, y, z);
    return pt2;
  };

  return mesh::impl::create_mesh(spec, func, comm, createTriangles);
}

std::shared_ptr<mesh::Mesh> make_annulus_mesh(const int elemR, const int nelemTheta, const double rIn,
                                              const double rOut, double dtheta, MPI_Comm comm = MPI_COMM_WORLD, bool createTriangles=false);

} // namespace impl
} // namespace middle_mesh
} // namespace stk
#endif
