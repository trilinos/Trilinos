#include "change_of_basis.hpp"

#include <cassert>
#include <cmath>
#include <stdexcept>

#include "plane_projection.hpp"

namespace stk {
namespace middle_mesh {
namespace utils {
namespace impl {

// returns an orthonormal basis given a normal vector
// The first two vectors lie in a plane defined by the normal vector,
// the third vector is normal to the plane (parallel to the input)
std::array<Point, 3> compute_basis(Point norm)
{
  const double eps = 1e-13;
  double len       = std::sqrt(dot(norm, norm));
  if (len < eps)
    throw std::invalid_argument("normal vector cannot be zero");

  norm = norm / len;

  const bool nxZero = std::abs(norm.x) < eps;
  const bool nyZero = std::abs(norm.y) < eps;
  const bool nzZero = std::abs(norm.z) < eps;
  int nzero         = 0;
  nzero += nxZero;
  nzero += nyZero;
  nzero += nzZero;

  if (nzero == 3)
    throw std::invalid_argument("normal vector cannot be zero");

  Point b1, b2;
  if (nzero == 2) // doubly degenerate case
  {
    std::array<Point, 3> pts;
    if (!nxZero)
      if (norm.x > 0)
        pts = get_plane_projection_points(PlaneProjection::YzPos);
      else
        pts = get_plane_projection_points(PlaneProjection::YzNeg);
    else if (!nyZero)
      if (norm.y > 0)
        pts = get_plane_projection_points(PlaneProjection::XzPos);
      else
        pts = get_plane_projection_points(PlaneProjection::XzNeg);
    else // !nz_zero
      if (norm.z > 0)
        pts = get_plane_projection_points(PlaneProjection::XyPos);
      else
        pts = get_plane_projection_points(PlaneProjection::XyNeg);

    b1 = pts[1] - pts[0];
    b2 = pts[2] - pts[0];
  } else if (nxZero) // handle singlely degenerate cases
  {
    // these are actually (y - y0)
    double y1 = -std::copysign(1, norm.z);
    double y2 = -y1;
    double z1 = -norm.y * y1 / norm.z;
    double z2 = -norm.y * y2 / norm.z;

    double x1 = 1;
    double x2 = -(y1 * y2 + z1 * z2);

    b1 = Point(x1, y1, z1);
    b2 = Point(x2, y2, z2);
  } else if (nyZero)
  {
    // these are actually (x - x0)
    double x1 = std::copysign(1, norm.z);
    double x2 = -x1;
    double z1 = -norm.x * x1 / norm.z;
    double z2 = -norm.x * x2 / norm.z;

    double y1 = 1;
    double y2 = -(x1 * x2 + z1 * z2);

    b1 = Point(x1, y1, z1);
    b2 = Point(x2, y2, z2);
  } else if (nzZero)
  {
    // these are actually (x - x0)
    double x1 = -std::copysign(-1, norm.y);
    double x2 = -x1;
    double y1 = -norm.x * x1 / norm.y;
    double y2 = -norm.x * x2 / norm.y;

    // chose z such that the two vectors are orthogonal
    double z1 = 1;
    double z2 = -(x1 * x2 + y1 * y2); // divided by z1, which is = 1
    // double fac = x2/x1 + 1;

    b1 = Point(x1, y1, z1);
    b2 = Point(x2, y2, z2);
  } else // general (non-degenerate) case
  {
    double x1 = -std::copysign(1, norm.z), y1 = 0;
    double z1 = -norm.x * x1 / norm.z;

    b1 = Point(x1, y1, z1);
    b2 = cross(norm, b1);
  }

  // orthonormalize
  b1 = b1 / std::sqrt(dot(b1, b1));
  // b2 = b2 - dot(b1, b2)*b1;
  b2 = b2 / std::sqrt(dot(b2, b2));

  //  std::cout << "dot1 = " << dot(b1, b2) << std::endl;
  //  std::cout << "dot2 = " << dot(b1, norm) << std::endl;
  //  std::cout << "dot3 = " << dot(b2, norm) << std::endl;
  // verify orthonormal
  assert(dot(b1, b2) < eps);
  assert(dot(b1, norm) < eps);
  assert(dot(b2, norm) < eps);

  return {b1, b2, norm};
}

} // namespace impl

} // namespace utils
} // namespace middle_mesh
} // namespace stk
