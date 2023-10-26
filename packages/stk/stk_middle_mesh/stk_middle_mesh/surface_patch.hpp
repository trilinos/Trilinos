#ifndef SURFACE_PATCH_H
#define SURFACE_PATCH_H

#include "matrix.hpp"
#include "projection.hpp"

namespace stk {
namespace middle_mesh {
namespace utils {
namespace impl {

class SurfacePatch
{
  public:
    virtual ~SurfacePatch() {}

    // constructs a patch centered at pts[0]
    void construct_patch(const std::vector<Point>& pts) { construct_patch(pts, pts[0]); }

    virtual void construct_patch(const std::vector<Point>& pts, const Point& pt0) = 0;

    // evaluates a point on the patch at (x, y) coordinates
    virtual Point eval_point(const double x, const double y) = 0;

    // evaluates  dz/dx and dz/dy.  derivs must be a 3 x 2 matrix
    virtual void eval_deriv(const double x, const double y, double derivs[2]) = 0;
};

} // namespace impl

} // namespace utils
} // namespace middle_mesh
} // namespace stk
#endif
