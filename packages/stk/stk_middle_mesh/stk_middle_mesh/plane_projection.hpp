#ifndef PLANE_PROJECTION_H
#define PLANE_PROJECTION_H

#include <array>
#include <ostream>

#include "projection.hpp"

namespace stk {
namespace middle_mesh {
namespace utils {
namespace impl {

// enum describes which plane has normal vector closest to the given
// normal vector.  The positive or negative modifier describes the
// sign of the third component of the normal vector that gives a positively
// oriented element
enum class PlaneProjection
{
  YzPos = 0,
  YzNeg = 1,
  XzPos = 2,
  XzNeg = 3,
  XyPos = 4,
  XyNeg = 5,
};

std::ostream& operator<<(std::ostream& os, PlaneProjection val);

PlaneProjection get_plane_projection_enum(const Point& pt1, const Point& pt2, const Point& pt3);

// given a vector that defines a line, gives the plane that the vector
// has the largest projection onto
PlaneProjection get_line_projection_enum(const Point& b);

std::array<Point, 3> get_plane_projection_points(PlaneProjection val, const double scaleX = 1, const double scaleY = 1,
                                                 const double scaleZ = 1);

Projection* get_plane_projection(std::array<Point, 3> pts);

Projection* get_plane_projection(const Point& pt1, const Point& pt2, const Point& pt3);

Projection* get_scaled_plane_projection(const Point& pt1, const Point& pt2, const Point& pt3);

Point apply_plane_projection(PlaneProjection val, const Point& pt);

Point apply_inverse_plane_projection(PlaneProjection val, const Point& pt);

} // namespace impl

} // namespace utils
} // namespace middle_mesh
} // namespace stk
#endif
