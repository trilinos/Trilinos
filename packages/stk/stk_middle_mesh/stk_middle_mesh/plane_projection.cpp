#include "plane_projection.hpp"

namespace stk {
namespace middle_mesh {
namespace utils {
namespace impl {

std::ostream& operator<<(std::ostream& os, PlaneProjection val)
{
  switch (val)
  {
    case PlaneProjection::YzPos: {
      os << "YzPos";
      break;
    }
    case PlaneProjection::YzNeg: {
      os << "YzNeg";
      break;
    }
    case PlaneProjection::XzPos: {
      os << "XzPos";
      break;
    }
    case PlaneProjection::XzNeg: {
      os << "XzNeg";
      break;
    }
    case PlaneProjection::XyPos: {
      os << "XyPos";
      break;
    }
    case PlaneProjection::XyNeg: {
      os << "XyNeg";
      break;
    }
  }

  return os;
}

PlaneProjection get_plane_projection_enum(const Point& pt1, const Point& pt2, const Point& pt3)
{
  // project onto plane whose normal vector is closest to the normal
  // vector of the plan el1 lies in
  auto b1  = pt2 - pt1;
  auto b2  = pt3 - pt1;
  auto nrm = cross(b1, b2);

  // normal is zero if b1 and b2 are parallel
  double nrmMag = std::sqrt(dot(nrm, nrm));
  if (nrmMag < 1e-13)
    return get_line_projection_enum(b1);

  if (std::abs(nrm.x) >= std::abs(nrm.y) && std::abs(nrm.x) >= std::abs(nrm.z))
  {
    // y-z plane
    if (nrm.x > 0)
      return PlaneProjection::YzPos;
    else
      return PlaneProjection::YzNeg;
  } else if (std::abs(nrm.y) >= std::abs(nrm.x) && std::abs(nrm.y) >= std::abs(nrm.z))
  {
    // x-z plane
    if (nrm.y > 0)
      return PlaneProjection::XzPos;
    else
      return PlaneProjection::XzNeg;
  } else if (std::abs(nrm.z) >= std::abs(nrm.x) && std::abs(nrm.z) >= std::abs(nrm.y))
  {
    // x-y plane
    if (nrm.z > 0)
      return PlaneProjection::XyPos;
    else
      return PlaneProjection::XyNeg;
  }

  throw std::invalid_argument("failed to find plane projection");
}

// given a vector that defines a line, gives the plane that the vector
// has the largest projection onto
PlaneProjection get_line_projection_enum(const Point& b)
{
  double xy = b.x * b.x + b.y * b.y;
  double yz = b.y * b.y + b.z * b.z;
  double xz = b.x * b.x + b.z * b.z;

  if (yz >= xy && yz >= xz)
    return PlaneProjection::YzPos;
  else if (xz >= xy && xz >= yz)
    return PlaneProjection::XzPos;
  else if (xy >= yz && xy >= xz)
    return PlaneProjection::XyPos;

  throw std::invalid_argument("failed to find line projection");
}

std::array<Point, 3> get_plane_projection_points(PlaneProjection val, const double scaleX, const double scaleY,
                                                 const double scaleZ)
{
  // project onto plane whose normal vector is closest to the normal
  // vector of the plan el1 lies in

  Point ptA, ptB;
  double ax = std::abs(scaleX);
  double ay = std::abs(scaleY);
  double az = std::abs(scaleZ);

  if (val == PlaneProjection::YzPos)
  {
    ptA = Point(0, ay, 0);
    ptB = Point(0, 0, az);
  } else if (val == PlaneProjection::YzNeg)
  {
    ptA = Point(0, 0, az);
    ptB = Point(0, ay, 0);
  } else if (val == PlaneProjection::XzPos)
  {
    ptA = Point(0, 0, az);
    ptB = Point(ax, 0, 0);
  } else if (val == PlaneProjection::XzNeg)
  {
    ptA = Point(ax, 0, 0);
    ptB = Point(0, 0, az);
  } else if (val == PlaneProjection::XyPos)
  {
    ptA = Point(ax, 0, 0);
    ptB = Point(0, ay, 0);
  } else if (val == PlaneProjection::XyNeg)
  {
    ptA = Point(0, ay, 0);
    ptB = Point(ax, 0, 0);
  }

  return {Point(0, 0, 0), ptA, ptB};
}

// given the result of getPlaneProjection, creates the projection operator
Projection* get_plane_projection(std::array<Point, 3> pts)
{
  return new Projection(pts[0], pts[1], pts[2]);
}

Projection* get_plane_projection(const Point& pt1, const Point& pt2, const Point& pt3)
{
  auto val = get_plane_projection_enum(pt1, pt2, pt3);
  auto pts = get_plane_projection_points(val);
  return get_plane_projection(pts);
}

Projection* get_scaled_plane_projection(const Point& pt1, const Point& pt2, const Point& pt3)
{
  auto val = get_plane_projection_enum(pt1, pt2, pt3);

  double ax = std::max(std::abs(pt1.x), std::abs(pt2.x));
  ax        = std::max(ax, std::abs(pt3.x));
  double ay = std::max(std::abs(pt1.y), std::abs(pt2.y));
  ay        = std::max(ay, std::abs(pt3.y));
  double az = std::max(std::abs(pt1.z), std::abs(pt2.z));
  az        = std::max(az, std::abs(pt3.z));

  auto pts = get_plane_projection_points(val, ax, ay, az);
  return get_plane_projection(pts);
}

Point apply_plane_projection(PlaneProjection val, const Point& pt)
{
  if (val == PlaneProjection::YzPos)
    return Point(pt.y, pt.z, pt.x);
  else if (val == PlaneProjection::YzNeg)
    return Point(pt.z, pt.y, -pt.x);
  else if (val == PlaneProjection::XzPos)
    return Point(pt.z, pt.x, pt.y);
  else if (val == PlaneProjection::XzNeg)
    return Point(pt.x, pt.z, -pt.y);
  else if (val == PlaneProjection::XyPos)
    return Point(pt.x, pt.y, pt.z);
  else if (val == PlaneProjection::XyNeg)
    return Point(pt.y, pt.x, -pt.z);

  throw std::invalid_argument("enum value not recognized");
}

// given the result of applyPlaneProjection, computes the corresponding
// input
Point apply_inverse_plane_projection(PlaneProjection val, const Point& pt)
{
  if (val == PlaneProjection::YzPos)
    return Point(pt.z, pt.x, pt.y);
  else if (val == PlaneProjection::YzNeg)
    return Point(-pt.z, pt.y, pt.x);
  else if (val == PlaneProjection::XzPos)
    return Point(pt.y, pt.z, pt.x);
  else if (val == PlaneProjection::XzNeg)
    return Point(pt.x, -pt.z, pt.y);
  else if (val == PlaneProjection::XyPos)
    return Point(pt.x, pt.y, pt.z);
  else if (val == PlaneProjection::XyNeg)
    return Point(pt.y, pt.x, -pt.z);

  throw std::invalid_argument("enum value not recognized");
}

} // namespace impl

} // namespace utils
} // namespace middle_mesh
} // namespace stk
