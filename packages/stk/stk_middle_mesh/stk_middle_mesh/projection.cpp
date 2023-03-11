#include "projection.hpp"

namespace stk {
namespace middle_mesh {
namespace utils {
namespace impl {

Projection::Projection(const Point& pt1, const Point& pt2, const Point& pt3)
{
  // compute orthonormal basis
  m_q1 = pt2 - pt1;
  m_q2 = pt3 - pt1;
  m_q2 = m_q2 - project(m_q2, m_q1);
  m_d1 = std::sqrt(dot(m_q1, m_q1));
  m_d2 = std::sqrt(dot(m_q2, m_q2));
  m_q1 = m_q1 / m_d1;
  m_q2 = m_q2 / m_d2;

  m_v0 = pt1;

  m_norm = cross(m_q1, m_q2);
}

Point Projection::project_to_plane(const Point& pt) const
{
  auto v = pt - m_v0;

  // project v onto n (n is unit length, so no need for denominator)
  auto vn = v - dot(v, m_norm) * m_norm;
  return m_v0 + vn;
}

Point Projection::project_to_plane_rev(const Point& pt, const Point& ptBar) const
{
  auto v = pt - m_v0;

  // project v onto n (n is unit length, so no need for denominator)
  // auto t1 = dot(v, m_norm);
  // auto vn = v - t1 * m_norm;
  // return m_v0 + vn;

  auto vnBar = ptBar;
  auto vBar  = vnBar;
  auto t1Bar = -m_norm.x * vnBar.x - m_norm.y * vnBar.y - m_norm.z * vnBar.z;

  Point vBarTmp, mNormBar;
  dot_rev(v, m_norm, vBarTmp, mNormBar, t1Bar);
  vBar += vBarTmp;

  return vBar;
}

// compute coordinates on basis [m_q1, m_q2]
// the point *must* lie in the plane
// The output point has z coordinate = 0
Point Projection::compute_plane_coords(const Point& pt) const
{
  // solve least squares problem [m_q1, m_q2] D xi = pt
  // [m_q1, m_q2] is unitary, so the solution is D^-1 [m_q1, m_q2]^T * pt
  auto xi1 = dot(m_q1, pt) / m_d1;
  auto xi2 = dot(m_q2, pt) / m_d2;

  return {xi1, xi2, 0};
}

Point Projection::compute_plane_coords_rev(const Point& pt, const Point& ptBar) const
{
  // auto xi1 = dot(m_q1, pt)/m_d1;
  // auto xi2 = dot(m_q2, pt)/m_d2;
  // return {xi1, xi2, 0};

  Point qBar, ptBarTmp, ptBarSum;
  dot_rev(m_q2, pt, qBar, ptBarTmp, ptBar.y / m_d2);
  ptBarSum += ptBarTmp;

  dot_rev(m_q1, pt, qBar, ptBarTmp, ptBar.x / m_d1);
  ptBarSum += ptBarTmp;

  return ptBarSum;
}

// combines the above into one function call
Point Projection::project_plane_coords(const Point& pt) const
{
  auto pt1 = project_to_plane(pt);
  auto pt2 = compute_plane_coords(pt1);
  return pt2;
}

Point Projection::project_plane_coords_rev(const Point& pt, const Point& pt2Bar) const
{
  auto pt1 = project_to_plane(pt);
  // auto pt2 = computePlaneCoords(pt1);

  auto pt1Bar = compute_plane_coords_rev(pt1, pt2Bar);
  auto ptBar  = project_to_plane_rev(pt, pt1Bar);

  return ptBar;
}

std::ostream& operator<<(std::ostream& os, const Projection& p)
{
  os << "projection with basis vectors " << p.m_q1 << ", " << p.m_q2 << ", m_normal = " << p.m_norm;
  return os;
}

} // namespace impl

} // namespace utils
} // namespace middle_mesh
} // namespace stk
