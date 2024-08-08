#include "SideGeometry.hpp"
#include "stk_util/util/ReportHandler.hpp"
#include <algorithm>

namespace stk { namespace math {

SideGeometry::SideGeometry(size_t numNodes)
  : m_numNodes(numNodes)
{}

double
SideGeometry::min_distance_to_point(const stk::math::Vector3d & point) const
{
  stk::math::Vector3d projection = closest_proj_on_face(point);
  return (point - projection).length();
}

bool
SideGeometry::are_nodes_close_to_side(const SideGeometry & otherSide, double tolerance)
{
  if (otherSide.min_distance_to_point(centroid()) < tolerance) {
    return true;
  }
  for (size_t i = 0; i < m_numNodes; ++i) {
    if (otherSide.min_distance_to_point(node(i)) < tolerance) {
      return true;
    }
  }
  return false;
}

PointGeometry::PointGeometry(const stk::math::Vector3d & n)
  : SideGeometry(1),
    m_nodeData(n)
{}

stk::math::Vector3d
PointGeometry::centroid() const
{
  return node(0);
}

stk::math::Vector3d
PointGeometry::closest_proj_on_face(const stk::math::Vector3d & p) const
{
  return node(0);
}

LineGeometry::LineGeometry(const stk::math::Vector3d & n0,
                           const stk::math::Vector3d & n1)
  : SideGeometry(2),
    m_nodeData{n0, n1}
{}
 
stk::math::Vector3d
LineGeometry::centroid() const
{
  stk::math::Vector3d cent;
  for (unsigned i = 0u; i < 3u; ++i) {
    cent[i] = (node(0)[i] + node(1)[i]) / 2.0;
  }
  return cent;
}

stk::math::Vector3d
LineGeometry::closest_proj_on_face(const stk::math::Vector3d & p) const
{
  const stk::math::Vector3d & a(node(0));
  const stk::math::Vector3d & b(node(1));

  stk::math::Vector3d ab(b - a);

  stk::math::Vector3d ap(p - a);
  stk::math::Vector3d bp(p - b);

  double alpha(stk::math::Dot(ab, ap));
  double beta(stk::math::Dot(-ab, bp));

  // Projects to node a
  if (alpha <= 0.0) {
    return a;
  }

  // Projects to node b
  if (beta <= 0.0) {
    return b;
  }

  // Projects to middle
  stk::math::Vector3d abHat = ab.unit_vector();
  double scalarProjection = stk::math::Dot(ap, abHat);
  return a + scalarProjection*abHat;
}

TriGeometry::TriGeometry(const stk::math::Vector3d & n0,
                             const stk::math::Vector3d & n1,
                             const stk::math::Vector3d & n2)
  : SideGeometry(3),
    m_nodeData{n0, n1, n2}
{}

stk::math::Vector3d
TriGeometry::centroid() const
{
  stk::math::Vector3d cent;
  for (unsigned i = 0u; i < 3u; ++i) {
    cent[i] = (node(0)[i] + node(1)[i] + node(2)[i]) / 3.0;
  }
  return cent;
}

stk::math::Vector3d
TriGeometry::closest_proj_on_face(const stk::math::Vector3d & p) const
{
  const stk::math::Vector3d & a(node(0));
  const stk::math::Vector3d & b(node(1));
  const stk::math::Vector3d & c(node(2));

  stk::math::Vector3d ab(b - a);
  stk::math::Vector3d ac(c - a);
  stk::math::Vector3d ap(p - a);

  double d1(stk::math::Dot(ab, ap));
  double d2(stk::math::Dot(ac, ap));

  if (d1 <= 0.0 && d2 <= 0.0) {
    return a;  // Node projection
  }

  stk::math::Vector3d bp(p - b);
  double d3(stk::math::Dot(ab, bp));
  double d4(stk::math::Dot(ac, bp));

  if (d3 >= 0.0 && d4 <= d3) {
    return b;  // Node projection
  }

  double vc(d1 * d4 - d3 * d2);
  if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
    double v = d1 / (d1 - d3);
    return (a + v * ab);  // Edge projection
  }

  stk::math::Vector3d cp(p - c);
  double d5(stk::math::Dot(ab, cp));
  double d6(stk::math::Dot(ac, cp));
  if (d6 >= 0.0 && d5 <= d6) {
    return c;  // Node projection
  }

  double vb(d5 * d2 - d1 * d6);
  if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
    double w(d2 / (d2 - d6));
    return (a + w * ac);  // Edge projection
  }

  double va(d3 * d6 - d5 * d4);
  if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
    double w((d4 - d3) / ((d4 - d3) + (d5 - d6)));
    return (b + w * (c - b));  // Edge projection
  }

  double denom(1.0 / (va + vb + vc));
  double v(vb * denom);
  double w(vc * denom);
  return (a + ab * v + ac * w); // Face projection
}


QuadGeometry::QuadGeometry(const stk::math::Vector3d & n0,
               const stk::math::Vector3d & n1,
               const stk::math::Vector3d & n2,
               const stk::math::Vector3d & n3)
  : SideGeometry(4),
    m_nodeData{n0, n1, n2, n3}
{
}

stk::math::Vector3d
QuadGeometry::centroid() const
{
  stk::math::Vector3d cent;
  for (unsigned i = 0u; i < 3u; ++i) {
    cent[i] = (node(0)[i] + node(1)[i] + node(2)[i] + node(3)[i]) / 4.0;
  }
  return cent;
}

stk::math::Vector3d
QuadGeometry::normal_dir() const
{
  return stk::math::Cross(node(2) - node(0), node(3) - node(1));
}

stk::math::Vector3d
QuadGeometry::closest_proj_on_face(const stk::math::Vector3d & point) const
{
  stk::math::Vector3d face_norm_dir = normal_dir();

  stk::math::Vector3d edgeA(node(0) - node(3));
  stk::math::Vector3d point3(point - node(3));
  stk::math::Vector3d point0(point - node(0));

  double L9 = -stk::math::Dot(stk::math::Cross(edgeA, face_norm_dir), point3);

  stk::math::Vector3d edgeC(node(1) - node(0));
  stk::math::Vector3d point1(point - node(1));

  double L10 = -stk::math::Dot(stk::math::Cross(edgeC, face_norm_dir), point0);

  stk::math::Vector3d edgeB(node(2) - node(1));
  stk::math::Vector3d point2(point - node(2));
  double L11 = -stk::math::Dot(stk::math::Cross(edgeB, face_norm_dir), point1);

  stk::math::Vector3d edgeD(node(3) - node(2));
  double L12 = -stk::math::Dot(stk::math::Cross(edgeD, face_norm_dir), point2);

  if (L9 >= 0.0 && L10 >= 0.0 && L11 >= 0.0 && L12 >= 0.0) {
    return (point - (stk::math::Dot(point3, face_norm_dir) * face_norm_dir) / face_norm_dir.length_squared());  // Face projection
  }

  double L1 = stk::math::Dot(edgeA, point3);
  double L2 = stk::math::Dot(edgeA, point0);
  double L3 = stk::math::Dot(edgeC, point0);
  double L4 = stk::math::Dot(edgeC, point1);
  double L5 = stk::math::Dot(edgeB, point1);
  double L6 = stk::math::Dot(edgeB, point2);
  double L7 = stk::math::Dot(edgeD, point2);
  double L8 = stk::math::Dot(edgeD, point3);

  if (L2 >= 0.0 && L3 <= 0.0) {
    return node(0);  // Node projection
  }
  if (L4 >= 0.0 && L5 <= 0.0) {
    return node(1);  // Node projection
  }
  if (L6 >= 0.0 && L7 <= 0.0) {
    return node(2);  // Node projection
  }
  if (L8 >= 0.0 && L1 <= 0.0) {
    return node(3);  // Node projection
  }

  if (L1 >= 0.0 && L9 <= 0.0 && L2 <= 0.0) {
    return (node(3) + (stk::math::Dot(edgeA, point3) * edgeA) / edgeA.length_squared());  // Edge projection
  }

  if (L3 >= 0.0 && L10 <= 0.0 && L4 <= 0.0) {
    return (node(0) + (stk::math::Dot(edgeC, point0) * edgeC) / edgeC.length_squared());  // Edge projection
  }

  if (L5 >= 0.0 && L11 <= 0.0 && L6 <= 0.0) {
    return (node(1) + (stk::math::Dot(edgeB, point1) * edgeB) / edgeB.length_squared());  // Edge projection
  }

  if (L7 >= 0.0 && L12 <= 0.0 && L8 <= 0.0) {
    return (node(2) + (stk::math::Dot(edgeD, point2) * edgeD) / edgeD.length_squared());  // Edge projection
  }

  return stk::math::Vector3d();
}


}}

