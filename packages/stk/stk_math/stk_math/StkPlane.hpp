/*--------------------------------------------------------------------*/
/*    Copyright 2008 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef STK_MATH_PLANE_HPP
#define STK_MATH_PLANE_HPP

#include <cmath>
#include <limits>
#include <stk_math/StkVector.hpp>

namespace stk {
namespace math {

///
/// The plane is represented as Dot(N,X) = c where N is a unit-length
/// normal vector, c is the plane constant, and X is any point on the
/// plane.
///

template<class REAL>
class Plane3 {
 public:
  typedef REAL Real;
  Plane3() { m_normal.set_invalid(); m_constant = std::nan(""); }

  /// N is specified, c = Dot(N,P) where P is on the plane, N need not be a unit-vector
  Plane3 (const Vec<REAL,3>& normal_dir, const Vec<REAL,3>& P) : m_normal(normal_dir) { m_normal.unitize(); m_constant = Dot(m_normal,P); }

  /// N = Cross(P1-P0,P2-P0)/Length(Cross(P1-P0,P2-P0)), c = Dot(N,P0) where
  /// P0, P1, P2 are points on the plane.
  Plane3 (const Vec<REAL,3>& P0, const Vec<REAL,3>& P1, const Vec<REAL,3>& P2) { m_normal = Cross(P1 - P0, P2 - P0); m_normal.unitize(); m_constant = Dot(m_normal,P0); }

  bool is_valid() const
  {
      return m_normal.is_valid() && !std::isnan(m_constant);
  }

  /// Compute d = Dot(N,P)-c where N is the plane normal and c is the plane
  /// constant.  This is a signed distance.  The sign of the return value is
  /// positive if the point is on the side in which the normal vector points.
  Real signed_distance (const Vec<REAL,3>& P) const { return Dot(m_normal,P) - m_constant; }

  const Vec<REAL,3> & normal() const { return m_normal; }
  Vec<REAL,3> & normal() { return m_normal; }

  const Real & constant() const { return const_cast< const Real &>(m_constant); }
  Real & constant() { return m_constant; }

  bool intersects_segment(const std::array<Vec<REAL,3>,2> & segment_nodes, REAL & location) const;

private:
  Vec<REAL,3>  m_normal;
  Real   m_constant;
};

template<class REAL>
std::ostream& operator<<( std::ostream& out, const Plane3<REAL>& plane )
{
  out << "Plane3: Normal=" << plane.normal()[0] << ", " << plane.normal()[1] << ", " << plane.normal()[2] << ", constant= " << plane.constant();
  return out;
}

template<class REAL>
inline bool Plane3<REAL>::intersects_segment(const std::array<Vec<REAL,3>,2> & segment_nodes, REAL & location) const
{
  const Real s[2] = { signed_distance(segment_nodes[0]), signed_distance(segment_nodes[1]) };
  //
  // Above and below the plane
  // or One end in plane and the other above the plane
  if ( (s[0] <= 0.0 && s[1] > 0.0) ||
       (s[0] > 0.0 && s[1] <= 0.0 ) ) {
    Real const tolerance = 10000.0 * std::numeric_limits<Real>::denorm_min();
    Real diff = s[1] - s[0];
    if ( diff < tolerance && diff > -tolerance ) {
      location = 0.5; // avoid NaN when edge is nearly in the plane (cut in the middle)
    } else {
      location = s[1] / diff;
    }
    return true;
  }
  // Both above or below the plane
  // One end in plane and the other below the plane
  // Both in the plane
  location = std::nan("");
  return false;
}

typedef Plane3<double> Plane3d;

}
}

#endif
