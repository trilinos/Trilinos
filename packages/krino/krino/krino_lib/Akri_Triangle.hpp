// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_Triangle_h
#define Akri_Triangle_h

#include <iostream>
#include <assert.h>
#include <cmath>
#include <array>

#include <Akri_Vec.hpp>

namespace krino {

enum class ProjectionType{NODE, EDGE, FACE, NUM_PROJ_TYPE};

template<class REAL>
class Triangle3 {
 public:
  typedef REAL Real;


  Triangle3() {} /// Default, all coords zero
  Triangle3(const Vec<REAL,3> &n0, const Vec<REAL,3> &n1, const Vec<REAL,3> &n2) : nodes{{n0, n1, n2}} {} /// Explicit set three coordinates
  explicit Triangle3(const MemberInit) {/* MemberInit::NONE */} /// Special non-initialize constructor (for performance)

  Vec<REAL,3> normal() const { return normal_dir().unit_vector(); } /// Unit vector
  Vec<REAL,3> normal_dir() const {return Cross(GetNode(1)-GetNode(0),GetNode(2)-GetNode(0)); } /// Non-unit normal (faster)
  Real area() const { return 0.5*normal_dir().length(); }

  Vec<REAL,3> ParametricToRealCoords(const Vec<REAL,2> & Param) const {
    // Avoids temporary constructions
    return Vec<REAL,3>(
        GetNode(0)[0]*(1.0-Param[0]-Param[1]) + GetNode(1)[0]*Param[0] + GetNode(2)[0]*Param[1],
        GetNode(0)[1]*(1.0-Param[0]-Param[1]) + GetNode(1)[1]*Param[0] + GetNode(2)[1]*Param[1],
        GetNode(0)[2]*(1.0-Param[0]-Param[1]) + GetNode(1)[2]*Param[0] + GetNode(2)[2]*Param[1]
        );
  }

  ///  Compute the closest point on the triangle to an input point.  Return the type of projection, i.e.,
  ///  is the projected point on a node, edge, surface of the triangle.  Optionally calculate parametric coordinates.
  template<class T = std::nullptr_t> ProjectionType ClosestPoint(const Vec<REAL,3>& p, Vec<REAL,3>& ClosestPt, T ParamPt = nullptr) const;

  Real DistanceSquared(const Vec<REAL,3>& P, Vec<REAL,2>& ParamPt) const;
  Real DistanceSquared(const Vec<REAL,3>& P) const;

  const Vec<REAL,3>& GetNode(const int index) const { assert(index >= 0 && index < 3); return nodes[index]; }
  Vec<REAL,3>& GetNode(const int index) { assert(index >= 0 && index < 3); return nodes[index]; }

  void SetNodes(const Vec<REAL,3>& n0, const Vec<REAL,3>& n1, const Vec<REAL,3>& n2) { nodes = {{n0, n1, n2}}; } // Change triangle coordinates
private:
  std::array<Vec<REAL,3>,3> nodes;
};

namespace detail {
template <class REAL, class T>
struct assign_parametric_coords {
  static void apply(const REAL & x, const REAL & y, T parametric_coords)
  {
    static_assert(std::is_pointer<T>::value, "Expecting pointer");
    parametric_coords[0] = x;
    parametric_coords[1] = y;
  }
};
template <class REAL>
struct assign_parametric_coords<REAL, std::nullptr_t>  {
  static void apply(const REAL, const REAL, std::nullptr_t) {}
};
}

template<class REAL>
inline REAL Triangle3<REAL>::DistanceSquared(const Vec<REAL,3>& P, Vec<REAL,2>& ParamPt) const
{
  Vec<REAL,3> ClosestPt(MemberInit::NONE);
  ClosestPoint(P, ClosestPt, ParamPt.data());
  return (P-ClosestPt).length_squared();
}

template<class REAL>
inline REAL Triangle3<REAL>::DistanceSquared(const Vec<REAL,3>& P) const
{
  Vec<REAL,3> ClosestPt(MemberInit::NONE);
  ClosestPoint(P, ClosestPt);
  return (P-ClosestPt).length_squared();
}

//  Adapted from closest face projection from "Real time Collision Detection" text, highly optimized
template<class REAL>
template<class T>
ProjectionType Triangle3<REAL>::ClosestPoint(const Vec<REAL,3>& p, Vec<REAL,3>& ClosestPt, T ParamPt) const
{
  const Vec<REAL,3>& a = GetNode(0);
  const Vec<REAL,3>& b = GetNode(1);
  const Vec<REAL,3>& c = GetNode(2);

  Vec<REAL,3> ab(b-a);
  Vec<REAL,3> ac(c-a);
  Vec<REAL,3> ap(p-a);

  Real d1(Dot(ab,ap));
  Real d2(Dot(ac,ap));

  if(d1 <= 0.0 && d2 <= 0.0) {
    ClosestPt = a;
    detail::assign_parametric_coords<REAL,T>::apply(0.0, 0.0, ParamPt);
    return ProjectionType::NODE;
   }

  Vec<REAL,3> bp(p-b);
  Real d3(Dot(ab,bp));
  Real d4(Dot(ac,bp));

  if(d3 >= 0.0 && d4  <= d3) {
    ClosestPt = b;
    detail::assign_parametric_coords<REAL,T>::apply(1.0, 0.0, ParamPt);
    return ProjectionType::NODE;
  }

  Real vc(d1*d4-d3*d2);
  if(vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
    Real v = d1/(d1-d3);
    ClosestPt = a + v * ab;
    detail::assign_parametric_coords<REAL,T>::apply(v, 0.0, ParamPt);
    return ProjectionType::EDGE;
  }

  Vec<REAL,3> cp(p-c);
  Real d5(Dot(ab,cp));
  Real d6(Dot(ac,cp));
  if(d6 >= 0.0 && d5 <= d6) {
    ClosestPt = c;
    detail::assign_parametric_coords<REAL,T>::apply(0.0, 1.0, ParamPt);
    return ProjectionType::NODE;
  }

  Real vb(d5*d2 - d1*d6);
  if(vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
    Real w(d2/(d2-d6));
    ClosestPt = a + w * ac;
    detail::assign_parametric_coords<REAL,T>::apply(0.0, w, ParamPt);
    return ProjectionType::EDGE;
  }

  Real va(d3*d6 - d5*d4);
  if(va <= 0.0 && (d4-d3) >= 0.0 && (d5-d6) >= 0.0) {
    Real w((d4-d3)/((d4-d3) + (d5-d6)));
    ClosestPt =  b + w * (c-b);
    detail::assign_parametric_coords<REAL,T>::apply(1.0-w, w, ParamPt);
    return ProjectionType::EDGE;
  }

  Real denom(1.0/(va+vb+vc));
  Real v(vb*denom);
  Real w(vc*denom);
  ClosestPt = a + ab*v + ac*w;
  detail::assign_parametric_coords<REAL,T>::apply(v, w, ParamPt);
  return ProjectionType::FACE;
}

typedef Triangle3<double> Triangle3d;

}
#endif // Akri_Triangle_h
