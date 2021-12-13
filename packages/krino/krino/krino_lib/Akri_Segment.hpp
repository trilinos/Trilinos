// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_Segment_h
#define Akri_Segment_h

#include <Akri_Vec.hpp>
#include <cmath>
#include <iostream>
#include <array>

namespace krino {

template<class REAL>
class Segment3 {
 public:
  typedef REAL Real;
  Segment3() : nodes{} {}  ///< Default all coords zero
  Segment3(const Vec<REAL,3>& n0, const Vec<REAL,3>& n1) : nodes{{n0, n1}} {} ///< Explicit set two end points
  Segment3(MemberInit type) {/* type == NONE */}

  Real Length() const { return (nodes[1]-nodes[0]).length(); }
  Real LengthSquared() const { return (nodes[1]-nodes[0]).length_squared(); }

  const Vec<REAL,3>& GetNode(const int index) const { return nodes[index]; }
  Vec<REAL,3>& GetNode(const int index) { return nodes[index]; }

  const std::array<Vec<REAL,3>,2>& GetNodes() const { return nodes; }
  std::array<Vec<REAL,3>,2>& GetNodes() { return nodes; }

  //
  //  Find the closest point projection between point and the face.
  //
  Real closest_projection(const Vec<REAL,3> &point) const {
    Vec<REAL,3> edge_dir(nodes[1]-nodes[0]);
    Real dotValA = Dot(edge_dir,(point-nodes[0]));
    if(dotValA <= 0.0) { return 0.0; }
    Real dotValB = Dot(edge_dir,(point-nodes[1]));
    if(dotValB >= 0.0) { return 1.0; }
    Real lenSquared = edge_dir.length_squared();
    if(lenSquared == 0.0) { return 0.0; }
    return dotValA / lenSquared;
  }
  Real closest_projection(const Vec<REAL,3> &point, Vec<REAL,3> &proj_point) const {
    const Real location = closest_projection(point);
    proj_point = (1.0-location)*GetNode(0) + location*GetNode(1);
    return location;
  }

  Real DistanceSquared(const Vec<REAL,3>& x, Real & parametric_coord) const {
    Vector3d closest_pt(MemberInit::NONE);
    parametric_coord = closest_projection(x, closest_pt);
    return (x - closest_pt).length_squared();
  }
  Real DistanceSquared(const Vec<REAL,3>& x) const {
    Vector3d closest_pt(MemberInit::NONE);
    closest_projection(x, closest_pt);
    return (x - closest_pt).length_squared();
  }

  friend std::ostream& operator<<( std::ostream& out, const Segment3<REAL>& seg )
  {
    out << "Segment3:"
         << " Node0= " << seg.GetNode(0)[0] << ", " << seg.GetNode(0)[1] << ", " << seg.GetNode(0)[2]
         << " Node1= " << seg.GetNode(1)[0] << ", " << seg.GetNode(1)[1] << ", " << seg.GetNode(1)[2]
         << "\n";
    return out;
  }

 private:
  std::array<Vec<REAL,3>,2> nodes;
};

typedef Segment3<double> Segment3d;

} // namespace krino

#endif // Akri_Segment_h
