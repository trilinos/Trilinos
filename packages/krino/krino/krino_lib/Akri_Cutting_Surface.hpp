// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_Cutting_Surface_h
#define Akri_Cutting_Surface_h

#include <stk_math/StkPlane.hpp>
#include <Akri_Segment.hpp>
#include <stk_math/StkVector.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <string>

namespace krino {

class Cutting_Surface
{
public:
  Cutting_Surface() {}
  virtual ~Cutting_Surface() {}
  virtual bool on_surface(const stk::math::Vector3d & p_coords, const double tol) const = 0;
  virtual int sign_at_position(const stk::math::Vector3d & p_coords) const = 0;
  virtual double interface_crossing_position(const std::array<stk::math::Vector3d,2> & edgeNodeCoords) const = 0;
  virtual std::string print() const = 0;
  virtual std::string print_visualization(std::array<int,4> & geomIds, const std::string & description) const = 0;
  virtual const stk::math::Plane3d & get_plane() const = 0;
};

class Plane_Cutting_Surface : public Cutting_Surface
{
public:
  // for 3D from 3 points on plane
  Plane_Cutting_Surface(const stk::math::Vector3d & p0, const stk::math::Vector3d & p1, const stk::math::Vector3d & p2) : Cutting_Surface() {my_plane.set_from_most_orthogonal_angle_of_triangle(p0,p1,p2);}
  // for 2D or 3D from direction and point on plane
  Plane_Cutting_Surface(const stk::math::Vector3d & direction, const stk::math::Vector3d & p0) : my_plane(direction,p0) {}
  virtual ~Plane_Cutting_Surface() {}
  virtual bool on_surface(const stk::math::Vector3d & p_coords, const double tol) const override;
  int sign_at_position(const stk::math::Vector3d & p_coords) const override;
  double interface_crossing_position(const std::array<stk::math::Vector3d,2> & edgeNodeCoords) const override;
  std::string print() const override;
  virtual std::string print_visualization(std::array<int,4> & geomIds, const std::string & description) const override;
  virtual const stk::math::Plane3d & get_plane() const override { return my_plane; }
private:
  stk::math::Plane3d my_plane;
};

class Plane_Cutting_Surface_2D : public Plane_Cutting_Surface
{
public:
  // for 2D from 2 points on line
  Plane_Cutting_Surface_2D(const stk::math::Vector3d & p0, const stk::math::Vector3d & p1) : Plane_Cutting_Surface(stk::math::Vector3d(p1[1]-p0[1],p0[0]-p1[0],0),p0) {}
  virtual ~Plane_Cutting_Surface_2D() {}
private:
};


class Intersecting_Planes_Cutting_Surface : public Cutting_Surface
{
public:
  Intersecting_Planes_Cutting_Surface(const stk::math::Vector3d & p0, const stk::math::Vector3d & p1, const stk::math::Vector3d & p2, const stk::math::Vector3d & p3); // 3D
  virtual ~Intersecting_Planes_Cutting_Surface() {}
  virtual bool on_surface(const stk::math::Vector3d & p_coords, const double tol) const override;
  int sign_at_position(const stk::math::Vector3d & p_coords) const override;
  double interface_crossing_position(const std::array<stk::math::Vector3d,2> & edgeNodeCoords) const override;
  std::string print() const override;
  virtual std::string print_visualization(std::array<int,4> & geomIds, const std::string & description) const override;
  virtual const stk::math::Plane3d & get_plane() const override { /* THIS IS A BAD IDEA. SHOULD WE TEST PLANARITY? */ return myTriangleForPlane0IsLarger ? my_plane0 : my_plane1; }
private:
  bool myTriangleForPlane0IsLarger;
  stk::math::Plane3d my_plane0;
  stk::math::Plane3d my_plane1;
  bool my_positive_dihedral;
};

} // namespace krino

#endif // Akri_Cutting_Surface_h
