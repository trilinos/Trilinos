// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_Cutting_Surface.hpp>
#include <Akri_DiagWriter.hpp>

#include <math.h>
#include <iostream>

namespace krino{

int
Plane_Cutting_Surface::sign_at_position(const stk::math::Vector3d & p_coords) const
{
  return my_plane.signed_distance(p_coords) < 0. ? -1 : 1;
}

bool
Plane_Cutting_Surface::on_surface(const stk::math::Vector3d & p_coords, const double tol) const
{
  return std::abs(my_plane.signed_distance(p_coords)) < tol;
}

double
Plane_Cutting_Surface::interface_crossing_position(const std::array<stk::math::Vector3d,2> & edgeNodeCoords) const
{
  const double signed_dist_node0 = my_plane.signed_distance(edgeNodeCoords[0]);
  const double signed_dist_node1 = my_plane.signed_distance(edgeNodeCoords[1]);

  if((signed_dist_node0 < 0) == (signed_dist_node1 < 0))
  {
    std::stringstream str;
    str << "Failed to find intersection of plane " << my_plane << " with segment with nodes " << edgeNodeCoords[0] << " and " << edgeNodeCoords[1];
    throw std::runtime_error(str.str());
  }

  const double pos = signed_dist_node0 / (signed_dist_node0-signed_dist_node1);

  return pos;
}

std::string Plane_Cutting_Surface::print() const
{
  std::ostringstream os;
  os << "Plane cutting surface, plane = " << my_plane;
  return os.str();
}

Intersecting_Planes_Cutting_Surface::Intersecting_Planes_Cutting_Surface(const stk::math::Vector3d & p0, const stk::math::Vector3d & p1, const stk::math::Vector3d & p2, const stk::math::Vector3d & p3)
: Cutting_Surface()
{
  // Points p0 and p2 lie on the line common to the two planes
  my_plane0.set_from_most_orthogonal_angle_of_triangle(p0,p1,p2);
  my_plane1.set_from_most_orthogonal_angle_of_triangle(p2,p3,p0);

  myTriangleForPlane0IsLarger = Cross(p2-p0, p1-p0).length_squared() > Cross(p2-p0, p3-p0).length_squared();

  // Need to keep track of whether the normals to the plane point "toward" each other (positive_dihedral = true)
  // or "away" from each other (positive_dihedral = false).
  // We check this by seeing if p1 and p2 are on the positive side of the other plane.
  const double dist1 = my_plane1.signed_distance(p1);
  const double dist3 = my_plane0.signed_distance(p3);
  my_positive_dihedral = (std::abs(dist1) > std::abs(dist3)) ? (dist1>0) : (dist3>0);
}

int
Intersecting_Planes_Cutting_Surface::sign_at_position(const stk::math::Vector3d & p_coords) const
{
  const double signed_dist0 = my_plane0.signed_distance(p_coords);
  const double signed_dist1 = my_plane1.signed_distance(p_coords);
  if (my_positive_dihedral)
  {
    const double min = std::min(signed_dist0,signed_dist1);
    return min < 0. ? -1 : 1;
  }
  const double max = std::max(signed_dist0,signed_dist1);
  return max < 0. ? -1 : 1;
}

bool
Intersecting_Planes_Cutting_Surface::on_surface(const stk::math::Vector3d & p_coords, const double tol) const
{
  const double signed_dist0 = my_plane0.signed_distance(p_coords);
  const double signed_dist1 = my_plane1.signed_distance(p_coords);
  if (my_positive_dihedral)
  {
    return (signed_dist0 > -tol && std::abs(signed_dist1) < tol) ||
        (signed_dist1 > -tol && std::abs(signed_dist0) < tol);
  }
  return (signed_dist0 < tol && std::abs(signed_dist1) < tol) ||
        (signed_dist1 < tol && std::abs(signed_dist0) < tol);
}

double
Intersecting_Planes_Cutting_Surface::interface_crossing_position(const std::array<stk::math::Vector3d,2> & edgeNodeCoords) const
{
  const double signed_dist0_node0 = my_plane0.signed_distance(edgeNodeCoords[0]);
  const double signed_dist0_node1 = my_plane0.signed_distance(edgeNodeCoords[1]);
  const double signed_dist1_node0 = my_plane1.signed_distance(edgeNodeCoords[0]);
  const double signed_dist1_node1 = my_plane1.signed_distance(edgeNodeCoords[1]);

  double pos0 = signed_dist0_node0 / (signed_dist0_node0-signed_dist0_node1);
  double pos1 = signed_dist1_node0 / (signed_dist1_node0-signed_dist1_node1);

  const int sign_case_id =
      (((signed_dist0_node0 < 0) == my_positive_dihedral) ? 0 : 1) +
      (((signed_dist1_node0 < 0) == my_positive_dihedral) ? 0 : 2) +
      (((signed_dist0_node1 < 0) == my_positive_dihedral) ? 0 : 4) +
      (((signed_dist1_node1 < 0) == my_positive_dihedral) ? 0 : 8);
  static const int case_id_from_sign_case_id [] = {-1,-1,-1, 2,-1,-1,-1, 1,-1,-1,-1, 0, 3, 1, 0,-1};

  const int case_id = case_id_from_sign_case_id[sign_case_id];
  switch (case_id)
  {
    case 0: return pos0;
    case 1: return pos1;
    case 2: return (pos0 < pos1) ? pos0 : pos1;
    case 3: return (pos0 < pos1) ? pos1 : pos0;
  }

  // Either 0 or 2 true crossings.  Not ok.
  std::stringstream str;
  str << "Failed to find intersection of Intersecting_Planes_Cutting_Surface with planes "
      << my_plane0 << " and " << my_plane1
      << " with segment with nodes " << edgeNodeCoords[0] << " and " << edgeNodeCoords[1]
      << " with crossings " << pos0 << " and " << pos1 << " with diff " << pos0-pos1 << " sign_case_id " << sign_case_id;
  throw std::runtime_error(str.str());
}

std::string Intersecting_Planes_Cutting_Surface::print() const
{
  std::ostringstream os;
  os << "Intersecting planes cutting surface. Plane 0 = " << my_plane0 << " Plane 1 = " << my_plane1
     << " positive dihedral = " << std::boolalpha << my_positive_dihedral;
  return os.str();
}

std::string print_plane_visualization(const stk::math::Plane3d & plane, std::array<int,4> & geomIds, const std::string & description)
{
  std::ostringstream os;
  const stk::math::Vector3d x0 = plane.constant() * plane.normal();
  const stk::math::Vector3d x1 = x0 + plane.normal();
  os << "Create vertex " << x0[0] << " " << x0[1] << " " << x0[2] << std::endl;
  os << "Create vertex " << x1[0] << " " << x1[1] << " " << x1[2] << std::endl;
  os << "Create curve vertex " << ++geomIds[0] << " " << ++geomIds[0] << std::endl;
  os << "Create planar surface curve " << ++geomIds[1] << " distance 0. # Plane " << ++geomIds[2] << ": " << description << std::endl;
  geomIds[0] += 4; // 4 vertices created for plane
  geomIds[1] += 4; // 4 curves created for plane
  return os.str();
}

std::string Plane_Cutting_Surface::print_visualization(std::array<int,4> & geomIds, const std::string & description) const
{
  return print_plane_visualization(my_plane, geomIds, description);
}

std::string Intersecting_Planes_Cutting_Surface::print_visualization(std::array<int,4> & geomIds, const std::string & description) const
{
  const std::string viz0 = print_plane_visualization(my_plane0, geomIds, description+", plane 1/2");
  const std::string viz1 = print_plane_visualization(my_plane1, geomIds, description+", plane 2/2");
  return viz0 + viz1;
}

} // namespace krino
