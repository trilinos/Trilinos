// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_GEOMETRY_AKRI_CUBITJOURNAL_HPP_
#define KRINO_KRINO_GEOMETRY_AKRI_CUBITJOURNAL_HPP_

#include <stk_math/StkPlane.hpp>
#include <stk_math/StkVector.hpp>
#include <array>
#include <functional>
#include <map>
#include <string>
#include <vector>

namespace krino {

// Collects points, segments, tris, tets and sections of planes for debugging, and writes them as a
// Cubit journal (play it back with "playback 'file.jou'", or "cubit -input file.jou").
// - Tets, tris, segments and polygons are mesh entities. Each group gets one block per element type,
//   named <group>_tet, <group>_tri, <group>_edge, so the groups survive "export mesh".
// - Points are geometry vertices, in the Cubit group <group>_points.
// - Nodes are shared by exact coordinates. Node and vertex ids are relative to the ids in use when the
//   journal is played, so several journals can be played into one session; give each its own
//   firstBlockId, as Cubit has no Id("block").
// - Coordinates are written with 17 significant digits. 2D coordinates have z = 0.
class CubitJournal
{
public:
  explicit CubitJournal(const std::string & title = "", const int firstBlockId = 1);

  // Later entities go in this group, made if new. A color is a Cubit color name, applied to the
  // group's blocks and points. Before any call, entities go in the group "default".
  void set_group(const std::string & name, const std::string & color = "");
  // A comment, written at the start of the current group. Each line is prefixed with "# ".
  void add_comment(const std::string & text);

  void add_point(const stk::math::Vector3d & x, const std::string & label = "");
  void add_segment(const stk::math::Vector3d & x0, const stk::math::Vector3d & x1);
  void add_tri(const std::array<stk::math::Vector3d,3> & x);
  void add_tet(const std::array<stk::math::Vector3d,4> & x);
  // A convex polygon with its vertices in order around it, as a fan of tris.
  void add_convex_polygon(const std::vector<stk::math::Vector3d> & x);
  // A section from the section_of_* functions below: a point, a segment or a convex polygon. Empty
  // sections add nothing.
  void add_section(const std::vector<stk::math::Vector3d> & section);

  void add_plane_in_tet(const stk::math::Plane3d & plane, const std::array<stk::math::Vector3d,4> & tetCoords);
  void add_plane_in_tri(const stk::math::Plane3d & plane, const std::array<stk::math::Vector3d,3> & triCoords);
  void add_plane_in_box(const stk::math::Plane3d & plane, const stk::math::Vector3d & boxMin, const stk::math::Vector3d & boxMax);

  std::string str() const;
  void write(const std::string & fileName) const;

private:
  struct Group
  {
    std::string name;
    std::string color;
    std::vector<std::string> comments;
    std::vector<std::pair<stk::math::Vector3d,std::string>> points; // location, label
    std::vector<std::array<int,2>> segments;
    std::vector<std::array<int,3>> tris;
    std::vector<std::array<int,4>> tets;
  };
  int node_index(const stk::math::Vector3d & x);
  Group & current_group() { return myGroups[myCurrentGroup]; }

  std::string myTitle;
  int myFirstBlockId;
  std::vector<stk::math::Vector3d> myNodes;
  std::map<std::array<double,3>,int> myNodeIndices;
  std::vector<Group> myGroups;
  size_t myCurrentGroup;
};

// The section of a tet, a tri or an axis-aligned box by a plane: its crossings of their edges, with
// the vertices on the plane, in order around the plane's normal. For a tet or a box it is a convex
// polygon, for a tri (in 2D, a plane with no z component is a line) a segment; it can also be a single
// point, or empty when the plane misses.
std::vector<stk::math::Vector3d> section_of_tet_by_plane(const stk::math::Plane3d & plane, const std::array<stk::math::Vector3d,4> & tetCoords);
std::vector<stk::math::Vector3d> section_of_tri_by_plane(const stk::math::Plane3d & plane, const std::array<stk::math::Vector3d,3> & triCoords);
std::vector<stk::math::Vector3d> section_of_box_by_plane(const stk::math::Plane3d & plane, const stk::math::Vector3d & boxMin, const stk::math::Vector3d & boxMax);

// The part of a section (point, segment or convex polygon in order) on the positive side of a plane,
// or on it to within roundoff of the section's size.
std::vector<stk::math::Vector3d> clip_section_to_positive_side_of_plane(const std::vector<stk::math::Vector3d> & section, const stk::math::Plane3d & plane);

// The zero contour of a function in a box, approximated on a grid of numCellsPerSide cells a side:
// each cell is split into tets (triangles if the box has no extent in z), and the contour of the
// function's linear interpolant in each is added as a polygon (a segment in 2D).  Crossings on edges
// shared by tets are computed the same way in each, so the polygons share their nodes.
void add_zero_contour_in_box(CubitJournal & journal,
    const std::function<double(const stk::math::Vector3d &)> & fn,
    const stk::math::Vector3d & boxMin,
    const stk::math::Vector3d & boxMax,
    const int numCellsPerSide);

}

#endif /* KRINO_KRINO_GEOMETRY_AKRI_CUBITJOURNAL_HPP_ */
