// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_CubitJournal.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace krino {

CubitJournal::CubitJournal(const std::string & title, const int firstBlockId)
: myTitle(title),
  myFirstBlockId(firstBlockId),
  myCurrentGroup(0)
{
  myGroups.push_back(Group{"default", "", {}, {}, {}, {}, {}});
}

void CubitJournal::set_group(const std::string & name, const std::string & color)
{
  auto iter = std::find_if(myGroups.begin(), myGroups.end(), [&](const Group & group) { return group.name == name; });
  if (iter == myGroups.end())
  {
    myGroups.push_back(Group{name, color, {}, {}, {}, {}, {}});
    myCurrentGroup = myGroups.size()-1;
    return;
  }
  myCurrentGroup = std::distance(myGroups.begin(), iter);
  if (!color.empty())
    iter->color = color;
}

void CubitJournal::add_comment(const std::string & text)
{
  current_group().comments.push_back(text);
}

int CubitJournal::node_index(const stk::math::Vector3d & x)
{
  const auto result = myNodeIndices.emplace(std::array<double,3>{x[0], x[1], x[2]}, myNodes.size());
  if (result.second)
    myNodes.push_back(x);
  return result.first->second;
}

void CubitJournal::add_point(const stk::math::Vector3d & x, const std::string & label)
{
  current_group().points.emplace_back(x, label);
}

void CubitJournal::add_segment(const stk::math::Vector3d & x0, const stk::math::Vector3d & x1)
{
  current_group().segments.push_back({{node_index(x0), node_index(x1)}});
}

void CubitJournal::add_tri(const std::array<stk::math::Vector3d,3> & x)
{
  current_group().tris.push_back({{node_index(x[0]), node_index(x[1]), node_index(x[2])}});
}

void CubitJournal::add_tet(const std::array<stk::math::Vector3d,4> & x)
{
  current_group().tets.push_back({{node_index(x[0]), node_index(x[1]), node_index(x[2]), node_index(x[3])}});
}

void CubitJournal::add_convex_polygon(const std::vector<stk::math::Vector3d> & x)
{
  for (size_t i=1; i+1<x.size(); ++i)
    add_tri({{x[0], x[i], x[i+1]}});
}

void CubitJournal::add_section(const std::vector<stk::math::Vector3d> & section)
{
  if (section.size() == 1)
    add_point(section[0]);
  else if (section.size() == 2)
    add_segment(section[0], section[1]);
  else
    add_convex_polygon(section);
}

void CubitJournal::add_plane_in_tet(const stk::math::Plane3d & plane, const std::array<stk::math::Vector3d,4> & tetCoords)
{
  add_section(section_of_tet_by_plane(plane, tetCoords));
}

void CubitJournal::add_plane_in_tri(const stk::math::Plane3d & plane, const std::array<stk::math::Vector3d,3> & triCoords)
{
  add_section(section_of_tri_by_plane(plane, triCoords));
}

void CubitJournal::add_plane_in_box(const stk::math::Plane3d & plane, const stk::math::Vector3d & boxMin, const stk::math::Vector3d & boxMax)
{
  add_section(section_of_box_by_plane(plane, boxMin, boxMax));
}

static void write_comment_lines(std::ostream & os, const std::string & text)
{
  std::istringstream lines(text);
  std::string line;
  while (std::getline(lines, line))
    os << "# " << line << "\n";
}

static std::string journal_name(const std::string & groupName, const std::string & suffix)
{
  std::string name = groupName + "_" + suffix;
  std::replace_if(name.begin(), name.end(), [](const char c) { return !(std::isalnum(static_cast<unsigned char>(c)) || c == '_'); }, '_');
  return name;
}

static std::string node_ref(const int index)
{
  return "{krinoNode0+" + std::to_string(index+1) + "}";
}

static std::string vertex_ref(const int index)
{
  return "{krinoVertex0+" + std::to_string(index+1) + "}";
}

template<size_t N>
static void write_block_of_elements(std::ostream & os,
    const std::vector<std::array<int,N>> & elements,
    const std::string & cubitType,
    const std::string & groupName,
    const std::string & color,
    int & blockId)
{
  if (elements.empty()) return;
  for (auto && element : elements)
  {
    os << "create " << cubitType << " node";
    for (int index : element)
      os << " " << node_ref(index);
    os << "\n";
  }
  const std::string idFn = "Id(\"" + cubitType + "\")";
  os << "block " << blockId << " " << cubitType << " ";
  if (elements.size() == 1)
    os << "{" << idFn << "}\n";
  else
    os << "{" << idFn << "-" << elements.size()-1 << "} to {" << idFn << "}\n";
  os << "block " << blockId << " name \"" << journal_name(groupName, cubitType) << "\"\n";
  if (!color.empty())
    os << "color block " << blockId << " " << color << "\n";
  ++blockId;
}

std::string CubitJournal::str() const
{
  std::ostringstream os;
  os << std::setprecision(17);
  os << "# Cubit journal written by krino";
  if (!myTitle.empty())
    os << ": " << myTitle;
  os << "\n";
  os << "#{krinoNode0=Id(\"node\")}\n";
  os << "#{krinoVertex0=Id(\"vertex\")}\n";

  for (auto && x : myNodes)
    os << "create node location " << x[0] << " " << x[1] << " " << x[2] << "\n";

  int blockId = myFirstBlockId;
  int numVerticesWritten = 0;
  for (auto && group : myGroups)
  {
    const bool isEmpty = group.comments.empty() && group.points.empty() && group.segments.empty() && group.tris.empty() && group.tets.empty();
    if (isEmpty) continue;
    os << "\n# Group " << group.name << "\n";
    for (auto && comment : group.comments)
      write_comment_lines(os, comment);

    // Edges first: Cubit makes edges for tris and tets, which would break up the range of new edges.
    write_block_of_elements(os, group.segments, "edge", group.name, group.color, blockId);
    write_block_of_elements(os, group.tris, "tri", group.name, group.color, blockId);
    write_block_of_elements(os, group.tets, "tet", group.name, group.color, blockId);

    if (!group.points.empty())
    {
      for (auto && [x, label] : group.points)
      {
        os << "create vertex " << x[0] << " " << x[1] << " " << x[2];
        if (!label.empty())
          os << " # " << label;
        os << "\n";
      }
      const std::string firstRef = vertex_ref(numVerticesWritten);
      numVerticesWritten += group.points.size();
      const std::string lastRef = vertex_ref(numVerticesWritten-1);
      const std::string range = (group.points.size() == 1) ? firstRef : (firstRef + " to " + lastRef);
      os << "group \"" << journal_name(group.name, "points") << "\" add vertex " << range << "\n";
      if (!group.color.empty())
        os << "color vertex " << range << " " << group.color << "\n";
    }
  }
  return os.str();
}

void CubitJournal::write(const std::string & fileName) const
{
  std::ofstream file(fileName);
  STK_ThrowRequireMsg(file.good(), "Could not open " << fileName << " to write a Cubit journal.");
  file << str();
}

static stk::math::Vector3d point_in_plane_orthogonal_to(const stk::math::Vector3d & normal)
{
  // The axis least aligned with the normal, made orthogonal to it.
  int axis = 0;
  for (int i=1; i<3; ++i)
    if (std::abs(normal[i]) < std::abs(normal[axis]))
      axis = i;
  stk::math::Vector3d direction(0.,0.,0.);
  direction[axis] = 1.;
  return (direction - Dot(direction, normal)*normal).unit_vector();
}

static void sort_around_normal(std::vector<stk::math::Vector3d> & points, const stk::math::Vector3d & normal)
{
  if (points.size() < 3) return;
  stk::math::Vector3d centroid(0.,0.,0.);
  for (auto && x : points)
    centroid += x;
  centroid /= points.size();
  const stk::math::Vector3d u = point_in_plane_orthogonal_to(normal);
  const stk::math::Vector3d v = Cross(normal, u);
  std::sort(points.begin(), points.end(), [&](const stk::math::Vector3d & a, const stk::math::Vector3d & b)
  {
    return std::atan2(Dot(a-centroid, v), Dot(a-centroid, u)) < std::atan2(Dot(b-centroid, v), Dot(b-centroid, u));
  });
}

template<size_t NVERTS, size_t NEDGES>
static std::vector<stk::math::Vector3d> section_of_convex_polytope_by_plane(const stk::math::Plane3d & plane,
    const std::array<stk::math::Vector3d,NVERTS> & vertices,
    const std::array<std::array<int,2>,NEDGES> & edges)
{
  std::array<double,NVERTS> dist;
  for (size_t i=0; i<NVERTS; ++i)
    dist[i] = plane.signed_distance(vertices[i]);

  std::vector<stk::math::Vector3d> section;
  for (size_t i=0; i<NVERTS; ++i)
    if (dist[i] == 0.)
      section.push_back(vertices[i]);
  for (auto && edge : edges)
  {
    const double d0 = dist[edge[0]];
    const double d1 = dist[edge[1]];
    if ((d0 < 0. && d1 > 0.) || (d0 > 0. && d1 < 0.))
    {
      const double pos = d0/(d0-d1);
      section.push_back((1.-pos)*vertices[edge[0]] + pos*vertices[edge[1]]);
    }
  }
  sort_around_normal(section, plane.normal());
  return section;
}

std::vector<stk::math::Vector3d> section_of_tet_by_plane(const stk::math::Plane3d & plane, const std::array<stk::math::Vector3d,4> & tetCoords)
{
  static const std::array<std::array<int,2>,6> tetEdges{{ {{0,1}}, {{1,2}}, {{2,0}}, {{0,3}}, {{1,3}}, {{2,3}} }};
  return section_of_convex_polytope_by_plane(plane, tetCoords, tetEdges);
}

std::vector<stk::math::Vector3d> section_of_tri_by_plane(const stk::math::Plane3d & plane, const std::array<stk::math::Vector3d,3> & triCoords)
{
  static const std::array<std::array<int,2>,3> triEdges{{ {{0,1}}, {{1,2}}, {{2,0}} }};
  return section_of_convex_polytope_by_plane(plane, triCoords, triEdges);
}

std::vector<stk::math::Vector3d> section_of_box_by_plane(const stk::math::Plane3d & plane, const stk::math::Vector3d & boxMin, const stk::math::Vector3d & boxMax)
{
  std::array<stk::math::Vector3d,8> corners;
  for (int i=0; i<8; ++i)
    corners[i] = stk::math::Vector3d((i&1) ? boxMax[0] : boxMin[0], (i&2) ? boxMax[1] : boxMin[1], (i&4) ? boxMax[2] : boxMin[2]);
  static const std::array<std::array<int,2>,12> boxEdges{{
    {{0,1}}, {{2,3}}, {{4,5}}, {{6,7}},
    {{0,2}}, {{1,3}}, {{4,6}}, {{5,7}},
    {{0,4}}, {{1,5}}, {{2,6}}, {{3,7}} }};
  return section_of_convex_polytope_by_plane(plane, corners, boxEdges);
}

static double roundoff_tolerance_for_section(const std::vector<stk::math::Vector3d> & section)
{
  double extent = 0.;
  for (auto && x : section)
    extent = std::max(extent, (x-section[0]).length());
  return 1.e-12*extent;
}

std::vector<stk::math::Vector3d> clip_section_to_positive_side_of_plane(const std::vector<stk::math::Vector3d> & section, const stk::math::Plane3d & plane)
{
  // Points within roundoff of the plane are on it, so that a section in the plane is kept.
  const double tol = roundoff_tolerance_for_section(section);
  std::vector<stk::math::Vector3d> clipped;
  const size_t numPoints = section.size();
  if (numPoints < 3)
  {
    // A point or a segment: keep the points on the positive side, and the crossing.
    for (auto && x : section)
      if (plane.signed_distance(x) >= -tol)
        clipped.push_back(x);
    if (numPoints == 2 && clipped.size() == 1)
    {
      const double d0 = plane.signed_distance(section[0]);
      const double d1 = plane.signed_distance(section[1]);
      if (std::abs(d0) > tol && std::abs(d1) > tol)
      {
        const double pos = d0/(d0-d1);
        clipped.push_back((1.-pos)*section[0] + pos*section[1]);
      }
    }
    return clipped;
  }

  // Sutherland-Hodgman, for one plane.
  for (size_t i=0; i<numPoints; ++i)
  {
    const stk::math::Vector3d & x0 = section[i];
    const stk::math::Vector3d & x1 = section[(i+1)%numPoints];
    const double d0 = plane.signed_distance(x0);
    const double d1 = plane.signed_distance(x1);
    if (d0 >= -tol)
      clipped.push_back(x0);
    if ((d0 < -tol && d1 > tol) || (d0 > tol && d1 < -tol))
    {
      const double pos = d0/(d0-d1);
      clipped.push_back((1.-pos)*x0 + pos*x1);
    }
  }
  return clipped;
}


static stk::math::Vector3d gradient_of_linear_interpolant_in_tet(const std::array<stk::math::Vector3d,4> & x, const std::array<double,4> & f)
{
  const stk::math::Vector3d e1 = x[1]-x[0];
  const stk::math::Vector3d e2 = x[2]-x[0];
  const stk::math::Vector3d e3 = x[3]-x[0];
  const double det = Dot(e1, Cross(e2, e3));
  if (det == 0.) return stk::math::Vector3d(0.,0.,0.);
  return ((f[1]-f[0])*Cross(e2, e3) + (f[2]-f[0])*Cross(e3, e1) + (f[3]-f[0])*Cross(e1, e2))/det;
}

template<size_t NVERTS>
static std::vector<stk::math::Vector3d> zero_contour_of_linear_interpolant_in_simplex(
    const std::array<size_t,NVERTS> & gridIndices,
    const std::array<stk::math::Vector3d,NVERTS> & x,
    const std::array<double,NVERTS> & f)
{
  std::vector<stk::math::Vector3d> contour;
  for (size_t i=0; i<NVERTS; ++i)
    if (f[i] == 0.)
      contour.push_back(x[i]);
  for (size_t i=0; i<NVERTS; ++i)
  {
    for (size_t j=i+1; j<NVERTS; ++j)
    {
      if (!((f[i] < 0. && f[j] > 0.) || (f[i] > 0. && f[j] < 0.))) continue;
      // From the lower grid index, so that tets sharing the edge find the same point.
      const size_t a = (gridIndices[i] < gridIndices[j]) ? i : j;
      const size_t b = (a == i) ? j : i;
      const double pos = f[a]/(f[a]-f[b]);
      contour.push_back((1.-pos)*x[a] + pos*x[b]);
    }
  }
  return contour;
}

void add_zero_contour_in_box(CubitJournal & journal,
    const std::function<double(const stk::math::Vector3d &)> & fn,
    const stk::math::Vector3d & boxMin,
    const stk::math::Vector3d & boxMax,
    const int numCellsPerSide)
{
  const bool is2D = (boxMin[2] == boxMax[2]);
  const int n = numCellsPerSide;
  const int nz = is2D ? 0 : n;
  const size_t numPointsPerSide = static_cast<size_t>(n) + 1;
  auto grid_index = [&](const int i, const int j, const int k)
  { return static_cast<size_t>(i) + numPointsPerSide*(static_cast<size_t>(j) + numPointsPerSide*static_cast<size_t>(k)); };
  std::vector<stk::math::Vector3d> x(numPointsPerSide*numPointsPerSide*(static_cast<size_t>(nz) + 1));
  std::vector<double> f(x.size());
  for (int k=0; k<=nz; ++k)
    for (int j=0; j<=n; ++j)
      for (int i=0; i<=n; ++i)
      {
        const size_t index = grid_index(i, j, k);
        x[index] = stk::math::Vector3d(
            boxMin[0] + (boxMax[0]-boxMin[0])*i/n,
            boxMin[1] + (boxMax[1]-boxMin[1])*j/n,
            is2D ? boxMin[2] : (boxMin[2] + (boxMax[2]-boxMin[2])*k/n));
        f[index] = fn(x[index]);
      }

  if (is2D)
  {
    // Each square as two triangles along its diagonal from corner 0 to corner 3.
    static const std::array<std::array<int,3>,2> squareTris{{ {{0,1,3}}, {{0,3,2}} }};
    for (int j=0; j<n; ++j)
      for (int i=0; i<n; ++i)
        for (auto && tri : squareTris)
        {
          std::array<size_t,3> indices;
          for (int v=0; v<3; ++v)
            indices[v] = grid_index(i + (tri[v]&1), j + ((tri[v]>>1)&1), 0);
          const std::vector<stk::math::Vector3d> contour = zero_contour_of_linear_interpolant_in_simplex<3>(
              indices, {{x[indices[0]], x[indices[1]], x[indices[2]]}}, {{f[indices[0]], f[indices[1]], f[indices[2]]}});
          if (contour.size() == 2)
            journal.add_segment(contour[0], contour[1]);
        }
    return;
  }

  // Each cube as the 6 Kuhn tets along its diagonal from corner 0 to corner 7 (corner bits are
  // x, y, z), which match across neighboring cubes.
  static const std::array<std::array<int,4>,6> cubeTets{{
    {{0,1,3,7}}, {{0,1,5,7}}, {{0,2,3,7}}, {{0,2,6,7}}, {{0,4,5,7}}, {{0,4,6,7}} }};
  for (int k=0; k<n; ++k)
    for (int j=0; j<n; ++j)
      for (int i=0; i<n; ++i)
        for (auto && tet : cubeTets)
        {
          std::array<size_t,4> indices;
          for (int v=0; v<4; ++v)
            indices[v] = grid_index(i + (tet[v]&1), j + ((tet[v]>>1)&1), k + ((tet[v]>>2)&1));
          const std::array<stk::math::Vector3d,4> tetX{{x[indices[0]], x[indices[1]], x[indices[2]], x[indices[3]]}};
          const std::array<double,4> tetF{{f[indices[0]], f[indices[1]], f[indices[2]], f[indices[3]]}};
          std::vector<stk::math::Vector3d> contour = zero_contour_of_linear_interpolant_in_simplex<4>(indices, tetX, tetF);
          if (contour.size() < 3) continue;
          const stk::math::Vector3d gradient = gradient_of_linear_interpolant_in_tet(tetX, tetF);
          if (gradient.length() == 0.) continue;
          sort_around_normal(contour, gradient.unit_vector());
          journal.add_convex_polygon(contour);
        }
}

}
