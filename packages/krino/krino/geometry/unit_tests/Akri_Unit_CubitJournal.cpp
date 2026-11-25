// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>
#include <Akri_CubitJournal.hpp>
#include <Akri_UnitCubitJournalUtils.hpp>
#include <stk_math/StkPlane.hpp>
#include <stk_math/StkVector.hpp>
#include <array>
#include <cmath>
#include <sstream>
#include <string>
#include <vector>

namespace krino {

static int count_journal_lines_starting_with(const std::string & journal, const std::string & start)
{
  std::istringstream lines(journal);
  std::string line;
  int count = 0;
  while (std::getline(lines, line))
    if (line.rfind(start, 0) == 0)
      ++count;
  return count;
}

static const std::array<stk::math::Vector3d,4> unitRightTet{{ {0.,0.,0.}, {1.,0.,0.}, {0.,1.,0.}, {0.,0.,1.} }};

static void expect_convex_section_ordered_around_normal(const std::vector<stk::math::Vector3d> & section, const stk::math::Plane3d & plane)
{
  for (size_t i=0; i<section.size(); ++i)
  {
    EXPECT_NEAR(0., plane.signed_distance(section[i]), 1.e-15);
    const stk::math::Vector3d & x0 = section[i];
    const stk::math::Vector3d & x1 = section[(i+1)%section.size()];
    const stk::math::Vector3d & x2 = section[(i+2)%section.size()];
    EXPECT_GT(Dot(Cross(x1-x0, x2-x1), plane.normal()), 0.);
  }
}

TEST(CubitJournal, twoTetsSharingAFace_fiveNodesAndTwoTetsInOneBlock)
{
  CubitJournal journal;
  journal.add_tet(unitRightTet);
  journal.add_tet({{ {1.,0.,0.}, {0.,1.,0.}, {0.,0.,1.}, {1.,1.,1.} }});
  const std::string str = journal.str();

  EXPECT_EQ(5, count_journal_lines_starting_with(str, "create node location "));
  EXPECT_EQ(2, count_journal_lines_starting_with(str, "create tet node "));
  EXPECT_EQ(1, count_journal_lines_starting_with(str, "block 1 tet {Id(\"tet\")-1} to {Id(\"tet\")}"));
  EXPECT_EQ(1, count_journal_lines_starting_with(str, "block 1 name \"default_tet\""));
  EXPECT_EQ(1, count_journal_lines_starting_with(str, "create tet node {krinoNode0+2} {krinoNode0+3} {krinoNode0+4} {krinoNode0+5}"));
}

TEST(CubitJournal, groups_oneBlockPerElementTypeWithNameAndColor)
{
  CubitJournal journal("groups", 7);
  journal.set_group("leaves", "blue");
  journal.add_tri({{ {0.,0.,0.}, {1.,0.,0.}, {0.,1.,0.} }});
  journal.add_tet(unitRightTet);
  journal.set_group("crossings", "red");
  journal.add_segment({0.,0.,0.}, {1.,1.,1.});
  journal.add_point({0.5,0.,0.}, "a crossing");
  journal.set_group("leaves");
  journal.add_comment("two\nlines");
  const std::string str = journal.str();

  EXPECT_EQ(1, count_journal_lines_starting_with(str, "block 7 tri {Id(\"tri\")}"));
  EXPECT_EQ(1, count_journal_lines_starting_with(str, "block 7 name \"leaves_tri\""));
  EXPECT_EQ(1, count_journal_lines_starting_with(str, "color block 7 blue"));
  EXPECT_EQ(1, count_journal_lines_starting_with(str, "block 8 name \"leaves_tet\""));
  EXPECT_EQ(1, count_journal_lines_starting_with(str, "block 9 name \"crossings_edge\""));
  EXPECT_EQ(1, count_journal_lines_starting_with(str, "color block 9 red"));
  EXPECT_EQ(1, count_journal_lines_starting_with(str, "create vertex 0.5 0 0 # a crossing"));
  EXPECT_EQ(1, count_journal_lines_starting_with(str, "group \"crossings_points\" add vertex {krinoVertex0+1}"));
  EXPECT_EQ(1, count_journal_lines_starting_with(str, "color vertex {krinoVertex0+1} red"));
  EXPECT_EQ(1, count_journal_lines_starting_with(str, "# two"));
  EXPECT_EQ(1, count_journal_lines_starting_with(str, "# lines"));
  EXPECT_EQ(0, count_journal_lines_starting_with(str, "# Group default"));
}

TEST(CubitJournal, coordinates_writtenExactly)
{
  const stk::math::Vector3d x(1./3., -2./7., 1.e-9/3.);
  CubitJournal journal;
  journal.add_point(x);
  const std::vector<stk::math::Vector3d> written = parse_cubit_journal(journal.str()).vertices;

  ASSERT_EQ(1u, written.size());
  for (int i=0; i<3; ++i)
  {
    EXPECT_EQ(x[i], written[0][i]);
  }
}

TEST(CubitJournal, planeNearTetCorner_triangleSection)
{
  const stk::math::Plane3d plane(stk::math::Vector3d(1.,0.,0.), stk::math::Vector3d(0.75,0.,0.));
  const std::vector<stk::math::Vector3d> section = section_of_tet_by_plane(plane, unitRightTet);

  EXPECT_EQ(3u, section.size());
  expect_convex_section_ordered_around_normal(section, plane);
}

TEST(CubitJournal, planeBetweenTwoPairsOfTetEdges_quadSection)
{
  const stk::math::Plane3d plane(stk::math::Vector3d(1.,1.,0.), stk::math::Vector3d(0.25,0.25,0.));
  const std::vector<stk::math::Vector3d> section = section_of_tet_by_plane(plane, unitRightTet);

  EXPECT_EQ(4u, section.size());
  expect_convex_section_ordered_around_normal(section, plane);

  CubitJournal journal;
  journal.add_plane_in_tet(plane, unitRightTet);
  EXPECT_EQ(2, count_journal_lines_starting_with(journal.str(), "create tri node "));
}

TEST(CubitJournal, planeMissingTet_nothingWritten)
{
  const stk::math::Plane3d plane(stk::math::Vector3d(1.,0.,0.), stk::math::Vector3d(2.,0.,0.));
  EXPECT_TRUE(section_of_tet_by_plane(plane, unitRightTet).empty());

  CubitJournal journal;
  journal.add_plane_in_tet(plane, unitRightTet);
  EXPECT_EQ(0, count_journal_lines_starting_with(journal.str(), "create "));
}

TEST(CubitJournal, planeThroughTetVertexAndOppositeEdge_triangleWithTheVertex)
{
  const stk::math::Plane3d plane(stk::math::Vector3d(1.,-1.,0.), stk::math::Vector3d(0.,0.,0.));
  const std::vector<stk::math::Vector3d> section = section_of_tet_by_plane(plane, unitRightTet);

  EXPECT_EQ(3u, section.size());
  expect_convex_section_ordered_around_normal(section, plane);
}

TEST(CubitJournal, planeInBox_squareSection)
{
  const stk::math::Plane3d plane(stk::math::Vector3d(0.,0.,1.), stk::math::Vector3d(0.,0.,0.5));
  const std::vector<stk::math::Vector3d> section = section_of_box_by_plane(plane, {0.,0.,0.}, {1.,2.,3.});

  EXPECT_EQ(4u, section.size());
  expect_convex_section_ordered_around_normal(section, plane);
}

TEST(CubitJournal, lineInTri_segmentSection)
{
  const stk::math::Plane3d line(stk::math::Vector3d(1.,0.,0.), stk::math::Vector3d(0.5,0.,0.));
  const std::vector<stk::math::Vector3d> section = section_of_tri_by_plane(line, {{ {0.,0.,0.}, {1.,0.,0.}, {0.,1.,0.} }});

  ASSERT_EQ(2u, section.size());
  for (auto && x : section)
  {
    EXPECT_DOUBLE_EQ(0.5, x[0]);
  }
}

TEST(CubitJournal, clipQuadByPlane_partOnPositiveSide)
{
  const std::vector<stk::math::Vector3d> square{{0.,0.,0.}, {1.,0.,0.}, {1.,1.,0.}, {0.,1.,0.}};
  const stk::math::Plane3d plane(stk::math::Vector3d(1.,0.,0.), stk::math::Vector3d(0.25,0.,0.));
  const std::vector<stk::math::Vector3d> clipped = clip_section_to_positive_side_of_plane(square, plane);

  ASSERT_EQ(4u, clipped.size());
  for (auto && x : clipped)
  {
    EXPECT_GE(x[0], 0.25);
  }
  EXPECT_TRUE(clip_section_to_positive_side_of_plane(square, stk::math::Plane3d(stk::math::Vector3d(1.,0.,0.), stk::math::Vector3d(2.,0.,0.))).empty());
}

TEST(CubitJournal, clipSectionInThePlaneWithRoundoff_keptWhole)
{
  const stk::math::Plane3d plane(stk::math::Vector3d(1.,1.,1.), stk::math::Vector3d(0.,0.,0.5));
  const std::vector<stk::math::Vector3d> section = section_of_tet_by_plane(plane, unitRightTet);
  const stk::math::Plane3d reversedPlane(-plane.normal(), stk::math::Vector3d(0.,0.,0.5));

  EXPECT_EQ(3u, clip_section_to_positive_side_of_plane(section, plane).size());
  EXPECT_EQ(3u, clip_section_to_positive_side_of_plane(section, reversedPlane).size());
}

TEST(CubitJournal, clipSegmentByPlane_endAtTheCrossing)
{
  const std::vector<stk::math::Vector3d> segment{{0.,0.,0.}, {1.,0.,0.}};
  const stk::math::Plane3d plane(stk::math::Vector3d(-1.,0.,0.), stk::math::Vector3d(0.25,0.,0.));
  const std::vector<stk::math::Vector3d> clipped = clip_section_to_positive_side_of_plane(segment, plane);

  ASSERT_EQ(2u, clipped.size());
  EXPECT_DOUBLE_EQ(0., clipped[0][0]);
  EXPECT_DOUBLE_EQ(0.25, clipped[1][0]);
}

TEST(CubitJournal, zeroContourOfLinearFunction_exactSectionOfTheBox)
{
  const stk::math::Plane3d plane(stk::math::Vector3d(1.,2.,3.), stk::math::Vector3d(0.,0.,0.3));
  CubitJournal journal;
  add_zero_contour_in_box(journal, [&](const stk::math::Vector3d & x) { return plane.signed_distance(x); }, {0.,0.,0.}, {1.,1.,1.}, 4);
  const ParsedCubitJournal parsed = parse_cubit_journal(journal.str());

  ASSERT_FALSE(parsed.tris.empty());
  for (auto && x : parsed.nodes)
  {
    EXPECT_NEAR(0., plane.signed_distance(x), 1.e-15);
  }
  const std::vector<stk::math::Vector3d> section = section_of_box_by_plane(plane, {0.,0.,0.}, {1.,1.,1.});
  double sectionArea = 0.;
  for (size_t i=1; i+1<section.size(); ++i)
    sectionArea += 0.5*Cross(section[i]-section[0], section[i+1]-section[0]).length();
  EXPECT_NEAR(sectionArea, total_area_of_tris(parsed), 1.e-14);
}

TEST(CubitJournal, zeroContourOfSphereDistance_closeToTheSphere)
{
  const double radius = 0.6;
  CubitJournal journal;
  add_zero_contour_in_box(journal, [&](const stk::math::Vector3d & x) { return x.length() - radius; }, {-1.,-1.,-1.}, {1.,1.,1.}, 16);
  const ParsedCubitJournal parsed = parse_cubit_journal(journal.str());

  // The crossings of the interpolant are within a small part of the grid spacing of the sphere.
  const double gridSpacing = 2./16.;
  for (auto && x : parsed.nodes)
  {
    EXPECT_NEAR(radius, x.length(), 0.1*gridSpacing);
  }
  const double sphereArea = 4.*std::acos(-1.)*radius*radius;
  EXPECT_NEAR(sphereArea, total_area_of_tris(parsed), 0.02*sphereArea);
}

TEST(CubitJournal, zeroContourIn2D_segmentsOnTheLine)
{
  CubitJournal journal;
  add_zero_contour_in_box(journal, [](const stk::math::Vector3d & x) { return x[0] - 0.3; }, {0.,0.,0.}, {1.,1.,0.}, 4);
  const ParsedCubitJournal parsed = parse_cubit_journal(journal.str());

  EXPECT_TRUE(parsed.tris.empty());
  for (auto && x : parsed.nodes)
  {
    EXPECT_DOUBLE_EQ(0.3, x[0]);
  }
  EXPECT_NEAR(1., total_length_of_edges(parsed), 1.e-15);
}

}
