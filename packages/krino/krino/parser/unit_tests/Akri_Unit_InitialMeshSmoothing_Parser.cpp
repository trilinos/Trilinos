// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <optional>
#include <string>

#include <gtest/gtest.h>

#include <Akri_InitialMeshSmoothing.hpp>
#include <Akri_Parser.hpp>
#include <Akri_Region_Parser.hpp>
#include <Akri_YAML.hpp>
#include <stk_util/environment/RuntimeDoomed.hpp>

namespace krino {
namespace {

#ifdef KRINO_HAVE_YAML

// Region_Parser::parse_initial_mesh_smoothing reads the region's smooth_initial_mesh.  RuntimeDoomedAdHoc
// increments a process-global counter rather than throwing: error cases check that it grew, the others
// that it did not.
InitialMeshSmoothing parse_smoothing_of_region(const std::string & regionYaml)
{
  return Region_Parser::parse_initial_mesh_smoothing(Parser::Node(YAML::Load(regionYaml)));
}

void expect_smoothing(const InitialMeshSmoothing & smoothing, const bool smooth, const bool boundaryNodesSlide,
    const std::optional<double> sharpFeatureAngleInDegrees)
{
  EXPECT_EQ(smooth, smoothing.smooth);
  EXPECT_EQ(boundaryNodesSlide, smoothing.boundaryNodesSlide);
  EXPECT_EQ(sharpFeatureAngleInDegrees, smoothing.sharpFeatureAngleInDegrees);
}

TEST(InitialMeshSmoothingParser, absent_noSmoothing)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  expect_smoothing(parse_smoothing_of_region("name: region1\n"), false, true, std::nullopt);
  EXPECT_EQ(oldDoomedCount, stk::get_doomed_count());
}

TEST(InitialMeshSmoothingParser, trueOrFalse_smoothingWithSlidingBoundaryNodesOrNone)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  expect_smoothing(parse_smoothing_of_region("smooth_initial_mesh: true\n"), true, true, std::nullopt);
  expect_smoothing(parse_smoothing_of_region("smooth_initial_mesh: false\n"), false, true, std::nullopt);
  EXPECT_EQ(oldDoomedCount, stk::get_doomed_count());
}

TEST(InitialMeshSmoothingParser, emptyBlock_smoothingWithSlidingBoundaryNodes)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  expect_smoothing(parse_smoothing_of_region("name: region1\nsmooth_initial_mesh:\n"), true, true, std::nullopt);
  EXPECT_EQ(oldDoomedCount, stk::get_doomed_count());
}

TEST(InitialMeshSmoothingParser, blockWithSharpFeatureAngle_slidingBoundaryNodesWithTheAngle)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  expect_smoothing(parse_smoothing_of_region(R"yaml(
smooth_initial_mesh:
  sharp_feature_angle: 135
)yaml"), true, true, 135.);
  expect_smoothing(parse_smoothing_of_region(R"yaml(
smooth_initial_mesh:
  boundary_nodes: slide
  sharp_feature_angle: 120
)yaml"), true, true, 120.);
  EXPECT_EQ(oldDoomedCount, stk::get_doomed_count());
}

TEST(InitialMeshSmoothingParser, blockWithFixedBoundaryNodes_smoothingWithFixedBoundaryNodes)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  expect_smoothing(parse_smoothing_of_region(R"yaml(
smooth_initial_mesh:
  boundary_nodes: fixed
)yaml"), true, false, std::nullopt);
  EXPECT_EQ(oldDoomedCount, stk::get_doomed_count());
}

TEST(InitialMeshSmoothingParser, unknownBoundaryNodes_error)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  parse_smoothing_of_region(R"yaml(
smooth_initial_mesh:
  boundary_nodes: float
)yaml");
  EXPECT_LT(oldDoomedCount, stk::get_doomed_count());
}

TEST(InitialMeshSmoothingParser, sharpFeatureAngleOutOfRange_error)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  parse_smoothing_of_region(R"yaml(
smooth_initial_mesh:
  sharp_feature_angle: 200
)yaml");
  EXPECT_LT(oldDoomedCount, stk::get_doomed_count());
}

TEST(InitialMeshSmoothingParser, sharpFeatureAngleWithFixedBoundaryNodes_error)
{
  const unsigned oldDoomedCount = stk::get_doomed_count();
  parse_smoothing_of_region(R"yaml(
smooth_initial_mesh:
  boundary_nodes: fixed
  sharp_feature_angle: 135
)yaml");
  EXPECT_LT(oldDoomedCount, stk::get_doomed_count());
}

#endif

} // namespace
} // namespace krino
