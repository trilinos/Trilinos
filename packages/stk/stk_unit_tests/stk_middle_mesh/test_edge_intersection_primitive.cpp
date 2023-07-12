#include "gtest/gtest.h"

#include "stk_middle_mesh/predicates/edge_intersection_primitive.hpp"

namespace stk {
namespace middle_mesh {
namespace impl {

using namespace predicates::impl;

TEST(EdgeIntersectionPrimitive, LinearbInteriorgInteriorIntersectionComputation)
{
  utils::Point b0(0.5, -0.5, 0), b1(0.5, 0.5, 0);
  utils::Point g0(0, 0, 0), g1(1, 0, 0);
  utils::Point d0(0, 0, 1), d1(0, 0, 1);

  EdgeIntersectionData edgeIntersectionData(b0, b1, g0, g1, d0, d1);
  edgeIntersectionData.compute_beta();

  EXPECT_TRUE(edgeIntersectionData.is_beta_1_valid());
  EXPECT_DOUBLE_EQ(edgeIntersectionData.get_beta_1_value(), 0.5);
  EXPECT_FALSE(edgeIntersectionData.is_beta_2_valid());

  edgeIntersectionData.compute_alpha();

  EXPECT_TRUE(edgeIntersectionData.is_alpha_1_valid());
  EXPECT_DOUBLE_EQ(edgeIntersectionData.get_alpha_1_value(), 0.5);
  EXPECT_FALSE(edgeIntersectionData.is_alpha_2_valid());
}

TEST(EdgeIntersectionPrimitive, LinearbInteriorgInteriorIntersectionResult)
{
  utils::Point b0(0.5, -0.5, 0), b1(0.5, 0.5, 0);
  utils::Point g0(0, 0, 0), g1(1, 0, 0);
  utils::Point d0(0, 0, 1), d1(0, 0, 1);

  EdgeIntersectionResult result = compute_edge_intersection(b0, b1, g0, g1, d0, d1);
  EXPECT_TRUE(result.intersection1_found());
  EXPECT_DOUBLE_EQ(result.get_beta1(), 0.5);
  EXPECT_DOUBLE_EQ(result.get_alpha1(), 0.5);
}

TEST(EdgeIntersectionPrimitive, LinearbVertex0gVertex0IntersectionComputation)
{
  utils::Point b0(0, 0, 0), b1(0, 1, 0);
  utils::Point g0(0, 0, 0), g1(1, 0, 0);
  utils::Point d0(0, 0, 1), d1(0, 0, 1);

  EdgeIntersectionData edgeIntersectionData(b0, b1, g0, g1, d0, d1);
  edgeIntersectionData.compute_beta();

  EXPECT_TRUE(edgeIntersectionData.is_beta_1_valid());
  EXPECT_DOUBLE_EQ(edgeIntersectionData.get_beta_1_value(), 0);
  EXPECT_FALSE(edgeIntersectionData.is_beta_2_valid());

  edgeIntersectionData.compute_alpha();

  EXPECT_TRUE(edgeIntersectionData.is_alpha_1_valid());
  EXPECT_DOUBLE_EQ(edgeIntersectionData.get_alpha_1_value(), 0);
  EXPECT_FALSE(edgeIntersectionData.is_alpha_2_valid());
}

TEST(EdgeIntersectionPrimitive, LinearbVertex0gVertex0IntersectionResult)
{
  utils::Point b0(0, 0, 0), b1(0, 1, 0);
  utils::Point g0(0, 0, 0), g1(1, 0, 0);
  utils::Point d0(0, 0, 1), d1(0, 0, 1);

  EdgeIntersectionResult result = compute_edge_intersection(b0, b1, g0, g1, d0, d1);
  EXPECT_TRUE(result.intersection1_found());
  EXPECT_DOUBLE_EQ(result.get_beta1(), 0);
  EXPECT_DOUBLE_EQ(result.get_alpha1(), 0);
}

TEST(EdgeIntersectionPrimitive, LinearbVertex1gVertex0IntersectionComputation)
{
  utils::Point b0(1, 0, 0), b1(1, 1, 0);
  utils::Point g0(0, 0, 0), g1(1, 0, 0);
  utils::Point d0(0, 0, 1), d1(0, 0, 1);

  EdgeIntersectionData edgeIntersectionData(b0, b1, g0, g1, d0, d1);
  edgeIntersectionData.compute_beta();

  EXPECT_TRUE(edgeIntersectionData.is_beta_1_valid());
  EXPECT_DOUBLE_EQ(edgeIntersectionData.get_beta_1_value(), 1);
  EXPECT_FALSE(edgeIntersectionData.is_beta_2_valid());

  edgeIntersectionData.compute_alpha();

  EXPECT_TRUE(edgeIntersectionData.is_alpha_1_valid());
  EXPECT_DOUBLE_EQ(edgeIntersectionData.get_alpha_1_value(), 0);
  EXPECT_FALSE(edgeIntersectionData.is_alpha_2_valid());
}

TEST(EdgeIntersectionPrimitive, LinearbVertex1gVertex0IntersectionResult)
{
  utils::Point b0(1, 0, 0), b1(1, 1, 0);
  utils::Point g0(0, 0, 0), g1(1, 0, 0);
  utils::Point d0(0, 0, 1), d1(0, 0, 1);

  EdgeIntersectionResult result = compute_edge_intersection(b0, b1, g0, g1, d0, d1);
  EXPECT_TRUE(result.intersection1_found());
  EXPECT_DOUBLE_EQ(result.get_beta1(), 1);
  EXPECT_DOUBLE_EQ(result.get_alpha1(), 0);
}

TEST(EdgeIntersectionPrimitive, LinearbVertex0gVertex1IntersectionComputation)
{
  utils::Point b0(0, -1, 0), b1(0, 0, 0);
  utils::Point g0(0, 0, 0), g1(1, 0, 0);
  utils::Point d0(0, 0, 1), d1(0, 0, 1);

  EdgeIntersectionData edgeIntersectionData(b0, b1, g0, g1, d0, d1);
  edgeIntersectionData.compute_beta();

  EXPECT_TRUE(edgeIntersectionData.is_beta_1_valid());
  EXPECT_DOUBLE_EQ(edgeIntersectionData.get_beta_1_value(), 0);
  EXPECT_FALSE(edgeIntersectionData.is_beta_2_valid());

  edgeIntersectionData.compute_alpha();

  EXPECT_TRUE(edgeIntersectionData.is_alpha_1_valid());
  EXPECT_DOUBLE_EQ(edgeIntersectionData.get_alpha_1_value(), 1);
  EXPECT_FALSE(edgeIntersectionData.is_alpha_2_valid());
}

TEST(EdgeIntersectionPrimitive, LinearbVertex0gVertex1IntersectionResult)
{
  utils::Point b0(0, -1, 0), b1(0, 0, 0);
  utils::Point g0(0, 0, 0), g1(1, 0, 0);
  utils::Point d0(0, 0, 1), d1(0, 0, 1);

  EdgeIntersectionResult result = compute_edge_intersection(b0, b1, g0, g1, d0, d1);
  EXPECT_TRUE(result.intersection1_found());
  EXPECT_DOUBLE_EQ(result.get_beta1(), 0);
  EXPECT_DOUBLE_EQ(result.get_alpha1(), 1);
}

TEST(EdgeIntersectionPrimitive, LinearbVertex1gVertex1IntersectionComputation)
{
  utils::Point b0(1, -1, 0), b1(1, 0, 0);
  utils::Point g0(0, 0, 0), g1(1, 0, 0);
  utils::Point d0(0, 0, 1), d1(0, 0, 1);

  EdgeIntersectionData edgeIntersectionData(b0, b1, g0, g1, d0, d1);
  edgeIntersectionData.compute_beta();

  EXPECT_TRUE(edgeIntersectionData.is_beta_1_valid());
  EXPECT_DOUBLE_EQ(edgeIntersectionData.get_beta_1_value(), 1);
  EXPECT_FALSE(edgeIntersectionData.is_beta_2_valid());

  edgeIntersectionData.compute_alpha();

  EXPECT_TRUE(edgeIntersectionData.is_alpha_1_valid());
  EXPECT_DOUBLE_EQ(edgeIntersectionData.get_alpha_1_value(), 1);
  EXPECT_FALSE(edgeIntersectionData.is_alpha_2_valid());
}

TEST(EdgeIntersectionPrimitive, LinearbVertex1gVertex1IntersectionResult)
{
  utils::Point b0(1, -1, 0), b1(1, 0, 0);
  utils::Point g0(0, 0, 0), g1(1, 0, 0);
  utils::Point d0(0, 0, 1), d1(0, 0, 1);

  EdgeIntersectionResult result = compute_edge_intersection(b0, b1, g0, g1, d0, d1);
  EXPECT_TRUE(result.intersection1_found());
  EXPECT_DOUBLE_EQ(result.get_beta1(), 1);
  EXPECT_DOUBLE_EQ(result.get_alpha1(), 1);
}

TEST(EdgeIntersectionPrimitive, LinearbInteriorgVertexIntersectionComputation)
{
  utils::Point b0(0, -0.5, 0), b1(0, 0.5, 0);
  utils::Point g0(0, 0, 0), g1(1, 0, 0);
  utils::Point d0(0, 0, 1), d1(0, 0, 1);

  EdgeIntersectionData edgeIntersectionData(b0, b1, g0, g1, d0, d1);
  edgeIntersectionData.compute_beta();

  EXPECT_TRUE(edgeIntersectionData.is_beta_1_valid());
  EXPECT_DOUBLE_EQ(edgeIntersectionData.get_beta_1_value(), 0);
  EXPECT_FALSE(edgeIntersectionData.is_beta_2_valid());

  edgeIntersectionData.compute_alpha();

  EXPECT_TRUE(edgeIntersectionData.is_alpha_1_valid());
  EXPECT_DOUBLE_EQ(edgeIntersectionData.get_alpha_1_value(), 0.5);
  EXPECT_FALSE(edgeIntersectionData.is_alpha_2_valid());
}

TEST(EdgeIntersectionPrimitive, LinearbInteriorgVertexIntersectionResult)
{
  utils::Point b0(0, -0.5, 0), b1(0, 0.5, 0);
  utils::Point g0(0, 0, 0), g1(1, 0, 0);
  utils::Point d0(0, 0, 1), d1(0, 0, 1);

  EdgeIntersectionResult result = compute_edge_intersection(b0, b1, g0, g1, d0, d1);
  EXPECT_TRUE(result.intersection1_found());
  EXPECT_DOUBLE_EQ(result.get_beta1(), 0);
  EXPECT_DOUBLE_EQ(result.get_alpha1(), 0.5);
}

TEST(EdgeIntersectionPrimitive, LinearbVertexgInteriorIntersectionComputation)
{
  utils::Point b0(0.5, 0, 0), b1(0.5, 1, 0);
  utils::Point g0(0, 0, 0), g1(1, 0, 0);
  utils::Point d0(0, 0, 1), d1(0, 0, 1);

  EdgeIntersectionData edgeIntersectionData(b0, b1, g0, g1, d0, d1);
  edgeIntersectionData.compute_beta();

  EXPECT_TRUE(edgeIntersectionData.is_beta_1_valid());
  EXPECT_DOUBLE_EQ(edgeIntersectionData.get_beta_1_value(), 0.5);
  EXPECT_FALSE(edgeIntersectionData.is_beta_2_valid());

  edgeIntersectionData.compute_alpha();

  EXPECT_TRUE(edgeIntersectionData.is_alpha_1_valid());
  EXPECT_DOUBLE_EQ(edgeIntersectionData.get_alpha_1_value(), 0);
  EXPECT_FALSE(edgeIntersectionData.is_alpha_2_valid());
}

TEST(EdgeIntersectionPrimitive, LinearbVertexgInteriorIntersectionResult)
{
  utils::Point b0(0.5, 0, 0), b1(0.5, 1, 0);
  utils::Point g0(0, 0, 0), g1(1, 0, 0);
  utils::Point d0(0, 0, 1), d1(0, 0, 1);

  EdgeIntersectionResult result = compute_edge_intersection(b0, b1, g0, g1, d0, d1);
  EXPECT_TRUE(result.intersection1_found());
  EXPECT_DOUBLE_EQ(result.get_beta1(), 0.5);
  EXPECT_DOUBLE_EQ(result.get_alpha1(), 0);
}

TEST(EdgeIntersectionPrimitive, LinearbgParallelIntersectionComputation1)
{
  utils::Point b0(0, 0, 0), b1(1, 0, 0);
  utils::Point g0(0, 1, 0), g1(1, 1, 0);
  utils::Point d0(0, 0, 1), d1(0, 0, 1);

  EdgeIntersectionData edgeIntersectionData(b0, b1, g0, g1, d0, d1);
  edgeIntersectionData.compute_beta();

  EXPECT_FALSE(edgeIntersectionData.is_beta_1_valid());
  EXPECT_FALSE(edgeIntersectionData.is_beta_2_valid());

  edgeIntersectionData.compute_alpha();

  EXPECT_FALSE(edgeIntersectionData.is_alpha_1_valid());
  EXPECT_FALSE(edgeIntersectionData.is_alpha_2_valid());
}

TEST(EdgeIntersectionPrimitive, LinearbgParallelIntersectionResult1)
{
  utils::Point b0(0, 0, 0), b1(1, 0, 0);
  utils::Point g0(0, 1, 0), g1(1, 1, 0);
  utils::Point d0(0, 0, 1), d1(0, 0, 1);

  EdgeIntersectionResult result = compute_edge_intersection(b0, b1, g0, g1, d0, d1);
  EXPECT_FALSE(result.intersection1_found());
}

TEST(EdgeIntersectionPrimitive, LinearbgParallelIntersectionComputation2)
{
  utils::Point b0(0, 0, 0), b1(1, 0, 0);
  utils::Point g0(0, 1, 0), g1(1, 1, 0);
  utils::Point d0(0, 0, 1), d1(1 / std::sqrt(3), 1 / std::sqrt(3), -1 / std::sqrt(3));

  EdgeIntersectionData edgeIntersectionData(b0, b1, g0, g1, d0, d1);
  edgeIntersectionData.compute_beta();

  EXPECT_TRUE(edgeIntersectionData.is_beta_1_valid());
  EXPECT_NEAR(edgeIntersectionData.get_beta_1_value(), 0.634, 1e-4);
  EXPECT_FALSE(edgeIntersectionData.is_beta_2_valid());

  edgeIntersectionData.compute_alpha();

  EXPECT_TRUE(edgeIntersectionData.is_alpha_1_valid());
  EXPECT_FALSE(edgeIntersectionData.is_alpha_2_valid());
}

TEST(EdgeIntersectionPrimitive, LinearbgParallelIntersectionResult2)
{
  utils::Point b0(0, 0, 0), b1(1, 0, 0);
  utils::Point g0(0, 1, 0), g1(1, 1, 0);
  utils::Point d0(0, 0, 1), d1(1 / std::sqrt(3), 1 / std::sqrt(3), -1 / std::sqrt(3));

  EdgeIntersectionResult result = compute_edge_intersection(b0, b1, g0, g1, d0, d1);
  EXPECT_TRUE(result.intersection1_found());
}

TEST(EdgeIntersectionPrimitive, QuadraticbInteriorgInteriorOneSolutionIntersectionComputation)
{
  utils::Point b0(0.5, -0.5, 0), b1(0.5, 0.5, 0);
  utils::Point g0(0, 0, 0), g1(1, 0, 0);
  utils::Point d0(0, 0, 1), d1(1 / std::sqrt(3), 1 / std::sqrt(3), -1 / std::sqrt(3));

  EdgeIntersectionData edgeIntersectionData(b0, b1, g0, g1, d0, d1);
  edgeIntersectionData.compute_beta();

  EXPECT_TRUE(edgeIntersectionData.is_beta_1_valid());
  EXPECT_NEAR(edgeIntersectionData.get_beta_1_value(), 0.5, 1e-13);
  EXPECT_TRUE(edgeIntersectionData.is_beta_2_valid());
  EXPECT_NEAR(edgeIntersectionData.get_beta_2_value(), 0.634, 1e-4);

  edgeIntersectionData.compute_alpha();

  EXPECT_TRUE(edgeIntersectionData.is_alpha_1_valid());
  EXPECT_DOUBLE_EQ(edgeIntersectionData.get_alpha_1_value(), 0.5);
  EXPECT_TRUE(edgeIntersectionData.is_alpha_2_valid());
  EXPECT_NEAR(edgeIntersectionData.get_alpha_2_value(), 0.366025, 1e-4);
}

TEST(EdgeIntersectionPrimitive, QuadraticbInteriorgInteriorOneSolutionIntersectionResult)
{
  utils::Point b0(0.5, -0.5, 0), b1(0.5, 0.5, 0);
  utils::Point g0(0, 0, 0), g1(1, 0, 0);
  utils::Point d0(0, 0, 1), d1(1 / std::sqrt(3), 1 / std::sqrt(3), -1 / std::sqrt(3));

  EdgeIntersectionResult result = compute_edge_intersection(b0, b1, g0, g1, d0, d1);
  EXPECT_TRUE(result.intersection1_found());
  EXPECT_NEAR(result.get_beta1(), 0.5, 1e-13);
  EXPECT_NEAR(result.get_alpha1(), 0.5, 1e-13);
}

TEST(EdgeIntersectionPrimitive, QuadraticbVertex0gVertex0OneSolutionIntersectionComputation)
{
  utils::Point b0(0, 0, 0), b1(0, 1, 0);
  utils::Point g0(0, 0, 0), g1(1, 0, 0);
  utils::Point d0(0, 0, 1), d1(1 / std::sqrt(3), 1 / std::sqrt(3), -1 / std::sqrt(3));

  EdgeIntersectionData edgeIntersectionData(b0, b1, g0, g1, d0, d1);
  edgeIntersectionData.compute_beta();

  EXPECT_TRUE(edgeIntersectionData.is_beta_1_valid());
  EXPECT_DOUBLE_EQ(edgeIntersectionData.get_beta_1_value(), 0);
  EXPECT_TRUE(edgeIntersectionData.is_beta_2_valid());
  EXPECT_NEAR(edgeIntersectionData.get_beta_2_value(), 0.634, 1e-4);

  edgeIntersectionData.compute_alpha();

  EXPECT_TRUE(edgeIntersectionData.is_alpha_1_valid());
  EXPECT_DOUBLE_EQ(edgeIntersectionData.get_alpha_1_value(), 0);
  EXPECT_TRUE(edgeIntersectionData.is_alpha_2_valid());
  EXPECT_NEAR(edgeIntersectionData.get_alpha_2_value(), -0.633975, 1e-4);
}

TEST(EdgeIntersectionPrimitive, QuadraticbVertex0gVertex0OneSolutionIntersectionResult)
{
  utils::Point b0(0, 0, 0), b1(0, 1, 0);
  utils::Point g0(0, 0, 0), g1(1, 0, 0);
  utils::Point d0(0, 0, 1), d1(1 / std::sqrt(3), 1 / std::sqrt(3), -1 / std::sqrt(3));

  EdgeIntersectionResult result = compute_edge_intersection(b0, b1, g0, g1, d0, d1);
  EXPECT_TRUE(result.intersection1_found());
  EXPECT_DOUBLE_EQ(result.get_beta1(), 0);
  EXPECT_DOUBLE_EQ(result.get_alpha1(), 0);
}

TEST(EdgeIntersectionPrimitive, QuadraticbVertex1gVertex0IntersectionOneSolutionComputation)
{
  utils::Point b0(1, 0, 0), b1(1, 1, 0);
  utils::Point g0(0, 0, 0), g1(1, 0, 0);
  utils::Point d0(0, 0, 1), d1(1 / std::sqrt(3), 1 / std::sqrt(3), -1 / std::sqrt(3));

  EdgeIntersectionData edgeIntersectionData(b0, b1, g0, g1, d0, d1);
  edgeIntersectionData.compute_beta();

  EXPECT_TRUE(edgeIntersectionData.is_beta_1_valid());
  EXPECT_NEAR(edgeIntersectionData.get_beta_1_value(), 0.634, 1e-4);
  EXPECT_TRUE(edgeIntersectionData.is_beta_2_valid());
  EXPECT_DOUBLE_EQ(edgeIntersectionData.get_beta_2_value(), 1);

  edgeIntersectionData.compute_alpha();

  EXPECT_TRUE(edgeIntersectionData.is_alpha_1_valid());
  EXPECT_NEAR(edgeIntersectionData.get_alpha_1_value(), 0.366025, 1e-4);
  EXPECT_TRUE(edgeIntersectionData.is_alpha_2_valid());
  EXPECT_NEAR(edgeIntersectionData.get_alpha_2_value(), 0, 1e-13);
}

TEST(EdgeIntersectionPrimitive, QuadraticbVertex1gVertex0OneSolutionIntersectionResult)
{
  utils::Point b0(1, 0, 0), b1(1, 1, 0);
  utils::Point g0(0, 0, 0), g1(1, 0, 0);
  utils::Point d0(0, 0, 1), d1(1 / std::sqrt(3), 1 / std::sqrt(3), -1 / std::sqrt(3));

  EdgeIntersectionResult result = compute_edge_intersection(b0, b1, g0, g1, d0, d1);
  EXPECT_TRUE(result.intersection1_found());
  EXPECT_TRUE(result.intersection2_found());

  EXPECT_NEAR(result.get_beta1(), 0.633975, 1e-4);
  EXPECT_NEAR(result.get_alpha1(), 0.366025, 1e-4);

  EXPECT_NEAR(result.get_beta2(), 1, 1e-13);
  EXPECT_NEAR(result.get_alpha2(), 0, 1e-13);
}

TEST(EdgeIntersectionPrimitive, QuadraticbVertex0gVertex1OneSolutionIntersectionComputation)
{
  utils::Point b0(0, -1, 0), b1(0, 0, 0);
  utils::Point g0(0, 0, 0), g1(1, 0, 0);
  utils::Point d0(0, 0, 1), d1(1 / std::sqrt(3), 1 / std::sqrt(3), -1 / std::sqrt(3));

  EdgeIntersectionData edgeIntersectionData(b0, b1, g0, g1, d0, d1);
  edgeIntersectionData.compute_beta();

  EXPECT_TRUE(edgeIntersectionData.is_beta_1_valid());
  EXPECT_DOUBLE_EQ(edgeIntersectionData.get_beta_1_value(), 0);
  EXPECT_TRUE(edgeIntersectionData.is_beta_2_valid());
  EXPECT_NEAR(edgeIntersectionData.get_beta_2_value(), 0.634, 1e-4);

  edgeIntersectionData.compute_alpha();

  EXPECT_TRUE(edgeIntersectionData.is_alpha_1_valid());
  EXPECT_DOUBLE_EQ(edgeIntersectionData.get_alpha_1_value(), 1);
  EXPECT_TRUE(edgeIntersectionData.is_alpha_2_valid());
  EXPECT_NEAR(edgeIntersectionData.get_alpha_2_value(), 0.366025, 1e-4);
}

TEST(EdgeIntersectionPrimitive, QuadraticbVertex0gVertex1OneSolutionIntersectionResult)
{
  utils::Point b0(0, -1, 0), b1(0, 0, 0);
  utils::Point g0(0, 0, 0), g1(1, 0, 0);
  utils::Point d0(0, 0, 1), d1(1 / std::sqrt(3), 1 / std::sqrt(3), -1 / std::sqrt(3));

  EdgeIntersectionResult result = compute_edge_intersection(b0, b1, g0, g1, d0, d1);
  EXPECT_TRUE(result.intersection1_found());
  EXPECT_DOUBLE_EQ(result.get_beta1(), 0);
  EXPECT_DOUBLE_EQ(result.get_alpha1(), 1);
}

TEST(EdgeIntersectionPrimitive, QuadraticbVertex1gVertex1OneSolutionIntersectionComputation)
{
  utils::Point b0(1, -1, 0), b1(1, 0, 0);
  utils::Point g0(0, 0, 0), g1(1, 0, 0);
  utils::Point d0(0, 0, 1), d1(1 / std::sqrt(3), 1 / std::sqrt(3), -1 / std::sqrt(3));

  EdgeIntersectionData edgeIntersectionData(b0, b1, g0, g1, d0, d1);
  edgeIntersectionData.compute_beta();

  EXPECT_TRUE(edgeIntersectionData.is_beta_1_valid());
  EXPECT_NEAR(edgeIntersectionData.get_beta_1_value(), 0.634, 1e-4);
  EXPECT_TRUE(edgeIntersectionData.is_beta_2_valid());
  EXPECT_DOUBLE_EQ(edgeIntersectionData.get_beta_2_value(), 1);

  edgeIntersectionData.compute_alpha();

  EXPECT_TRUE(edgeIntersectionData.is_alpha_1_valid());
  EXPECT_NEAR(edgeIntersectionData.get_alpha_1_value(), 1.36603, 1e-4);
  EXPECT_TRUE(edgeIntersectionData.is_alpha_2_valid());
  EXPECT_DOUBLE_EQ(edgeIntersectionData.get_alpha_2_value(), 1);
}

TEST(EdgeIntersectionPrimitive, QuadraticbVertex1gVertex1OneSolutionIntersectionResult)
{
  utils::Point b0(1, -1, 0), b1(1, 0, 0);
  utils::Point g0(0, 0, 0), g1(1, 0, 0);
  utils::Point d0(0, 0, 1), d1(1 / std::sqrt(3), 1 / std::sqrt(3), -1 / std::sqrt(3));

  EdgeIntersectionResult result = compute_edge_intersection(b0, b1, g0, g1, d0, d1);
  EXPECT_TRUE(result.intersection1_found());
  EXPECT_TRUE(result.intersection2_found());

  EXPECT_NEAR(result.get_beta1(), 0.633975, 1e-4);
  EXPECT_NEAR(result.get_alpha1(), 1.36603, 1e-4);

  EXPECT_DOUBLE_EQ(result.get_beta2(), 1);
  EXPECT_DOUBLE_EQ(result.get_alpha2(), 1);
}

TEST(EdgeIntersectionPrimitive, QuadraticbInteriorgVertexOneSolutionIntersectionComputation)
{
  utils::Point b0(0, -0.5, 0), b1(0, 0.5, 0);
  utils::Point g0(0, 0, 0), g1(1, 0, 0);
  utils::Point d0(0, 0, 1), d1(1 / std::sqrt(3), 1 / std::sqrt(3), -1 / std::sqrt(3));

  EdgeIntersectionData edgeIntersectionData(b0, b1, g0, g1, d0, d1);
  edgeIntersectionData.compute_beta();

  EXPECT_TRUE(edgeIntersectionData.is_beta_1_valid());
  EXPECT_DOUBLE_EQ(edgeIntersectionData.get_beta_1_value(), 0);
  EXPECT_TRUE(edgeIntersectionData.is_beta_2_valid());
  EXPECT_NEAR(edgeIntersectionData.get_beta_2_value(), 0.634, 1e-4);

  edgeIntersectionData.compute_alpha();

  EXPECT_TRUE(edgeIntersectionData.is_alpha_1_valid());
  EXPECT_DOUBLE_EQ(edgeIntersectionData.get_alpha_1_value(), 0.5);
  EXPECT_TRUE(edgeIntersectionData.is_alpha_2_valid());
  EXPECT_NEAR(edgeIntersectionData.get_alpha_2_value(), -0.133975, 1e-4);
}

TEST(EdgeIntersectionPrimitive, QuadraticbInteriorgVertexOneSolutionIntersectionResult)
{
  utils::Point b0(0, -0.5, 0), b1(0, 0.5, 0);
  utils::Point g0(0, 0, 0), g1(1, 0, 0);
  utils::Point d0(0, 0, 1), d1(1 / std::sqrt(3), 1 / std::sqrt(3), -1 / std::sqrt(3));

  EdgeIntersectionResult result = compute_edge_intersection(b0, b1, g0, g1, d0, d1);
  EXPECT_TRUE(result.intersection1_found());
  EXPECT_DOUBLE_EQ(result.get_beta1(), 0);
  EXPECT_DOUBLE_EQ(result.get_alpha1(), 0.5);
}

TEST(EdgeIntersectionPrimitive, QuadraticbVertexgInteriorOneSolutionIntersectionComputation)
{
  utils::Point b0(0.5, 0, 0), b1(0.5, 1, 0);
  utils::Point g0(0, 0, 0), g1(1, 0, 0);
  utils::Point d0(0, 0, 1), d1(1 / std::sqrt(3), 1 / std::sqrt(3), -1 / std::sqrt(3));

  EdgeIntersectionData edgeIntersectionData(b0, b1, g0, g1, d0, d1);
  edgeIntersectionData.compute_beta();

  EXPECT_TRUE(edgeIntersectionData.is_beta_1_valid());
  EXPECT_NEAR(edgeIntersectionData.get_beta_1_value(), 0.5, 1e-13);
  EXPECT_TRUE(edgeIntersectionData.is_beta_2_valid());
  EXPECT_NEAR(edgeIntersectionData.get_beta_2_value(), 0.634, 1e-4);

  edgeIntersectionData.compute_alpha();

  EXPECT_TRUE(edgeIntersectionData.is_alpha_1_valid());
  EXPECT_NEAR(edgeIntersectionData.get_alpha_1_value(), 0, 1e-13);
  EXPECT_TRUE(edgeIntersectionData.is_alpha_2_valid());
  EXPECT_NEAR(edgeIntersectionData.get_alpha_2_value(), -0.133975, 1e-4);
}

TEST(EdgeIntersectionPrimitive, QuadraticbVertexgInteriorOneSolutionIntersectionResult)
{
  utils::Point b0(0.5, 0, 0), b1(0.5, 1, 0);
  utils::Point g0(0, 0, 0), g1(1, 0, 0);
  utils::Point d0(0, 0, 1), d1(1 / std::sqrt(3), 1 / std::sqrt(3), -1 / std::sqrt(3));

  EdgeIntersectionResult result = compute_edge_intersection(b0, b1, g0, g1, d0, d1);
  EXPECT_TRUE(result.intersection1_found());
  EXPECT_NEAR(result.get_beta1(), 0.5, 1e-13);
  EXPECT_NEAR(result.get_alpha1(), 0, 1e-13);
}

TEST(EdgeIntersectionPrimitive, Regression)
{
  // from James Overfelt, 5/25/2023
  const stk::middle_mesh::utils::Point edge1Pt1  (1.0000000000000002, 0, 0.39999999999999997);
  const stk::middle_mesh::utils::Point edge1Pt2  (1, 0.20000000000000107, 0.40000000000000019);
  const stk::middle_mesh::utils::Point edge2Pt1  (1, -3.388131789017199e-18, 0.44999999999999996);
  const stk::middle_mesh::utils::Point edge2Pt2  (1.0000000000000002, -1.3552527156068807e-18, 0.41999999999999998);
  const stk::middle_mesh::utils::Point pt1Normal (0.031079502397333692, -2.091728074970464e-16, 1.1501726385137537e-16);
  const stk::middle_mesh::utils::Point pt2Normal (0.02493346738921879, -0.016448179238063677, 1.0961226633804594e-16);
  EdgeIntersectionResult result = stk::middle_mesh::predicates::impl::compute_edge_intersection(edge1Pt1, edge1Pt2, edge2Pt1, edge2Pt2,
                                                                                                pt1Normal, pt2Normal, 1e-12, 0.25);

  EXPECT_TRUE(result.intersection1_found());
  EXPECT_FALSE(result.intersection2_found());
}

} // namespace impl
} // namespace middle_mesh
} // namespace stk
