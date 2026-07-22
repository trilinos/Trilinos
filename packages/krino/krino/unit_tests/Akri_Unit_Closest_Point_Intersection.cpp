// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "gtest/gtest.h"
#include <algorithm>
#include <random>
#include "Akri_Intersection_Points.hpp"
#include "Akri_Unit_Single_Element_Fixtures.hpp"
#include "Akri_ElementCutterUtils.hpp"
#include "krino/math_utils/Akri_MathUtil.hpp"

namespace 
{
  void compare_intersections(const std::vector<krino::IntersectionPoint> & intersectionPoints, unsigned startIdx, 
    const std::vector<stk::mesh::Entity> & nodeOrder, const std::vector<std::vector<double>> & expWeights)
  {
    constexpr double tol = 1e-8;
    ASSERT_EQ(expWeights.size(), intersectionPoints.size() - startIdx);
    for(unsigned w=0; w<expWeights.size(); w++)
    {
      ASSERT_EQ(expWeights[w].size(), nodeOrder.size());
    }
    for(unsigned p=startIdx; p<intersectionPoints.size(); p++)
    {
      ASSERT_EQ(intersectionPoints[p].get_nodes().size(), nodeOrder.size());
      std::vector<int> indexVec(nodeOrder.size());
      for(unsigned n=0; n<nodeOrder.size(); n++)
      {
        EXPECT_EQ(nodeOrder[n], intersectionPoints[p].get_nodes()[n]);
        indexVec[n] = n;
      }

      bool foundMatch = false;
      const auto & actWeightVec = intersectionPoints[p].get_weights();
      for(unsigned w=0; w<expWeights.size(); w++)
      {
        const auto & expWeightVec = expWeights[w];
        if(std::all_of(indexVec.begin(), indexVec.end(), [&expWeightVec, &actWeightVec] (int idx) 
          {return std::fabs(expWeightVec[idx] - actWeightVec[idx]) < tol;})) 
        {
          foundMatch = true;
          break;
        }
      }
      EXPECT_TRUE(foundMatch);
    }
  }
}

namespace krino
{

std::mt19937 mt(std::mt19937::default_seed);
using ClosestPointFixture = SingleElementFixture;

TEST(ClosestPointFixture, TriIntersection)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) GTEST_SKIP();
  stk::topology tet4 = stk::topology::TETRAHEDRON_4;
  SingleElementFixture test_fixture(tet4);

  test_fixture.generate_mesh();
  auto & bulk = test_fixture.stk_fixture.bulk_data();

  auto elem = test_fixture.my_elem;
  std::vector<stk::mesh::Entity> myNodes;
  for(auto node = bulk.begin_nodes(elem); node < bulk.end_nodes(elem); node++)
  {
    myNodes.push_back(*node);
  }

  std::vector<IntersectionPoint> intersectionPoints;

  for(unsigned i=0; i<3; i++)
  {
    IntersectionPoint intPt(true, {myNodes[0], myNodes[i+1]}, {0.5, 0.5}, {});
    intersectionPoints.push_back(intPt);
  }

  auto endInputPts = intersectionPoints.size();

  append_closest_point_intersections(intersectionPoints,
    bulk,
    test_fixture.stk_fixture.meta_data().universal_part(),
    test_fixture.coord_field,
    keep_all_intersection_points_filter());

  compare_intersections(intersectionPoints, endInputPts, myNodes, {{0.5, 0.5/3., 0.5/3., 0.5/3.}});
}

TEST(ClosestPointFixture, QuadIntersection)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) GTEST_SKIP();
  stk::topology tet4 = stk::topology::TETRAHEDRON_4;
  SingleElementFixture test_fixture(tet4);

  test_fixture.generate_mesh();
  auto & bulk = test_fixture.stk_fixture.bulk_data();

  auto elem = test_fixture.my_elem;
  std::vector<stk::mesh::Entity> myNodes;
  for(auto node = bulk.begin_nodes(elem); node < bulk.end_nodes(elem); node++)
  {
    myNodes.push_back(*node);
  }

  std::vector<IntersectionPoint> intersectionPoints;

  for(unsigned i=2; i<4; i++)
  {
    {
      IntersectionPoint intPt(true, {myNodes[0], myNodes[i]}, {0.5, 0.5}, {});
      intersectionPoints.push_back(intPt);
    }
    {
      IntersectionPoint intPt(true, {myNodes[1], myNodes[i]}, {0.5, 0.5}, {});
      intersectionPoints.push_back(intPt);
    }
  }

  auto endInputPts = intersectionPoints.size();

  append_closest_point_intersections(intersectionPoints,
    bulk,
    test_fixture.stk_fixture.meta_data().universal_part(),
    test_fixture.coord_field,
    keep_all_intersection_points_filter());

  compare_intersections(intersectionPoints, endInputPts, myNodes, {{1./3., 1./6., 1./3., 1./6.},
    {0., 1./2., 1./4., 1./4.}, {1./2., 0., 1./4., 1./4.}});
}

TEST(ClosestPointFixture, OnlyTwoIntersections)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) GTEST_SKIP();
  stk::topology tet4 = stk::topology::TETRAHEDRON_4;
  SingleElementFixture test_fixture(tet4);

  test_fixture.generate_mesh();
  auto & bulk = test_fixture.stk_fixture.bulk_data();

  auto elem = test_fixture.my_elem;
  std::vector<stk::mesh::Entity> myNodes;
  for(auto node = bulk.begin_nodes(elem); node < bulk.end_nodes(elem); node++)
  {
    myNodes.push_back(*node);
  }

  std::vector<IntersectionPoint> intersectionPoints;

  for(unsigned i=2; i<4; i++)
  {
    IntersectionPoint intPt(true, {myNodes[1], myNodes[i]}, {0.5, 0.5}, {});
    intersectionPoints.push_back(intPt);
  }

  auto endInputPts = intersectionPoints.size();

  append_closest_point_intersections(intersectionPoints,
    bulk,
    test_fixture.stk_fixture.meta_data().universal_part(),
    test_fixture.coord_field,
    keep_all_intersection_points_filter());

  EXPECT_EQ(intersectionPoints.size(), endInputPts);
}

TEST(ClosestPointFixture, TriFuzzTest)
{
  std::uniform_real_distribution<double> realDist(0.01, 0.99);
  std::uniform_int_distribution<> intDist(0, 3);
  constexpr int dim = 3;
  constexpr double tol = 1e-12;

  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) GTEST_SKIP();
  stk::topology tet4 = stk::topology::TETRAHEDRON_4;
  SingleElementFixture test_fixture(tet4);

  test_fixture.generate_mesh();
  auto & bulk = test_fixture.stk_fixture.bulk_data();

  auto elem = test_fixture.my_elem;
  std::vector<stk::mesh::Entity> myNodes;
  std::vector<stk::math::Vector3d> myVerts;
  for(auto node = bulk.begin_nodes(elem); node < bulk.end_nodes(elem); node++)
  {
    myNodes.push_back(*node);
    myVerts.push_back(stk::math::Vector3d(field_data<double>(test_fixture.coord_field, *node), dim));
  }

  for(unsigned t=0; t<1000; t++)
  {
    std::vector<IntersectionPoint> intersectionPoints;
    int vert = intDist(mt);

    std::vector<stk::math::Vector3d> edgeInters;

    for(int i=0; i<4; i++)
    {
      if(i==vert) continue;
      double wt = realDist(mt);
      IntersectionPoint intPt(true, {myNodes[vert], myNodes[i]}, {wt, 1.-wt}, {});
      intersectionPoints.push_back(intPt);

      edgeInters.push_back(wt * myVerts[vert] + (1.-wt) * myVerts[i]);
    }

    auto endInputPts = intersectionPoints.size();

    append_closest_point_intersections(intersectionPoints,
      bulk,
      test_fixture.stk_fixture.meta_data().universal_part(),
      test_fixture.coord_field,
      keep_all_intersection_points_filter());

    for(unsigned i=endInputPts; i<intersectionPoints.size(); i++)
    {
      //Expectation is that each new intersection point lies in the triangle formed by the
      //intersected edges and is closer to a tet node than any of the intersected edges
      const auto & intPt = intersectionPoints[i];
      auto intLoc = stk::math::Vector3d::ZERO;
      for(unsigned w=0; w<intPt.get_weights().size(); w++) intLoc += myVerts[w] * intPt.get_weights()[w];
      auto parCoords = get_parametric_coordinates_of_point(edgeInters, intLoc);
      EXPECT_DOUBLE_EQ(parCoords[2], 0.);
      EXPECT_GT(parCoords[0], -tol);
      EXPECT_GT(parCoords[1], -tol);
      EXPECT_GT(1.-parCoords[0]-parCoords[1], -tol);
      
      bool isValid = false;
      for(unsigned n=0; n<myNodes.size(); n++)
      {
        double minLen = std::numeric_limits<double>::max();
        for(auto && edgeInter : edgeInters) minLen = std::fmin(minLen, (myVerts[n] - edgeInter).length());
        if((intLoc-myVerts[n]).length() < minLen)
        {
          isValid = true;
          break;
        }
      }
      EXPECT_TRUE(isValid);
    }
  }
}

TEST(ClosestPointFixture, QuadFuzzTest)
{
  std::uniform_real_distribution<double> realDist(0.01, 0.99);
  std::uniform_int_distribution<> intDist(0, 3);
  constexpr int dim = 3;
  constexpr double tol = 1e-12;

  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) GTEST_SKIP();
  stk::topology tet4 = stk::topology::TETRAHEDRON_4;
  SingleElementFixture test_fixture(tet4);

  test_fixture.generate_mesh();
  auto & bulk = test_fixture.stk_fixture.bulk_data();

  auto elem = test_fixture.my_elem;
  std::vector<stk::mesh::Entity> myNodes;
  std::vector<stk::math::Vector3d> myVerts;
  for(auto node = bulk.begin_nodes(elem); node < bulk.end_nodes(elem); node++)
  {
    myNodes.push_back(*node);
    myVerts.push_back(stk::math::Vector3d(field_data<double>(test_fixture.coord_field, *node), dim));
  }

  for(unsigned t=0; t<1000; t++)
  {
    std::vector<IntersectionPoint> intersectionPoints;
    int vert0 = intDist(mt);
    int vert1 = intDist(mt);
    while(vert1 == vert0) vert1 = intDist(mt);

    std::vector<stk::math::Vector3d> edgeInters;

    for(int i=0; i<4; i++)
    {
      if(i==vert0 || i==vert1) continue;
      {
        double wt = realDist(mt);
        IntersectionPoint intPt(true, {myNodes[vert0], myNodes[i]}, {wt, 1.-wt}, {});
        intersectionPoints.push_back(intPt);
        edgeInters.push_back(wt * myVerts[vert0] + (1.-wt) * myVerts[i]);
      }
      {
        double wt = realDist(mt);
        IntersectionPoint intPt(true, {myNodes[vert1], myNodes[i]}, {wt, 1.-wt}, {});
        intersectionPoints.push_back(intPt);
        edgeInters.push_back(wt * myVerts[vert1] + (1.-wt) * myVerts[i]);
      }
    }

    const std::vector<std::array<stk::math::Vector3d, 3>> tris = make_tris_from_intersections(edgeInters);
    ASSERT_EQ(tris.size(), 2u);

    auto endInputPts = intersectionPoints.size();

    append_closest_point_intersections(intersectionPoints,
      bulk,
      test_fixture.stk_fixture.meta_data().universal_part(),
      test_fixture.coord_field,
      keep_all_intersection_points_filter());

    for(unsigned i=endInputPts; i<intersectionPoints.size(); i++)
    {
      //Expectation is that each new intersection point lies in one of the triangles formed by the
      //intersected edges and is closer to a tet node than any of the intersected edges in that tri
      const auto & intPt = intersectionPoints[i];
      auto intLoc = stk::math::Vector3d::ZERO;
      for(unsigned w=0; w<intPt.get_weights().size(); w++) intLoc += myVerts[w] * intPt.get_weights()[w];
      bool liesInTri = false;
      for(auto && tri : tris)
      {
        auto parCoords = get_parametric_coordinates_of_point(tri, intLoc);
        EXPECT_DOUBLE_EQ(parCoords[2], 0.);
        if(parCoords[0] > -tol && parCoords[1] > -tol && 1.-parCoords[0]-parCoords[1] > -tol)
        {
          liesInTri = true;
          break;
        }
      }
      EXPECT_TRUE(liesInTri);
      
      bool isValid = false;
      for(auto && tri : tris)
      {
        for(unsigned n=0; n<myNodes.size(); n++)
        {
          double minLen = std::numeric_limits<double>::max();
          for(auto && edgeInter : tri) minLen = std::fmin(minLen, (myVerts[n] - edgeInter).length());
          if((intLoc-myVerts[n]).length() < minLen)
          {
            isValid = true;
            break;
          }
        }
      }
      EXPECT_TRUE(isValid);
    }
  }
}

TEST(ClosestPointFixture, NearZeroArea)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) GTEST_SKIP();
  stk::topology tet4 = stk::topology::TETRAHEDRON_4;
  SingleElementFixture test_fixture(tet4);

  test_fixture.generate_mesh();
  auto & bulk = test_fixture.stk_fixture.bulk_data();

  auto elem = test_fixture.my_elem;
  std::vector<stk::mesh::Entity> myNodes;
  for(auto node = bulk.begin_nodes(elem); node < bulk.end_nodes(elem); node++)
  {
    myNodes.push_back(*node);
  }

  std::vector<IntersectionPoint> intersectionPoints;

  for(unsigned i=0; i<3; i++)
  {
    IntersectionPoint intPt(true, {myNodes[0], myNodes[i+1]}, {1.0-1e-12, 1e-12}, {});
    intersectionPoints.push_back(intPt);
  }

  auto endInputPts = intersectionPoints.size();

  append_closest_point_intersections(intersectionPoints,
    bulk,
    test_fixture.stk_fixture.meta_data().universal_part(),
    test_fixture.coord_field,
    keep_all_intersection_points_filter());
  
  EXPECT_EQ(intersectionPoints.size(), endInputPts);
}

}
