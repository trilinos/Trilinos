// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <stk_math/StkVector.hpp>
#include <Akri_CurvatureLeastSquares.hpp>
#include "Akri_UnitTestUtils.hpp"

namespace krino
{

class NodePatchInterface
{
public:
  virtual ~NodePatchInterface() = default;
  virtual const std::vector<stk::math::Vector3d>  & get_halo_node_locations() const = 0;
  virtual const std::vector<std::array<int,2>>& get_halo_segments() const = 0;
};

void test_rotation(const stk::math::Vector3d & normal)
{
  static const stk::math::Vector3d zDir(0.,0.,1.);

  std::array<std::array<double,3>,3> rotationMatrix;
  set_rotation_matrix_for_rotating_normal_to_zDir(rotationMatrix, normal);

  const stk::math::Vector3d rotatedNormal = rotate_3d_vector(rotationMatrix, normal);

  expect_eq(zDir, rotatedNormal);

  const stk::math::Vector3d reverseRotatedZdir = reverse_rotate_3d_vector(rotationMatrix, zDir);

  expect_eq(normal, reverseRotatedZdir);
}

TEST(CurvatureLeastSquares, rotationTest)
{
  test_rotation(stk::math::Vector3d(0.,0.,1.));
  test_rotation(stk::math::Vector3d(1.,0.,0.));
  test_rotation(stk::math::Vector3d(0.,1.,0.));

  const double cos45 = std::sqrt(2.)/2.;
  test_rotation(stk::math::Vector3d(cos45, cos45, 0.));
  test_rotation(stk::math::Vector3d(cos45, 0., cos45));
  test_rotation(stk::math::Vector3d(0., cos45, cos45));
}

class PolygonalPatchOnSphere : public NodePatchInterface
{
public:
  PolygonalPatchOnSphere(const stk::math::Vector3d & normalDir, const double halfCurvature, const int numHaloPts)
  {
    const double radius = 2.0/halfCurvature;

    std::array<std::array<double,3>,3> rotationMatrix;
    set_rotation_matrix_for_rotating_normal_to_zDir(rotationMatrix, normalDir);

    const stk::math::Vector3d unrotatedNodeLoc(0.,0.,radius);
    myNodeLoc = reverse_rotate_3d_vector(rotationMatrix, unrotatedNodeLoc);

    const double phi = 15.*M_PI/180.;
    myHaloNodeLocs.reserve(numHaloPts);
    myHaloSegments.reserve(numHaloPts);
    const double dTheta = 2.*M_PI/numHaloPts;
    for (int i=0; i<numHaloPts; ++i)
    {
      const stk::math::Vector3d unrotatedPtLoc(radius*std::cos(i*dTheta)*std::sin(phi),radius*std::sin(i*dTheta)*std::sin(phi),radius*std::cos(phi));
      myHaloNodeLocs.push_back(reverse_rotate_3d_vector(rotationMatrix, unrotatedPtLoc) - myNodeLoc);
      myHaloSegments.push_back({{i,(i+1)%numHaloPts}});
    }
  }
  virtual const std::vector<stk::math::Vector3d> & get_halo_node_locations() const override { return myHaloNodeLocs; }
  virtual const std::vector<std::array<int,2>> & get_halo_segments() const override { return myHaloSegments; }

private:
  stk::math::Vector3d myNodeLoc;
  std::vector<stk::math::Vector3d> myHaloNodeLocs;
  std::vector<std::array<int,2>> myHaloSegments;
};

class PolygonalPatchOnPlane : public NodePatchInterface
{
public:
  PolygonalPatchOnPlane(const stk::math::Vector3d & normalDir, const int numHaloPts)
  {
    std::array<std::array<double,3>,3> rotationMatrix;
    set_rotation_matrix_for_rotating_normal_to_zDir(rotationMatrix, normalDir);

    myHaloNodeLocs.reserve(numHaloPts);
    myHaloSegments.reserve(numHaloPts);
    const double dTheta = 2.*M_PI/numHaloPts;
    for (int i=0; i<numHaloPts; ++i)
    {
      const stk::math::Vector3d unrotatedPtLoc(std::cos(i*dTheta),std::sin(i*dTheta),0.);
      myHaloNodeLocs.push_back(reverse_rotate_3d_vector(rotationMatrix, unrotatedPtLoc));
      myHaloSegments.push_back({{i,(i+1)%numHaloPts}});
    }
  }
  virtual const std::vector<stk::math::Vector3d> & get_halo_node_locations() const override { return myHaloNodeLocs; }
  virtual const std::vector<std::array<int,2>> & get_halo_segments() const override { return myHaloSegments; }

private:
  std::vector<stk::math::Vector3d> myHaloNodeLocs;
  std::vector<std::array<int,2>> myHaloSegments;
};

void test_flat_triangle_patch_with_normal_gives_zero_normalCurvature(const stk::math::Vector3d & normal)
{
  PolygonalPatchOnPlane patch(normal, 3);

  const stk::math::Vector3d normalCurvature = compute_least_squares_curvature_times_normal(patch.get_halo_node_locations(), patch.get_halo_segments());
  expect_eq_absolute(stk::math::Vector3d::ZERO, normalCurvature, 1.e-6);
}

TEST(CurvatureLeastSquares, Flat3TrianglePatches_zeroNormalCurvature)
{
  test_flat_triangle_patch_with_normal_gives_zero_normalCurvature(stk::math::Vector3d(0.,0.,1.));
  test_flat_triangle_patch_with_normal_gives_zero_normalCurvature(stk::math::Vector3d(1.,0.,0.));
  test_flat_triangle_patch_with_normal_gives_zero_normalCurvature(stk::math::Vector3d(0.,1.,0.));

  const double cos45 = std::sqrt(2.)/2.;
  test_flat_triangle_patch_with_normal_gives_zero_normalCurvature(stk::math::Vector3d(cos45, cos45, 0.));
  test_flat_triangle_patch_with_normal_gives_zero_normalCurvature(stk::math::Vector3d(cos45, 0., cos45));
}

void test_normalCurvature_for_curved_patch(const stk::math::Vector3d & normalDir, const double curvature, const int numHaloNodes)
{
  const stk::math::Vector3d normal = normalDir.unit_vector();
  PolygonalPatchOnSphere patch(normal, curvature, numHaloNodes);
  const stk::math::Vector3d goldNormalCurvature = curvature*normal;

  const stk::math::Vector3d normalCurvature = compute_least_squares_curvature_times_normal(patch.get_halo_node_locations(), patch.get_halo_segments());
  expect_eq(goldNormalCurvature, normalCurvature, 1.e-2);
}

void test_normalCurvature_for_curved_patches(const int numHaloNodes)
{
  const double curvature = 0.1;
  test_normalCurvature_for_curved_patch(stk::math::Vector3d(0.,0.,1.), curvature, numHaloNodes);
  test_normalCurvature_for_curved_patch(stk::math::Vector3d(1.,0.,0.), curvature, numHaloNodes);
  test_normalCurvature_for_curved_patch(stk::math::Vector3d(0.,1.,0.), curvature, numHaloNodes);

  const double cos45 = std::sqrt(2.)/2.;
  test_normalCurvature_for_curved_patch(stk::math::Vector3d(cos45, cos45, 0.), curvature, numHaloNodes);
  test_normalCurvature_for_curved_patch(stk::math::Vector3d(cos45, 0., cos45), curvature, numHaloNodes);
}

TEST(CurvatureLeastSquares, CurvedPatchesOfVariousSizes_correctNormalCurvature)
{
  // curvature only fit
  test_normalCurvature_for_curved_patches(3);
  test_normalCurvature_for_curved_patches(4);

  // curvature and normal fit
  test_normalCurvature_for_curved_patches(5);
  test_normalCurvature_for_curved_patches(7);
}

} // namespace krino
