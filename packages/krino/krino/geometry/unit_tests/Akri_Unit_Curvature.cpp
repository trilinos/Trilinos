/*
 * Akri_Unit_Curvature.cpp
 *
 *  Created on: Jan 27, 2026
 *      Author: drnoble
 */

#include <Akri_AnalyticSurf.hpp>
#include <gtest/gtest.h>

#include <stk_math/StkVector.hpp>
#include <Akri_MeshSpecs.hpp>
#include <Akri_SimplexCurvature.hpp>
#include <Akri_StkMeshFixture.hpp>

namespace krino
{

TEST(LaplacianCurvature, patchOfTetsAroundNode_getFiniteElementCurvature)
{
  OctahedralPatchAroundNode meshSpec;
  const double radius = 1.;
  const std::vector<double> phi = {0., radius, radius, radius, radius, radius, radius};
  const double curvature = laplacian_of_field_at_node_of_patch_of_tetrahedra(0, meshSpec.allElementConn, meshSpec.nodeLocs, phi);
  EXPECT_NEAR(12., curvature, 1.e-3);
}

TEST(LaplacianCurvature, patchOfTrisAroundNode_getFiniteElementCurvature)
{
  PatchOfRegularTrisAroundNode meshSpec;
  const double radius = 1.;
  const std::vector<double> phi = {0., radius, radius, radius, radius, radius, radius};
  std::vector<stk::math::Vector3d> nodeLocs;
  for (auto & nodeLoc2d : meshSpec.nodeLocs)
    nodeLocs.emplace_back(nodeLoc2d[0], nodeLoc2d[1], 0.);
  const double curvature = laplacian_of_field_at_node_of_patch_of_triangles(0, meshSpec.allElementConn, nodeLocs, phi);
  EXPECT_NEAR(4., curvature, 1.e-3);
}

double finite_difference_curvature_2d(const Surface & surf, const stk::math::Vector3d & x0, const double delta)
{
  auto f = [&surf](const stk::math::Vector3d & x) { return surf.point_signed_distance(x); };
  const stk::math::Vector3d dx(delta,0,0);
  const stk::math::Vector3d dy(0,delta,0);
  const double normalDivergence = finite_difference_divergence_of_normal_2d(delta, f(x0),
      f(x0-dx), f(x0+dx), f(x0-dy), f(x0+dy),
      f(x0-dx-dy), f(x0-dx+dy), f(x0+dx-dy), f(x0+dx+dy));
  return normalDivergence;
}

double finite_difference_mean_curvature_3d(const Surface & surf, const stk::math::Vector3d & x0, const double delta)
{
  auto f = [&surf](const stk::math::Vector3d & x) { return surf.point_signed_distance(x); };
  const stk::math::Vector3d dx(delta,0,0);
  const stk::math::Vector3d dy(0,delta,0);
  const stk::math::Vector3d dz(0,0,delta);
  const double normalDivergence = finite_difference_divergence_of_normal_3d(delta, f(x0),
      f(x0-dx), f(x0+dx), f(x0-dy), f(x0+dy), f(x0-dz), f(x0+dz),
      f(x0-dx-dy), f(x0-dx+dy), f(x0+dx-dy), f(x0+dx+dy),
      f(x0-dx-dz), f(x0-dx+dz), f(x0+dx-dz), f(x0+dx+dz),
      f(x0-dy-dz), f(x0-dy+dz), f(x0+dy-dz), f(x0+dy+dz));
  return 0.5*normalDivergence;
}

std::vector<stk::math::Vector3d> get_scaled_node_locations(const std::vector<stk::math::Vector3d> & origNodeLocations, const double scale)
{
  std::vector<stk::math::Vector3d> nodeLocs;
  nodeLocs.reserve(origNodeLocations.size());
  for (auto & nodeLoc : origNodeLocations)
    nodeLocs.push_back(scale*nodeLoc);
  return nodeLocs;
}

std::vector<stk::math::Vector3d> get_closest_point_normals(const Surface & surface, const std::vector<stk::math::Vector3d> & nodeLocations)
{
  std::vector<stk::math::Vector3d> nodeNormals;
  nodeNormals.reserve(nodeLocations.size());
  for (auto & nodeLoc : nodeLocations)
    nodeNormals.push_back(surface.closest_point_normal(nodeLoc));
  return nodeNormals;
}

double curvature_from_node_normals_2d(const Surface & surface, const std::vector<std::array<unsigned,3>> & triConn, const std::vector<stk::math::Vector3d> & nodeLocs)
{
  const std::vector<stk::math::Vector3d> nodeNormals = get_closest_point_normals(surface, nodeLocs);
  return divergence_of_node_normals_on_patch_of_elements(triConn, nodeLocs, nodeNormals);
}

double mean_curvature_from_node_normals_3d(const Surface & surface, const std::vector<std::array<unsigned,4>> & tetConn, const std::vector<stk::math::Vector3d> & nodeLocs)
{
  const std::vector<stk::math::Vector3d> nodeNormals = get_closest_point_normals(surface, nodeLocs);
  return 0.5*divergence_of_node_normals_on_patch_of_elements(tetConn, nodeLocs, nodeNormals);
}

void expect_FD_and_FEM_sphere_curvature(const std::vector<std::array<unsigned,4>> & tetConn, const std::vector<stk::math::Vector3d> & nodeLocs, const unsigned centroidNodeIndex, const double delta, const double errTol)
{
  const auto scaledNodeLocations = get_scaled_node_locations(nodeLocs, delta);
  const double sphereRadius = 3.;
  const stk::math::Vector3d dir{1.,2.,3.};
  const stk::math::Vector3d sphereCenter = sphereRadius * dir/dir.length();
  const Sphere sphere(sphereCenter, sphereRadius);

  const double finiteDiffCurvature = finite_difference_mean_curvature_3d(sphere, scaledNodeLocations[centroidNodeIndex], delta);
  EXPECT_NEAR(1./sphereRadius, finiteDiffCurvature, errTol);

  const double nodeNormalCurvature = mean_curvature_from_node_normals_3d(sphere, tetConn, scaledNodeLocations);
  EXPECT_NEAR(1./sphereRadius, nodeNormalCurvature, errTol);
}

TEST(LaplacianCurvature, patchOfTetsAroundNode_expectFiniteDifferenceAndFiniteElementSphereCurvature)
{
  const double delta = 0.2;
  const double errTol = 2.e-4;
  {
    OctahedralPatchAroundNode meshSpec;
    expect_FD_and_FEM_sphere_curvature(meshSpec.allElementConn, meshSpec.nodeLocs, 0, delta, errTol);
    expect_FD_and_FEM_sphere_curvature(meshSpec.allElementConn, meshSpec.nodeLocs, 0, 0.5*delta, 0.25*errTol);
  }

  {
    PatchOf32Tets meshSpec;
    expect_FD_and_FEM_sphere_curvature(meshSpec.allElementConn, meshSpec.nodeLocs, 0, delta, errTol);
    expect_FD_and_FEM_sphere_curvature(meshSpec.allElementConn, meshSpec.nodeLocs, 0, 0.5*delta, 0.25*errTol);
  }

  {
    CubeOf24Tets meshSpec;
    expect_FD_and_FEM_sphere_curvature(meshSpec.allElementConn, meshSpec.nodeLocs, 14, delta, errTol);
    expect_FD_and_FEM_sphere_curvature(meshSpec.allElementConn, meshSpec.nodeLocs, 14, 0.5*delta, 0.25*errTol);
  }
}

void expect_FD_and_FEM_circle_curvature(const std::vector<std::array<unsigned,3>> & triConn, const std::vector<stk::math::Vector2d> & nodeLocs2d, const unsigned centroidNodeIndex, const double delta, const double errTol)
{
  std::vector<stk::math::Vector3d> nodeLocs;
  for (auto & nodeLoc2d : nodeLocs2d)
    nodeLocs.emplace_back(nodeLoc2d[0], nodeLoc2d[1], 0.);
  const auto scaledNodeLocations = get_scaled_node_locations(nodeLocs, delta);
  const double circleRadius = 3.;
  const stk::math::Vector3d dir{1.,2.,0.};
  const stk::math::Vector3d sphereCenter = circleRadius * dir/dir.length();
  const Sphere sphere(sphereCenter, circleRadius);

  const double finiteDiffCurvature = finite_difference_curvature_2d(sphere, scaledNodeLocations[centroidNodeIndex], delta);
  EXPECT_NEAR(1./circleRadius, finiteDiffCurvature, errTol);

  const double nodeNormalCurvature = curvature_from_node_normals_2d(sphere, triConn, scaledNodeLocations);
  EXPECT_NEAR(1./circleRadius, nodeNormalCurvature, errTol);
}

TEST(LaplacianCurvature, patchOfTrisAroundNode_expectFiniteDifferenceAndFiniteElementCircleCurvature)
{
  const double delta = 0.2;
  const double errTol = 2.e-4;
  PatchOfRegularTrisAroundNode meshSpec;
  expect_FD_and_FEM_circle_curvature(meshSpec.allElementConn, meshSpec.nodeLocs, 0, delta, errTol);
  expect_FD_and_FEM_circle_curvature(meshSpec.allElementConn, meshSpec.nodeLocs, 0, 0.5*delta, 0.25*errTol);
}


}


