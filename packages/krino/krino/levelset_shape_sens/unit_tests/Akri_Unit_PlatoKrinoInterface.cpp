#include <gtest/gtest.h>

#include <cmath>
#include <string>

// Deliberately the ONLY krino header this test includes, plus the shared unit
// test helpers that write and remove a temp mesh file. This exercises (and
// documents) the contract that a consumer needs nothing but Akri_PlatoKrinoInterface.hpp.
#include <Akri_PlatoKrinoInterface.hpp>

#include <Akri_ScopedTestFile.hpp>
#include <Akri_UnitMeshUtils.hpp>
#include <stk_topology/topology.hpp>

namespace
{
using namespace PlatoKrinoInterface;

LevelSetMesh setup_sphere_mesh(const std::string & meshFile, const double meshSize, const double radius,
    const bool doSnapping, const bool doClosestPoint)
{
  krino::generate_and_write_bounding_box_mesh(stk::topology::TETRAHEDRON_4, {-1.,-1.,-1.}, {1.,1.,1.}, meshSize, meshFile);

  // generate_and_write_bounding_box_mesh already writes pre-spread per-rank
  // files, so read them back with the identity decomposition ("") rather than
  // re-decomposing a single serial file (the "rib" default of the facade).
  LevelSetMesh mesh = read_and_setup_mesh_for_levelset_decomposition(meshFile, MPI_COMM_WORLD, 1, doSnapping, doClosestPoint, "");

  Geometry geometry;
  geometry.spheres.push_back(SphereSpec{{{0.,0.,0.}}, radius});
  initialize_levelset_field_from_geometry(mesh, geometry);

  decompose_mesh(mesh);
  return mesh;
}

TEST(PlatoKrinoInterface, setupDecomposeAndGetSensitivitiesForSphere_returnsSensitivitiesWithTwoParents)
{
  const std::string meshFile = "plato_interface_sphere.g";
  const krino::ScopedTestFile cleanup(meshFile);

  const double radius = 0.5;
  LevelSetMesh mesh = setup_sphere_mesh(meshFile, 0.4, radius, false, false);

  EXPECT_EQ(3u, mesh.spatial_dimension());
  EXPECT_EQ(1u, mesh.num_level_sets());

  const std::vector<ShapeSensitivity> sensitivities = get_shape_sensitivities(mesh);
  EXPECT_FALSE(sensitivities.empty());

  for (const auto & sens : sensitivities)
  {
    EXPECT_EQ(2u, sens.parentNodeIds.size());
    EXPECT_EQ(2u, sens.dCoordsdParentLevelSets.size());
    EXPECT_NE(0u, sens.interfaceNodeId);
  }
}

TEST(PlatoKrinoInterface, closestPointSensitivitiesForSphere_summedSensitivityIsInwardRadialUnit)
{
  const std::string meshFile = "plato_interface_sphere_cp.g";
  const krino::ScopedTestFile cleanup(meshFile);

  const double radius = 0.5;
  LevelSetMesh mesh = setup_sphere_mesh(meshFile, 0.2, radius, false, true);

  const std::vector<ShapeSensitivity> sensitivities = get_shape_sensitivities(mesh);
  EXPECT_FALSE(sensitivities.empty());

  // Summed closest-point sensitivity should be the inward unit normal (magnitude 1).
  for (const auto & sens : sensitivities)
  {
    std::array<double,3> sum{{0.,0.,0.}};
    for (const auto & d : sens.dCoordsdParentLevelSets)
      for (int i = 0; i < 3; ++i)
        sum[i] += d[i];
    const double mag = std::sqrt(sum[0]*sum[0] + sum[1]*sum[1] + sum[2]*sum[2]);
    EXPECT_NEAR(1.0, mag, 0.1);
  }
}

TEST(PlatoKrinoInterface, initializeFromFunction_matchesGeometryInitialization)
{
  const std::string meshFile = "plato_interface_fn.g";
  const krino::ScopedTestFile cleanup(meshFile);

  const double radius = 0.5;
  krino::generate_and_write_bounding_box_mesh(stk::topology::TETRAHEDRON_4, {-1.,-1.,-1.}, {1.,1.,1.}, 0.4, meshFile);

  LevelSetMesh mesh = read_and_setup_mesh_for_levelset_decomposition(meshFile, MPI_COMM_WORLD, 1, false, false, "");

  // Signed distance to a sphere of the given radius, stamped through the facade's
  // function initializer (an optimizer would supply its own design values here).
  const LevelSetValueFunction sphereDistance = [radius](const std::array<double,3> & x)
  {
    return std::sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]) - radius;
  };
  initialize_levelset_field_from_function(mesh, sphereDistance);

  decompose_mesh(mesh);

  const std::vector<ShapeSensitivity> sensitivities = get_shape_sensitivities(mesh);
  EXPECT_FALSE(sensitivities.empty());
  for (const auto & sens : sensitivities)
    EXPECT_EQ(2u, sens.parentNodeIds.size());
}

}
