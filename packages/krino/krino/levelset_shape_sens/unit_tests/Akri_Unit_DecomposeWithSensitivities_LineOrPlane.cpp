#include <gtest/gtest.h>

#include <Akri_Unit_DecomposeWithSensitivities_TestUtils.hpp>

#include <Akri_BoundingBoxMesh.hpp>
#include <Akri_LevelSetPolicy.hpp>
#include <Akri_LevelSetShapeSensitivities.hpp>
#include <Akri_NodalSurfaceDistance.hpp>
#include <Akri_UnitMeshUtils.hpp>
#include <Akri_UnitTestUtils.hpp>
#include <stk_topology/topology.hpp>

constexpr auto EDGE_DIRECTION_SENS = krino::ShapeSensitivityMethod::EDGE_DIRECTION_USING_EDGE_GRADIENT;
constexpr auto SURFACE_NORMAL_SENS = krino::ShapeSensitivityMethod::SURFACE_NORMAL_USING_EDGE_GRADIENT;

class DecomposeMeshAndComputeSensitivitiesForLineOrPlane : public ::testing::Test
{
protected:
  void build_line_or_plane_conforming_mesh_and_test_sensitivity(const int dim,
      const double meshSize,
      const double planeOffset,
      const krino::ShapeSensitivityMethod method,
      const bool doSnapping,
      const bool doWriteMesh)
  {
    const stk::topology elemTopo = (dim == 2) ? stk::topology::TRIANGLE_3_2D : stk::topology::TETRAHEDRON_4;
    krino::BoundingBoxMesh bboxMesh(elemTopo, MPI_COMM_WORLD);

    const std::vector<krino::LS_Field> lsFields = krino::LSPerInterfacePolicy::setup_levelsets_on_all_blocks_with_void_phase_for_any_negative_levelset(bboxMesh.meta_data(), 1);
    krino::ConformingDecompositionWithSensitivities decomposition(
        bboxMesh.meta_data(), lsFields, doSnapping);

    krino::populate_bounding_box_mesh_and_activate(bboxMesh, {0.,0.,0.}, {1.,1.,1.}, meshSize);

    compute_nodal_distance_from_plane(bboxMesh.bulk_data(), bboxMesh.meta_data().coordinate_field(), lsFields[0].isovar, {1.,0.,0.}, planeOffset);

    decomposition.decompose();

    std::vector<krino::LevelSetShapeSensitivity> sensitivities =
        decomposition.sensitivities(0, method);

    if (method == krino::ShapeSensitivityMethod::SURFACE_NORMAL_USING_EDGE_GRADIENT)
    {
      for (auto & sens : sensitivities)
        test_normal_sensitivity_for_plane(sens, {-1., 0., 0.});
    }
    else
    {
      for (auto & sens : sensitivities)
        test_sensitivity_for_plane(sens, 0);
    }

    if (doWriteMesh)
    {
      const std::string filename = std::string(dim == 2 ? "line" : "plane") + std::string(doSnapping ? "_snapped.e" : "_cut.e");
      output_full_mesh(bboxMesh.bulk_data(), filename, 1, 0.0);
    }
  }

  void test_sensitivity_for_plane(const krino::LevelSetShapeSensitivity & sens, const int iPlaneCoordDir)
  {
    double sum = 0;
    for (auto & dCoordsdParentLevelSet : sens.dCoordsdParentLevelSets)
      sum += dCoordsdParentLevelSet[iPlaneCoordDir];

    EXPECT_NEAR(-1., sum, 1.e-4);
  }

  // The surface normal of a plane is its normal, so each node moves along it at unit speed.
  void test_normal_sensitivity_for_plane(
      const krino::LevelSetShapeSensitivity & sens, const stk::math::Vector3d & goldSens)
  {
    stk::math::Vector3d sumSens = stk::math::Vector3d::ZERO;
    for (auto & dCoordsdParentLevelSet : sens.dCoordsdParentLevelSets)
      sumSens += dCoordsdParentLevelSet;

    krino::expect_near(goldSens, sumSens);
  }
};

TEST_F(DecomposeMeshAndComputeSensitivitiesForLineOrPlane, createDecomposedMeshForPlaneNotThroughAnyBackgroundNodes_testSensitivities)
{
  build_line_or_plane_conforming_mesh_and_test_sensitivity(
      3, 0.3333, -0.2, EDGE_DIRECTION_SENS, false, false);
}

TEST_F(DecomposeMeshAndComputeSensitivitiesForLineOrPlane, createSnappedMeshForPlaneNotThroughAnyBackgroundNodes_testSensitivities)
{
  build_line_or_plane_conforming_mesh_and_test_sensitivity(
      3, 0.3333, -0.2, EDGE_DIRECTION_SENS, true, false);
}

TEST_F(DecomposeMeshAndComputeSensitivitiesForLineOrPlane, createDecomposedMeshForPlaneThroughSomeBackgroundNodes_testSensitivities)
{
  build_line_or_plane_conforming_mesh_and_test_sensitivity(
      3, 0.5, -0.5, EDGE_DIRECTION_SENS, false, false);
}

TEST_F(DecomposeMeshAndComputeSensitivitiesForLineOrPlane, createSnappedMeshForPlaneThroughSomeBackgroundNodes_testSensitivities)
{
  build_line_or_plane_conforming_mesh_and_test_sensitivity(
      3, 0.5, -0.5, EDGE_DIRECTION_SENS, true, false);
}

TEST_F(DecomposeMeshAndComputeSensitivitiesForLineOrPlane, createDecomposedMeshForLineNotThroughAnyBackgroundNodes_testSensitivities)
{
  build_line_or_plane_conforming_mesh_and_test_sensitivity(
      2, 0.3333, -0.2, EDGE_DIRECTION_SENS, false, false);
}

TEST_F(DecomposeMeshAndComputeSensitivitiesForLineOrPlane, createSnappedMeshForLineNotThroughAnyBackgroundNodes_testSensitivities)
{
  build_line_or_plane_conforming_mesh_and_test_sensitivity(
      2, 0.3333, -0.2, EDGE_DIRECTION_SENS, true, false);
}

TEST_F(DecomposeMeshAndComputeSensitivitiesForLineOrPlane, createDecomposedMeshForLineThroughSomeBackgroundNodes_testSensitivities)
{
  build_line_or_plane_conforming_mesh_and_test_sensitivity(
      2, 0.5, -0.5, EDGE_DIRECTION_SENS, false, false);
}

TEST_F(DecomposeMeshAndComputeSensitivitiesForLineOrPlane, createSnappedMeshForLineThroughSomeBackgroundNodes_testSensitivities)
{
  build_line_or_plane_conforming_mesh_and_test_sensitivity(
      2, 0.5, -0.5, EDGE_DIRECTION_SENS, true, false);
}

TEST_F(DecomposeMeshAndComputeSensitivitiesForLineOrPlane,
    createDecomposedMeshForPlaneNotThroughAnyBackgroundNodes_testSurfaceNormalSensitivities)
{
  build_line_or_plane_conforming_mesh_and_test_sensitivity(
      3, 0.3333, -0.2, SURFACE_NORMAL_SENS, false, false);
}

TEST_F(DecomposeMeshAndComputeSensitivitiesForLineOrPlane,
    createSnappedMeshForPlaneNotThroughAnyBackgroundNodes_testSurfaceNormalSensitivities)
{
  build_line_or_plane_conforming_mesh_and_test_sensitivity(
      3, 0.3333, -0.2, SURFACE_NORMAL_SENS, true, false);
}

TEST_F(DecomposeMeshAndComputeSensitivitiesForLineOrPlane,
    createDecomposedMeshForPlaneThroughSomeBackgroundNodes_testSurfaceNormalSensitivities)
{
  build_line_or_plane_conforming_mesh_and_test_sensitivity(
      3, 0.5, -0.5, SURFACE_NORMAL_SENS, false, false);
}

TEST_F(DecomposeMeshAndComputeSensitivitiesForLineOrPlane,
    createSnappedMeshForPlaneThroughSomeBackgroundNodes_testSurfaceNormalSensitivities)
{
  build_line_or_plane_conforming_mesh_and_test_sensitivity(
      3, 0.5, -0.5, SURFACE_NORMAL_SENS, true, false);
}

TEST_F(DecomposeMeshAndComputeSensitivitiesForLineOrPlane,
    createDecomposedMeshForLineNotThroughAnyBackgroundNodes_testSurfaceNormalSensitivities)
{
  build_line_or_plane_conforming_mesh_and_test_sensitivity(
      2, 0.3333, -0.2, SURFACE_NORMAL_SENS, false, false);
}

TEST_F(DecomposeMeshAndComputeSensitivitiesForLineOrPlane,
    createSnappedMeshForLineNotThroughAnyBackgroundNodes_testSurfaceNormalSensitivities)
{
  build_line_or_plane_conforming_mesh_and_test_sensitivity(
      2, 0.3333, -0.2, SURFACE_NORMAL_SENS, true, false);
}

TEST_F(DecomposeMeshAndComputeSensitivitiesForLineOrPlane,
    createDecomposedMeshForLineThroughSomeBackgroundNodes_testSurfaceNormalSensitivities)
{
  build_line_or_plane_conforming_mesh_and_test_sensitivity(
      2, 0.5, -0.5, SURFACE_NORMAL_SENS, false, false);
}

TEST_F(DecomposeMeshAndComputeSensitivitiesForLineOrPlane,
    createSnappedMeshForLineThroughSomeBackgroundNodes_testSurfaceNormalSensitivities)
{
  build_line_or_plane_conforming_mesh_and_test_sensitivity(
      2, 0.5, -0.5, SURFACE_NORMAL_SENS, true, false);
}
