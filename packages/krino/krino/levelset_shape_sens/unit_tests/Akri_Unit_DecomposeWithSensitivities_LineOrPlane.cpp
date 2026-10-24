#include <gtest/gtest.h>

#include <Akri_Unit_DecomposeWithSensitivities_TestUtils.hpp>

#include <Akri_BoundingBoxMesh.hpp>
#include <Akri_LevelSetPolicy.hpp>
#include <Akri_LevelSetShapeSensitivities.hpp>
#include <Akri_NodalSurfaceDistance.hpp>
#include <Akri_UnitMeshUtils.hpp>
#include <Akri_UnitTestUtils.hpp>
#include <stk_topology/topology.hpp>

class DecomposeMeshAndComputeSensitivitiesForLineOrPlane : public ::testing::Test
{
protected:
  void build_line_or_plane_conforming_mesh_and_test_sensitivity(const int dim, const double meshSize, const double planeOffset, const bool doComputeClosestPointSensitivities, const bool doSnapping, const bool doWriteMesh)
  {
    const stk::topology elemTopo = (dim == 2) ? stk::topology::TRIANGLE_3_2D : stk::topology::TETRAHEDRON_4;
    krino::BoundingBoxMesh bboxMesh(elemTopo, MPI_COMM_WORLD);

    const std::vector<krino::LS_Field> lsFields = krino::LSPerInterfacePolicy::setup_levelsets_on_all_blocks_with_void_phase_for_any_negative_levelset(bboxMesh.meta_data(), 1);
    setup_fields_for_conforming_decomposition(bboxMesh.meta_data(), lsFields, doSnapping, doComputeClosestPointSensitivities);

    krino::populate_bounding_box_mesh_and_activate(bboxMesh, {0.,0.,0.}, {1.,1.,1.}, meshSize);

    compute_nodal_distance_from_plane(bboxMesh.bulk_data(), bboxMesh.meta_data().coordinate_field(), lsFields[0].isovar, {1.,0.,0.}, planeOffset);

    decompose_mesh_to_conform_to_levelsets(bboxMesh.bulk_data(), lsFields, doSnapping, doComputeClosestPointSensitivities);

    std::vector<krino::LevelSetShapeSensitivity> sensitivities = get_levelset_shape_sensitivities(bboxMesh.bulk_data(), lsFields[0], doComputeClosestPointSensitivities);

    if (doComputeClosestPointSensitivities)
    {
      for (auto & sens : sensitivities)
        test_closest_point_sensitivity_for_plane(sens, {-1.,0.,0.});
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

  void test_closest_point_sensitivity_for_plane(const krino::LevelSetShapeSensitivity & sens, const stk::math::Vector3d & goldSens)
  {
    stk::math::Vector3d sumSens = stk::math::Vector3d::ZERO;
    for (auto & dCoordsdParentLevelSet : sens.dCoordsdParentLevelSets)
      sumSens += dCoordsdParentLevelSet;

    krino::expect_near(goldSens, sumSens);
  }
};

TEST_F(DecomposeMeshAndComputeSensitivitiesForLineOrPlane, createDecomposedMeshForPlaneNotThroughAnyBackgroundNodes_testSensitivities)
{
  build_line_or_plane_conforming_mesh_and_test_sensitivity(3, 0.3333, -0.2, false, false, false);
}

TEST_F(DecomposeMeshAndComputeSensitivitiesForLineOrPlane, createSnappedMeshForPlaneNotThroughAnyBackgroundNodes_testSensitivities)
{
  build_line_or_plane_conforming_mesh_and_test_sensitivity(3, 0.3333, -0.2, false, true, false);
}

TEST_F(DecomposeMeshAndComputeSensitivitiesForLineOrPlane, createDecomposedMeshForPlaneThroughSomeBackgroundNodes_testSensitivities)
{
  build_line_or_plane_conforming_mesh_and_test_sensitivity(3, 0.5, -0.5, false, false, false);
}

TEST_F(DecomposeMeshAndComputeSensitivitiesForLineOrPlane, createSnappedMeshForPlaneThroughSomeBackgroundNodes_testSensitivities)
{
  build_line_or_plane_conforming_mesh_and_test_sensitivity(3, 0.5, -0.5, false, true, false);
}

TEST_F(DecomposeMeshAndComputeSensitivitiesForLineOrPlane, createDecomposedMeshForLineNotThroughAnyBackgroundNodes_testSensitivities)
{
  build_line_or_plane_conforming_mesh_and_test_sensitivity(2, 0.3333, -0.2, false, false, false);
}

TEST_F(DecomposeMeshAndComputeSensitivitiesForLineOrPlane, createSnappedMeshForLineNotThroughAnyBackgroundNodes_testSensitivities)
{
  build_line_or_plane_conforming_mesh_and_test_sensitivity(2, 0.3333, -0.2, false, true, false);
}

TEST_F(DecomposeMeshAndComputeSensitivitiesForLineOrPlane, createDecomposedMeshForLineThroughSomeBackgroundNodes_testSensitivities)
{
  build_line_or_plane_conforming_mesh_and_test_sensitivity(2, 0.5, -0.5, false, false, false);
}

TEST_F(DecomposeMeshAndComputeSensitivitiesForLineOrPlane, createSnappedMeshForLineThroughSomeBackgroundNodes_testSensitivities)
{
  build_line_or_plane_conforming_mesh_and_test_sensitivity(2, 0.5, -0.5, false, true, false);
}

// closest point version of plane tests

TEST_F(DecomposeMeshAndComputeSensitivitiesForLineOrPlane, createDecomposedMeshForPlaneNotThroughAnyBackgroundNodes_testClosestPointSensitivities)
{
  build_line_or_plane_conforming_mesh_and_test_sensitivity(3, 0.3333, -0.2, true, false, false);
}

TEST_F(DecomposeMeshAndComputeSensitivitiesForLineOrPlane, createSnappedMeshForPlaneNotThroughAnyBackgroundNodes_testClosestPointSensitivities)
{
  build_line_or_plane_conforming_mesh_and_test_sensitivity(3, 0.3333, -0.2, true, true, false);
}

TEST_F(DecomposeMeshAndComputeSensitivitiesForLineOrPlane, createDecomposedMeshForPlaneThroughSomeBackgroundNodes_testClosestPointSensitivities)
{
  build_line_or_plane_conforming_mesh_and_test_sensitivity(3, 0.5, -0.5, true, false, false);
}

TEST_F(DecomposeMeshAndComputeSensitivitiesForLineOrPlane, createSnappedMeshForPlaneThroughSomeBackgroundNodes_testClosestPointSensitivities)
{
  build_line_or_plane_conforming_mesh_and_test_sensitivity(3, 0.5, -0.5, true, true, false);
}

TEST_F(DecomposeMeshAndComputeSensitivitiesForLineOrPlane, createDecomposedMeshForLineNotThroughAnyBackgroundNodes_testClosestPointSensitivities)
{
  build_line_or_plane_conforming_mesh_and_test_sensitivity(2, 0.3333, -0.2, true, false, false);
}

TEST_F(DecomposeMeshAndComputeSensitivitiesForLineOrPlane, createSnappedMeshForLineNotThroughAnyBackgroundNodes_testClosestPointSensitivities)
{
  build_line_or_plane_conforming_mesh_and_test_sensitivity(2, 0.3333, -0.2, true, true, false);
}

TEST_F(DecomposeMeshAndComputeSensitivitiesForLineOrPlane, createDecomposedMeshForLineThroughSomeBackgroundNodes_testClosestPointSensitivities)
{
  build_line_or_plane_conforming_mesh_and_test_sensitivity(2, 0.5, -0.5, true, false, false);
}

TEST_F(DecomposeMeshAndComputeSensitivitiesForLineOrPlane, createSnappedMeshForLineThroughSomeBackgroundNodes_testClosestPointSensitivities)
{
  build_line_or_plane_conforming_mesh_and_test_sensitivity(2, 0.5, -0.5, true, true, false);
}
