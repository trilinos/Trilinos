#include <gtest/gtest.h>

#include <Akri_AllReduce.hpp>
#include <Akri_BoundingBoxMesh.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_LevelSetPolicy.hpp>
#include <Akri_LevelSetShapeSensitivities.hpp>
#include <Akri_NodalSurfaceDistance.hpp>
#include <Akri_UnitMeshUtils.hpp>
#include <Akri_UnitTestUtils.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_topology/topology.hpp>

#include <memory>
#include <vector>

// Tests of the SURFACE_NORMAL_USING_EDGE_GRADIENT shape sensitivities on a plane, where the surface
// normal is the plane's normal: every node moves along the normal, at the edge direction form's
// normal speed.  Volume sensitivities against the edge direction form are checked in the finite
// difference tests.

namespace
{

constexpr auto SURFACE_NORMAL_EDGE_GRADIENT = krino::ShapeSensitivityMethod::SURFACE_NORMAL_USING_EDGE_GRADIENT;

void test_plane(const unsigned dim, const double offset, const bool doSnapping)
{
  const stk::topology elemTopo =
      (dim == 2) ? stk::topology::TRIANGLE_3_2D : stk::topology::TETRAHEDRON_4;
  auto bboxMesh = std::make_unique<krino::BoundingBoxMesh>(elemTopo, MPI_COMM_WORLD);
  const std::vector<krino::LS_Field> lsFields = krino::LSPerInterfacePolicy::
      setup_levelsets_on_all_blocks_with_void_phase_for_any_negative_levelset(
          bboxMesh->meta_data(), 1);
  krino::ConformingDecompositionWithSensitivities decomposition(
      bboxMesh->meta_data(), lsFields, doSnapping);
  krino::populate_bounding_box_mesh_and_activate(*bboxMesh, {-1., -1., -1.}, {1., 1., 1.}, 0.25);
  stk::mesh::BulkData & mesh = bboxMesh->bulk_data();

  const stk::math::Vector3d normal =
      (dim == 2) ? stk::math::Vector3d(0.6, 0.8, 0.) : stk::math::Vector3d(0.48, 0.6, 0.64);
  compute_nodal_distance_from_plane(
      mesh, mesh.mesh_meta_data().coordinate_field(), lsFields[0].isovar, normal, offset);
  decomposition.decompose();

  const std::vector<krino::LevelSetShapeSensitivity> sensitivities =
      decomposition.sensitivities(0, SURFACE_NORMAL_EDGE_GRADIENT);

  size_t numSensitivities = sensitivities.size();
  krino::all_reduce_sum(mesh.parallel(), numSensitivities);
  EXPECT_GT(numSensitivities, 0u);

  for (auto & sens : sensitivities)
  {
    stk::math::Vector3d sum = stk::math::Vector3d::ZERO;
    for (auto & d : sens.dCoordsdParentLevelSets)
    {
      krino::expect_near_absolute(Dot(d, normal) * normal, d, 1.e-12);
      sum += d;
    }
    krino::expect_near_absolute(-normal, sum, 1.e-12);
  }
}

} // namespace

TEST(SurfaceNormalSensitivities, plane2d_offsetFromNodes) { test_plane(2, 0.1, false); }

TEST(SurfaceNormalSensitivities, plane2d_offsetFromNodes_snapping) { test_plane(2, 0.1, true); }

TEST(SurfaceNormalSensitivities, plane2d_nearNodes_snapping) { test_plane(2, 0.0, true); }

TEST(SurfaceNormalSensitivities, plane3d_offsetFromNodes) { test_plane(3, 0.1, false); }

TEST(SurfaceNormalSensitivities, plane3d_offsetFromNodes_snapping) { test_plane(3, 0.1, true); }

TEST(SurfaceNormalSensitivities, plane3d_nearNodes_snapping) { test_plane(3, 0.0, true); }
