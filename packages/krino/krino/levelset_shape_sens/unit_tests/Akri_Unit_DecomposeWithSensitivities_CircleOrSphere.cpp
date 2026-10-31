#include <Akri_AllReduce.hpp>
#include <gtest/gtest.h>

#include <Akri_Unit_DecomposeWithSensitivities_TestUtils.hpp>

#include <Akri_BoundingBoxMesh.hpp>
#include <Akri_CDFEM_Support.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_LevelSetShapeSensitivities.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_UnitTestUtils.hpp>
#include <stk_util/environment/EnvData.hpp>

static stk::math::Vector3d get_snap_displacements(const int dim, const stk::mesh::BulkData &bulk, const stk::mesh::EntityId node)
{
  krino::CDFEM_Support & cdfemSupport = krino::CDFEM_Support::get(bulk.mesh_meta_data());
  krino::FieldRef disp = cdfemSupport.get_cdfem_snap_displacements_field();
  if (disp.valid())
    return get_vector_field(bulk, disp, bulk.get_entity(stk::topology::NODE_RANK,node), dim);
  return stk::math::Vector3d::ZERO;
}

static stk::math::Vector3d get_node_coordinates(const int dim, const stk::mesh::BulkData &bulk, const stk::mesh::EntityId node)
{
  krino::FieldRef coords(bulk.mesh_meta_data().coordinate_field());
  return get_vector_field(bulk, coords, bulk.get_entity(stk::topology::NODE_RANK,node), dim);
}

static std::vector<stk::math::Vector3d> get_parent_node_coordinates_for_sensitivitiy(const int dim, const krino::LevelSetShapeSensitivity & sens,
  stk::mesh::BulkData & bulk)
{
  std::vector<stk::math::Vector3d> result;
  for(auto && nodeId : sens.parentNodeIds)
  {
    result.push_back(get_node_coordinates(dim, bulk, nodeId) - get_snap_displacements(dim, bulk, nodeId));
  }
  return result;
}

class DecomposeMeshAndComputeSensitivitiesForCircleOrSphere : public ::testing::Test
{
protected:
  void build_circle_or_sphere_conforming_mesh_and_test_sensitivity(const int dim, const double meshSize, const double radius,
    const bool doComputeClosestPointSensitivities, const bool doSnapping, const bool doWriteMesh)
  {
    std::vector<krino::LS_Field> lsFields;
    std::unique_ptr<krino::BoundingBoxMesh> mesh = build_circle_or_sphere_conforming_bounding_box_mesh(dim, {-1.,-1.,-1.}, {1.,1.,1.}, meshSize, {0.0,0.0,0.0}, radius, doComputeClosestPointSensitivities, doSnapping, doWriteMesh, lsFields);

    std::vector<krino::LevelSetShapeSensitivity> sensitivities = get_levelset_shape_sensitivities(mesh->bulk_data(), lsFields[0], doComputeClosestPointSensitivities);
    if (doComputeClosestPointSensitivities)
    {
      for (auto & sens : sensitivities)
        test_closest_point_sensitivity_for_circle_or_sphere_at_origin(dim, radius, sens, mesh->bulk_data());
    }
    else
    {
      for (auto & sens : sensitivities)
        test_sensitivity_for_circle_or_sphere_at_origin(dim, sens, mesh->bulk_data());
    }
  }

  void test_sensitivity_for_circle_or_sphere_at_origin(const int dim, const krino::LevelSetShapeSensitivity & sens, stk::mesh::BulkData & bulk)
  {
    ASSERT_EQ(sens.parentNodeIds.size(), 2u);

    const auto parent_coords = get_parent_node_coordinates_for_sensitivitiy(dim, sens, bulk);
    double radiusDiff = parent_coords[1].length() - parent_coords[0].length();

    for(int d=0; d<dim; d++)
    {
      double sum = 0;
      for (auto & dCoordsdParentLevelSet : sens.dCoordsdParentLevelSets)
        sum += dCoordsdParentLevelSet[d];
      EXPECT_DOUBLE_EQ(sum, (parent_coords[0][d]-parent_coords[1][d])/radiusDiff);
    }
  }

  void test_closest_point_sensitivity_for_circle_or_sphere_at_origin(const int dim, const double radius, const krino::LevelSetShapeSensitivity & sens, stk::mesh::BulkData & bulk)
  {
    const stk::math::Vector3d nodeCoords = get_node_coordinates(dim, bulk, sens.interfaceNodeId);
    const stk::math::Vector3d goldSens = -nodeCoords/radius;
    for(int d=0; d<dim; d++)
    {
      stk::math::Vector3d sumSens = stk::math::Vector3d::ZERO;
      for (auto & dCoordsdParentLevelSet : sens.dCoordsdParentLevelSets)
        sumSens += dCoordsdParentLevelSet;

      krino::expect_near(goldSens, sumSens, 0.05);
    }
  }
};

TEST_F(DecomposeMeshAndComputeSensitivitiesForCircleOrSphere, createDecomposedMeshForSphereNotThroughAnyBackgroundNodes_testSensitivities)
{
  build_circle_or_sphere_conforming_mesh_and_test_sensitivity(3, 0.4, 0.5, false, false, false);
}

TEST_F(DecomposeMeshAndComputeSensitivitiesForCircleOrSphere, createSnappedMeshForSphereNotThroughAnyBackgroundNodes_testSensitivities)
{
  build_circle_or_sphere_conforming_mesh_and_test_sensitivity(3, 0.4, 0.5, false, true, false);
}

TEST_F(DecomposeMeshAndComputeSensitivitiesForCircleOrSphere, createDecomposedMeshForSphereThroughSomeBackgroundNodes_testSensitivities)
{
  build_circle_or_sphere_conforming_mesh_and_test_sensitivity(3, 0.2, 0.6, false, false, false);
}

TEST_F(DecomposeMeshAndComputeSensitivitiesForCircleOrSphere, createSnappedMeshForSphereThroughSomeBackgroundNodes_testSensitivities)
{
  build_circle_or_sphere_conforming_mesh_and_test_sensitivity(3, 0.2, 0.6, false, true, false);
}

TEST_F(DecomposeMeshAndComputeSensitivitiesForCircleOrSphere, createDecomposedMeshForCircleNotThroughAnyBackgroundNodes_testSensitivities)
{
  build_circle_or_sphere_conforming_mesh_and_test_sensitivity(2, 0.4, 0.5, false, false, false);
}

TEST_F(DecomposeMeshAndComputeSensitivitiesForCircleOrSphere, createSnappedMeshForCircleNotThroughAnyBackgroundNodes_testSensitivities)
{
  build_circle_or_sphere_conforming_mesh_and_test_sensitivity(2, 0.4, 0.5, false, true, false);
}

TEST_F(DecomposeMeshAndComputeSensitivitiesForCircleOrSphere, createDecomposedMeshForCircleThroughSomeBackgroundNodes_testSensitivities)
{
  build_circle_or_sphere_conforming_mesh_and_test_sensitivity(2, 0.2, 0.6, false, false, false);
}

TEST_F(DecomposeMeshAndComputeSensitivitiesForCircleOrSphere, createSnappedMeshForCircleThroughSomeBackgroundNodes_testSensitivities)
{
  build_circle_or_sphere_conforming_mesh_and_test_sensitivity(2, 0.2, 0.6, false, true, false);
}

TEST_F(DecomposeMeshAndComputeSensitivitiesForCircleOrSphere, createDecomposedMeshForSphereNotThroughAnyBackgroundNodes_testClosestPointSensitivities)
{
  build_circle_or_sphere_conforming_mesh_and_test_sensitivity(3, 0.2, 0.5, true, false, false);
}

TEST_F(DecomposeMeshAndComputeSensitivitiesForCircleOrSphere, createSnappedMeshForSphereNotThroughAnyBackgroundNodes_testClosestPointSensitivities)
{
  build_circle_or_sphere_conforming_mesh_and_test_sensitivity(3, 0.2, 0.5, true, true, false);
}

TEST_F(DecomposeMeshAndComputeSensitivitiesForCircleOrSphere, createDecomposedMeshForSphereThroughSomeBackgroundNodes_testClosestPointSensitivities)
{
  build_circle_or_sphere_conforming_mesh_and_test_sensitivity(3, 0.2, 0.6, true, false, false);
}

TEST_F(DecomposeMeshAndComputeSensitivitiesForCircleOrSphere, createSnappedMeshForSphereThroughSomeBackgroundNodes_testClosestPointSensitivities)
{
  build_circle_or_sphere_conforming_mesh_and_test_sensitivity(3, 0.2, 0.6, true, true, false);
}

TEST_F(DecomposeMeshAndComputeSensitivitiesForCircleOrSphere, createDecomposedMeshForCircleNotThroughAnyBackgroundNodes_testClosestPointSensitivities)
{
  build_circle_or_sphere_conforming_mesh_and_test_sensitivity(2, 0.2, 0.5, true, false, false);
}

TEST_F(DecomposeMeshAndComputeSensitivitiesForCircleOrSphere, createSnappedMeshForCircleNotThroughAnyBackgroundNodes_testClosestPointSensitivities)
{
  build_circle_or_sphere_conforming_mesh_and_test_sensitivity(2, 0.2, 0.5, true, true, false);
}

TEST_F(DecomposeMeshAndComputeSensitivitiesForCircleOrSphere, createDecomposedMeshForCircleThroughSomeBackgroundNodes_testClosestPointSensitivities)
{
  build_circle_or_sphere_conforming_mesh_and_test_sensitivity(2, 0.2, 0.6, true, false, false);
}

TEST_F(DecomposeMeshAndComputeSensitivitiesForCircleOrSphere, createSnappedMeshForCircleThroughSomeBackgroundNodes_testClosestPointSensitivities)
{
  build_circle_or_sphere_conforming_mesh_and_test_sensitivity(2, 0.2, 0.6, true, true, false);
}

class DecomposeMeshAndCheckConvergenceForSensitivitiesForSphere : public krino::OptimizedOnlyTest
{
protected:
  double build_sphere_conforming_mesh_and_compute_sensitivity_error(const double meshSize, const double radius, const bool doComputeClosestPointSensitivities, const bool doSnapping, const bool doWriteMesh)
  {
    std::vector<krino::LS_Field> lsFields;
    std::unique_ptr<krino::BoundingBoxMesh> mesh = build_circle_or_sphere_conforming_bounding_box_mesh(3, {-1.,-1.,-0.1}, {1.,1.,0.1}, meshSize, {0.0,0.0,0.0}, radius, doComputeClosestPointSensitivities, doSnapping, doWriteMesh, lsFields);

    std::vector<krino::LevelSetShapeSensitivity> sensitivities = get_levelset_shape_sensitivities(mesh->bulk_data(), lsFields[0], doComputeClosestPointSensitivities);
    double errorSumSquare = 0.;
    size_t sensCount = 0;
    double nodeError = 0.;
    for (auto & sens : sensitivities)
    {
      if(mesh->bulk_data().bucket(mesh->bulk_data().get_entity(stk::topology::NODE_RANK,sens.interfaceNodeId)).owned())
      {
        if (doComputeClosestPointSensitivities)
          nodeError = difference_against_analytical_closest_point_sensitivity_for_sphere_at_origin(sens, mesh->bulk_data(), radius);
        else
          nodeError = difference_against_analytical_sensitivity_for_sphere_at_origin(sens, mesh->bulk_data(), radius);

        errorSumSquare += nodeError;
        sensCount += 1;
      }
    }

    double globalErrorSumSquare;
    size_t globalSensCount;
    stk::all_reduce_sum(mesh->bulk_data().parallel(), &errorSumSquare, &globalErrorSumSquare, 1);
    stk::all_reduce_sum(mesh->bulk_data().parallel(), &sensCount, &globalSensCount, 1);
    return std::sqrt(globalErrorSumSquare)/globalSensCount;
  }

  void check_convergence_of_sensitivities_for_sphere_conforming_meshes(const double radius, const bool doComputeClosestPointSensitivities, const bool doSnapping, const bool doWriteMesh)
  {
    double meshSize = 0.2;
    constexpr unsigned nRefine = 3;
    std::array<double, nRefine> errorNorms;

    for(unsigned r=0; r<nRefine; r++)
    {
      errorNorms[r] = build_sphere_conforming_mesh_and_compute_sensitivity_error(meshSize, radius, doComputeClosestPointSensitivities, doSnapping, doWriteMesh);
      if (0 == stk::EnvData::parallel_rank())
        std::cout << "Convergence " << meshSize << " " << errorNorms[r] << std::endl;
      meshSize /= 2.;
    }
    EXPECT_LT(errorNorms[nRefine-1], errorNorms[nRefine-2]/2.);
  }

  double calc_weight_of_exact_sphere_intersection_point(const std::vector<stk::math::Vector3d> & parent_coords, double radius)
  {
    STK_ThrowRequire(parent_coords.size() == 2u);
    auto deltaX = parent_coords[1] - parent_coords[0];
    double a = deltaX.length_squared();
    double b = 2.*stk::math::Dot(deltaX, parent_coords[0]);
    double c = parent_coords[0].length_squared() - radius*radius;

    double det = b*b - 4*a*c;
    STK_ThrowRequire(det >= 0);
    double t0 = (-b - sqrt(det))/2./a;
    double t1 = (-b + sqrt(det))/2./a;
    STK_ThrowRequire((t0 < 0. || t0 > 1.) || (t1 < 0. || t1 > 1.));
    STK_ThrowRequire((t0 >= 0. && t0 <= 1.) || (t1 >= 0. || t1 <= 1.));
    return (t0 >= 0. && t0 <= 1.) ? t0 : t1;
  }

  double difference_against_analytical_sensitivity_for_sphere_at_origin(const krino::LevelSetShapeSensitivity & sens, stk::mesh::BulkData & bulk, double radius)
  {
    static constexpr int DIM = 3;
    static constexpr double tol = 1e-8;
    STK_ThrowRequire(sens.parentNodeIds.size() == 2u);

    const auto parent_coords = get_parent_node_coordinates_for_sensitivitiy(DIM, sens, bulk);
    const auto deltaX = parent_coords[1] - parent_coords[0];
    const double t = calc_weight_of_exact_sphere_intersection_point(parent_coords, radius);

    const auto xIntersection = parent_coords[0] + t * deltaX;
    STK_ThrowRequire(std::fabs(xIntersection.length()-radius) <= tol*radius);

    const double cosTheta = stk::math::Dot(xIntersection.unit_vector(), deltaX.unit_vector());


    stk::math::Vector3d expVal;
    stk::math::Vector3d result;
    for(int d=0; d<DIM; d++)
    {
      double sum = 0;
      for (auto & dCoordsdParentLevelSet : sens.dCoordsdParentLevelSets)
      {
        sum += dCoordsdParentLevelSet[d];
      }
      if(std::fabs(deltaX[d]) < tol)
      {
        expVal[d] = 0.;
        result[d] = sum;
      }
      else
      {
        expVal[d] = cosTheta;
        result[d] = -deltaX[d]/sum/deltaX.length();
      }
    }
    return (result-expVal).length_squared();
  }

  double difference_against_analytical_closest_point_sensitivity_for_sphere_at_origin(const krino::LevelSetShapeSensitivity & sens, stk::mesh::BulkData & bulk, double radius)
  {
    static constexpr int DIM = 3;
    const stk::math::Vector3d nodeCoords = get_node_coordinates(DIM, bulk, sens.interfaceNodeId);

    stk::math::Vector3d goldVal = -nodeCoords / nodeCoords.length();
    stk::math::Vector3d result = stk::math::Vector3d::ZERO;
    for (auto & dCoordsdParentLevelSet : sens.dCoordsdParentLevelSets)
      result += dCoordsdParentLevelSet;
    return (result-goldVal).length_squared();

  }
};

TEST_F(DecomposeMeshAndCheckConvergenceForSensitivitiesForSphere, createDecomposedMeshForSphereNotThroughBackgroundNodes_testAnalyticalSensitivities)
{
  check_convergence_of_sensitivities_for_sphere_conforming_meshes(0.869, false, false, false);
}

TEST_F(DecomposeMeshAndCheckConvergenceForSensitivitiesForSphere, createSnappedMeshForSphereNotThroughBackgroundNodes_testAnalyticalSensitivities)
{
  check_convergence_of_sensitivities_for_sphere_conforming_meshes(0.869, false, true, false);
}

TEST_F(DecomposeMeshAndCheckConvergenceForSensitivitiesForSphere, createDecomposedMeshForSphereNotThroughBackgroundNodes_testAnalyticalClosestPointSensitivities)
{
  check_convergence_of_sensitivities_for_sphere_conforming_meshes(0.869, true, false, false);
}

TEST_F(DecomposeMeshAndCheckConvergenceForSensitivitiesForSphere, createSnappedMeshForSphereNotThroughBackgroundNodes_testAnalyticalClosestPointSensitivities)
{
  check_convergence_of_sensitivities_for_sphere_conforming_meshes(0.869, true, true, false);
}
