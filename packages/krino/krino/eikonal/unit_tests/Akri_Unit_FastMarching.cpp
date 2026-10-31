#include <gtest/gtest.h>

#include <Akri_AuxMetaData.hpp>
#include <Akri_BoundingBoxMesh.hpp>
#include <Akri_Fast_Marching.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_NodalSurfaceDistance.hpp>
#include <Akri_OutputUtils.hpp>
#include <Akri_PostProcess.hpp>
#include <stk_math/StkPlane.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_topology/topology.hpp>
#include <stk_util/environment/EnvData.hpp>

namespace krino {

static FieldRef register_distance_field(stk::mesh::MetaData & meta)
{
  FieldRef distanceField = AuxMetaData::get(meta).register_field("distance", FieldType::REAL, stk::topology::NODE_RANK, 2, 1, meta.universal_part());
  return distanceField;
}

static void generate_bounding_box_mesh(krino::BoundingBoxMesh & bboxMesh, const stk::math::Vector3d & minCorner, const stk::math::Vector3d & maxCorner, const double meshSize)
{
  bboxMesh.set_domain(krino::BoundingBoxMesh::BoundingBoxType(minCorner, maxCorner), meshSize);
  bboxMesh.populate_mesh();
  stk::mesh::BulkData & mesh = bboxMesh.bulk_data();
  krino::activate_all_entities(mesh, krino::AuxMetaData::get(mesh.mesh_meta_data()).active_part());
  populate_stk_local_ids(bboxMesh.bulk_data());
}

static double redistance_and_compute_error(const stk::mesh::BulkData & mesh, const FieldRef coordsField, const FieldRef distanceField, const std::function<double(const stk::math::Vector3d &)> & analytic_fn)
{
  std::function<double(ParallelErrorMessage& err, stk::mesh::Entity)> get_interface_speed;

  Fast_Marching fm(mesh,
        mesh.mesh_meta_data().universal_part(),
        coordsField,
        distanceField,
        get_interface_speed,
        sierra::Diag::sierraTimer());
  fm.redistance();

  const bool doOutput = false;
  if (doOutput)
    output_composed_mesh_with_fields(mesh, mesh.mesh_meta_data().universal_part(), "fastMarching.e", 1, 0.);

  const double err = compute_relative_nodal_RMS_error(mesh, coordsField, distanceField, analytic_fn);
  return err;
}

static double build_mesh_initialize_sphere_redistance_and_compute_error(const stk::topology elemTopology, const double meshSize, stk::diag::Timer & testTimer)
{
  stk::diag::TimeBlock timer__(testTimer);
  krino::BoundingBoxMesh bboxMesh(elemTopology, MPI_COMM_WORLD);
  FieldRef distanceField = register_distance_field(bboxMesh.meta_data());
  generate_bounding_box_mesh(bboxMesh, stk::math::Vector3d{0,0,0}, stk::math::Vector3d{1,1,1}, meshSize);

  const stk::math::Vector3d center{-0.05,-0.05,0};
  const double radius = 0.7;
  const std::vector<std::pair<stk::math::Vector3d,double>> sphere { { center, radius } };
  auto analytic_fn = [&center, radius](const stk::math::Vector3d &x) { return (x-center).length() - radius; };

  const FieldRef coordsField = bboxMesh.meta_data().coordinate_field();
  compute_nodal_distance_from_spheres(bboxMesh.bulk_data(), coordsField, distanceField, sphere);

  return redistance_and_compute_error(bboxMesh.bulk_data(), coordsField, distanceField, analytic_fn);
}

static void test_fast_marching_error_for_circle_or_sphere(const stk::topology elemTopology, const double coarseMeshSize, const std::vector<double> & goldErrors, stk::diag::Timer & testTimer)
{
  double meshSize = coarseMeshSize;
  for(size_t i=0; i<goldErrors.size(); ++i)
  {
    const double error = build_mesh_initialize_sphere_redistance_and_compute_error(elemTopology, meshSize, testTimer);
    if (0 == stk::EnvData::parallel_rank())
      std::cout << "Error for " << elemTopology.name() << " at size " << meshSize << " = " << error << ", CPU time = " << testTimer.getMetric<stk::diag::CPUTime>().getLap() << std::endl;
    const double tol = 0.05*(error + goldErrors[i]);
    EXPECT_NEAR(goldErrors[i], error, tol);
    meshSize /= 2.;
  }
}

static double build_mesh_initialize_plane_redistance_and_compute_error(const stk::topology elemTopology, const double meshSize, stk::diag::Timer & testTimer)
{
  stk::diag::TimeBlock timer__(testTimer);
  krino::BoundingBoxMesh bboxMesh(elemTopology, MPI_COMM_WORLD);
  FieldRef distanceField = register_distance_field(bboxMesh.meta_data());
  generate_bounding_box_mesh(bboxMesh, stk::math::Vector3d{0,0,0}, stk::math::Vector3d{1,1,1}, meshSize);

  const stk::math::Vector3d normal = stk::math::Vector3d{-std::sqrt(2.)/2., std::sqrt(2.)/2., 0.};
  const double offset = 0.;
  auto analytic_fn = [&normal, offset](const stk::math::Vector3d &x) { return Dot(x,normal) + offset; };

  const FieldRef coordsField = bboxMesh.meta_data().coordinate_field();
  compute_nodal_distance_from_plane(bboxMesh.bulk_data(), coordsField, distanceField, normal, offset);

  return redistance_and_compute_error(bboxMesh.bulk_data(), coordsField, distanceField, analytic_fn);
}

static void test_fast_marching_error_for_plane(const stk::topology elemTopology, const double meshSize, const double tol, stk::diag::Timer & testTimer)
{
  const double error = build_mesh_initialize_plane_redistance_and_compute_error(elemTopology, meshSize, testTimer);
  if (0 == stk::EnvData::parallel_rank())
    std::cout << "Error for " << elemTopology.name() << " = " << error << ", CPU time = " << testTimer.getMetric<stk::diag::CPUTime>().getLap() << std::endl;
  EXPECT_NEAR(0., error, tol);
}

TEST(FastMarching, redistanceForCircleOrSphere_firstOrderAccurateDistance)
{
  const double coarseMeshSize = 0.5;
  const std::vector<double> goldHexErrors{0.076, 0.035, 0.023, 0.013}; // First order convergence at finer resolutions than shown here
  const std::vector<double> goldTetErrors{0.068, 0.027, 0.015, 0.0089}; // First order convergence at finer resolutions than shown here
  const std::vector<double> goldTriErrors{0.060, 0.029, 0.017, 0.0085, 0.0044, 0.0023}; // First order convergence as shown here
  const std::vector<double> goldQuadErrors{0.048, 0.048, 0.028, 0.015, 0.0083, 0.0043}; // First order convergence as shown here

  stk::diag::Timer testTimer("FastMarching Sphere/Circle Test", sierra::Diag::sierraTimer());

  test_fast_marching_error_for_circle_or_sphere(stk::topology::TRIANGLE_3_2D, coarseMeshSize, goldTriErrors, testTimer);
  test_fast_marching_error_for_circle_or_sphere(stk::topology::QUADRILATERAL_4_2D, coarseMeshSize, goldQuadErrors, testTimer);
  test_fast_marching_error_for_circle_or_sphere(stk::topology::TETRAHEDRON_4, coarseMeshSize, goldTetErrors, testTimer);
  test_fast_marching_error_for_circle_or_sphere(stk::topology::HEXAHEDRON_8, coarseMeshSize, goldHexErrors, testTimer);
}

TEST(FastMarching, redistanceForPlaneWhereClosestPointEveryWhereIsInDomain_machinePrecisionDistance)
{
  const double meshSize = 0.2;
  const double tol = 1.e-12;

  stk::diag::Timer testTimer("FastMarching Plane Test", sierra::Diag::sierraTimer());

  test_fast_marching_error_for_plane(stk::topology::TRIANGLE_3_2D, meshSize, tol, testTimer);
  test_fast_marching_error_for_plane(stk::topology::QUADRILATERAL_4_2D, meshSize, tol, testTimer);
  test_fast_marching_error_for_plane(stk::topology::TETRAHEDRON_4, meshSize, tol, testTimer);
  test_fast_marching_error_for_plane(stk::topology::HEXAHEDRON_8, meshSize, tol, testTimer);
}


}
