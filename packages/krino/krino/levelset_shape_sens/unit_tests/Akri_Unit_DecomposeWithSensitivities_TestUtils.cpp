#include <Akri_Unit_DecomposeWithSensitivities_TestUtils.hpp>

#include <Akri_AuxMetaData.hpp>
#include <Akri_BoundingBoxMesh.hpp>
#include <Akri_LevelSetPolicy.hpp>
#include <Akri_LevelSetShapeSensitivities.hpp>
#include <Akri_NodalSurfaceDistance.hpp>
#include <Akri_OutputUtils.hpp>
#include <Akri_UnitMeshUtils.hpp>
#include <stk_io/IossBridge.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_topology/topology.hpp>

stk::mesh::Selector build_output_selector(const stk::mesh::MetaData & meta, const stk::mesh::Part & activePart)
{
  stk::mesh::PartVector outputParts;
  for (auto * part : meta.get_parts())
    if (stk::io::is_part_io_part(*part) && part->name().find("_void") == std::string::npos)
      outputParts.push_back(part);
  return activePart & stk::mesh::selectUnion(outputParts);
}

void output_nonvoid_mesh(stk::mesh::BulkData & mesh, const std::string & fileName, int step, double time)
{
  stk::mesh::Selector outputSelector = build_output_selector(mesh.mesh_meta_data(), krino::AuxMetaData::get(mesh.mesh_meta_data()).active_part());
  krino::fix_ownership_and_output_composed_mesh_with_fields(mesh, outputSelector, fileName, step, time);
}

void output_full_mesh(stk::mesh::BulkData & mesh, const std::string & fileName, int step, double time)
{
  krino::fix_ownership_and_output_composed_mesh_with_fields(mesh, krino::AuxMetaData::get(mesh.mesh_meta_data()).active_part(), fileName, step, time);
}

std::unique_ptr<krino::BoundingBoxMesh> build_circle_or_sphere_conforming_bounding_box_mesh(
    const int dim,
    const stk::math::Vector3d & minCorner,
    const stk::math::Vector3d & maxCorner,
    const double meshSize,
    const stk::math::Vector3d & sphereCenter,
    const double sphereRadius,
    const bool doComputeClosestPointSensitivities,
    const bool doSnapping,
    const bool doWriteMesh,
    std::vector<krino::LS_Field> & lsFields)
{
  const stk::topology elemTopo = (dim == 2) ? stk::topology::TRIANGLE_3_2D : stk::topology::TETRAHEDRON_4;
  std::unique_ptr<krino::BoundingBoxMesh> bboxMesh = std::make_unique<krino::BoundingBoxMesh>(elemTopo, MPI_COMM_WORLD);

  lsFields = krino::LSPerInterfacePolicy::setup_levelsets_on_all_blocks_with_void_phase_for_any_negative_levelset(bboxMesh->meta_data(), 1);
  setup_fields_for_conforming_decomposition(bboxMesh->meta_data(), lsFields, doSnapping, doComputeClosestPointSensitivities);

  krino::populate_bounding_box_mesh_and_activate(*bboxMesh, minCorner, maxCorner, meshSize);

  const std::vector<std::pair<stk::math::Vector3d,double>> spheres{ { sphereCenter, sphereRadius } };

  compute_nodal_distance_from_spheres(bboxMesh->bulk_data(), bboxMesh->meta_data().coordinate_field(), lsFields[0].isovar, spheres);

  decompose_mesh_to_conform_to_levelsets(bboxMesh->bulk_data(), lsFields, doSnapping, doComputeClosestPointSensitivities);

  if (doWriteMesh)
  {
    const std::string filename = std::string(dim == 2 ? "circle" : "sphere") + std::string(doSnapping ? "_snapped.e" : "_cut.e");
    output_full_mesh(bboxMesh->bulk_data(), filename, 1, 0.0);
  }

  return bboxMesh;
}
