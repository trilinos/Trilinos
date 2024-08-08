#include <gtest/gtest.h>

#include <Akri_BoundingBoxMesh.hpp>
#include <Akri_CDFEM_Support.hpp>
#include <Akri_Composite_Surface.hpp>
#include <Akri_AnalyticSurf.hpp>
#include <Akri_CDMesh.hpp>
#include <Akri_CreateInterfaceGeometry.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_LevelSetPolicy.hpp>
#include <Akri_LevelSetShapeSensitivities.hpp>
#include <Akri_MeshFromFile.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_NodalSurfaceDistance.hpp>
#include <Akri_OutputUtils.hpp>
#include <Akri_Phase_Support.hpp>
#include <Akri_SideAttachedElements.hpp>
#include <Akri_Unit_LogRedirecter.hpp>
#include <stk_io/IossBridge.hpp>
#include <stk_topology/topology.hpp>
#include <stk_mesh/base/MetaData.hpp>

static void setup_fields_for_conforming_decomposition(const stk::mesh::MetaData & meta)
{
  krino::CDFEM_Support & cdfemSupport = krino::CDFEM_Support::get(meta);
  const krino::FieldRef coordsField = meta.coordinate_field();

  cdfemSupport.set_coords_field(coordsField);
  cdfemSupport.add_edge_interpolation_field(coordsField);
  cdfemSupport.register_parent_node_ids_field();
}

static void decompose_mesh_to_conform_to_levelsets(stk::mesh::BulkData & mesh, const std::vector<krino::LS_Field> & lsFields)
{
  stk::mesh::MetaData & meta = mesh.mesh_meta_data();
  krino::AuxMetaData & auxMeta = krino::AuxMetaData::get(meta);
  krino::CDFEM_Support & cdfemSupport = krino::CDFEM_Support::get(meta);
  krino::Phase_Support & phaseSupport = krino::Phase_Support::get(meta);

  std::unique_ptr<krino::InterfaceGeometry> interfaceGeometry = krino::create_levelset_geometry(meta.spatial_dimension(), auxMeta.active_part(), cdfemSupport, phaseSupport, lsFields);
  auxMeta.clear_force_64bit_flag();
  krino::CDMesh::decompose_mesh(mesh, *interfaceGeometry);
}

static void generate_bounding_box_mesh(krino::BoundingBoxMesh & bboxMesh, const stk::math::Vector3d & minCorner, const stk::math::Vector3d & maxCorner, const double meshSize)
{
  bboxMesh.set_domain(krino::BoundingBoxMesh::BoundingBoxType(minCorner, maxCorner), meshSize);
  bboxMesh.set_mesh_structure_type(krino::FLAT_WALLED_BCC_BOUNDING_BOX_MESH);
  bboxMesh.populate_mesh();
  stk::mesh::BulkData & mesh = bboxMesh.bulk_data();
  krino::activate_all_entities(mesh, krino::AuxMetaData::get(mesh.mesh_meta_data()).active_part());
}

static void generate_and_write_bounding_box_mesh(const stk::topology elemTopology, const stk::math::Vector3d & minCorner, const stk::math::Vector3d & maxCorner, const double meshSize, const std::string & filename)
{
  krino::BoundingBoxMesh bboxMesh(elemTopology, MPI_COMM_WORLD);
  generate_bounding_box_mesh(bboxMesh, minCorner, maxCorner, meshSize);
  krino::output_mesh_with_fields(bboxMesh.bulk_data(), krino::AuxMetaData::get(bboxMesh.meta_data()).active_part(), filename, 1, 0.0);
}

std::map<std::string,size_t> get_part_counts(const stk::mesh::BulkData & mesh)
{
  std::map<std::string,size_t> partNameToCount;
  for (auto & partPtr : mesh.mesh_meta_data().get_mesh_parts())
    partNameToCount[partPtr->name()] = krino::get_global_num_entities(mesh, *partPtr);
  return partNameToCount;
}

void expect_part_counts(const stk::mesh::BulkData & mesh, const std::map<std::string,size_t> & goldPartCounts)
{
  for (auto & partPtr : mesh.mesh_meta_data().get_mesh_parts())
  {
    const auto iter = goldPartCounts.find(partPtr->name());
    EXPECT_FALSE(iter == goldPartCounts.end());
    if (iter != goldPartCounts.end())
    {
      EXPECT_EQ(iter->second, krino::get_global_num_entities(mesh, *partPtr)) << "Mismatch for part count for part " << partPtr->name();
    }
  }
}

void test_sensitivity_for_plane(const krino::LevelSetShapeSensitivity & sens, const int iPlaneCoord)
{
  double sum = 0;
  for (auto & dCoordsdParentLevelSet : sens.dCoordsdParentLevelSets)
    sum += dCoordsdParentLevelSet[iPlaneCoord];

  EXPECT_NEAR(-1., sum, 1.e-4);
}

stk::mesh::Selector build_output_selector(const stk::mesh::MetaData & meta, const stk::mesh::Part & activePart)
{
  stk::mesh::PartVector outputParts;
  for (auto * part : meta.get_parts())
    if (stk::io::is_part_io_part(*part) && part->name().find("_void") == std::string::npos)
      outputParts.push_back(part);
  return activePart & stk::mesh::selectUnion(outputParts);
}

void output_mesh(const stk::mesh::BulkData & mesh, const std::string & fileName, int step, double time)
{
  stk::mesh::Selector outputSelector = build_output_selector(mesh.mesh_meta_data(), krino::AuxMetaData::get(mesh.mesh_meta_data()).active_part());
  krino::output_composed_mesh_with_fields(mesh, outputSelector, fileName, step, time);
}

TEST(DecomposeMeshAndComputeSensitivities, createDecomposedMeshForPlaneNotThroughAnyBackgroundNodes_testSensitivities)
{
  krino::BoundingBoxMesh bboxMesh(stk::topology::TETRAHEDRON_4, MPI_COMM_WORLD);

  const std::vector<krino::LS_Field> lsFields = krino::LSPerInterfacePolicy::setup_levelsets_on_all_blocks_with_void_phase_for_any_negative_levelset(bboxMesh.meta_data(), 1);
  setup_fields_for_conforming_decomposition(bboxMesh.meta_data());

  generate_bounding_box_mesh(bboxMesh, {0.,0.,0.}, {1.,1.,1.}, 0.3333);

  compute_nodal_distance_from_plane(bboxMesh.bulk_data(), bboxMesh.meta_data().coordinate_field(), lsFields[0].isovar, {1.,0.,0.}, -0.2);

  decompose_mesh_to_conform_to_levelsets(bboxMesh.bulk_data(), lsFields);

  std::vector<krino::LevelSetShapeSensitivity> sensitivities = get_levelset_shape_sensitivities(bboxMesh.bulk_data(), lsFields);
  for (auto & sens : sensitivities)
    test_sensitivity_for_plane(sens, 0);

  const bool doWriteMesh = false;
  if (doWriteMesh)
    output_mesh(bboxMesh.bulk_data(), "output.e", 1, 0.0);
}

TEST(DecomposeMeshAndComputeSensitivities, createDecomposedMeshForPlaneThrowSomeBackgroundNodes_testSensitivities)
{
  krino::BoundingBoxMesh bboxMesh(stk::topology::TETRAHEDRON_4, MPI_COMM_WORLD);

  const std::vector<krino::LS_Field> lsFields = krino::LSPerInterfacePolicy::setup_levelsets_on_all_blocks_with_void_phase_for_any_negative_levelset(bboxMesh.meta_data(), 1);
  setup_fields_for_conforming_decomposition(bboxMesh.meta_data());

  generate_bounding_box_mesh(bboxMesh, {0.,0.,0.}, {1.,1.,1.}, 1.0);

  compute_nodal_distance_from_plane(bboxMesh.bulk_data(), bboxMesh.meta_data().coordinate_field(), lsFields[0].isovar, {1.,0.,0.}, -0.5);

  decompose_mesh_to_conform_to_levelsets(bboxMesh.bulk_data(), lsFields);

  std::vector<krino::LevelSetShapeSensitivity> sensitivities = get_levelset_shape_sensitivities(bboxMesh.bulk_data(), lsFields);
  for (auto & sens : sensitivities)
    test_sensitivity_for_plane(sens, 0);

  const bool doWriteMesh = false;
  if (doWriteMesh)
    output_mesh(bboxMesh.bulk_data(), "output.e", 1, 0.0);
}

stk::mesh::PartVector get_nonvoid_parts_of_rank(const stk::mesh::MetaData & meta, const stk::mesh::EntityRank entityRank)
{
  stk::mesh::PartVector parts;
  for (auto * part : meta.get_parts())
    if (stk::io::is_part_io_part(*part) && part->primary_entity_rank() == entityRank && part->name().find("_void") == std::string::npos)
      parts.push_back(part);
  return parts;
}

stk::mesh::Selector get_nonvoid_selector(const stk::mesh::MetaData & meta, const stk::mesh::EntityRank entityRank)
{
  return stk::mesh::selectUnion(get_nonvoid_parts_of_rank(meta, entityRank));
}

std::map<int,int> get_nonvoid_to_void_part_ordinal_map(const stk::mesh::MetaData & meta)
{
  const krino::Phase_Support & phaseSupport = krino::Phase_Support::get(meta);

  std::map<int,int> partOrdinalMapping;
  for (auto * part : meta.get_parts())
  {
    if ((part->primary_entity_rank() == stk::topology::ELEMENT_RANK || part->primary_entity_rank() == meta.side_rank()) &&
        stk::io::is_part_io_part(*part))
    {
      if (phaseSupport.is_interface(part))
      {
        partOrdinalMapping[part->mesh_meta_data_ordinal()] = -1; // negative value used to indicate that we will remove interface parts
      }
      else if (part->name().find("_void") == std::string::npos)
      {
        stk::mesh::Part * voidPart = meta.get_part(part->name() + "_void");
        if (voidPart)
          partOrdinalMapping[part->mesh_meta_data_ordinal()] = voidPart->mesh_meta_data_ordinal();
      }
    }
  }
  return partOrdinalMapping;
}

void move_disconnected_elements_to_void(stk::mesh::BulkData & mesh)
{
  const stk::mesh::Selector sideSelector = get_nonvoid_selector(mesh.mesh_meta_data(), mesh.mesh_meta_data().side_rank());
  const stk::mesh::Selector elementSelector = get_nonvoid_selector(mesh.mesh_meta_data(), stk::topology::ELEMENT_RANK);
  const std::vector<stk::mesh::Entity> unattachedElems = krino::get_selected_owned_side_unattached_elements(mesh, elementSelector, sideSelector);

  const std::map<int,int> partOrdinalMapping = get_nonvoid_to_void_part_ordinal_map(mesh.mesh_meta_data());
  krino::batch_convert_elements_and_their_sides(mesh, partOrdinalMapping, unattachedElems);
}

void move_elements_not_in_largest_group_to_void(stk::mesh::BulkData & mesh)
{
  const stk::mesh::Selector elementSelector = get_nonvoid_selector(mesh.mesh_meta_data(), stk::topology::ELEMENT_RANK);
  const std::vector<stk::mesh::Entity> elemsNotInLargestGroup = krino::find_owned_elements_that_are_not_in_the_largest_group_of_selected_side_attached_elements(mesh, elementSelector);

  const std::map<int,int> partOrdinalMapping = get_nonvoid_to_void_part_ordinal_map(mesh.mesh_meta_data());
  krino::batch_convert_elements_and_their_sides(mesh, partOrdinalMapping, elemsNotInLargestGroup);
}

TEST(DecomposeMeshAndComputeSensitivities, readMeshInitializeDecomposeResetInitializeDecompose)
{
  const std::string initialMeshName = "mesh.g";
  generate_and_write_bounding_box_mesh(stk::topology::TETRAHEDRON_4, {0.,0.,0.}, {1.,1.,1.}, 0.3333, initialMeshName);

  krino::MeshFromFile meshFromFile(initialMeshName, MPI_COMM_WORLD);

  const std::vector<krino::LS_Field> lsFields = krino::LSPerInterfacePolicy::setup_levelsets_on_all_blocks_with_void_phase_for_any_negative_levelset(meshFromFile.meta_data(), 1);
  setup_fields_for_conforming_decomposition(meshFromFile.meta_data());

  meshFromFile.populate_mesh();
  krino::activate_all_entities(meshFromFile.bulk_data(), krino::AuxMetaData::get(meshFromFile.meta_data()).active_part());

  const std::map<std::string,size_t> origPartCounts = get_part_counts(meshFromFile.bulk_data());

  const std::vector<std::pair<stk::math::Vector3d,double>> spheres
    {
      { {0.2,0.2,0.2}, 0.25 },
      { {0.2,0.2,0.7}, 0.25 },
      { {0.2,0.7,0.2}, 0.25 },
      { {0.2,0.7,0.7}, 0.25 },
      { {0.7,0.2,0.2}, 0.25 },
      { {0.7,0.2,0.7}, 0.25 },
      { {0.7,0.7,0.2}, 0.25 },
      { {0.7,0.7,0.7}, 0.25 },
    };

  compute_nodal_distance_from_spheres(meshFromFile.bulk_data(), meshFromFile.meta_data().coordinate_field(), lsFields[0].isovar, spheres);

  decompose_mesh_to_conform_to_levelsets(meshFromFile.bulk_data(), lsFields);

  const bool doWriteMesh = false;
  if (doWriteMesh)
    output_mesh(meshFromFile.bulk_data(), "output1.e", 1, 0.0);

  const std::map<std::string,size_t> decompPartCounts = get_part_counts(meshFromFile.bulk_data());

  krino::CDMesh::reset_mesh_to_original_undecomposed_state(meshFromFile.bulk_data());

  expect_part_counts(meshFromFile.bulk_data(), origPartCounts);

  if (doWriteMesh)
    output_mesh(meshFromFile.bulk_data(), "reset.e", 1, 1.0);

  compute_nodal_distance_from_spheres(meshFromFile.bulk_data(), meshFromFile.meta_data().coordinate_field(), lsFields[0].isovar, spheres );

  decompose_mesh_to_conform_to_levelsets(meshFromFile.bulk_data(), lsFields);

  if (doWriteMesh)
    output_mesh(meshFromFile.bulk_data(), "output2.e", 1, 1.0);

  expect_part_counts(meshFromFile.bulk_data(), decompPartCounts);
}

void test_moving_islands_to_void(const std::function<void(stk::mesh::BulkData&)> & island_removal_method)
{
  const std::string initialMeshName = "mesh.g";
  generate_and_write_bounding_box_mesh(stk::topology::TETRAHEDRON_4, {0.,0.,0.}, {1.,1.,1.}, 0.3333, initialMeshName);

  krino::MeshFromFile meshFromFile(initialMeshName, MPI_COMM_WORLD);

  const std::vector<krino::LS_Field> lsFields = krino::LSPerInterfacePolicy::setup_levelsets_on_all_blocks_with_void_phase_for_any_negative_levelset(meshFromFile.meta_data(), 1);
  setup_fields_for_conforming_decomposition(meshFromFile.meta_data());

  meshFromFile.populate_mesh();
  krino::activate_all_entities(meshFromFile.bulk_data(), krino::AuxMetaData::get(meshFromFile.meta_data()).active_part());

  const std::vector<std::pair<stk::math::Vector3d,double>> oneSphere
    {
      { {0.0,0.0,0.0}, 0.7 },
    };

  compute_nodal_distance_from_spheres(meshFromFile.bulk_data(), meshFromFile.meta_data().coordinate_field(), lsFields[0].isovar, oneSphere, -1);

  decompose_mesh_to_conform_to_levelsets(meshFromFile.bulk_data(), lsFields);

  const bool doWriteMesh = false;
  if (doWriteMesh)
    output_mesh(meshFromFile.bulk_data(), "output1.e", 1, 0.0);

  const stk::mesh::Part & block_1 = *meshFromFile.meta_data().get_part("block_1");
  const size_t numElementsInBlockFromOneSphere = krino::get_global_num_entities(meshFromFile.bulk_data(), block_1);

  krino::CDMesh::reset_mesh_to_original_undecomposed_state(meshFromFile.bulk_data());

  const std::vector<std::pair<stk::math::Vector3d,double>> oneConnectedAndOneDisconnectedSphere
    {
      { {0.0,0.0,0.0}, 0.7 },
      { {0.7,0.7,0.7}, 0.25 }
    };

  compute_nodal_distance_from_spheres(meshFromFile.bulk_data(), meshFromFile.meta_data().coordinate_field(), lsFields[0].isovar, oneConnectedAndOneDisconnectedSphere, -1);

  decompose_mesh_to_conform_to_levelsets(meshFromFile.bulk_data(), lsFields);

  island_removal_method(meshFromFile.bulk_data());

  if (doWriteMesh)
    output_mesh(meshFromFile.bulk_data(), "output2.e", 1, 1.0);

  const size_t numElementsInBlockAfterDeletingSecondDisconnectedSphere = krino::get_global_num_entities(meshFromFile.bulk_data(), block_1);

  EXPECT_EQ(numElementsInBlockFromOneSphere, numElementsInBlockAfterDeletingSecondDisconnectedSphere);
}

TEST(DecomposeMeshAndComputeSensitivities, readMeshInitializeDecomposeAndRemoveVolNotConnectedToSide)
{
  auto island_removal_method = [](stk::mesh::BulkData& mesh){ move_disconnected_elements_to_void(mesh); };
  test_moving_islands_to_void(island_removal_method);
}

TEST(DecomposeMeshAndComputeSensitivities, readMeshInitializeDecomposeAndRemoveVolNotConnectedToLargestVol)
{
  auto island_removal_method = [](stk::mesh::BulkData& mesh){ move_elements_not_in_largest_group_to_void(mesh); };
  test_moving_islands_to_void(island_removal_method);
}
