#include <Akri_AllReduce.hpp>
#include <gtest/gtest.h>

#include <string>

#include <Akri_Unit_DecomposeWithSensitivities_TestUtils.hpp>

#include <Akri_AuxMetaData.hpp>
#include <Akri_BoundingBoxMesh.hpp>
#include <Akri_CDFEM_Support.hpp>
#include <Akri_CDMesh.hpp>
#include <Akri_LevelSetPolicy.hpp>
#include <Akri_LevelSetShapeSensitivities.hpp>
#include <Akri_MeshFromFile.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_MeshSpecs.hpp>
#include <Akri_NodalSurfaceDistance.hpp>
#include <Akri_Phase_Support.hpp>
#include <Akri_ScopedTestFile.hpp>
#include <Akri_SegmentWithSensitivities.hpp>
#include <Akri_SideAttachedElements.hpp>
#include <Akri_StkMeshFixture.hpp>
#include <Akri_TriangleWithSensitivities.hpp>
#include <Akri_UnitMeshUtils.hpp>
#include <Akri_UnitTestUtils.hpp>
#include <stk_io/IossBridge.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_topology/topology.hpp>
#include <stk_util/environment/EnvData.hpp>

static void create_extra_phase_per_block(stk::mesh::MetaData & meta, const std::string & phaseName)
{
  krino::AuxMetaData & auxMeta = krino::AuxMetaData::get(meta);
  krino::Phase_Support & phaseSupport = krino::Phase_Support::get(meta);
  for(auto * part : meta.get_parts())
  {
    if (part->primary_entity_rank() == stk::topology::ELEMENT_RANK &&
        stk::io::is_part_io_part(*part))
    {
      const stk::mesh::Part & originalPart = phaseSupport.find_original_part(*part);
      if (originalPart.name() == part->name())
      {
        const std::string phasePartName = part->name() + "_" + phaseName;
        if (!auxMeta.has_part(phasePartName))
          auxMeta.declare_io_part_with_topology(phasePartName, part->topology());
      }
    }
  }
}

static std::map<std::string,size_t> get_part_counts(const stk::mesh::BulkData & mesh)
{
  std::map<std::string,size_t> partNameToCount;
  for (auto & partPtr : mesh.mesh_meta_data().get_mesh_parts())
    partNameToCount[partPtr->name()] = krino::get_global_num_entities(mesh, *partPtr);
  return partNameToCount;
}

static void expect_part_counts(const stk::mesh::BulkData & mesh, const std::map<std::string,size_t> & goldPartCounts)
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

TEST(CopiedLevelSetValuesForSnapping, checkThatOriginalLevelsetIsNotModifiedByDecomposition)
{
  const stk::topology elemTopo = stk::topology::TETRAHEDRON_4;
  std::unique_ptr<krino::BoundingBoxMesh> bboxMesh = std::make_unique<krino::BoundingBoxMesh>(elemTopo, MPI_COMM_WORLD);

  auto lsFields = krino::LSPerInterfacePolicy::setup_levelsets_on_all_blocks_with_void_phase_for_any_negative_levelset(bboxMesh->meta_data(), 1);
  const bool doComputeClosestPointSensitivities = false;
  setup_fields_for_conforming_decomposition(bboxMesh->meta_data(), lsFields, true, doComputeClosestPointSensitivities);

  krino::populate_bounding_box_mesh_and_activate(*bboxMesh, {-1.,-1.,-0.1}, {1.,1.,0.1}, 0.1);

  const std::vector<std::pair<stk::math::Vector3d,double>> spheres{ { {0, 0, 0}, 0.869 } };

  compute_nodal_distance_from_spheres(bboxMesh->bulk_data(), bboxMesh->meta_data().coordinate_field(), lsFields[0].isovar, spheres);

  std::map<stk::mesh::EntityId, double> fieldValuesByEntityId;

  std::vector<stk::mesh::Entity> nodes;
  stk::mesh::get_entities(bboxMesh->bulk_data(), stk::topology::NODE_RANK, nodes);
  for (stk::mesh::Entity node : nodes)
  {
    stk::mesh::EntityId entityId = bboxMesh->bulk_data().identifier(node);
    fieldValuesByEntityId[entityId] = lsFields[0].isovar.field().data<double>().entity_values(node)();
  }

  decompose_mesh_to_conform_to_levelsets(bboxMesh->bulk_data(), lsFields, true, doComputeClosestPointSensitivities);

  // Confirm that snapping/decomposition has not modified the field
  for(auto && entry : fieldValuesByEntityId)
  {
    auto n = bboxMesh->bulk_data().get_entity(stk::topology::NODE_RANK, entry.first);
    EXPECT_DOUBLE_EQ(entry.second,
      lsFields[0].isovar.field().data<double>().entity_values(n)());
  }
}

static stk::mesh::PartVector get_nonvoid_parts_of_rank(const stk::mesh::MetaData & meta, const stk::mesh::EntityRank entityRank)
{
  stk::mesh::PartVector parts;
  for (auto * part : meta.get_parts())
    if (stk::io::is_part_io_part(*part) && part->primary_entity_rank() == entityRank && part->name().find("_void") == std::string::npos)
      parts.push_back(part);
  return parts;
}

static stk::mesh::Selector get_nonvoid_selector(const stk::mesh::MetaData & meta, const stk::mesh::EntityRank entityRank)
{
  return stk::mesh::selectUnion(get_nonvoid_parts_of_rank(meta, entityRank));
}

static std::map<int,int> get_nonvoid_to_phase_part_ordinal_map(const stk::mesh::MetaData & meta, const std::string & phaseName)
{
  const krino::Phase_Support & phaseSupport = krino::Phase_Support::get(meta);
  const std::string phaseSuffix = "_" + phaseName;

  std::map<int,int> partOrdinalMapping;
  for (auto * part : meta.get_mesh_parts())
  {
    if ((part->primary_entity_rank() == stk::topology::ELEMENT_RANK || part->primary_entity_rank() == meta.side_rank()))
    {
      if (phaseSupport.is_interface(*part))
      {
        partOrdinalMapping[part->mesh_meta_data_ordinal()] = -1; // negative value used to indicate that we will remove interface parts
      }
      else if (part->name().find(phaseSuffix) == std::string::npos)
      {
        stk::mesh::Part * phasePart = meta.get_part(part->name() + phaseSuffix);
        if (phasePart)
          partOrdinalMapping[part->mesh_meta_data_ordinal()] = phasePart->mesh_meta_data_ordinal();
      }
    }
  }
  return partOrdinalMapping;
}

static void move_disconnected_elements_to_phase(stk::mesh::BulkData & mesh, const std::string & phaseName)
{
  const stk::mesh::Selector sideSelector = get_nonvoid_selector(mesh.mesh_meta_data(), mesh.mesh_meta_data().side_rank());
  const stk::mesh::Selector elementSelector = get_nonvoid_selector(mesh.mesh_meta_data(), stk::topology::ELEMENT_RANK);
  const std::vector<stk::mesh::Entity> unattachedElems = krino::get_selected_owned_side_unattached_elements(mesh, elementSelector, sideSelector);

  const std::map<int,int> partOrdinalMapping = get_nonvoid_to_phase_part_ordinal_map(mesh.mesh_meta_data(), phaseName);
  krino::batch_convert_elements_and_their_sides(mesh, partOrdinalMapping, unattachedElems);
}

static void move_elements_not_in_largest_group_to_phase(stk::mesh::BulkData & mesh, const std::string & phaseName)
{
  const stk::mesh::Selector elementSelector = get_nonvoid_selector(mesh.mesh_meta_data(), stk::topology::ELEMENT_RANK);
  const std::vector<stk::mesh::Entity> elemsNotInLargestGroup = krino::find_owned_elements_that_are_not_in_the_largest_group_of_selected_side_attached_elements(mesh, elementSelector);

  const std::map<int,int> partOrdinalMapping = get_nonvoid_to_phase_part_ordinal_map(mesh.mesh_meta_data(), phaseName);
  krino::batch_convert_elements_and_their_sides(mesh, partOrdinalMapping, elemsNotInLargestGroup);
}

TEST(DecomposeMeshAndComputeSensitivities, readMeshInitializeDecomposeResetInitializeDecompose)
{
  const std::string initialMeshName = "mesh.g";
  krino::generate_and_write_bounding_box_mesh(stk::topology::TETRAHEDRON_4, {0.,0.,0.}, {1.,1.,1.}, 0.3333, initialMeshName);
  const krino::ScopedTestFile generatedMeshCleanup(initialMeshName);

  krino::MeshFromFile meshFromFile(initialMeshName, MPI_COMM_WORLD);

  const std::vector<krino::LS_Field> lsFields = krino::LSPerInterfacePolicy::setup_levelsets_on_all_blocks_with_void_phase_for_any_negative_levelset(meshFromFile.meta_data(), 1);
  const bool doComputeClosestPointSensitivities = false;
  setup_fields_for_conforming_decomposition(meshFromFile.meta_data(), lsFields, false, doComputeClosestPointSensitivities);

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

  decompose_mesh_to_conform_to_levelsets(meshFromFile.bulk_data(), lsFields, false, doComputeClosestPointSensitivities);

  const bool doWriteMesh = false;
  if (doWriteMesh)
    output_full_mesh(meshFromFile.bulk_data(), "output1.e", 1, 0.0);

  const std::map<std::string,size_t> decompPartCounts = get_part_counts(meshFromFile.bulk_data());

  krino::CDMesh::reset_mesh_to_original_undecomposed_state(meshFromFile.bulk_data());

  expect_part_counts(meshFromFile.bulk_data(), origPartCounts);

  if (doWriteMesh)
    output_full_mesh(meshFromFile.bulk_data(), "reset.e", 1, 1.0);

  compute_nodal_distance_from_spheres(meshFromFile.bulk_data(), meshFromFile.meta_data().coordinate_field(), lsFields[0].isovar, spheres );

  decompose_mesh_to_conform_to_levelsets(meshFromFile.bulk_data(), lsFields, false, doComputeClosestPointSensitivities);

  if (doWriteMesh)
    output_full_mesh(meshFromFile.bulk_data(), "output2.e", 1, 1.0);

  expect_part_counts(meshFromFile.bulk_data(), decompPartCounts);
}

static void test_moving_islands_to_separate_phase(const std::function<void(stk::mesh::BulkData&)> & island_removal_method, const std::string & islandPhaseName)
{
  const std::string initialMeshName = "mesh.g";
  krino::generate_and_write_bounding_box_mesh(stk::topology::TETRAHEDRON_4, {0.,0.,0.}, {1.,1.,1.}, 0.3333, initialMeshName);
  const krino::ScopedTestFile generatedMeshCleanup(initialMeshName);

  krino::MeshFromFile meshFromFile(initialMeshName, MPI_COMM_WORLD);

  const std::vector<krino::LS_Field> lsFields = krino::LSPerInterfacePolicy::setup_levelsets_on_all_blocks_with_void_phase_for_any_negative_levelset(meshFromFile.meta_data(), 1);
  create_extra_phase_per_block(meshFromFile.meta_data(), islandPhaseName);
  const bool doComputeClosestPointSensitivities = false;
  setup_fields_for_conforming_decomposition(meshFromFile.meta_data(), lsFields, false, doComputeClosestPointSensitivities);

  meshFromFile.populate_mesh();
  krino::activate_all_entities(meshFromFile.bulk_data(), krino::AuxMetaData::get(meshFromFile.meta_data()).active_part());

  const std::vector<std::pair<stk::math::Vector3d,double>> oneSphere
    {
      { {0.0,0.0,0.0}, 0.7 },
    };

  compute_nodal_distance_from_spheres(meshFromFile.bulk_data(), meshFromFile.meta_data().coordinate_field(), lsFields[0].isovar, oneSphere, -1);

  decompose_mesh_to_conform_to_levelsets(meshFromFile.bulk_data(), lsFields, false, doComputeClosestPointSensitivities);

  const bool doWriteMesh = false;
  if (doWriteMesh)
    output_full_mesh(meshFromFile.bulk_data(), "output1.e", 1, 0.0);

  const stk::mesh::Part & block_1 = *meshFromFile.meta_data().get_part("block_1");
  const size_t numElementsInBlockFromOneSphere = krino::get_global_num_entities(meshFromFile.bulk_data(), block_1);

  krino::CDMesh::reset_mesh_to_original_undecomposed_state(meshFromFile.bulk_data());

  const std::vector<std::pair<stk::math::Vector3d,double>> oneConnectedAndOneDisconnectedSphere
    {
      { {0.0,0.0,0.0}, 0.7 },
      { {0.7,0.7,0.7}, 0.25 }
    };

  compute_nodal_distance_from_spheres(meshFromFile.bulk_data(), meshFromFile.meta_data().coordinate_field(), lsFields[0].isovar, oneConnectedAndOneDisconnectedSphere, -1);

  decompose_mesh_to_conform_to_levelsets(meshFromFile.bulk_data(), lsFields, false, doComputeClosestPointSensitivities);

  island_removal_method(meshFromFile.bulk_data());

  if (doWriteMesh)
    output_full_mesh(meshFromFile.bulk_data(), "output2.e", 1, 1.0);

  const size_t numElementsInBlockAfterDeletingSecondDisconnectedSphere = krino::get_global_num_entities(meshFromFile.bulk_data(), block_1);

  EXPECT_EQ(numElementsInBlockFromOneSphere, numElementsInBlockAfterDeletingSecondDisconnectedSphere);
}

TEST(DecomposeMeshAndComputeSensitivities, readMeshInitializeDecomposeAndMoveVolNotConnectedToSideToVoid)
{
  const std::string islandPhaseName = "void";
  auto island_removal_method = [=](stk::mesh::BulkData& mesh){ move_disconnected_elements_to_phase(mesh, islandPhaseName); };
  test_moving_islands_to_separate_phase(island_removal_method, islandPhaseName);
}

TEST(DecomposeMeshAndComputeSensitivities, readMeshInitializeDecomposeAndMoveVolNotConnectedToLargestVolToDisconnected)
{
  const std::string islandPhaseName = "disconnected";
  auto island_removal_method = [=](stk::mesh::BulkData& mesh){ move_elements_not_in_largest_group_to_phase(mesh, islandPhaseName); };
  test_moving_islands_to_separate_phase(island_removal_method, islandPhaseName);
}

class DecomposeMeshAndCheckConvergenceOfAreaAndSensitivity : public krino::OptimizedOnlyTest
{
protected:

  double facet_area_and_sensitivities(const std::array<stk::math::Vector3d,2> & facetCoords, double * dArea)
  {
    return krino::SegmentWithSens::length_and_optional_sensitivities(facetCoords[0], facetCoords[1], dArea);
  }

  double facet_area_and_sensitivities(const std::array<stk::math::Vector3d,3> & facetCoords, double * dArea)
  {
    return krino::TriangleWithSens::area_and_optional_sensitivities(facetCoords[0], facetCoords[1], facetCoords[2], dArea);
  }

  template <unsigned DIM>
  std::pair<double,double> build_circle_or_sphere_conforming_mesh_and_compute_area_and_sensitivity(const double meshSize, const double radius, const bool doComputeClosestPointSensitivities, const bool doSnapping, const bool doWriteMesh)
  {
    constexpr unsigned numFacetNodes = DIM;
    std::vector<krino::LS_Field> lsFields;

    std::unique_ptr<krino::BoundingBoxMesh> mesh = build_circle_or_sphere_conforming_bounding_box_mesh(DIM, {-1.,-1.,-1.}, {1.,1.,1.}, meshSize, {0.0,0.0,0.0}, radius, doComputeClosestPointSensitivities, doSnapping, doWriteMesh, lsFields);

    std::vector<krino::LevelSetShapeSensitivity> sensitivities;
    std::vector<std::array<size_t,numFacetNodes>> facetsSensitivityIndices;
    krino::fill_levelset_facets_and_shape_sensitivities(mesh->bulk_data(), lsFields[0], doComputeClosestPointSensitivities, facetsSensitivityIndices, sensitivities);

    double area = 0.;
    double dArea_dLs = 0.;

    double dArea[3*numFacetNodes];
    for (auto & facetSensIndices : facetsSensitivityIndices)
    {
      const std::array<const krino::LevelSetShapeSensitivity*,numFacetNodes> facetNodeSensitivities = get_facet_node_sensitivities(sensitivities, facetSensIndices);
      const std::array<stk::mesh::Entity,numFacetNodes> facetNodes = get_facet_nodes(mesh->bulk_data(), facetNodeSensitivities);
      const std::array<stk::math::Vector3d,numFacetNodes> facetCoords = krino::get_nodes_vector(mesh->bulk_data(), mesh->meta_data().coordinate_field(), facetNodes, DIM);
      const double facetArea = facet_area_and_sensitivities(facetCoords, dArea);

      double dFacetArea_dLs = 0.;
      for (unsigned n=0; n<numFacetNodes; ++n)
      {
        const stk::math::Vector3d dArea_dX(dArea+3*n);
        for (auto & dX_dLs : facetNodeSensitivities[n]->dCoordsdParentLevelSets)
          dFacetArea_dLs += Dot(dArea_dX, dX_dLs);
      }

      area += facetArea;
      dArea_dLs += dFacetArea_dLs;
    }

    krino::all_reduce_sum(mesh->bulk_data().parallel(), area, dArea_dLs);
    return std::make_pair(area, dArea_dLs);
  }

  std::pair<double,double> compute_gold_area_and_sensitivity(const unsigned dim, const double radius)
  {
    if (2 == dim)
      return std::make_pair(2.*M_PI*radius, -2.*M_PI);
    return std::make_pair(4*M_PI*radius*radius, -8*M_PI*radius);
  }

  template <unsigned DIM>
  void check_convergence_of_area_and_sensitivity_for_circle_or_sphere_conforming_meshes(const unsigned nRefine, const bool doComputeClosestPointSensitivities, const bool doSnapping, const bool doWriteMesh)
  {
    const double radius = 0.76;
    const auto & [goldArea, gold_dArea_dLS] = compute_gold_area_and_sensitivity(DIM, radius);

    for(unsigned r=0; r<nRefine; r++)
    {
      const double meshSize = 0.33 / std::pow(2., r);
      const auto & [area, dArea_dLs] = build_circle_or_sphere_conforming_mesh_and_compute_area_and_sensitivity<DIM>(meshSize, radius, doComputeClosestPointSensitivities, doSnapping, doWriteMesh);

      const double errArea = std::abs(area-goldArea);
      const double errdArea = std::abs(dArea_dLs-gold_dArea_dLS);

      if (0 == stk::EnvData::parallel_rank())
        std::cout << "Convergence " << meshSize << " " << errArea << " " << errdArea << std::endl;

      EXPECT_LT(errArea, 2.5*meshSize*meshSize);
      EXPECT_LT(errdArea, 15.*errArea); // also appears to converge quadratically, but with larger magnitude and much more variation
    }
  }
};

TEST_F(DecomposeMeshAndCheckConvergenceOfAreaAndSensitivity, checkSphereAreaAndSensitivity)
{
  check_convergence_of_area_and_sensitivity_for_circle_or_sphere_conforming_meshes<3>(3, false, false, false);
}

TEST_F(DecomposeMeshAndCheckConvergenceOfAreaAndSensitivity, checkSnappedSphereAreaAndSensitivity)
{
  check_convergence_of_area_and_sensitivity_for_circle_or_sphere_conforming_meshes<3>(3, false, true, false);
}

TEST_F(DecomposeMeshAndCheckConvergenceOfAreaAndSensitivity, checkCircleAreaAndSensitivity)
{
  check_convergence_of_area_and_sensitivity_for_circle_or_sphere_conforming_meshes<2>(5, false, false, false);
}

TEST_F(DecomposeMeshAndCheckConvergenceOfAreaAndSensitivity, checkSnappedCircleAreaAndSensitivity)
{
  check_convergence_of_area_and_sensitivity_for_circle_or_sphere_conforming_meshes<2>(5, false, true, false);
}

TEST_F(DecomposeMeshAndCheckConvergenceOfAreaAndSensitivity, checkSphereAreaAndClosestPointSensitivity)
{
  check_convergence_of_area_and_sensitivity_for_circle_or_sphere_conforming_meshes<3>(3, true, false, false);
}

TEST_F(DecomposeMeshAndCheckConvergenceOfAreaAndSensitivity, checkSnappedSphereAreaAndClosestPointSensitivity)
{
  check_convergence_of_area_and_sensitivity_for_circle_or_sphere_conforming_meshes<3>(3, true, true, false);
}

TEST_F(DecomposeMeshAndCheckConvergenceOfAreaAndSensitivity, checkCircleAreaAndClosestPointSensitivity)
{
  check_convergence_of_area_and_sensitivity_for_circle_or_sphere_conforming_meshes<2>(5, true, false, false);
}

TEST_F(DecomposeMeshAndCheckConvergenceOfAreaAndSensitivity, checkSnappedCircleAreaAndClosestPointSensitivity)
{
  check_convergence_of_area_and_sensitivity_for_circle_or_sphere_conforming_meshes<2>(5, true, true, false);
}

static void expect_planar_sensitivity(const stk::mesh::BulkData & mesh,
    const stk::math::Vector3d & planeNormal,
    const stk::mesh::Part & childNodePart,
    const stk::mesh::Selector & lsInterfaceSurf,
    const bool doSnapping,
    const bool doComputeClosestPointSensitivities,
    const bool doOffset,
    const krino::LevelSetShapeSensitivity & sens)
{
  const stk::mesh::Entity interfaceNode = mesh.get_entity(stk::topology::NODE_RANK, sens.interfaceNodeId);
  EXPECT_TRUE(lsInterfaceSurf(mesh.bucket(interfaceNode)));

  for (auto parentNodeId : sens.parentNodeIds)
  {
    const stk::mesh::Entity parentNode = mesh.get_entity(stk::topology::NODE_RANK, parentNodeId);
    EXPECT_FALSE(mesh.bucket(parentNode).member(childNodePart));
  }

  stk::math::Vector3d sensTot = stk::math::Vector3d::ZERO;
  for (auto & dCoordsdParentLevelSet : sens.dCoordsdParentLevelSets)
  {
    EXPECT_TRUE(Dot(-planeNormal, dCoordsdParentLevelSet) >= 0.);
    sensTot += dCoordsdParentLevelSet;
  }

  if (doComputeClosestPointSensitivities)
  {
    krino::expect_near(-planeNormal, sensTot);
  }
  else
  {
    const double tol = 1.e-6;
    EXPECT_NEAR(1., Dot(-planeNormal, sensTot), tol);
  }
}

static void expect_correct_planar_sensitivities(const stk::mesh::BulkData & mesh,
    const std::vector<krino::LS_Field> & lsFields,
    const std::vector<stk::math::Vector3d> & planeNormals,
    const bool doSnapping,
    const bool doComputeClosestPointSensitivities,
    const bool doOffset,
    const std::vector<std::vector<krino::LevelSetShapeSensitivity>> & lsSens)
{
  STK_ThrowRequire(lsFields.size() == planeNormals.size() && lsFields.size() == lsSens.size());
  const krino::CDFEM_Support & cdfemSupport = krino::CDFEM_Support::get(mesh.mesh_meta_data());
  auto & phaseSupport = krino::Phase_Support::get(mesh.mesh_meta_data());

  for (unsigned ls=0; ls<lsFields.size(); ++ls)
  {
    const stk::mesh::Selector lsInterfaceSurf = phaseSupport.get_negative_levelset_interface_selector(lsFields[ls].identifier);

    for(auto && sens : lsSens[ls])
      expect_planar_sensitivity(mesh,
          planeNormals[ls],
          cdfemSupport.get_child_node_part(),
          lsInterfaceSurf,
          doSnapping,
          doComputeClosestPointSensitivities,
          doOffset,
          sens);
  }
}

static void test_correct_planar_sensitivities(const bool doSnapping, const bool doComputeClosestPointSensitivities, const bool doOffset, const bool doOutput)
{
  std::cout << "Do snapping = " << doSnapping
      << ", do compute closest point sensitivities = " << doComputeClosestPointSensitivities
      << ", do offset from mesh nodes = " << doOffset << "\n";
  const double meshSize = 1.;

  const stk::topology elemTopo = stk::topology::TETRAHEDRON_4;
  krino::BoundingBoxMesh bboxMesh(elemTopo, MPI_COMM_WORLD);

  const std::vector<krino::LS_Field> lsFields = krino::LSPerInterfacePolicy::setup_levelsets_on_all_blocks(bboxMesh.meta_data(), 2);
  setup_fields_for_conforming_decomposition(bboxMesh.meta_data(), lsFields, doSnapping, doComputeClosestPointSensitivities);

  krino::populate_bounding_box_mesh_and_activate(bboxMesh, {-0.5,-0.5,-0.5}, {1.5,1.5,1.5}, meshSize);

  const std::vector<stk::math::Vector3d> planeNormals = { {1.,0.,0.}, {0.,1.,0.} };
  const double offset = doOffset ? 0.1 : 0.0;
  compute_nodal_distance_from_plane(bboxMesh.bulk_data(), bboxMesh.meta_data().coordinate_field(), lsFields[0].isovar, planeNormals[0], offset);
  compute_nodal_distance_from_plane(bboxMesh.bulk_data(), bboxMesh.meta_data().coordinate_field(), lsFields[1].isovar, planeNormals[1], offset);

  decompose_mesh_to_conform_to_levelsets(bboxMesh.bulk_data(), lsFields, doSnapping, doComputeClosestPointSensitivities);

  const auto lsSens = krino::get_levelset_shape_sensitivities(bboxMesh.bulk_data(), lsFields, doComputeClosestPointSensitivities);
  expect_correct_planar_sensitivities(bboxMesh.bulk_data(), lsFields, planeNormals, doSnapping, doComputeClosestPointSensitivities, doOffset, lsSens);

  if (doOutput)
  {
    const std::string filename = "plane_" + std::string(doSnapping ? "snapped" : "cut") +
        std::string(doComputeClosestPointSensitivities ? "_closestPtSens" : "_cutEdgeSens") +
        std::string(doOffset ? "_offset.e" : "_throughNodes.e");
    output_full_mesh(bboxMesh.bulk_data(), filename, 1, 0.0);
  }
}

TEST(DecomposeMeshMultipleLS, TwoLevelSets)
{
  const bool doOutput = false;
  test_correct_planar_sensitivities(/*doSnapping=*/ true, /*doClosePt=*/ true, /*doOffset=*/false, doOutput);
  test_correct_planar_sensitivities(/*doSnapping=*/false, /*doClosePt=*/ true, /*doOffset=*/false, doOutput);

  test_correct_planar_sensitivities(/*doSnapping=*/ true, /*doClosePt=*/ true, /*doOffset=*/ true, doOutput);
  test_correct_planar_sensitivities(/*doSnapping=*/false, /*doClosePt=*/ true, /*doOffset=*/ true, doOutput);

  test_correct_planar_sensitivities(/*doSnapping=*/ true, /*doClosePt=*/false, /*doOffset=*/false, doOutput);
  test_correct_planar_sensitivities(/*doSnapping=*/false, /*doClosePt=*/false, /*doOffset=*/false, doOutput);

  test_correct_planar_sensitivities(/*doSnapping=*/ true, /*doClosePt=*/false, /*doOffset=*/ true, doOutput);
  test_correct_planar_sensitivities(/*doSnapping=*/false, /*doClosePt=*/false, /*doOffset=*/ true, doOutput);
}

class SensitivityRightTrisFixture : public krino::StkMeshTriFixture
{
public:
  SensitivityRightTrisFixture()
  {
    set_valid_proc_sizes_for_test({1,2,4});
  }
  void populate_mesh()
  {
    const std::vector<unsigned> elementBlockIDs = {1,2,2,1};
    this->build_mesh_from_spec_in_blocks(meshSpec, elementBlockIDs);
  }
  using krino::StkMeshTriFixture::mMesh;
  using krino::StkMeshTriFixture::mBuilder;

  void setup_levelsets_and_populate_mesh(const bool doSnapping, const bool doComputeClosestPointSensitivities)
  {
    mBuilder.create_block_parts({1,2});
    lsFields = krino::LSPerInterfacePolicy::setup_levelsets_on_all_blocks(mMesh.mesh_meta_data(), 1);
    setup_fields_for_conforming_decomposition(mMesh.mesh_meta_data(), lsFields, doSnapping, doComputeClosestPointSensitivities);

    populate_mesh();
    krino::activate_all_entities(mMesh, krino::AuxMetaData::get(mMesh.mesh_meta_data()).active_part());
  }

  void test_correct_sensitivities_for_plane_with_tiny_offset(const bool doSnapping, const bool doComputeClosestPointSensitivities, const bool doOutput)
  {
    std::cout << "Do snapping = " << doSnapping << ", do compute closest point sensitivities = " << doComputeClosestPointSensitivities << "\n";
    const bool doLargeOffsetFromBackgroundNodes = false;

    setup_levelsets_and_populate_mesh(doSnapping, doComputeClosestPointSensitivities);

    const krino::CDFEM_Support & cdfemSupport = krino::CDFEM_Support::get(mMesh.mesh_meta_data());
    const double snapTol = cdfemSupport.get_snapper().get_edge_tolerance();

    const stk::math::Vector3d planeNormal{0.,1.,0.};
    compute_nodal_distance_from_plane(mMesh, mMesh.mesh_meta_data().coordinate_field(), lsFields[0].isovar, planeNormal, -0.5*snapTol);

    decompose_mesh_to_conform_to_levelsets(mMesh, lsFields, doSnapping, doComputeClosestPointSensitivities);

    const auto LSsens = krino::get_levelset_shape_sensitivities(mMesh, lsFields[0], doComputeClosestPointSensitivities);


    const krino::Phase_Support & phaseSupport = krino::Phase_Support::get(mMesh.mesh_meta_data());
    const stk::mesh::Selector lsInterfaceSurf = phaseSupport.get_negative_levelset_interface_selector(lsFields[0].identifier);
    for(auto && sens : LSsens)
    {
      expect_planar_sensitivity(mMesh,
                planeNormal,
                cdfemSupport.get_child_node_part(),
                lsInterfaceSurf,
                doSnapping,
                doComputeClosestPointSensitivities,
                doLargeOffsetFromBackgroundNodes,
                sens);
    }

    if (doOutput)
    {
      const std::string filename = "decomposed_four_tri" + std::string(doSnapping ? "_snapped.e" : "_cut.e");
      output_full_mesh(mMesh, filename, 1, 0.0);
    }
  }
protected:
  krino::FourRightTris meshSpec;
  std::vector<krino::LS_Field> lsFields;
};

// NOTE:
// 06/2026: Would it be better to be able to run all of these scenarios in one test?
// If so, we need to have the setup not depend on doSnapping and doComputeClosestPointSensitivities.
// As it is currently, you are locked in after calling setup_fields_for_conforming_decomposition.
// That seems fragile
TEST_F(SensitivityRightTrisFixture, interfaceSnapToNodeOnBlockBoundary_doSnapping_doClosestPt)
{
  test_correct_sensitivities_for_plane_with_tiny_offset(/*doSnapping=*/ true, /*doClosePt=*/ true, /*doOutput=*/false);
}

TEST_F(SensitivityRightTrisFixture, interfaceSnapToNodeOnBlockBoundary_noSnapping_doClosestPt)
{
  test_correct_sensitivities_for_plane_with_tiny_offset(/*doSnapping=*/false, /*doClosePt=*/ true, /*doOutput=*/false);
}

TEST_F(SensitivityRightTrisFixture, interfaceSnapToNodeOnBlockBoundary_doSnapping_noClosestPt)
{
  test_correct_sensitivities_for_plane_with_tiny_offset(/*doSnapping=*/ true, /*doClosePt=*/false, /*doOutput=*/false);
}

TEST_F(SensitivityRightTrisFixture, interfaceSnapToNodeOnBlockBoundary_noSnapping_noClosestPt)
{
  test_correct_sensitivities_for_plane_with_tiny_offset(/*doSnapping=*/false, /*doClosePt=*/false, /*doOutput=*/false);
}

static void test_find_edge_location_for_node_along_edge(const stk::math::Vector3d & edgeCoord0,
  const stk::math::Vector3d & edgeCoord1,
  const double edgeLoc)
{
  const stk::math::Vector3d snapDisp = edgeLoc*(edgeCoord1-edgeCoord0);
  const auto & [foundEdgeLoc, edgeLocErr] =  krino::snapped_node_edge_position_and_square_error(3, snapDisp, edgeCoord0, edgeCoord1);
  const double tol = std::numeric_limits<double>::epsilon();
  EXPECT_LT(edgeLocErr, tol);
  EXPECT_LT(std::abs(edgeLoc-foundEdgeLoc), tol);
}

static void test_find_edge_location_for_node_NOT_along_edge(const stk::math::Vector3d & edgeCoord0,
  const stk::math::Vector3d & edgeCoord1,
  const stk::math::Vector3d & snapDisp)
{
  const auto & [foundEdgeLoc, edgeLocErr] =  krino::snapped_node_edge_position_and_square_error(3, snapDisp, edgeCoord0, edgeCoord1);
  const double tol = std::numeric_limits<double>::epsilon();
  EXPECT_GT(edgeLocErr, tol);
}

TEST(findSnapLocation, nodeAlongEdgeOrOffEdge_computeLocationAndError)
{
  test_find_edge_location_for_node_along_edge(stk::math::Vector3d{0,0,0.1}, stk::math::Vector3d{1,2,1}, 1.e-1);
  test_find_edge_location_for_node_along_edge(stk::math::Vector3d{0.1,0,0}, stk::math::Vector3d{1,1,2}, 1.e-6);
  test_find_edge_location_for_node_along_edge(stk::math::Vector3d{0,0.1,0}, stk::math::Vector3d{2,1,1}, 1.e-9);
  test_find_edge_location_for_node_along_edge(stk::math::Vector3d{0,0,0}, stk::math::Vector3d{1,1,1}, 1.e-12);
  test_find_edge_location_for_node_along_edge(stk::math::Vector3d{0,0,0}, stk::math::Vector3d{1,1,1}, 1.e-15);
  test_find_edge_location_for_node_along_edge(stk::math::Vector3d{0,0,0}, stk::math::Vector3d{1,1,1}, 0);

  test_find_edge_location_for_node_NOT_along_edge(stk::math::Vector3d{0,0,0}, stk::math::Vector3d{1,1,1}, stk::math::Vector3d{0.1,0.100001,0.1});
  test_find_edge_location_for_node_NOT_along_edge(stk::math::Vector3d{0,0,0}, stk::math::Vector3d{1,1,1}, stk::math::Vector3d{-0.1,-0.1,-0.1});
}
