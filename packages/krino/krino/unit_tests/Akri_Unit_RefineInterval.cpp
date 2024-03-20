#include <Akri_AuxMetaData.hpp>
#include <Akri_MeshHelpers.hpp>
#include <gtest/gtest.h>
#include <Akri_MeshSpecs.hpp>
#include <Akri_StkMeshBuilder.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <Akri_RefineNearLevelSets.hpp>
#include <Akri_StkMeshFixture.hpp>
#include <Akri_UnitMeshUtils.hpp>


stk::mesh::Field<double> & create_levelset_field(stk::mesh::BulkData & mesh, const std::string & levelsetName)
{
  stk::mesh::Field<double> & lsField = mesh.mesh_meta_data().declare_field<double>(stk::topology::NODE_RANK, levelsetName);
  stk::mesh::put_field_on_mesh(lsField, mesh.mesh_meta_data().universal_part(), 1, 1, nullptr);
  return lsField;
}

std::vector<stk::mesh::Field<double>*> create_levelset_fields(stk::mesh::BulkData & mesh)
{
  mesh.mesh_meta_data().enable_late_fields();
  std::vector<stk::mesh::Field<double>*> levelSetFields;
  levelSetFields.push_back(&create_levelset_field(mesh,"LS1"));
  levelSetFields.push_back(&create_levelset_field(mesh,"LS2"));
  return levelSetFields;
}

void initialize_circle_step_function(const stk::mesh::BulkData & mesh, const stk::mesh::Field<double> & lsField, const std::array<double,2> & circleCenter, const double circleRadius)
{
  const auto & coordsField = static_cast<const stk::mesh::Field<double>&>(*mesh.mesh_meta_data().coordinate_field());
  for (auto & bucket : mesh.buckets(stk::topology::NODE_RANK))
  {
    for (auto node : *bucket)
    {
      double * coords = stk::mesh::field_data(coordsField, node);
      double * ls = stk::mesh::field_data(lsField, node);
      const double circleDist = std::sqrt((coords[0] - circleCenter[0])*(coords[0] - circleCenter[0]) + (coords[1] - circleCenter[1])*(coords[1] - circleCenter[1]));
      *ls = ((circleDist - circleRadius) < 0. ? -1. : 1.);
    }
  }
}

std::function<void()> build_levelset_initialization_function(const stk::mesh::BulkData & mesh, const std::vector<stk::mesh::Field<double>*> & levelSetFields)
{
  std::function<void()> initialize_levelsets =
    [&mesh, &levelSetFields]()
    {
      initialize_circle_step_function(mesh, *levelSetFields[0], {{0.15, -0.2}}, 0.25);
      initialize_circle_step_function(mesh, *levelSetFields[1], {{-0.2, 0.15}}, 0.3);
    };
  return initialize_levelsets;
}

void test_refinement_level_for_elements_attached_to_node(const stk::mesh::BulkData & mesh, const stk::mesh::Part & activePart, const stk::mesh::Entity node, const int goldRefinementLevel)
{
  const double goldSize = 1./(1<<goldRefinementLevel);
  if (mesh.is_valid(node))
  {
    const double maxElemSizeForNode = krino::compute_maximum_size_of_selected_elements_using_node(mesh, activePart, node);
    EXPECT_NEAR(goldSize, maxElemSizeForNode, 1.e-3*goldSize);
  }
}

TEST(RefineWithinDistanceOfLevelSets, twoIntersectingLevelSets_correctElementsGetRefined)
{
  krino::StkTriMeshAndBuilder meshAndBuilder;
  krino::QuadSplit4Tri meshSpec;
  meshAndBuilder.build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,1,1,1});

  stk::mesh::BulkData & mesh = meshAndBuilder.get_mesh();
  stk::mesh::Part & activePart = meshAndBuilder.get_aux_meta().active_part();

  const std::vector<stk::mesh::Field<double>*> levelSetFields = create_levelset_fields(mesh);
  const auto initialize_levelsets = build_levelset_initialization_function(mesh, levelSetFields);

  const std::array<double,2> refinementDistanceInterval{-0.02,0.05};
  const unsigned numRefinementLevels = 6;

  for (auto & levelSetField : levelSetFields)
    krino::refine_elements_that_intersect_distance_interval_from_levelset(mesh, activePart, *levelSetField, initialize_levelsets, refinementDistanceInterval, numRefinementLevels);

  const bool doWriteMesh = false;
  if (doWriteMesh)
    meshAndBuilder.write_mesh("test.e");

  test_refinement_level_for_elements_attached_to_node(mesh, activePart, meshAndBuilder.get_assigned_node_for_index(2), 2);
  test_refinement_level_for_elements_attached_to_node(mesh, activePart, meshAndBuilder.get_assigned_node_for_index(4), numRefinementLevels);

  if (mesh.parallel_rank() == 0) // Whole mesh on 0
  {
    stk::mesh::Entity nodeOffCenter = krino::find_local_node_closest_to_location(mesh, mesh.mesh_meta_data().universal_part(), meshAndBuilder.get_coordinates_field(), stk::math::Vector3d(0.05,-0.05,0));
    std::cout << "Checking element size for elements using node " << mesh.identifier(nodeOffCenter) << std::endl;
    test_refinement_level_for_elements_attached_to_node(mesh, activePart, nodeOffCenter, numRefinementLevels);
  }
}

