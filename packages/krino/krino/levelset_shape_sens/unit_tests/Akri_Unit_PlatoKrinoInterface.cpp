#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>

// Deliberately the ONLY krino header this test includes, plus the shared unit
// test helpers that write and remove a temp mesh file. This exercises (and
// documents) the contract that a consumer needs nothing but Akri_PlatoKrinoInterface.hpp
// (and the krino value types it includes) besides STK.
#include <Akri_PlatoKrinoInterface.hpp>

#include <Akri_ScopedTestFile.hpp>
#include <Akri_UnitMeshUtils.hpp>
#include <Ioss_NodeBlock.h>
#include <Ioss_Region.h>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_topology/topology.hpp>
#include <stk_util/parallel/Parallel.hpp>

namespace
{
using namespace PlatoKrinoInterface;

constexpr auto EDGE_DIRECTION = krino::ShapeSensitivityMethod::EDGE_DIRECTION_USING_EDGE_GRADIENT;
constexpr auto SURFACE_NORMAL = krino::ShapeSensitivityMethod::SURFACE_NORMAL_USING_EDGE_GRADIENT;

stk::math::Vector3d total_sensitivity(const krino::LevelSetShapeSensitivity & sens)
{
  stk::math::Vector3d sum = stk::math::Vector3d::ZERO;
  for (const auto & d : sens.dCoordsdParentLevelSets)
    sum += d;
  return sum;
}

double distance_from_origin(const std::array<double, 3> & x)
{
  return std::sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
}

stk::math::Vector3d node_coordinates(const stk::mesh::BulkData & mesh, const stk::mesh::Entity node)
{
  const auto coords = mesh.mesh_meta_data().coordinate_field()->data<double>().entity_values(node);
  return stk::math::Vector3d(coords(stk::mesh::ComponentIdx(0)),
      coords(stk::mesh::ComponentIdx(1)),
      coords(stk::mesh::ComponentIdx(2)));
}

std::array<stk::math::Vector3d, 3> facet_coordinates(const stk::mesh::BulkData & mesh,
    const krino::InterfaceFacetsWithSensitivities<3> & facets,
    const size_t facet)
{
  const std::array<stk::mesh::Entity, 3> nodes = facets.nodes(mesh, facet);
  return {{node_coordinates(mesh, nodes[0]),
      node_coordinates(mesh, nodes[1]),
      node_coordinates(mesh, nodes[2])}};
}

// The global number of sensitivities and the sum of the magnitudes of each node's total
// sensitivity, to compare the sensitivities of two decompositions.
std::pair<size_t, double>
global_count_and_checksum(const std::vector<krino::LevelSetShapeSensitivity> & sensitivities)
{
  unsigned long count = sensitivities.size();
  double checksum = 0.;
  for (const auto & sens : sensitivities)
    checksum += total_sensitivity(sens).length();
  MPI_Allreduce(MPI_IN_PLACE, &count, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &checksum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return {count, checksum};
}

// The global area of the interface facets of level set lsIndex of a 3D conforming mesh.
double
global_area(const ConformingMeshWithSensitivities & conformingMesh, const unsigned lsIndex = 0)
{
  const krino::InterfaceFacetsWithSensitivities<3> facets =
      conformingMesh.interface_facets<3>(lsIndex);
  double area = 0.;
  for (size_t facet = 0; facet < facets.num_facets(); ++facet)
  {
    const std::array<stk::math::Vector3d, 3> x =
        facet_coordinates(conformingMesh.bulk_data(), facets, facet);
    area += 0.5 * Cross(x[1] - x[0], x[2] - x[0]).length();
  }
  MPI_Allreduce(MPI_IN_PLACE, &area, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return area;
}

// The global number of sides in the part partName, 0 if there is none.
size_t global_num_sides_in_part(const stk::mesh::BulkData & mesh, const std::string & partName)
{
  const stk::mesh::MetaData & meta = mesh.mesh_meta_data();
  const stk::mesh::Part * part = meta.get_part(partName);
  unsigned long count = (part == nullptr) ?
      0 :
      stk::mesh::count_selected_entities(
          *part & meta.locally_owned_part(), mesh.buckets(meta.side_rank()));
  MPI_Allreduce(MPI_IN_PLACE, &count, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  return count;
}

// The global volume of the tetrahedra of a mesh.
double global_volume(const stk::mesh::BulkData & mesh, const stk::mesh::Selector & elementSelector)
{
  const stk::mesh::MetaData & meta = mesh.mesh_meta_data();
  double volume = 0.;
  for (auto * bucketPtr :
      mesh.get_buckets(stk::topology::ELEMENT_RANK, elementSelector & meta.locally_owned_part()))
  {
    for (auto elem : *bucketPtr)
    {
      const stk::mesh::Entity * nodes = mesh.begin_nodes(elem);
      const stk::math::Vector3d x0 = node_coordinates(mesh, nodes[0]);
      volume +=
          Dot(node_coordinates(mesh, nodes[1]) - x0,
              Cross(node_coordinates(mesh, nodes[2]) - x0, node_coordinates(mesh, nodes[3]) - x0)) /
          6.;
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &volume, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return volume;
}

double global_volume(const stk::mesh::BulkData & mesh)
{
  return global_volume(mesh, mesh.mesh_meta_data().universal_part());
}

// The global number of elements and the global area of the (triangular) sides of the selectors.
size_t global_num_elements(const stk::mesh::BulkData & mesh, const stk::mesh::Selector & selector)
{
  unsigned long count =
      stk::mesh::count_selected_entities(selector & mesh.mesh_meta_data().locally_owned_part(),
          mesh.buckets(stk::topology::ELEMENT_RANK));
  MPI_Allreduce(MPI_IN_PLACE, &count, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  return count;
}

double global_side_area(const stk::mesh::BulkData & mesh, const stk::mesh::Selector & selector)
{
  const stk::mesh::MetaData & meta = mesh.mesh_meta_data();
  double area = 0.;
  for (auto * bucketPtr : mesh.get_buckets(meta.side_rank(), selector & meta.locally_owned_part()))
  {
    for (auto side : *bucketPtr)
    {
      const stk::mesh::Entity * nodes = mesh.begin_nodes(side);
      const stk::math::Vector3d x0 = node_coordinates(mesh, nodes[0]);
      area += 0.5 *
          Cross(node_coordinates(mesh, nodes[1]) - x0, node_coordinates(mesh, nodes[2]) - x0)
              .length();
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &area, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return area;
}

std::set<std::string>
names_of_blocks_touching_surface(const stk::mesh::MetaData & meta, const std::string & surfaceName)
{
  std::set<std::string> blockNames;
  for (auto * block : meta.get_blocks_touching_surface(meta.get_part(surfaceName)))
    blockNames.insert(block->name());
  return blockNames;
}

// Writes, on processor 0, an ASCII STL of the closed box [-halfWidth, halfWidth]^3, with each
// triangle counterclockwise about its outward normal.
void write_box_stl(const std::string & fileName, const double halfWidth)
{
  if (stk::parallel_machine_rank(MPI_COMM_WORLD) == 0)
  {
    std::ofstream out(fileName);
    out << "solid box\n";
    for (int axis = 0; axis < 3; ++axis)
    {
      for (const double side : {-1., 1.})
      {
        const int b = (axis + 1) % 3;
        const int c = (axis + 2) % 3;
        std::array<stk::math::Vector3d, 4> corners;
        const std::array<std::array<double, 2>, 4> bc{
            {{{-1., -1.}}, {{1., -1.}}, {{1., 1.}}, {{-1., 1.}}}};
        for (int i = 0; i < 4; ++i)
        {
          corners[i][axis] = side * halfWidth;
          corners[i][b] = bc[i][0] * halfWidth;
          corners[i][c] = bc[i][1] * halfWidth;
        }
        // The corners go counterclockwise about +axis, so reverse them for the -axis side.
        if (side < 0.) std::swap(corners[1], corners[3]);
        stk::math::Vector3d normal = stk::math::Vector3d::ZERO;
        normal[axis] = side;
        for (const std::array<int, 3> & tri :
            {std::array<int, 3>{{0, 1, 2}}, std::array<int, 3>{{0, 2, 3}}})
        {
          out << " facet normal " << normal[0] << " " << normal[1] << " " << normal[2] << "\n";
          out << "  outer loop\n";
          for (const int i : tri)
            out << "   vertex " << corners[i][0] << " " << corners[i][1] << " " << corners[i][2]
                << "\n";
          out << "  endloop\n";
          out << " endfacet\n";
        }
      }
    }
    out << "endsolid box\n";
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

// The signed distance to the box [-halfWidth, halfWidth]^3, positive inside.
double box_level_set(const std::array<double, 3> & x, const double halfWidth)
{
  double maxAbs = 0.;
  double outsideSquared = 0.;
  for (int d = 0; d < 3; ++d)
  {
    maxAbs = std::max(maxAbs, std::abs(x[d]));
    const double excess = std::abs(x[d]) - halfWidth;
    if (excess > 0.) outsideSquared += excess * excess;
  }
  return (maxAbs <= halfWidth) ? halfWidth - maxAbs : -std::sqrt(outsideSquared);
}

// Expects level set fields 0 and 1 of the background mesh to agree on every owned node.
void expect_same_level_sets(const BackgroundMesh & background)
{
  const stk::mesh::BulkData & mesh = background.bulk_data();
  const auto ls0 = background.level_set_field(0).data<double>();
  const auto ls1 = background.level_set_field(1).data<double>();
  for (auto * bucketPtr :
      mesh.get_buckets(stk::topology::NODE_RANK, mesh.mesh_meta_data().locally_owned_part()))
    for (auto node : *bucketPtr)
      EXPECT_NEAR(ls1.entity_values(node)(), ls0.entity_values(node)(), 1.e-12)
          << node_coordinates(mesh, node);
}

BackgroundMesh box_mesh()
{
  return BackgroundMesh::bounding_box(3, {{-1., -1., -1.}}, {{1., 1., 1.}}, 0.4, MPI_COMM_WORLD);
}

Geometry sphere_at_origin(const double radius)
{
  Geometry geometry;
  geometry.spheres.push_back(SphereSpec{{{0., 0., 0.}}, radius});
  return geometry;
}

// The material half-space x < 0.1, whose interface x = 0.1 crosses the box [-1,1]^3 between
// the nodes of the meshes used here, with area 4.
Geometry half_space_below_x_equal_point_one()
{
  Geometry geometry;
  geometry.planes.push_back(PlaneSpec{{{1., 0., 0.}}, -0.1});
  return geometry;
}

// A box mesh, written to meshFile and read back, whose sides are in the sidesets Xminus, Xplus, ...
BackgroundMesh read_box_mesh_with_sidesets(const std::string & meshFile)
{
  krino::generate_and_write_bounding_box_mesh(
      stk::topology::TETRAHEDRON_4, {-1., -1., -1.}, {1., 1., 1.}, 0.4, meshFile, true);
  return BackgroundMesh::read(meshFile, MPI_COMM_WORLD, "");
}

ConformingMeshWithSensitivities decompose(const BackgroundMesh & background,
    const bool doSnapping,
    const krino::ShapeSensitivityMethod method)
{
  return ConformingMeshWithSensitivities(background.bulk_data(),
      {&background.level_set_field()},
      {.doSnapping = doSnapping, .sensitivityMethod = method});
}

ConformingMeshWithSensitivities decompose_geometry(BackgroundMesh & background,
    const Geometry & geometry,
    const bool doSnapping,
    const krino::ShapeSensitivityMethod method)
{
  initialize_field_from_geometry(background.bulk_data(), background.level_set_field(), geometry);
  return decompose(background, doSnapping, method);
}

void expect_same_sensitivities(const ConformingMeshWithSensitivities & conformingMesh,
    const ConformingMeshWithSensitivities & goldConformingMesh)
{
  const auto [count, checksum] = global_count_and_checksum(conformingMesh.shape_sensitivities());
  const auto [goldCount, goldChecksum] =
      global_count_and_checksum(goldConformingMesh.shape_sensitivities());
  EXPECT_EQ(goldCount, count);
  EXPECT_NEAR(goldChecksum, checksum, 1.e-10 * goldChecksum);
}

TEST(PlatoKrinoInterface, sphereOnBoundingBox_sensitivitiesHaveTwoParentsThatAreBackgroundNodes)
{
  BackgroundMesh background = box_mesh();
  EXPECT_EQ(1u, background.num_level_sets());
  const ConformingMeshWithSensitivities conformingMesh =
      decompose_geometry(background, sphere_at_origin(0.5), false, EDGE_DIRECTION);
  EXPECT_EQ(1u, conformingMesh.num_level_sets());

  const std::vector<krino::LevelSetShapeSensitivity> sensitivities =
      conformingMesh.shape_sensitivities();
  EXPECT_LT(0u, global_count_and_checksum(sensitivities).first);

  for (const auto & sens : sensitivities)
  {
    EXPECT_EQ(2u, sens.parentNodeIds.size());
    EXPECT_EQ(2u, sens.dCoordsdParentLevelSets.size());
    for (const auto parentId : sens.parentNodeIds)
      EXPECT_TRUE(background.bulk_data().is_valid(
          background.bulk_data().get_entity(stk::topology::NODE_RANK, parentId)));
  }
}

TEST(PlatoKrinoInterface, surfaceNormalSensitivitiesForPlane_summedSensitivityIsUnitNormalIntoVoid)
{
  BackgroundMesh background = box_mesh();
  const ConformingMeshWithSensitivities conformingMesh =
      decompose_geometry(background, half_space_below_x_equal_point_one(), false, SURFACE_NORMAL);

  const std::vector<krino::LevelSetShapeSensitivity> sensitivities =
      conformingMesh.shape_sensitivities();
  EXPECT_LT(0u, global_count_and_checksum(sensitivities).first);

  // The level set 0.1 - x is a distance function, so increasing it moves each node along +x, into
  // the void, at unit speed.
  for (const auto & sens : sensitivities)
  {
    const stk::math::Vector3d sum = total_sensitivity(sens);
    EXPECT_NEAR(1., sum[0], 1.e-12);
    EXPECT_NEAR(0., sum[1], 1.e-12);
    EXPECT_NEAR(0., sum[2], 1.e-12);
  }
}

// The material is inside the sphere: increasing the level set grows it, so the interface nodes
// move outward, and the facets are oriented out of the void phase, toward the center.
TEST(PlatoKrinoInterface, sphereGeometry_materialInsideSoNodesMoveOutwardAndFacetsPointInward)
{
  BackgroundMesh background = box_mesh();
  const ConformingMeshWithSensitivities conformingMesh =
      decompose_geometry(background, sphere_at_origin(0.5), false, EDGE_DIRECTION);
  const stk::mesh::BulkData & mesh = conformingMesh.bulk_data();

  const krino::InterfaceFacetsWithSensitivities<3> facets = conformingMesh.interface_facets<3>();
  for (const auto & sens : facets.sensitivities)
  {
    const stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, sens.interfaceNodeId);
    EXPECT_GT(Dot(total_sensitivity(sens), node_coordinates(mesh, node)), 0.);
  }

  for (size_t facet = 0; facet < facets.num_facets(); ++facet)
  {
    const std::array<stk::math::Vector3d, 3> x = facet_coordinates(mesh, facets, facet);
    const stk::math::Vector3d centroid = (x[0] + x[1] + x[2]) / 3.;
    EXPECT_LT(Dot(Cross(x[1] - x[0], x[2] - x[0]), centroid), 0.);
  }
}

TEST(PlatoKrinoInterface, planeInterfaceFacets_areaIsCrossSectionOfBox)
{
  BackgroundMesh background = box_mesh();
  const ConformingMeshWithSensitivities conformingMesh =
      decompose_geometry(background, half_space_below_x_equal_point_one(), false, EDGE_DIRECTION);

  const krino::InterfaceFacetsWithSensitivities<3> facets = conformingMesh.interface_facets<3>();
  for (size_t facet = 0; facet < facets.num_facets(); ++facet)
    for (const auto & x : facet_coordinates(conformingMesh.bulk_data(), facets, facet))
      EXPECT_NEAR(0.1, x[0], 1.e-12);
  EXPECT_NEAR(4., global_area(conformingMesh), 1.e-12);
}

TEST(PlatoKrinoInterface, initializeFromFunction_matchesGeometryInitialization)
{
  const double radius = 0.5;
  BackgroundMesh background = box_mesh();
  const ConformingMeshWithSensitivities geometryMesh =
      decompose_geometry(background, sphere_at_origin(radius), false, EDGE_DIRECTION);

  // Signed distance to the sphere, positive inside, stamped through the function initializer
  // (an optimizer would supply its own design values here).
  initialize_field_from_function(background.bulk_data(),
      background.level_set_field(),
      [radius](const std::array<double, 3> & x) { return radius - distance_from_origin(x); });
  const ConformingMeshWithSensitivities functionMesh = decompose(background, false, EDGE_DIRECTION);

  expect_same_sensitivities(functionMesh, geometryMesh);
}

TEST(PlatoKrinoInterface, setFieldValuesOfOwnedNodes_matchesFunctionInitialization)
{
  const double radius = 0.5;
  const LevelSetValueFunction sphereDistance = [radius](const std::array<double, 3> & x)
  { return radius - distance_from_origin(x); };

  BackgroundMesh background = box_mesh();
  initialize_field_from_function(
      background.bulk_data(), background.level_set_field(), sphereDistance);
  const ConformingMeshWithSensitivities functionMesh = decompose(background, true, EDGE_DIRECTION);

  const stk::mesh::BulkData & mesh = background.bulk_data();
  std::vector<stk::mesh::EntityId> nodeIds;
  std::vector<double> values;
  for (auto * bucketPtr :
      mesh.get_buckets(stk::topology::NODE_RANK, mesh.mesh_meta_data().locally_owned_part()))
  {
    for (auto node : *bucketPtr)
    {
      const stk::math::Vector3d x = node_coordinates(mesh, node);
      nodeIds.push_back(mesh.identifier(node));
      values.push_back(sphereDistance({{x[0], x[1], x[2]}}));
    }
  }
  initialize_field_from_function(background.bulk_data(),
      background.level_set_field(),
      [](const std::array<double, 3> &) { return 1.; });
  set_field_values(background.bulk_data(), background.level_set_field(), nodeIds, values);
  const ConformingMeshWithSensitivities valuesMesh = decompose(background, true, EDGE_DIRECTION);

  expect_same_sensitivities(valuesMesh, functionMesh);
}

TEST(PlatoKrinoInterface, decomposeWithSnapping_backgroundMeshUnchanged)
{
  BackgroundMesh background = box_mesh();
  const stk::mesh::BulkData & mesh = background.bulk_data();
  initialize_field_from_geometry(
      background.bulk_data(), background.level_set_field(), sphere_at_origin(0.5));

  std::vector<size_t> origCounts;
  stk::mesh::count_entities(mesh.mesh_meta_data().universal_part(), mesh, origCounts);
  std::map<stk::mesh::EntityId, stk::math::Vector3d> origCoords;
  for (auto * bucketPtr : mesh.buckets(stk::topology::NODE_RANK))
    for (auto node : *bucketPtr)
      origCoords[mesh.identifier(node)] = node_coordinates(mesh, node);

  {
    const ConformingMeshWithSensitivities conformingMesh =
        decompose(background, true, EDGE_DIRECTION);
    EXPECT_LT(0u, global_count_and_checksum(conformingMesh.shape_sensitivities()).first);
  }

  std::vector<size_t> counts;
  stk::mesh::count_entities(mesh.mesh_meta_data().universal_part(), mesh, counts);
  EXPECT_EQ(origCounts, counts);
  for (auto * bucketPtr : mesh.buckets(stk::topology::NODE_RANK))
    for (auto node : *bucketPtr)
      EXPECT_EQ(0., (node_coordinates(mesh, node) - origCoords[mesh.identifier(node)]).length());
}

void test_second_design_matches_fresh_background(const bool doSnapping)
{
  const double secondRadius = 0.6;

  BackgroundMesh background = box_mesh();
  {
    const ConformingMeshWithSensitivities firstDesign =
        decompose_geometry(background, sphere_at_origin(0.5), doSnapping, EDGE_DIRECTION);
  }
  const ConformingMeshWithSensitivities secondDesign =
      decompose_geometry(background, sphere_at_origin(secondRadius), doSnapping, EDGE_DIRECTION);

  BackgroundMesh freshBackground = box_mesh();
  const ConformingMeshWithSensitivities freshDesign = decompose_geometry(
      freshBackground, sphere_at_origin(secondRadius), doSnapping, EDGE_DIRECTION);

  expect_same_sensitivities(secondDesign, freshDesign);
  const double freshArea = global_area(freshDesign);
  EXPECT_NEAR(freshArea, global_area(secondDesign), 1.e-10 * freshArea);
}

TEST(PlatoKrinoInterface, secondDesignWithoutSnapping_matchesFreshBackground)
{
  test_second_design_matches_fresh_background(false);
}

TEST(PlatoKrinoInterface, secondDesignWithSnapping_matchesFreshBackground)
{
  test_second_design_matches_fresh_background(true);
}

TEST(PlatoKrinoInterface, meshReadFromFile_matchesBoundingBoxMesh)
{
  const std::string meshFile = "plato_interface_read.g";
  const krino::ScopedTestFile cleanup(meshFile);
  krino::generate_and_write_bounding_box_mesh(
      stk::topology::TETRAHEDRON_4, {-1., -1., -1.}, {1., 1., 1.}, 0.4, meshFile, true);

  // generate_and_write_bounding_box_mesh already writes pre-spread per-rank files, so read them
  // back with the identity decomposition ("").
  BackgroundMesh fileBackground = BackgroundMesh::read(meshFile, MPI_COMM_WORLD, "");
  const ConformingMeshWithSensitivities fileMesh =
      decompose_geometry(fileBackground, sphere_at_origin(0.5), false, EDGE_DIRECTION);

  BackgroundMesh boxBackground = box_mesh();
  const ConformingMeshWithSensitivities boxMesh =
      decompose_geometry(boxBackground, sphere_at_origin(0.5), false, EDGE_DIRECTION);

  expect_same_sensitivities(fileMesh, boxMesh);
}

// A consumer's own STK mesh, with a level set field of its own name.
TEST(PlatoKrinoInterface, separatelyGeneratedMeshAndField_matchesBoundingBoxMesh)
{
  const std::string meshFile = "plato_interface_consumer.g";
  const krino::ScopedTestFile cleanup(meshFile);
  krino::generate_and_write_bounding_box_mesh(
      stk::topology::TETRAHEDRON_4, {-1., -1., -1.}, {1., 1., 1.}, 0.4, meshFile, true);

  stk::io::StkMeshIoBroker ioBroker(MPI_COMM_WORLD);
  ioBroker.add_mesh_database(meshFile, stk::io::READ_MESH);
  ioBroker.create_input_mesh();
  stk::mesh::MetaData & meta = ioBroker.meta_data();
  auto & designField = meta.declare_field<double>(stk::topology::NODE_RANK, "designField");
  stk::mesh::put_field_on_mesh(designField, meta.universal_part(), nullptr);
  ioBroker.populate_bulk_data();

  initialize_field_from_geometry(ioBroker.bulk_data(), designField, sphere_at_origin(0.5));
  const ConformingMeshWithSensitivities consumerMesh(ioBroker.bulk_data(), {&designField});

  BackgroundMesh boxBackground = box_mesh();
  const ConformingMeshWithSensitivities boxMesh =
      decompose_geometry(boxBackground, sphere_at_origin(0.5), false, EDGE_DIRECTION);

  expect_same_sensitivities(consumerMesh, boxMesh);
}

// The sidesets are decomposed with the blocks: the sides of Xplus, all in the void of the
// half-space x < 0.1, are in Xplus_void, and those of Xminus, all in the material, are not in
// Xminus_void.
TEST(PlatoKrinoInterface, meshWithSidesets_sidesetsDecomposedIntoPhases)
{
  const std::string meshFile = "plato_interface_sidesets.g";
  const krino::ScopedTestFile cleanup(meshFile);
  BackgroundMesh background = read_box_mesh_with_sidesets(meshFile);
  const ConformingMeshWithSensitivities conformingMesh =
      decompose_geometry(background, half_space_below_x_equal_point_one(), false, EDGE_DIRECTION);

  const stk::mesh::MetaData & backgroundMeta = background.bulk_data().mesh_meta_data();
  const stk::mesh::MetaData & conformingMeta = conformingMesh.bulk_data().mesh_meta_data();
  const std::vector<const stk::mesh::Part *> surfaces =
      backgroundMeta.get_surfaces_in_surface_to_block_map();
  EXPECT_LT(0u, surfaces.size());
  for (auto * surface : surfaces)
    EXPECT_EQ(names_of_blocks_touching_surface(backgroundMeta, surface->name()),
        names_of_blocks_touching_surface(conformingMeta, surface->name()))
        << surface->name();

  const size_t numXplusSides = global_num_sides_in_part(background.bulk_data(), "Xplus");
  EXPECT_LT(0u, numXplusSides);
  EXPECT_EQ(numXplusSides, global_num_sides_in_part(conformingMesh.bulk_data(), "Xplus_void"));
  EXPECT_LT(0u, global_num_sides_in_part(conformingMesh.bulk_data(), "Xminus"));
  EXPECT_EQ(0u, global_num_sides_in_part(conformingMesh.bulk_data(), "Xminus_void"));
}

// The written material, x < 0.1, has the volume 1.1 * 2 * 2, and the level set is written as a
// transient field but the coordinates are not.
TEST(PlatoKrinoInterface,
    writeDecomposedMesh_readsBackMaterialAndLevelSetButNotCoordinatesAsTransient)
{
  const std::string meshFile = "plato_interface_write_input.g";
  const krino::ScopedTestFile cleanupInput(meshFile);
  BackgroundMesh background = read_box_mesh_with_sidesets(meshFile);
  ConformingMeshWithSensitivities conformingMesh =
      decompose_geometry(background, half_space_below_x_equal_point_one(), false, EDGE_DIRECTION);

  const std::string outputFile = "plato_interface_write_output.g";
  const krino::ScopedTestFile cleanupOutput(outputFile);
  conformingMesh.write_decomposed_mesh(outputFile, false);

  // The output is one composed file if parallel IO is enabled, and otherwise one file per
  // processor.
  stk::io::StkMeshIoBroker reader(MPI_COMM_WORLD);
  if (std::ifstream(outputFile).good())
    reader.property_add(Ioss::Property("DECOMPOSITION_METHOD", "rib"));
  reader.add_mesh_database(outputFile, stk::io::READ_MESH);
  reader.create_input_mesh();
  reader.populate_bulk_data();

  EXPECT_NEAR(4.4, global_volume(reader.bulk_data()), 1.e-12);

  const Ioss::NodeBlock * nodeBlock = reader.get_input_ioss_region()->get_node_blocks()[0];
  EXPECT_TRUE(nodeBlock->field_exists("LS"));
  EXPECT_FALSE(nodeBlock->field_exists("coordinates"));
}

TEST(PlatoKrinoInterface, levelSetFieldOfAnotherMesh_initializeAndDecomposeThrow)
{
  BackgroundMesh background = box_mesh();
  BackgroundMesh otherBackground = box_mesh();
  stk::mesh::FieldBase & otherField = otherBackground.level_set_field();
  EXPECT_ANY_THROW(
      initialize_field_from_geometry(background.bulk_data(), otherField, sphere_at_origin(0.5)));
  EXPECT_ANY_THROW(initialize_field_from_function(background.bulk_data(),
      otherField,
      [](const std::array<double, 3> & x) { return 0.1 - x[0]; }));
  EXPECT_ANY_THROW(set_field_values(background.bulk_data(), otherField, {}, {}));

  initialize_field_from_geometry(otherBackground.bulk_data(), otherField, sphere_at_origin(0.5));
  EXPECT_ANY_THROW((ConformingMeshWithSensitivities(background.bulk_data(), {&otherField})));
}

// A consumer's mesh read from a box mesh file with sidesets, with two fields that can't be level
// sets: the nodal double field "sidesetField" on the sideset Xminus only, and the nodal int field
// "intField" on every node.
std::unique_ptr<stk::io::StkMeshIoBroker>
read_box_mesh_with_invalid_level_set_fields(const std::string & meshFile)
{
  krino::generate_and_write_bounding_box_mesh(
      stk::topology::TETRAHEDRON_4, {-1., -1., -1.}, {1., 1., 1.}, 0.4, meshFile, true);
  auto ioBroker = std::make_unique<stk::io::StkMeshIoBroker>(MPI_COMM_WORLD);
  ioBroker->add_mesh_database(meshFile, stk::io::READ_MESH);
  ioBroker->create_input_mesh();
  stk::mesh::MetaData & meta = ioBroker->meta_data();
  auto & sidesetField = meta.declare_field<double>(stk::topology::NODE_RANK, "sidesetField");
  stk::mesh::put_field_on_mesh(sidesetField, *meta.get_part("Xminus"), nullptr);
  auto & intField = meta.declare_field<int>(stk::topology::NODE_RANK, "intField");
  stk::mesh::put_field_on_mesh(intField, meta.universal_part(), nullptr);
  ioBroker->populate_bulk_data();
  return ioBroker;
}

TEST(PlatoKrinoInterface, levelSetFieldOnSidesetOrInt_initializeAndDecomposeThrow)
{
  const std::string meshFile = "plato_interface_invalid_fields.g";
  const krino::ScopedTestFile cleanup(meshFile);
  const std::unique_ptr<stk::io::StkMeshIoBroker> ioBroker =
      read_box_mesh_with_invalid_level_set_fields(meshFile);
  const stk::mesh::MetaData & meta = ioBroker->meta_data();

  for (const std::string fieldName : {"sidesetField", "intField"})
  {
    stk::mesh::FieldBase & field = *meta.get_field(stk::topology::NODE_RANK, fieldName);
    EXPECT_ANY_THROW(initialize_field_from_function(
        ioBroker->bulk_data(), field, [](const std::array<double, 3> & x) { return 0.1 - x[0]; }));
    EXPECT_ANY_THROW((ConformingMeshWithSensitivities(ioBroker->bulk_data(), {&field})));
  }
}

// krino's planes unitize their normals, but the half-space is normal.x + offset < 0 for any normal.
TEST(PlatoKrinoInterface, nonUnitPlaneNormal_sameInterfaceAsUnitNormal)
{
  Geometry geometry;
  geometry.planes.push_back(PlaneSpec{{{2., 0., 0.}}, -0.2});
  BackgroundMesh background = box_mesh();
  const ConformingMeshWithSensitivities conformingMesh =
      decompose_geometry(background, geometry, false, EDGE_DIRECTION);

  const krino::InterfaceFacetsWithSensitivities<3> facets = conformingMesh.interface_facets<3>();
  for (size_t facet = 0; facet < facets.num_facets(); ++facet)
    for (const auto & x : facet_coordinates(conformingMesh.bulk_data(), facets, facet))
      EXPECT_NEAR(0.1, x[0], 1.e-12);
  EXPECT_NEAR(4., global_area(conformingMesh), 1.e-12);
}

TEST(PlatoKrinoInterface, emptyGeometryZeroPlaneNormalOrNonpositiveRadius_throws)
{
  BackgroundMesh background = box_mesh();
  Geometry zeroNormal;
  zeroNormal.planes.push_back(PlaneSpec{{{0., 0., 0.}}, 0.});
  EXPECT_ANY_THROW(initialize_field_from_geometry(
      background.bulk_data(), background.level_set_field(), Geometry()));
  EXPECT_ANY_THROW(initialize_field_from_geometry(
      background.bulk_data(), background.level_set_field(), zeroNormal));
  EXPECT_ANY_THROW(initialize_field_from_geometry(
      background.bulk_data(), background.level_set_field(), sphere_at_origin(0.)));
}

// Invalid arguments on processor 0 alone throw on every processor, rather than leaving the others
// waiting in the communication of the field values.
TEST(PlatoKrinoInterface, setFieldValuesInvalidOnlyOnOneProc_throwsOnEveryProc)
{
  BackgroundMesh background = box_mesh();
  const bool isProc0 = background.bulk_data().parallel_rank() == 0;

  const std::vector<stk::mesh::EntityId> missingNodeIds =
      isProc0 ? std::vector<stk::mesh::EntityId>{999999999} : std::vector<stk::mesh::EntityId>{};
  const std::vector<double> values = isProc0 ? std::vector<double>{1.} : std::vector<double>{};
  EXPECT_ANY_THROW(set_field_values(
      background.bulk_data(), background.level_set_field(), missingNodeIds, values));
  EXPECT_ANY_THROW(
      set_field_values(background.bulk_data(), background.level_set_field(), {}, values));
}

// The material is x < 0.1 and y < 0.1.  Each level set's facets and sensitivities are on the
// material boundary, x = 0.1 or y = 0.1 with the other coordinate below 0.1 and area 1.1 * 2, and
// not where it separates two void phases.
TEST(PlatoKrinoInterface, twoLevelSets_facetsAndSensitivitiesOfEachLevelSetOnMaterialBoundary)
{
  BackgroundMesh background = BackgroundMesh::bounding_box(
      3, {{-1., -1., -1.}}, {{1., 1., 1.}}, 0.4, MPI_COMM_WORLD, {.numLevelSets = 2});
  for (unsigned ls = 0; ls < 2; ++ls)
    initialize_field_from_function(background.bulk_data(),
        background.level_set_field(ls),
        [ls](const std::array<double, 3> & x) { return 0.1 - x[ls]; });
  const ConformingMeshWithSensitivities conformingMesh(
      background.bulk_data(), {&background.level_set_field(0), &background.level_set_field(1)});
  EXPECT_EQ(2u, conformingMesh.num_level_sets());
  const stk::mesh::BulkData & mesh = conformingMesh.bulk_data();
  const stk::mesh::Selector material = conformingMesh.material_selector();

  for (unsigned ls = 0; ls < 2; ++ls)
  {
    const unsigned otherLs = 1 - ls;
    const krino::InterfaceFacetsWithSensitivities<3> facets =
        conformingMesh.interface_facets<3>(ls);
    for (size_t facet = 0; facet < facets.num_facets(); ++facet)
    {
      for (const auto & x : facet_coordinates(mesh, facets, facet))
      {
        EXPECT_NEAR(0.1, x[ls], 1.e-12);
        EXPECT_LE(x[otherLs], 0.1 + 1.e-12);
      }
    }
    EXPECT_NEAR(2.2, global_area(conformingMesh, ls), 1.e-12);

    for (const auto & sens : conformingMesh.shape_sensitivities(ls))
      EXPECT_TRUE(
          material(mesh.bucket(mesh.get_entity(stk::topology::NODE_RANK, sens.interfaceNodeId))))
          << "node " << sens.interfaceNodeId;
  }
}

// The material is the half-space x < 0.1 of the box [-1,1]^3: its volume is 1.1 * 2 * 2 and the
// interface's area 4.
TEST(PlatoKrinoInterface, materialVoidAndInterfaceSelectors_halfSpaceVolumesAndInterfaceArea)
{
  BackgroundMesh background = box_mesh();
  const ConformingMeshWithSensitivities conformingMesh =
      decompose_geometry(background, half_space_below_x_equal_point_one(), false, EDGE_DIRECTION);
  const stk::mesh::BulkData & mesh = conformingMesh.bulk_data();

  EXPECT_NEAR(4.4, global_volume(mesh, conformingMesh.material_selector()), 1.e-12);
  EXPECT_NEAR(3.6, global_volume(mesh, conformingMesh.void_selector()), 1.e-12);
  EXPECT_EQ(0u,
      global_num_elements(
          mesh, conformingMesh.material_selector() & conformingMesh.void_selector()));
  EXPECT_NEAR(4., global_side_area(mesh, conformingMesh.interface_selector()), 1.e-12);
}

// The material is x < 0.1 and y < 0.1, with volume 1.1 * 1.1 * 2.  The interface sides of each
// level set are on the material boundary, with area 1.1 * 2, as are its facets.
TEST(PlatoKrinoInterface, twoLevelSetSelectors_materialIsBothHalfSpacesAndInterfacesBoundIt)
{
  BackgroundMesh background = BackgroundMesh::bounding_box(
      3, {{-1., -1., -1.}}, {{1., 1., 1.}}, 0.4, MPI_COMM_WORLD, {.numLevelSets = 2});
  for (unsigned ls = 0; ls < 2; ++ls)
    initialize_field_from_function(background.bulk_data(),
        background.level_set_field(ls),
        [ls](const std::array<double, 3> & x) { return 0.1 - x[ls]; });
  const ConformingMeshWithSensitivities conformingMesh(
      background.bulk_data(), {&background.level_set_field(0), &background.level_set_field(1)});
  const stk::mesh::BulkData & mesh = conformingMesh.bulk_data();

  EXPECT_NEAR(2.42, global_volume(mesh, conformingMesh.material_selector()), 1.e-12);
  EXPECT_NEAR(8. - 2.42, global_volume(mesh, conformingMesh.void_selector()), 1.e-12);
  EXPECT_EQ(0u,
      global_num_elements(
          mesh, conformingMesh.material_selector() & conformingMesh.void_selector()));
  for (unsigned ls = 0; ls < 2; ++ls)
    EXPECT_NEAR(2.2, global_side_area(mesh, conformingMesh.interface_selector(ls)), 1.e-12);
}

// The bounding box's walls are flat, and the sides of each are in a sideset touching the block.
TEST(PlatoKrinoInterface, boundingBox_sidesOfEachFlatWallInSidesetTouchingBlock)
{
  BackgroundMesh background = box_mesh();
  const stk::mesh::BulkData & mesh = background.bulk_data();
  const stk::mesh::MetaData & meta = mesh.mesh_meta_data();
  for (const std::string wallName : {"Xminus", "Xplus", "Yminus", "Yplus", "Zminus", "Zplus"})
  {
    const stk::mesh::Part * wall = meta.get_part(wallName);
    ASSERT_NE(nullptr, wall) << wallName;
    EXPECT_NEAR(4., global_side_area(mesh, *wall), 1.e-12) << wallName;
    EXPECT_EQ(std::set<std::string>{"block_1"}, names_of_blocks_touching_surface(meta, wallName))
        << wallName;
  }
}

BackgroundMesh box_mesh_with_two_level_sets()
{
  return BackgroundMesh::bounding_box(
      3, {{-1., -1., -1.}}, {{1., 1., 1.}}, 0.4, MPI_COMM_WORLD, {.numLevelSets = 2});
}

// The level set from a closed box's STL is the exact signed distance to the box.
TEST(PlatoKrinoInterface, facetFileBox_levelSetIsSignedDistanceToBox)
{
  const std::string stlFile = "plato_interface_box.stl";
  const krino::ScopedTestFile cleanup(stlFile);
  write_box_stl(stlFile, 0.5);

  BackgroundMesh background = box_mesh_with_two_level_sets();
  initialize_field_from_geometry(
      background.bulk_data(), background.level_set_field(0), Geometry{.facetFiles = {{stlFile}}});
  initialize_field_from_function(background.bulk_data(),
      background.level_set_field(1),
      [](const std::array<double, 3> & x) { return box_level_set(x, 0.5); });
  expect_same_level_sets(background);
}

// A scaled box united with a sphere is the larger of their level sets, and the outside of the box
// is the negated box level set.
TEST(PlatoKrinoInterface, facetFileUnitedWithSphereAndScaledAndOutside_matchesFunction)
{
  const std::string stlFile = "plato_interface_small_box.stl";
  const krino::ScopedTestFile cleanup(stlFile);
  write_box_stl(stlFile, 0.25);
  const FacetFileSpec scaledBox{.fileName = stlFile, .scale = {{2., 2., 2.}}};
  const stk::math::Vector3d sphereCenter(0.7, 0., 0.);
  const double sphereRadius = 0.3;

  BackgroundMesh background = box_mesh_with_two_level_sets();
  Geometry boxAndSphere{.facetFiles = {scaledBox}};
  boxAndSphere.spheres.push_back(SphereSpec{{{0.7, 0., 0.}}, sphereRadius});
  initialize_field_from_geometry(
      background.bulk_data(), background.level_set_field(0), boxAndSphere);
  initialize_field_from_function(background.bulk_data(),
      background.level_set_field(1),
      [&](const std::array<double, 3> & x)
      {
        const double sphere =
            sphereRadius - (stk::math::Vector3d(x[0], x[1], x[2]) - sphereCenter).length();
        return std::max(box_level_set(x, 0.5), sphere);
      });
  expect_same_level_sets(background);

  FacetFileSpec outsideOfBox = scaledBox;
  outsideOfBox.materialInside = false;
  initialize_field_from_geometry(background.bulk_data(),
      background.level_set_field(0),
      Geometry{.facetFiles = {outsideOfBox}});
  initialize_field_from_function(background.bulk_data(),
      background.level_set_field(1),
      [](const std::array<double, 3> & x) { return -box_level_set(x, 0.5); });
  expect_same_level_sets(background);
}

TEST(PlatoKrinoInterface, facetFileUnsupportedExtensionMissingFileOr2D_throws)
{
  BackgroundMesh background = box_mesh();
  EXPECT_ANY_THROW(initialize_field_from_geometry(background.bulk_data(),
      background.level_set_field(),
      Geometry{.facetFiles = {{"plato_interface_box.xyz"}}}));
  EXPECT_ANY_THROW(initialize_field_from_geometry(background.bulk_data(),
      background.level_set_field(),
      Geometry{.facetFiles = {{"plato_interface_no_such_file.stl"}}}));

  BackgroundMesh background2D =
      BackgroundMesh::bounding_box(2, {{-1., -1., 0.}}, {{1., 1., 0.}}, 0.4, MPI_COMM_WORLD);
  EXPECT_ANY_THROW(initialize_field_from_geometry(background2D.bulk_data(),
      background2D.level_set_field(),
      Geometry{.facetFiles = {{"plato_interface_box.stl"}}}));
}
}
