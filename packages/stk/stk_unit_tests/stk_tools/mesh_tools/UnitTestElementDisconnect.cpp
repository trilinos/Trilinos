#include <gtest/gtest.h>

#include <algorithm>
#include <sstream>
#include <stk_io/FillMesh.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_math/StkVector.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_unit_test_utils/BuildMesh.hpp>
#include <stk_unit_test_utils/ConstructedMesh.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <string>
#include <type_traits>

#include "stk_mesh/base/Bucket.hpp"
#include "stk_mesh/base/CreateFaces.hpp"
#include "stk_mesh/base/Entity.hpp"
#include "stk_mesh/base/FieldBase.hpp"
#include "stk_mesh/base/Selector.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_tools/mesh_tools/DisjointSet.hpp"
#include "stk_tools/mesh_tools/EntityDisconnectTool.hpp"
#include "stk_topology/topology.hpp"
#include "stk_util/parallel/GenerateParallelConsistentIDs.hpp"
#include "stk_util/parallel/Parallel.hpp"
#include "stk_util/parallel/ParallelReduce.hpp"
#include "stk_util/util/ReportHandler.hpp"
#include "stk_util/util/SortAndUnique.hpp"

namespace stk::experimental
{

class TestElementDisconnect : public stk::unit_test_util::MeshFixture
{
 public:
  TestElementDisconnect() : stk::unit_test_util::MeshFixture(3)
  {
    setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  }

  void create_mesh_with_faces(unsigned nx_, unsigned ny_, unsigned nz_)
  {
    nx = nx_;
    ny = ny_;
    nz = nz_;
    stk::io::fill_mesh(get_generated_mesh_description(), get_bulk());
    stk::mesh::create_faces(get_bulk());
  }

  template <typename FacePred>
  auto get_faces(FacePred predicate)
  {
    const auto& bulk = get_bulk();
    auto predBound = [&bulk, predicate](const stk::mesh::Entity& face) { return predicate(bulk, face); };
    auto facesFiltered = stk::mesh::get_entities(bulk, stk::topology::FACE_RANK);
    auto endFiltered = std::stable_partition(facesFiltered.begin(), facesFiltered.end(), predBound);
    facesFiltered.erase(endFiltered, facesFiltered.end());
    return EntityContainer(facesFiltered.begin(), facesFiltered.end());
  }

  static auto is_exterior(const stk::mesh::BulkData& bulk, const stk::mesh::Entity& face)
  {
    return bulk.get_connected_entities(face, stk::topology::ELEMENT_RANK).size() == 1U;
  }

  static auto is_interior(const stk::mesh::BulkData& bulk, const stk::mesh::Entity& face)
  {
    return bulk.get_connected_entities(face, stk::topology::ELEMENT_RANK).size() == 2U;
  }

  auto get_interior_faces() { return get_faces(is_interior); }
  auto get_exterior_faces() { return get_faces(is_exterior); }

  auto get_elements()
  {
    auto& bulk = get_bulk();
    auto& meta = get_meta();
    std::vector<stk::mesh::Entity> elems;
    bulk.get_entities(stk::topology::ELEMENT_RANK, meta.universal_part(), elems);
    return elems;
  }

  std::string get_generated_mesh_description() const
  {
    return "generated:" + std::to_string(nx) + "x" + std::to_string(ny) + "x" + std::to_string(nz);
  }

  auto original_num_exterior_faces() const { return 2 * (nx * ny + nx * nz + ny * nz); }
  auto original_num_interior_faces() const { return 3 * nx * ny * nz - (nx * ny + nx * nz + ny * nz); }
  auto original_num_faces() const { return original_num_exterior_faces() + original_num_interior_faces(); }
  auto original_num_nodes() const { return (nx + 1) * (ny + 1) * (nz + 1); }

  static auto compute_face_centroid(const stk::mesh::BulkData& bulk, const stk::mesh::Entity& entity)
  {
    const auto& coordField = *bulk.mesh_meta_data().coordinate_field();
    auto faceNodes = bulk.get_connected_entities(entity, stk::topology::NODE_RANK);
    stk::math::Vec<double, 3> centroid;
    const auto coordData = coordField.data<double, stk::mesh::ReadWrite>();
    for (auto n = 0U; n < faceNodes.size(); ++n) {
      const auto& node = faceNodes[n];
      const auto nodeCoords = coordData.entity_values(node);
      centroid[0] += nodeCoords(0_comp);
      centroid[1] += nodeCoords(1_comp);
      centroid[2] += nodeCoords(2_comp);
    }
    return 0.25 * centroid;
  }

  static auto num_entities(const stk::mesh::BulkData& bulk, const stk::topology::rank_t& rank)
  {
    const auto& meta = bulk.mesh_meta_data();
    stk::mesh::Selector owned = meta.universal_part() & meta.locally_owned_part();
    std::size_t localCount = 0;
    const auto buckets = bulk.get_buckets(rank, owned);
    for (const auto* bucket : buckets) {
      localCount += bucket->size();
    }
    return stk::get_global_sum(bulk.parallel(), localCount);
  }

  void check_face_pairs_are_different(const stk::experimental::EntityDisconnectTool& disconnecter)
  {
    // Old faces have identifiers ending in [7, 8, 9, 0]
    // New faces have identifiers ending in [1, 2, 3, 4, 5, 6] (faceId = 10 * elemId + sideIdx + 1)
    for (const auto& [faceOrig, faceNew] : disconnecter.get_elem_side_pairs()) {
      EXPECT_GE((get_bulk().identifier(faceOrig.face) - 1) % 10, 6U);
      EXPECT_LT((get_bulk().identifier(faceNew.face) - 1) % 10, 6U);
    }
  }

  using RankT = stk::topology::rank_t;
  static constexpr auto NodeRank = RankT::NODE_RANK;
  static constexpr auto FaceRank = RankT::FACE_RANK;
  static constexpr auto ElemRank = RankT::ELEM_RANK;

  unsigned nx = 2;
  unsigned ny = 2;
  unsigned nz = 2;
};

class ElementDisconnectAdjacencyTests : public TestElementDisconnect
{
 public:
  static auto face_is_contained_in_list(const stk::mesh::Entity& face, const EntityContainer& dcFaces)
  {
    return dcFaces.find(face) != dcFaces.end();
  }

  auto is_entity_node_adjacent_to_list(const stk::mesh::Entity& face, const EntityContainer& dcFaces)
  {
    auto adjFaces = get_adjacent_entities<NodeRank, FaceRank>(get_bulk(), {face},
        [&dcFaces](const stk::mesh::Entity& nodeAdjFace) { return face_is_contained_in_list(nodeAdjFace, dcFaces); });
    return adjFaces.size() != 0U;
  }

  static auto face_is_element_connected(const stk::mesh::BulkData& mesh, const stk::mesh::Entity& face)
  {
    return mesh.num_connectivity(face, stk::topology::ELEMENT_RANK) >= 1;
  }

  void check_faces_are_node_adjacent(
      const EntityContainer& dcFaces, const stk::experimental::EntityContainer& adjFaces)
  {
    for (const auto& face : adjFaces) {
      EXPECT_FALSE(face_is_contained_in_list(face, dcFaces))
          << "Node-adjacent face " << get_bulk().identifier(face) << " appears in the disconnect list!";
      EXPECT_TRUE(is_entity_node_adjacent_to_list(face, dcFaces))
          << "Face " << get_bulk().identifier(face)
          << " is not actually node-adjacent to any faces in disconnect list!";
    }
  }

  void check_adj_faces_are_element_connected(const EntityContainer& adjFaces)
  {
    // Not really checking anything, so mostly for completeness.  Unless something very wrong happens in the unit test
    // fixture, faces from `create_faces` should be guaranteed to be connnected to at least one element.
    for (const auto& face : adjFaces) {
      EXPECT_TRUE(face_is_element_connected(get_bulk(), face))
          << "Node-adjacent face " << get_bulk().identifier(face) << " is not connected to any elements!";
    }
  }

  void check_elements_are_adjacent(const EntityContainer& dcFaces, const EntityContainer& adjElems)
  {
    for (const auto& elem : adjElems) {
      EXPECT_TRUE((is_entity_node_adjacent_to_list(elem, dcFaces)))
          << "Elem " << get_bulk().identifier(elem) << " is not actually adjacent to any faces in the disconnect list!";
    }
  }
};

TEST_F(ElementDisconnectAdjacencyTests, emptyDisconnectList_gen2x2x2)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 1) GTEST_SKIP() << "Test only available in serial!";
  create_mesh_with_faces(2, 2, 2);

  EntityDisconnectTool disconnecter(get_bulk(), EntityContainer{});
  const auto& adjacentFaces = disconnecter.get_retained_faces();
  EXPECT_EQ(adjacentFaces.size(), 0U);
  EXPECT_EQ(disconnecter.get_adjacent_elements().size(), 0U);
}

TEST_F(ElementDisconnectAdjacencyTests, exteriorDisconnectList_gen2x2x2)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 1) GTEST_SKIP() << "Test only available in serial!";
  create_mesh_with_faces(2, 2, 2);

  auto disconnectFaces = get_exterior_faces();
  auto interiorFaces = get_interior_faces();

  EntityDisconnectTool disconnecter(get_bulk(), disconnectFaces);
  const auto& adjacentFaces = disconnecter.get_retained_faces();
  EXPECT_EQ(adjacentFaces.size(), 0U);
  check_faces_are_node_adjacent(disconnectFaces, adjacentFaces);
  check_elements_are_adjacent(disconnectFaces, disconnecter.get_adjacent_elements());
  check_adj_faces_are_element_connected(adjacentFaces);
}

TEST_F(ElementDisconnectAdjacencyTests, interiorDisconnectList_gen2x2x2)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 1) GTEST_SKIP() << "Test only available in serial!";
  create_mesh_with_faces(2, 2, 2);

  auto exteriorFaces = get_exterior_faces();
  auto disconnectFaces = get_interior_faces();

  EntityDisconnectTool disconnecter(get_bulk(), disconnectFaces);
  const auto& adjacentFaces = disconnecter.get_retained_faces();
  EXPECT_EQ(adjacentFaces.size(), exteriorFaces.size());
  check_faces_are_node_adjacent(disconnectFaces, adjacentFaces);
  check_elements_are_adjacent(disconnectFaces, disconnecter.get_adjacent_elements());
  check_adj_faces_are_element_connected(adjacentFaces);
}

TEST_F(ElementDisconnectAdjacencyTests, disconnectXPlane_gen2x2x2)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 1) GTEST_SKIP() << "Test only available in serial!";
  create_mesh_with_faces(2, 2, 2);

  auto on_x_plane = [](const stk::mesh::BulkData& bulk, const stk::mesh::Entity& entity) {
    auto centroid = compute_face_centroid(bulk, entity);
    constexpr double x_plane = 1.0;
    return (std::abs(centroid[0] - x_plane) <= 1.0e-10);
  };
  auto disconnectFaces = get_faces(on_x_plane);

  EntityDisconnectTool disconnecter(get_bulk(), disconnectFaces);
  const auto& adjacentFaces = disconnecter.get_retained_faces();
  EXPECT_EQ(adjacentFaces.size(), 24U);
  check_faces_are_node_adjacent(disconnectFaces, adjacentFaces);
  check_elements_are_adjacent(disconnectFaces, disconnecter.get_adjacent_elements());
  check_adj_faces_are_element_connected(adjacentFaces);
}

TEST_F(ElementDisconnectAdjacencyTests, disconnectXPlane_gen4x2x2)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 1) GTEST_SKIP() << "Test only available in serial!";
  create_mesh_with_faces(4, 2, 2);

  auto on_x_plane = [](const stk::mesh::BulkData& bulk, const stk::mesh::Entity& entity) {
    auto centroid = compute_face_centroid(bulk, entity);
    constexpr double x_plane = 2.0;
    return (std::abs(centroid[0] - x_plane) <= 1.0e-10);
  };
  auto disconnectFaces = get_faces(on_x_plane);

  EntityDisconnectTool disconnecter(get_bulk(), disconnectFaces);
  const auto& adjacentFaces = disconnecter.get_retained_faces();
  EXPECT_EQ(adjacentFaces.size(), 24U);
  check_faces_are_node_adjacent(disconnectFaces, adjacentFaces);
  check_elements_are_adjacent(disconnectFaces, disconnecter.get_adjacent_elements());
  check_adj_faces_are_element_connected(adjacentFaces);
}

TEST_F(ElementDisconnectAdjacencyTests, disconnectYPlane_gen2x2x2)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 1) GTEST_SKIP() << "Test only available in serial!";
  create_mesh_with_faces(2, 2, 2);

  auto on_y_plane = [](const stk::mesh::BulkData& bulk, const stk::mesh::Entity& entity) {
    auto centroid = compute_face_centroid(bulk, entity);
    constexpr double y_plane = 1.0;
    return (std::abs(centroid[1] - y_plane) <= 1.0e-10);
  };
  auto disconnectFaces = get_faces(on_y_plane);

  EntityDisconnectTool disconnecter(get_bulk(), disconnectFaces);
  const auto& adjacentFaces = disconnecter.get_retained_faces();
  EXPECT_EQ(adjacentFaces.size(), 24U);
  check_faces_are_node_adjacent(disconnectFaces, adjacentFaces);
  check_elements_are_adjacent(disconnectFaces, disconnecter.get_adjacent_elements());
  check_adj_faces_are_element_connected(adjacentFaces);
}

TEST_F(ElementDisconnectAdjacencyTests, disconnectYPlane_gen4x2x2)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 1) GTEST_SKIP() << "Test only available in serial!";
  create_mesh_with_faces(4, 2, 2);

  auto on_y_plane = [](const stk::mesh::BulkData& bulk, const stk::mesh::Entity& entity) {
    auto centroid = compute_face_centroid(bulk, entity);
    constexpr double y_plane = 1.0;
    return (std::abs(centroid[1] - y_plane) <= 1.0e-10);
  };
  auto disconnectFaces = get_faces(on_y_plane);

  EntityDisconnectTool disconnecter(get_bulk(), disconnectFaces);
  const auto& adjacentFaces = disconnecter.get_retained_faces();
  EXPECT_EQ(adjacentFaces.size(), 44U);
  check_faces_are_node_adjacent(disconnectFaces, adjacentFaces);
  check_elements_are_adjacent(disconnectFaces, disconnecter.get_adjacent_elements());
  check_adj_faces_are_element_connected(adjacentFaces);
}

TEST_F(TestElementDisconnect, disconnect_empty_list_noops)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 1) GTEST_SKIP() << "Test only available in serial!";
  create_mesh_with_faces(2, 2, 2);

  stk::mesh::EntityVector disconnectFaces;
  EntityDisconnectTool disconnecter(get_bulk(), disconnectFaces);

  auto newNodeIds = disconnecter.determine_new_nodes();
  EXPECT_EQ(newNodeIds.size(), 0U);
}

TEST_F(TestElementDisconnect, disconnect_exterior_faces_noops)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 1) GTEST_SKIP() << "Test only available in serial!";
  create_mesh_with_faces(2, 2, 2);

  auto disconnectFaces = get_faces(is_exterior);
  ASSERT_EQ(disconnectFaces.size(), original_num_exterior_faces());

  EntityDisconnectTool disconnecter(get_bulk(), disconnectFaces);
  auto newNodeIds = disconnecter.determine_new_nodes();
  EXPECT_EQ(newNodeIds.size(), 0U);
}

TEST_F(TestElementDisconnect, disconnect_interior_faces_gen2x1x1)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 1) GTEST_SKIP() << "Test only available in serial!";
  create_mesh_with_faces(2, 1, 1);

  auto disconnectFaces = get_faces(is_interior);
  ASSERT_EQ(disconnectFaces.size(), original_num_interior_faces());

  EntityDisconnectTool disconnecter(get_bulk(), disconnectFaces);
  auto newNodeIds = disconnecter.determine_new_nodes();
  EXPECT_EQ(newNodeIds.size(), 4U);

  disconnecter.modify_mesh();
  auto numNewNodes = num_entities(get_bulk(), stk::topology::NODE_RANK) - original_num_nodes();
  auto numNewFaces = num_entities(get_bulk(), stk::topology::FACE_RANK) - original_num_faces();
  EXPECT_EQ(numNewNodes, 4U);
  EXPECT_EQ(numNewFaces, 1U);

  check_face_pairs_are_different(disconnecter);
}

TEST_F(TestElementDisconnect, disconnect_interior_faces_gen3x1x1)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 1) GTEST_SKIP() << "Test only available in serial!";
  create_mesh_with_faces(3, 1, 1);

  auto disconnectFaces = get_faces(is_interior);
  ASSERT_EQ(disconnectFaces.size(), original_num_interior_faces());

  EntityDisconnectTool disconnecter(get_bulk(), disconnectFaces);
  auto newNodeIds = disconnecter.determine_new_nodes();
  EXPECT_EQ(newNodeIds.size(), 8U);

  disconnecter.modify_mesh();
  auto numNewNodes = num_entities(get_bulk(), stk::topology::NODE_RANK) - original_num_nodes();
  auto numNewFaces = num_entities(get_bulk(), stk::topology::FACE_RANK) - original_num_faces();
  EXPECT_EQ(numNewNodes, 8U);
  EXPECT_EQ(numNewFaces, 2U);

  check_face_pairs_are_different(disconnecter);
}

TEST_F(TestElementDisconnect, disconnect_all_faces_gen2x1x1)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 1) GTEST_SKIP() << "Test only available in serial!";
  auto& bulk = get_bulk();
  create_mesh_with_faces(2, 1, 1);

  auto disconnectFaces = stk::mesh::get_entities(bulk, stk::topology::FACE_RANK);
  EntityDisconnectTool disconnecter(get_bulk(), disconnectFaces);
  auto newNodeIds = disconnecter.determine_new_nodes();
  EXPECT_EQ(newNodeIds.size(), 4U);

  disconnecter.modify_mesh();
  auto numNewNodes = num_entities(get_bulk(), stk::topology::NODE_RANK) - original_num_nodes();
  auto numNewFaces = num_entities(get_bulk(), stk::topology::FACE_RANK) - original_num_faces();
  EXPECT_EQ(numNewNodes, 4U);
  EXPECT_EQ(numNewFaces, 1U);

  check_face_pairs_are_different(disconnecter);
}

TEST_F(TestElementDisconnect, disconnect_all_faces_gen3x1x1)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 1) GTEST_SKIP() << "Test only available in serial!";
  auto& bulk = get_bulk();
  create_mesh_with_faces(3, 1, 1);

  auto disconnectFaces = stk::mesh::get_entities(bulk, stk::topology::FACE_RANK);
  EntityDisconnectTool disconnecter(bulk, disconnectFaces);
  auto newNodeIds = disconnecter.determine_new_nodes();
  EXPECT_EQ(newNodeIds.size(), 8U);

  disconnecter.modify_mesh();
  auto numNewNodes = num_entities(get_bulk(), stk::topology::NODE_RANK) - original_num_nodes();
  auto numNewFaces = num_entities(get_bulk(), stk::topology::FACE_RANK) - original_num_faces();
  EXPECT_EQ(numNewNodes, 8U);
  EXPECT_EQ(numNewFaces, 2U);

  check_face_pairs_are_different(disconnecter);
}

TEST_F(TestElementDisconnect, disconnect_plane_x_eq_2_gen3x2x2)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 1) GTEST_SKIP() << "Test only available in serial!";
  auto& bulk = get_bulk();
  create_mesh_with_faces(3, 2, 2);

  constexpr double x_plane = 2.0;
  auto on_x_plane = [](const stk::mesh::BulkData& bulk_, const stk::mesh::Entity& entity) {
    auto centroid = compute_face_centroid(bulk_, entity);
    return (std::abs(centroid[0] - x_plane) <= 1.0e-10);
  };

  auto disconnectFaces = get_faces(on_x_plane);
  EntityDisconnectTool disconnecter(bulk, disconnectFaces);
  auto newNodeIds = disconnecter.determine_new_nodes();
  EXPECT_EQ(newNodeIds.size(), 9U);

  disconnecter.modify_mesh();
  auto numNewNodes = num_entities(get_bulk(), stk::topology::NODE_RANK) - original_num_nodes();
  auto numNewFaces = num_entities(get_bulk(), stk::topology::FACE_RANK) - original_num_faces();
  EXPECT_EQ(numNewNodes, 9U);
  EXPECT_EQ(numNewFaces, 4U);

  check_face_pairs_are_different(disconnecter);
}

TEST_F(TestElementDisconnect, disconnect_corner_elem_gen2x2x2)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 1) GTEST_SKIP() << "Test only available in serial!";
  auto& bulk = get_bulk();
  auto& meta = get_meta();
  create_mesh_with_faces(2, 2, 2);
  std::vector<stk::mesh::Entity> elems;
  bulk.get_entities(stk::topology::ELEMENT_RANK, meta.universal_part(), elems);

  std::vector<stk::mesh::Entity> disconnectFaces;
  auto elemFaces = bulk.get_connected_entities(elems[0], stk::topology::FACE_RANK);
  for (auto f = 0U; f < elemFaces.size(); ++f) disconnectFaces.push_back(elemFaces[f]);

  EntityDisconnectTool disconnecter(bulk, disconnectFaces);
  auto newNodeIds = disconnecter.determine_new_nodes();
  EXPECT_EQ(newNodeIds.size(), 7U);

  disconnecter.modify_mesh();
  auto numNewNodes = num_entities(get_bulk(), stk::topology::NODE_RANK) - original_num_nodes();
  auto numNewFaces = num_entities(get_bulk(), stk::topology::FACE_RANK) - original_num_faces();
  EXPECT_EQ(numNewNodes, 7U);
  EXPECT_EQ(numNewFaces, 3U);

  check_face_pairs_are_different(disconnecter);
}

TEST_F(TestElementDisconnect, disconnect_edge_elems_gen2x2x2)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 1) GTEST_SKIP() << "Test only available in serial!";
  auto& bulk = get_bulk();
  create_mesh_with_faces(2, 2, 2);

  constexpr double x_plane = 1.0;
  constexpr double y_plane = 1.0;
  auto is_lower_left_face = [](const stk::mesh::BulkData& bulk_, const stk::mesh::Entity& entity) {
    auto centroid = compute_face_centroid(bulk_, entity);
    auto is_lower_xplane = (std::abs(centroid[0] - x_plane) <= 1.0e-10 && centroid[1] < y_plane);
    auto is_lower_yplane = (std::abs(centroid[1] - y_plane) <= 1.0e-10 && centroid[0] < x_plane);
    return is_lower_xplane || is_lower_yplane;
  };

  auto disconnectFaces = get_faces(is_lower_left_face);
  ASSERT_EQ(disconnectFaces.size(), 4U);
  EntityDisconnectTool disconnecter(bulk, disconnectFaces);
  auto newNodeIds = disconnecter.determine_new_nodes();
  EXPECT_EQ(newNodeIds.size(), 9U);

  disconnecter.modify_mesh();
  auto numNewNodes = num_entities(get_bulk(), stk::topology::NODE_RANK) - original_num_nodes();
  auto numNewFaces = num_entities(get_bulk(), stk::topology::FACE_RANK) - original_num_faces();
  EXPECT_EQ(numNewNodes, 9U);
  EXPECT_EQ(numNewFaces, 4U);

  check_face_pairs_are_different(disconnecter);
}

auto elems_are_face_adjacent(const stk::mesh::ConnectedEntities& faces0, const stk::mesh::ConnectedEntities& faces1)
{
  std::vector<stk::mesh::Entity> faceList0(faces0.data(), faces0.data() + faces0.size());
  std::vector<stk::mesh::Entity> faceList1(faces1.data(), faces1.data() + faces1.size());
  stk::util::sort_and_unique(faceList0);
  stk::util::sort_and_unique(faceList1);
  std::vector<stk::mesh::Entity> sharedFaces;
  std::set_intersection(
      faceList0.begin(), faceList0.end(), faceList1.begin(), faceList1.end(), std::back_inserter(sharedFaces));
  return !sharedFaces.empty();
}

TEST_F(TestElementDisconnect, disconnect_two_adjacent_elems_gen2x2x2)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 1) GTEST_SKIP() << "Test only available in serial!";

  auto& bulk = get_bulk();
  create_mesh_with_faces(2, 2, 2);

  auto elems = get_elements();
  ASSERT_GE(elems.size(), 2U);
  auto elem0 = elems[0];
  auto elem1 = elems[1];
  auto faces0 = bulk.get_connected_entities(elem0, stk::topology::FACE_RANK);
  auto faces1 = bulk.get_connected_entities(elem1, stk::topology::FACE_RANK);
  ASSERT_TRUE(elems_are_face_adjacent(faces0, faces1));  // How do you do this with the ElemElemGraph?

  ASSERT_EQ(faces0.size(), faces1.size());
  std::vector<stk::mesh::Entity> disconnectFaces;
  for (auto n = 0U; n < faces0.size(); ++n) {
    disconnectFaces.push_back(faces0[n]);
    disconnectFaces.push_back(faces1[n]);
  }
  stk::util::sort_and_unique(disconnectFaces);

  // Each "separated" element will have 8 nodes apiece, and the original block will have 24 nodes remaining
  // 24 + 8 + 8 - 27 (original) = 13
  EntityDisconnectTool disconnecter(bulk, disconnectFaces);
  auto newNodeIds = disconnecter.determine_new_nodes();
  EXPECT_EQ(newNodeIds.size(), 13U);

  disconnecter.modify_mesh();
  auto numNewNodes = num_entities(get_bulk(), stk::topology::NODE_RANK) - original_num_nodes();
  auto numNewFaces = num_entities(get_bulk(), stk::topology::FACE_RANK) - original_num_faces();
  EXPECT_EQ(numNewNodes, 13U);
  EXPECT_EQ(numNewFaces, 5U);

  check_face_pairs_are_different(disconnecter);
}

TEST_F(TestElementDisconnect, disconnect_x_halfplane_gen2x2x2)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 1) GTEST_SKIP() << "Test only available in serial!";
  auto& bulk = get_bulk();
  create_mesh_with_faces(2, 2, 2);

  constexpr double x_plane = 1.0;
  constexpr double y_plane = 1.0;
  auto is_lower_left_face = [](const stk::mesh::BulkData& bulk_, const stk::mesh::Entity& entity) {
    auto centroid = compute_face_centroid(bulk_, entity);
    auto is_lower_xplane = (std::abs(centroid[0] - x_plane) <= 1.0e-10 && centroid[1] < y_plane);
    return is_lower_xplane;
  };

  auto disconnectFaces = get_faces(is_lower_left_face);
  ASSERT_EQ(disconnectFaces.size(), 2U);

  EntityDisconnectTool disconnecter(bulk, disconnectFaces);
  auto newNodeIds = disconnecter.determine_new_nodes();
  EXPECT_EQ(newNodeIds.size(), 3U);

  disconnecter.modify_mesh();

  for (const auto& face : disconnectFaces) {
    const auto faceElems = bulk.get_connected_entities(face, stk::topology::ELEM_RANK);
    EXPECT_EQ(1U, faceElems.size()) << "Face " << bulk.identifier(face)
                                    << " is still related to two elements after disconnect!";
  }
  auto numNewNodes = num_entities(get_bulk(), stk::topology::NODE_RANK) - original_num_nodes();
  auto numNewFaces = num_entities(get_bulk(), stk::topology::FACE_RANK) - original_num_faces();
  EXPECT_EQ(numNewNodes, 3U);
  EXPECT_EQ(numNewFaces, 2U);

  check_face_pairs_are_different(disconnecter);
}

using MeshParam = std::tuple<unsigned, unsigned, unsigned>;
std::string print_mesh_dimensions(const ::testing::TestParamInfo<MeshParam>& p)
{
  std::stringstream ss;
  const auto& [num_x, num_y, num_z] = p.param;
  ss << num_x << "x" << num_y << "x" << num_z;
  return ss.str();
}

class TestElementDisconnectSuite : public TestElementDisconnect, public ::testing::WithParamInterface<MeshParam>
{
 protected:
  TestElementDisconnectSuite() = default;
  auto get_mesh_dimensions() const { return GetParam(); }
};

auto generate_meshes()
{
  auto dimensions = ::testing::Values(1U, 2U, 3U, 4U, 5U);
  return ::testing::Combine(dimensions, dimensions, dimensions);
}

INSTANTIATE_TEST_SUITE_P(testDisconnectSuite, TestElementDisconnectSuite, generate_meshes(), print_mesh_dimensions);

TEST_P(TestElementDisconnectSuite, disconnect_all_elements)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) {
    GTEST_SKIP() << "Test only runs in serial";
  }

  auto [numElemsX, numElemsY, numElemsZ] = get_mesh_dimensions();
  create_mesh_with_faces(numElemsX, numElemsY, numElemsZ);

  auto interiorFaces = get_interior_faces();
  auto exteriorFaces = get_exterior_faces();

  unsigned expectedNumInteriorFaces = original_num_interior_faces();
  unsigned expectedNumExteriorFaces = original_num_exterior_faces();
  EXPECT_EQ(interiorFaces.size(), expectedNumInteriorFaces);
  EXPECT_EQ(exteriorFaces.size(), expectedNumExteriorFaces);

  stk::experimental::EntityDisconnectTool disconnecter(get_bulk(), interiorFaces);
  disconnecter.determine_new_nodes();
  disconnecter.modify_mesh();

  // Each interior face should have been disconnected, creating two new exterior faces for each original interior face
  interiorFaces = get_interior_faces();
  exteriorFaces = get_exterior_faces();
  EXPECT_EQ(interiorFaces.size(), 0U);
  EXPECT_EQ(exteriorFaces.size(), expectedNumExteriorFaces + expectedNumInteriorFaces * 2U);
}

}  // namespace stk::experimental
