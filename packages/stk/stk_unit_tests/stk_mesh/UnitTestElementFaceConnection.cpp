#include <gtest/gtest.h>

#include <vector>
#include <algorithm>
#include <functional>
#include <numeric>
#include <stdlib.h>

#include <stk_topology/topology.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/baseImpl/MeshImplUtils.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/SideSetEntry.hpp>
#include <stk_mesh/base/SidesetUpdater.hpp>
#include <stk_mesh/base/SkinMesh.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_mesh/base/CreateFaces.hpp>
#include <stk_mesh/base/SkinMeshUtil.hpp>
#include "stk_mesh/base/DestroyElements.hpp"  // for destroy_elements_no_mod...
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp>
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraphImpl.hpp>

#include <stk_util/environment/memory_util.hpp>
#include <stk_util/util/ReportHandler.hpp>

#include <stk_io/IossBridge.hpp>
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker

#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_unit_test_utils/getOption.h>

#include "stk_unit_test_utils/ElemGraphTestUtils.hpp"
#include "stk_unit_test_utils/unittestMeshUtils.hpp"
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/BuildMesh.hpp>

namespace {

using namespace stk::mesh::impl;
using namespace stk::mesh;
using stk::unit_test_util::build_mesh;

struct ElementConnectivityInfo {
  ElementConnectivityInfo()
    : element(Entity()), ordinal(INVALID_CONNECTIVITY_ORDINAL), permutation(INVALID_PERMUTATION)
  {  };
  ElementConnectivityInfo(Entity element_)
    : element(element_), ordinal(INVALID_CONNECTIVITY_ORDINAL), permutation(INVALID_PERMUTATION)
  {  }
  ElementConnectivityInfo(Entity element_, ConnectivityOrdinal ordinal_)
    : element(element_), ordinal(ordinal_), permutation(INVALID_PERMUTATION)
  {  }
  ElementConnectivityInfo(Entity element_, int ordinal_)
    : ElementConnectivityInfo(element_, static_cast<stk::mesh::ConnectivityOrdinal>(ordinal_))
  {  }
  ElementConnectivityInfo(Entity element_, stk::mesh::ConnectivityOrdinal ordinal_, stk::mesh::Permutation permutation_)
    : element(element_), ordinal(ordinal_), permutation(permutation_)
  {  }
  ElementConnectivityInfo(Entity element_, stk::mesh::ConnectivityOrdinal ordinal_, int permutation_)
    : ElementConnectivityInfo(element_, ordinal_, static_cast<stk::mesh::Permutation>(permutation_))
  {  }
  ElementConnectivityInfo(Entity element_, int ordinal_, int permutation_)
    : ElementConnectivityInfo(element_, static_cast<stk::mesh::ConnectivityOrdinal>(ordinal_), static_cast<stk::mesh::Permutation>(permutation_))
  {  }

  stk::mesh::Entity element;
  stk::mesh::ConnectivityOrdinal ordinal;
  stk::mesh::Permutation permutation;
};

using ConnectionTesterBeforeElementCreation = std::function<void(const stk::mesh::BulkData& bulk,
                                                                 const ElementConnectivityInfo& eci,
                                                                 stk::mesh::Entity face)>;
using ConnectionTesterAfterElementCreation = std::function<void(const stk::mesh::BulkData& bulk,
                                                                const ElementConnectivityInfo& eci1,
                                                                const ElementConnectivityInfo& eci2,
                                                                stk::mesh::Entity face)>;


std::shared_ptr<stk::mesh::BulkData>
setup_bulk_for_testing_connectivity_after_new_element_creation(unsigned spatialDimension,
                                                               stk::ParallelMachine comm,
                                                               stk::mesh::BulkData::AutomaticAuraOption auraOption = stk::mesh::BulkData::AUTO_AURA)
{
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = stk::unit_test_util::build_mesh(spatialDimension, comm, auraOption);
  stk::mesh::BulkData& bulk = *bulkPtr;
  stk::mesh::MetaData& meta = bulk.mesh_meta_data();

  if (!bulk.has_observer_type<stk::mesh::SidesetUpdater>()) {
    stk::mesh::Selector activeSelector = meta.universal_part();

    bulk.register_observer(std::make_shared<stk::mesh::IncrementalSidesetUpdater>(bulk, activeSelector));
  }

  return bulkPtr;
}

TEST(ElementFaceConnectivity, test_get_ordinals_without_connected_sides)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  if(stk::parallel_machine_size(comm) != 1) {
    GTEST_SKIP();
  }

  const unsigned spatialDimension = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = setup_bulk_for_testing_connectivity_after_new_element_creation(spatialDimension, comm);
  stk::mesh::BulkData& bulk = *bulkPtr;

  bulk.initialize_face_adjacent_element_graph();

  stk::topology hex = stk::topology::HEX_8;
  stk::mesh::Part& block_1 = bulk.mesh_meta_data().declare_part_with_topology("block_1", hex);

  bulk.modification_begin();
  stk::mesh::EntityIdVector elem1Nodes {1, 2, 3, 4, 5, 6, 7, 8};
  stk::mesh::EntityId elem1Id = 1;
  stk::mesh::Entity elem1 = stk::mesh::declare_element(bulk, block_1, elem1Id, elem1Nodes);
  bulk.modification_end();

  bulk.modification_begin();
  stk::mesh::Entity face = bulk.declare_element_side(elem1, 5, stk::mesh::PartVector{});
  bulk.modification_end();

  EXPECT_TRUE(bulk.is_valid(face));

  std::vector<stk::mesh::ConnectivityOrdinal> goldOrdinals{0, 1, 2, 3, 4};
  auto ordinals = stk::mesh::impl::get_ordinals_without_connected_sides(bulk, elem1);

  EXPECT_EQ(goldOrdinals, ordinals);
}

void test_side_connectivity_before_new_element_creation(const stk::mesh::BulkData& mesh,
                                                        const ElementConnectivityInfo& eci,
                                                        stk::mesh::Entity face)
{
  const unsigned numSideElements = mesh.num_elements(face);
  const stk::mesh::Entity* sideElements = mesh.begin_elements(face);
  const stk::mesh::ConnectivityOrdinal * sideOrdinals = mesh.begin_element_ordinals(face);
  const stk::mesh::Permutation * sidePermutations = mesh.begin_element_permutations(face);

  EXPECT_EQ(1u, numSideElements);
  EXPECT_EQ(eci.element, sideElements[0]);

  stk::topology elemTopo = mesh.bucket(sideElements[0]).topology();
  stk::mesh::EntityRank rank = mesh.entity_rank(face);

  unsigned ordinal = elemTopo.side_ordinal(sideOrdinals[0], rank);

  EXPECT_EQ(eci.ordinal, ordinal);
  EXPECT_EQ(eci.permutation, sidePermutations[0]);
}

void test_side_connectivity_after_new_element_creation(const stk::mesh::BulkData& mesh,
                                                       const ElementConnectivityInfo& eci1,
                                                       const ElementConnectivityInfo& eci2,
                                                       stk::mesh::Entity face)
{
  const unsigned numSideElements = mesh.num_elements(face);
  const stk::mesh::Entity* sideElements = mesh.begin_elements(face);
  const stk::mesh::ConnectivityOrdinal * sideOrdinals = mesh.begin_element_ordinals(face);
  const stk::mesh::Permutation * sidePermutations = mesh.begin_element_permutations(face);

  EXPECT_EQ(2u, numSideElements);

  EXPECT_TRUE(eci1.element == sideElements[0] || eci2.element == sideElements[0]);
  EXPECT_TRUE(eci1.element == sideElements[1] || eci2.element == sideElements[1]);

  unsigned index1 = 0;
  unsigned index2 = 1;

  if(eci1.element == sideElements[1]) {
    index1 = 1;
    index2 = 0;
  }

  EXPECT_EQ(eci1.element, sideElements[index1]);
  EXPECT_EQ(eci2.element, sideElements[index2]);

  stk::topology elemTopo1 = mesh.bucket(sideElements[index1]).topology();
  stk::topology elemTopo2 = mesh.bucket(sideElements[index2]).topology();

  stk::mesh::EntityRank rank = mesh.entity_rank(face);

  unsigned ordinal1 = elemTopo1.side_ordinal(sideOrdinals[index1], rank);
  unsigned ordinal2 = elemTopo2.side_ordinal(sideOrdinals[index2], rank);

  EXPECT_EQ(eci1.ordinal, ordinal1);
  EXPECT_EQ(eci2.ordinal, ordinal2);

  EXPECT_EQ(eci1.permutation, sidePermutations[index1]);
  EXPECT_EQ(eci2.permutation, sidePermutations[index2]);
}

void run_hex_face_connectivity_after_new_element_creation(stk::mesh::BulkData& bulk,
                                                          ConnectionTesterBeforeElementCreation testBeforeElementCreation,
                                                          ConnectionTesterAfterElementCreation testAfterElementCreation)
{
  stk::topology hex = stk::topology::HEX_8;
  stk::mesh::Part& block_1 = bulk.mesh_meta_data().declare_part_with_topology("block_1", hex);

  bulk.modification_begin();
  stk::mesh::EntityIdVector elem1Nodes {1, 2, 3, 4, 5, 6, 7, 8};
  stk::mesh::EntityId elem1Id = 1;
  stk::mesh::Entity elem1 = stk::mesh::declare_element(bulk, block_1, elem1Id, elem1Nodes);
  bulk.modification_end();

  bulk.modification_begin();
  stk::mesh::Entity face = bulk.declare_element_side(elem1, 5, stk::mesh::PartVector{});
  bulk.modification_end();

  EXPECT_EQ(stk::topology::FACE_RANK, bulk.entity_rank(face));

  int ordinal1 = 5;
  int permutation1 = 0;
  ElementConnectivityInfo eci1(elem1, ordinal1, permutation1);
  testBeforeElementCreation(bulk, eci1, face);

  bulk.modification_begin();
  stk::mesh::EntityIdVector elem2Nodes {5, 6, 7, 8, 9, 10, 11, 12};
  stk::mesh::EntityId elem2Id = 2;
  stk::mesh::Entity elem2 = stk::mesh::declare_element(bulk, block_1, elem2Id, elem2Nodes);
  bulk.modification_end();

  bulk.modification_begin();
  // Do this in different mod cycle so that the element has chance to be in the graph
  stk::mesh::impl::connect_element_to_existing_sides(bulk, elem2);
  bulk.modification_end();

  int ordinal2 = 4;
  int permutation2 = 4;
  ElementConnectivityInfo eci2(elem2, ordinal2, permutation2);
  testAfterElementCreation(bulk, eci1, eci2, face);
}

TEST(ElementFaceConnectivity, test_hex_face_connectivity_after_new_element_creation_with_graph)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  if(stk::parallel_machine_size(comm) != 1) {
    GTEST_SKIP();
  }

  const unsigned spatialDimension = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = setup_bulk_for_testing_connectivity_after_new_element_creation(spatialDimension, comm);
  stk::mesh::BulkData& bulk = *bulkPtr;

  bulk.initialize_face_adjacent_element_graph();
  const stk::mesh::ElemElemGraph& graph = bulk.get_face_adjacent_element_graph();

  ConnectionTesterBeforeElementCreation testBeforeElementCreation =
      [&graph](const stk::mesh::BulkData& mesh, const ElementConnectivityInfo& eci, stk::mesh::Entity face) {
    test_side_connectivity_before_new_element_creation(mesh, eci, face);
    ASSERT_EQ(1u, graph.size());
    ASSERT_EQ(0u, graph.get_num_connected_elems(eci.element));
  };

  ConnectionTesterAfterElementCreation testAfterElementCreation =
      [&graph](const stk::mesh::BulkData& mesh, const ElementConnectivityInfo& eci1, const ElementConnectivityInfo& eci2, stk::mesh::Entity face) {
    test_side_connectivity_after_new_element_creation(mesh, eci1, eci2, face);
    ASSERT_EQ(2u, graph.size());
    ASSERT_EQ(1u, graph.get_num_connected_elems(eci1.element));
    ASSERT_EQ(1u, graph.get_num_connected_elems(eci2.element));
  };

  run_hex_face_connectivity_after_new_element_creation(bulk, testBeforeElementCreation, testAfterElementCreation);
}

TEST(ElementFaceConnectivity, test_hex_face_connectivity_after_new_element_creation_no_graph)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  if(stk::parallel_machine_size(comm) != 1) {
    GTEST_SKIP();
  }

  const unsigned spatialDimension = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = setup_bulk_for_testing_connectivity_after_new_element_creation(spatialDimension, comm);
  stk::mesh::BulkData& bulk = *bulkPtr;

  ConnectionTesterBeforeElementCreation testBeforeElementCreation =
      [](const stk::mesh::BulkData& mesh, const ElementConnectivityInfo& eci, stk::mesh::Entity face) {
    test_side_connectivity_before_new_element_creation(mesh, eci, face);
  };

  ConnectionTesterAfterElementCreation testAfterElementCreation =
      [](const stk::mesh::BulkData& mesh, const ElementConnectivityInfo& eci1, const ElementConnectivityInfo& eci2, stk::mesh::Entity face) {
    test_side_connectivity_after_new_element_creation(mesh, eci1, eci2, face);
  };

  run_hex_face_connectivity_after_new_element_creation(bulk, testBeforeElementCreation, testAfterElementCreation);
}


void run_shell_face_connectivity_after_new_element_creation(stk::mesh::BulkData& bulk,
                                                            ConnectionTesterBeforeElementCreation testBeforeElementCreation,
                                                            ConnectionTesterAfterElementCreation testAfterElementCreation)
{
  /* shell-quad-4 mesh: */
  /*                    */
  /*                    */
  /*  4*-----*3----*6   */
  /*   |     |     |    */
  /*   |     |     |    */
  /*  1*-----*2----*5   */
  /*                    */

  stk::topology shell4 = stk::topology::SHELL_QUAD_4_ALL_FACE_SIDES;
  stk::mesh::Part& block_1 = bulk.mesh_meta_data().declare_part_with_topology("block_1", shell4);

  bulk.modification_begin();
  stk::mesh::EntityIdVector elem1Nodes {1, 2, 3, 4};
  stk::mesh::EntityId elem1Id = 1;
  stk::mesh::Entity elem1 = stk::mesh::declare_element(bulk, block_1, elem1Id, elem1Nodes);
  bulk.modification_end();

  bulk.modification_begin();
  stk::mesh::Entity face = bulk.declare_element_side(elem1, 3, stk::mesh::PartVector{});
  bulk.modification_end();

  EXPECT_EQ(stk::topology::FACE_RANK, bulk.entity_rank(face));

  int ordinal1 = 3;
  int permutation1 = 0;
  ElementConnectivityInfo eci1(elem1, ordinal1, permutation1);
  testBeforeElementCreation(bulk, eci1, face);

  bulk.modification_begin();
  stk::mesh::EntityIdVector elem2Nodes {2, 5, 6, 3};
  stk::mesh::EntityId elem2Id = 2;
  stk::mesh::Entity elem2 = stk::mesh::declare_element(bulk, block_1, elem2Id, elem2Nodes);
  bulk.modification_end();

  bulk.modification_begin();
  // Do this in different mod cycle so that the element has chance to be in the graph
  stk::mesh::impl::connect_element_to_existing_sides(bulk, elem2);
  bulk.modification_end();

  int ordinal2 = 5;
  int permutation2 = 1;
  ElementConnectivityInfo eci2(elem2, ordinal2, permutation2);
  testAfterElementCreation(bulk, eci1, eci2, face);
}

TEST(ElementFaceConnectivity, test_paved_shell_face_connectivity_after_new_element_creation_no_graph)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  if(stk::parallel_machine_size(comm) != 1) {
    GTEST_SKIP();
  }

  const unsigned spatialDimension = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = setup_bulk_for_testing_connectivity_after_new_element_creation(spatialDimension, comm);
  stk::mesh::BulkData& bulk = *bulkPtr;

  ConnectionTesterBeforeElementCreation testBeforeElementCreation =
      [](const stk::mesh::BulkData& mesh, const ElementConnectivityInfo& eci, stk::mesh::Entity face) {
    test_side_connectivity_before_new_element_creation(mesh, eci, face);
  };

  ConnectionTesterAfterElementCreation testAfterElementCreation =
      [](const stk::mesh::BulkData& mesh, const ElementConnectivityInfo& eci1, const ElementConnectivityInfo& eci2, stk::mesh::Entity face) {
    test_side_connectivity_after_new_element_creation(mesh, eci1, eci2, face);
  };

  run_shell_face_connectivity_after_new_element_creation(bulk, testBeforeElementCreation, testAfterElementCreation);
}

void run_shell_edge_connectivity_after_new_element_creation(stk::mesh::BulkData& bulk,
                                                            ConnectionTesterBeforeElementCreation testBeforeElementCreation,
                                                            ConnectionTesterAfterElementCreation testAfterElementCreation)
{
  /* shell-quad-4 mesh: */
  /*                    */
  /*                    */
  /*  4*-----*3----*6   */
  /*   |     |     |    */
  /*   |     |     |    */
  /*  1*-----*2----*5   */
  /*                    */

  stk::topology shell4 = stk::topology::SHELL_QUAD_4;
  stk::mesh::Part& block_1 = bulk.mesh_meta_data().declare_part_with_topology("block_1", shell4);

  bulk.modification_begin();
  stk::mesh::EntityIdVector elem1Nodes {1, 2, 3, 4};
  stk::mesh::EntityId elem1Id = 1;
  stk::mesh::Entity elem1 = stk::mesh::declare_element(bulk, block_1, elem1Id, elem1Nodes);
  bulk.modification_end();

  bulk.modification_begin();
  stk::mesh::Entity edge = bulk.declare_element_side(elem1, 3, stk::mesh::PartVector{});
  bulk.modification_end();

  EXPECT_EQ(stk::topology::EDGE_RANK, bulk.entity_rank(edge));

  int ordinal1 = 3;
  int permutation1 = 0;
  ElementConnectivityInfo eci1(elem1, ordinal1, permutation1);
  testBeforeElementCreation(bulk, eci1, edge);

  bulk.modification_begin();
  stk::mesh::EntityIdVector elem2Nodes {2, 5, 6, 3};
  stk::mesh::EntityId elem2Id = 2;
  stk::mesh::Entity elem2 = stk::mesh::declare_element(bulk, block_1, elem2Id, elem2Nodes);
  bulk.modification_end();

  bulk.modification_begin();
  // Do this in different mod cycle so that the element has chance to be in the graph
  stk::mesh::impl::connect_element_to_existing_sides(bulk, elem2);
  bulk.modification_end();

  int ordinal2 = 5;
  int permutation2 = 1;
  ElementConnectivityInfo eci2(elem2, ordinal2, permutation2);
  testAfterElementCreation(bulk, eci1, eci2, edge);
}

TEST(ElementFaceConnectivity, test_paved_shell_edge_connectivity_after_new_element_creation_no_graph)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  if(stk::parallel_machine_size(comm) != 1) {
    GTEST_SKIP();
  }

  const unsigned spatialDimension = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = setup_bulk_for_testing_connectivity_after_new_element_creation(spatialDimension, comm);
  stk::mesh::BulkData& bulk = *bulkPtr;

  ConnectionTesterBeforeElementCreation testBeforeElementCreation =
      [](const stk::mesh::BulkData& mesh, const ElementConnectivityInfo& eci, stk::mesh::Entity side) {
    test_side_connectivity_before_new_element_creation(mesh, eci, side);
  };

  ConnectionTesterAfterElementCreation testAfterElementCreation =
      [](const stk::mesh::BulkData& mesh, const ElementConnectivityInfo& eci1, const ElementConnectivityInfo& eci2, stk::mesh::Entity side) {
    test_side_connectivity_after_new_element_creation(mesh, eci1, eci2, side);
  };

  run_shell_edge_connectivity_after_new_element_creation(bulk, testBeforeElementCreation, testAfterElementCreation);
}

TEST(ElementFaceConnectivity, test_paved_shell_edge_connectivity_after_new_element_creation_with_graph)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  if(stk::parallel_machine_size(comm) != 1) {
    GTEST_SKIP();
  }

  const unsigned spatialDimension = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = setup_bulk_for_testing_connectivity_after_new_element_creation(spatialDimension, comm);
  stk::mesh::BulkData& bulk = *bulkPtr;

  bulk.initialize_face_adjacent_element_graph();
  const stk::mesh::ElemElemGraph& graph = bulk.get_face_adjacent_element_graph();

  ConnectionTesterBeforeElementCreation testBeforeElementCreation =
      [&graph](const stk::mesh::BulkData& mesh, const ElementConnectivityInfo& eci, stk::mesh::Entity side) {
    test_side_connectivity_before_new_element_creation(mesh, eci, side);
    ASSERT_EQ(1u, graph.size());
    ASSERT_EQ(0u, graph.get_num_connected_elems(eci.element));
  };

  ConnectionTesterAfterElementCreation testAfterElementCreation =
      [&graph](const stk::mesh::BulkData& mesh, const ElementConnectivityInfo& eci1, const ElementConnectivityInfo& eci2, stk::mesh::Entity side) {
    test_side_connectivity_after_new_element_creation(mesh, eci1, eci2, side);
    ASSERT_EQ(2u, graph.size());
    ASSERT_EQ(1u, graph.get_num_connected_elems(eci1.element));
    ASSERT_EQ(1u, graph.get_num_connected_elems(eci2.element));
  };

  run_shell_edge_connectivity_after_new_element_creation(bulk, testBeforeElementCreation, testAfterElementCreation);
}

bool element_external_faces(const stk::mesh::BulkData& mesh, const stk::mesh::Selector& selector,
                            const stk::mesh::Entity elem, std::vector<bool>& faceIsExternal)
{
  // SM code to figure out external faces to create shell contact elements on
  faceIsExternal.clear();
  const stk::topology elemTopology = mesh.bucket(elem).topology();
  const unsigned num_node = elemTopology.num_nodes();
  const unsigned num_face = elemTopology.num_faces();
  stk::mesh::Entity const* nodes = mesh.begin_nodes(elem);
  std::vector<std::set<stk::mesh::EntityId> > connElems(num_node);
  bool hasExternal = false;

  for (unsigned n = 0; n < num_node; ++n) {
    stk::mesh::Entity node = nodes[n];
    stk::mesh::Entity const* conn_elems = mesh.begin_elements(node);
    unsigned num_conn_elems = mesh.num_elements(node);
    std::set<stk::mesh::EntityId> conn_elem_ids;
    for (unsigned ie = 0; ie < num_conn_elems; ++ie) {
      if (selector(mesh.bucket(elem))) {
        conn_elem_ids.insert(mesh.identifier(conn_elems[ie]));
      }
    }
    connElems[n] = conn_elem_ids;
  }
  for (unsigned iface = 0; iface < num_face; ++iface) {
    stk::topology face_topology = elemTopology.face_topology(iface);
    const unsigned num_nodes_per_face = face_topology.num_nodes();
    STK_ThrowAssert(num_nodes_per_face <= 4);
    unsigned faceNodeIndex[4];
    elemTopology.face_node_ordinals(iface, faceNodeIndex);
    std::map<stk::mesh::EntityId, unsigned> elemIdFoundCount;
    for (unsigned inode = 0; inode < num_nodes_per_face; ++inode) {
      std::set<stk::mesh::EntityId>& connElemIDs = connElems[faceNodeIndex[inode]];
      for (stk::mesh::EntityId gid : connElemIDs) {
        if (elemIdFoundCount.find(gid) != elemIdFoundCount.end()) {
          elemIdFoundCount[gid] += 1;
        } else {
          elemIdFoundCount[gid] = 1;
        }
      }
    }
    unsigned numCommon = 0;
    for (auto& m : elemIdFoundCount) {
      if (m.second == num_nodes_per_face) {
        numCommon += 1;
      }
    }
    if (numCommon == 1) {
      faceIsExternal.push_back(true);
      hasExternal = true;
    } else {
      faceIsExternal.push_back(false);
    }
  }
  return hasExternal;
}

void delete_created_external_shells(stk::mesh::BulkData& bulk, std::vector<stk::mesh::Entity>& createdElems)
{
  stk::mesh::Selector orphansToDelete;
  bulk.modification_begin();
  stk::mesh::destroy_elements_no_mod_cycle(bulk, createdElems, orphansToDelete);
  bulk.modification_end();
}

// map the nodes of a tri3 to degenerated quad4
constexpr unsigned subQuadToTriNodeMap[] = {0, 1, 2, 2};

void create_external_shell(stk::mesh::BulkData& bulk,
                           const stk::mesh::Entity elem,
                           const std::vector<bool>& faceExternal,
                           stk::mesh::Part& blockPart,
                           const std::vector<unsigned>& shellElemIds,
                           std::vector<stk::mesh::Entity>& createdElems)
{
  stk::topology const elemTopology = bulk.bucket(elem).topology();

  std::vector<stk::mesh::Relation> nodeRelations;
  stk::mesh::EntityIdVector nodeIds;

  const stk::mesh::Entity* const elemNodes = bulk.begin_nodes(elem);
  const unsigned num_face = elemTopology.num_faces();
  for (unsigned iface = 0; iface < num_face; ++iface) {
    stk::topology faceTopology = elemTopology.face_topology(iface);
    STK_ThrowRequire(faceTopology == stk::topology::TRI_3 || faceTopology == stk::topology::QUAD_4);

    const unsigned numContSurfNodes = faceTopology.num_nodes();
    nodeRelations.reserve(numContSurfNodes);

    unsigned faceNodeIndex[4];
    elemTopology.face_node_ordinals(iface, faceNodeIndex);
    if (faceExternal[iface]) {
      nodeRelations.clear();
      nodeIds.clear();
      if (faceTopology == stk::topology::TRI_3) {
        for (unsigned inod = 0; inod < numContSurfNodes; ++inod) {
          auto node = elemNodes[faceNodeIndex[subQuadToTriNodeMap[inod]]];
          nodeRelations.emplace_back(stk::topology::NODE_RANK, node, stk::mesh::RelationType::USES, inod);
          nodeIds.push_back(bulk.identifier(node));
        }
      } else {
        for (unsigned inod = 0; inod < numContSurfNodes; ++inod) {
          auto node = elemNodes[faceNodeIndex[inod]];
          nodeRelations.emplace_back(stk::topology::NODE_RANK, node, stk::mesh::RelationType::USES, inod);
          nodeIds.push_back(bulk.identifier(node));
        }
      }
      std::cout << "P" << bulk.parallel_rank() << ": Creating external shell with id: "
                << shellElemIds[iface] << " on face: " << iface << " of element: " << bulk.identifier(elem) << std::endl;

      auto element = stk::mesh::declare_element(bulk, blockPart, shellElemIds[iface], nodeIds);
      createdElems.push_back(element);
    }
  }
}

TEST(ElementFaceConnectivity, test_SM_create_shell_and_delete_on_external_side)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  if(stk::parallel_machine_size(comm) != 2) {
    GTEST_SKIP();
  }

  const unsigned spatialDimension = 3;
  auto auraOption = stk::mesh::BulkData::NO_AUTO_AURA;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = setup_bulk_for_testing_connectivity_after_new_element_creation(spatialDimension, comm, auraOption);
  stk::mesh::BulkData& bulk = *bulkPtr;

  stk::mesh::Part& block_2 = bulk.mesh_meta_data().declare_part_with_topology("block_2", stk::topology::SHELL_QUAD_4);
  stk::mesh::Part& block_3 = bulk.mesh_meta_data().declare_part_with_topology("block_3", stk::topology::SHELL_QUAD_4);
  stk::io::put_io_part_attribute(block_2);
  stk::io::put_io_part_attribute(block_3);

  stk::mesh::PartVector blocks{&block_2, &block_3};

  stk::io::fill_mesh("textmesh:0,1,HEX_8, 1, 2, 3, 4, 5, 6, 7, 8, block_1\n\
                               1,2,HEX_8, 5, 6, 7, 8, 9,10,11,12, block_1\n\
                               |coordinates:0,0,0,  1,0,0,   1,1,0,   0,1,0,\n\
                                            0,0,1,  1,0,1,   1,1,1,   0,1,1,\n\
                                            0,0,2,  1,0,2,   1,1,2,   0,1,2", bulk);

  bulk.initialize_face_adjacent_element_graph();

  stk::mesh::EntityVector elems, createdElems;
  stk::mesh::get_entities(bulk, stk::topology::ELEM_RANK, elems);

  unsigned numElem = elems.size();
  EXPECT_EQ(1u, numElem);

  std::vector<unsigned> shellElemIds(6);
  std::iota(shellElemIds.begin(), shellElemIds.end(), 3 + bulk.parallel_rank()*6*numElem);

  stk::mesh::Selector selector = bulk.mesh_meta_data().universal_part();
  for (unsigned ielem = 0; ielem < numElem; ielem++) {
    stk::mesh::Entity const elem = elems[ielem];
    std::vector<bool> faceExternal;
    bool hasExternal = element_external_faces(bulk, selector, elem, faceExternal);
    EXPECT_TRUE(hasExternal);
    EXPECT_EQ(6u, faceExternal.size());
    for(bool faceStatus : faceExternal) {
      EXPECT_TRUE(faceStatus);
    }
  }

  bulk.modification_begin();
  for (unsigned ielem = 0; ielem < numElem; ielem++) {
    stk::mesh::Entity const elem = elems[ielem];

    // Side 0 only to create parallel paved configuration
    std::vector<bool> faceExternal(6, false);
    faceExternal[0] = true;

    create_external_shell(bulk, elem, faceExternal, *blocks[bulk.parallel_rank()], shellElemIds, createdElems);
  }
  bulk.modification_end();

  EXPECT_NO_THROW(delete_created_external_shells(bulk, createdElems));
}


} // namespace
