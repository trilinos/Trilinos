#include "gtest/gtest.h"
#include "mpi.h"
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp>
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraphUpdater.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/BuildMesh.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>

namespace
{
using stk::unit_test_util::build_mesh;

void setup_node_sharing(stk::mesh::BulkData &mesh, const std::vector< std::vector<unsigned> > & shared_nodeIDs_and_procs )
{
  const unsigned p_rank = mesh.parallel_rank();

  for (size_t nodeIdx = 0, end = shared_nodeIDs_and_procs.size(); nodeIdx < end; ++nodeIdx) {
    if (p_rank == shared_nodeIDs_and_procs[nodeIdx][0]) {
      stk::mesh::EntityId nodeID = shared_nodeIDs_and_procs[nodeIdx][1];
      int sharingProc = shared_nodeIDs_and_procs[nodeIdx][2];
      stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, nodeID);
      mesh.add_node_sharing(node, sharingProc);
    }
  }
}

class HexShellShell : public stk::unit_test_util::MeshFixture
{
protected:
  HexShellShell()
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  }

  void setup_hex_shell_shell_on_procs(std::vector<int> owningProcs)
  {
    stk::mesh::Part* hexPart = &get_meta().declare_part_with_topology("hex_part", stk::topology::HEX_8);
    stk::mesh::Part* shellPart = &get_meta().declare_part_with_topology("shell_part", stk::topology::SHELL_QUAD_4);
    stk::mesh::PartVector parts = {hexPart, shellPart, shellPart};
    declare_elements_on_procs_and_setup_node_sharing(owningProcs, parts);
  }

private:
  void declare_elements_on_procs_and_setup_node_sharing(const std::vector<int>& owningProcs, const stk::mesh::PartVector& parts)
  {
    get_bulk().modification_begin();
    declare_elements_on_procs(owningProcs, parts);
    if(get_bulk().parallel_size()==2)
    {
      setup_node_sharing(get_bulk(), shared_nodeIDs_and_procs2);
    }
    else
    {
      setup_node_sharing(get_bulk(), shared_nodeIDs_and_procs3);
    }
    get_bulk().modification_end();
  }

  void declare_elements_on_procs(const std::vector<int>& owningProcs, const stk::mesh::PartVector& parts)
  {
    for(size_t i = 0; i < nodeIDsPerElement.size(); ++i)
      if(owningProcs[i] == stk::parallel_machine_rank(get_comm()))
        stk::mesh::declare_element(get_bulk(), *parts[i], elemIDs[i], nodeIDsPerElement[i]);
  }

  std::vector<stk::mesh::EntityIdVector> nodeIDsPerElement { {1, 2, 3, 4, 5, 6, 7, 8}, {5, 6, 7, 8}, {5, 6, 7, 8}};
  stk::mesh::EntityIdVector elemIDs = {1, 2, 3};
  // list of triplets: (owner-proc, shared-nodeID, sharing-proc)
  std::vector<std::vector<unsigned> > shared_nodeIDs_and_procs2 = {
    {0, 5, 1}, {0, 6, 1}, {0, 7, 1}, {0, 8, 1},  // proc 0
    {1, 5, 0}, {1, 6, 0}, {1, 7, 0}, {1, 8, 0}}; // proc 1
  std::vector<std::vector<unsigned> > shared_nodeIDs_and_procs3 = {
    {0, 5, 1}, {0, 6, 1}, {0, 7, 1}, {0, 8, 1},  // proc 0
    {0, 5, 2}, {0, 6, 2}, {0, 7, 2}, {0, 8, 2},  // proc 0
    {1, 5, 0}, {1, 6, 0}, {1, 7, 0}, {1, 8, 0}, // proc 1
    {1, 5, 2}, {1, 6, 2}, {1, 7, 2}, {1, 8, 2}, // proc 1
    {2, 5, 0}, {2, 6, 0}, {2, 7, 0}, {2, 8, 0}, // proc 1
    {2, 5, 1}, {2, 6, 1}, {2, 7, 1}, {2, 8, 1}}; // proc 1
};


TEST_F(HexShellShell, Hex0Shell1Shell1Parallel)
{
  //  ID.proc
  //
  //          3.0------------7.0
  //          /|             /|
  //         / |            / |
  //        /  |           /  |
  //      4.0------------8.0  |
  //       |   |          |   |
  //       |   |   1.0    |2.1|
  //       |   |          |3.1|
  //       |  2.0---------|--6.0
  //       |  /           |  /
  //       | /            | /
  //       |/             |/
  //      1.0------------5.0
  //                      ^
  //                      |
  //                       ---- Two stacked shells

  if(stk::parallel_machine_size(get_comm()) == 2u)
  {
    setup_hex_shell_shell_on_procs({0, 1, 1});

    stk::mesh::ElemElemGraph elemElemGraph(get_bulk());

    if(stk::parallel_machine_rank(get_comm()) == 0)
    {
      const stk::mesh::Entity hex1 = get_bulk().get_entity(stk::topology::ELEM_RANK, 1);
      ASSERT_EQ(2u, elemElemGraph.get_num_connected_elems(hex1));
      EXPECT_EQ(2u, elemElemGraph.get_connected_remote_id_and_via_side(hex1, 0).id);
      EXPECT_EQ(5, elemElemGraph.get_connected_remote_id_and_via_side(hex1, 0).side);
      EXPECT_EQ(3u, elemElemGraph.get_connected_remote_id_and_via_side(hex1, 1).id);
      EXPECT_EQ(5, elemElemGraph.get_connected_remote_id_and_via_side(hex1, 1).side);
      EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));
      EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex1, 1));
      EXPECT_EQ(2u, elemElemGraph.num_edges());
      EXPECT_EQ(2u, elemElemGraph.num_parallel_edges());
    }
    else
    {
      const stk::mesh::Entity shell2 = get_bulk().get_entity(stk::topology::ELEM_RANK, 2);
      ASSERT_EQ(1u, elemElemGraph.get_num_connected_elems(shell2));
      EXPECT_EQ(1, elemElemGraph.get_connected_remote_id_and_via_side(shell2, 0).side);
      EXPECT_EQ(1u, elemElemGraph.get_connected_remote_id_and_via_side(shell2, 0).id);
      EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell2, 0));

      const stk::mesh::Entity shell3 = get_bulk().get_entity(stk::topology::ELEM_RANK, 3);
      ASSERT_EQ(1u, elemElemGraph.get_num_connected_elems(shell3));
      EXPECT_EQ(1, elemElemGraph.get_connected_remote_id_and_via_side(shell3, 0).side);
      EXPECT_EQ(1u, elemElemGraph.get_connected_remote_id_and_via_side(shell3, 0).id);
      EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell3, 0));
    }
  }
}

// disabled due to split coincident elements
TEST_F(HexShellShell, DISABLED_Hex0Shell0Shell1Parallel )
{
  //  ID.proc
  //
  //          3.0------------7.0
  //          /|             /|
  //         / |            / |
  //        /  |           /  |
  //      4.0------------8.0  |
  //       |   |          |   |
  //       |   |   1.0    |2.0|
  //       |   |          |3.1|
  //       |  2.0---------|--6.0
  //       |  /           |  /
  //       | /            | /
  //       |/             |/
  //      1.0------------5.0
  //                      ^
  //                      |
  //                       ---- Two stacked shells

  if(stk::parallel_machine_size(get_comm()) == 2u)
  {
    setup_hex_shell_shell_on_procs({0, 0, 1});

    stk::mesh::ElemElemGraph elemElemGraph(get_bulk());

    if(stk::parallel_machine_rank(get_comm()) == 0)
    {
      const stk::mesh::Entity hex1 = get_bulk().get_entity(stk::topology::ELEM_RANK, 1);
      const stk::mesh::Entity shell2 = get_bulk().get_entity(stk::topology::ELEM_RANK, 2);
      ASSERT_EQ(2u, elemElemGraph.get_num_connected_elems(hex1));
      ASSERT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));
      EXPECT_EQ(5, elemElemGraph.get_connected_element_and_via_side(hex1, 0).side);
      EXPECT_EQ(shell2, elemElemGraph.get_connected_element_and_via_side(hex1, 0).element);
      ASSERT_TRUE(!elemElemGraph.is_connected_elem_locally_owned(hex1, 1));
      EXPECT_EQ(3u, elemElemGraph.get_connected_remote_id_and_via_side(hex1, 1).id);
      EXPECT_EQ(5, elemElemGraph.get_connected_remote_id_and_via_side(hex1, 1).side);

      ASSERT_EQ(1u, elemElemGraph.get_num_connected_elems(shell2));
      EXPECT_EQ(1, elemElemGraph.get_connected_element_and_via_side(shell2, 0).side);
      EXPECT_EQ(hex1, elemElemGraph.get_connected_element_and_via_side(shell2, 0).element);
      EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell2, 0));
      EXPECT_EQ(3u, elemElemGraph.num_edges());
      EXPECT_EQ(1u, elemElemGraph.num_parallel_edges());
    }
    else
    {
      const stk::mesh::Entity shell3 = get_bulk().get_entity(stk::topology::ELEM_RANK, 3);
      ASSERT_EQ(1u, elemElemGraph.get_num_connected_elems(shell3));
      EXPECT_EQ(1, elemElemGraph.get_connected_remote_id_and_via_side(shell3, 0).side);
      EXPECT_EQ(1u, elemElemGraph.get_connected_remote_id_and_via_side(shell3, 0).id);
      EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell3, 0));
      EXPECT_EQ(1u, elemElemGraph.num_edges());
      EXPECT_EQ(1u, elemElemGraph.num_parallel_edges());
    }
  }
}

// disabled due to split coincident elements
TEST_F(HexShellShell, DISABLED_Skin)
{
  if(stk::parallel_machine_size(get_comm()) == 2u)
  {
    setup_hex_shell_shell_on_procs({0, 1, 0});

    stk::mesh::create_exposed_block_boundary_sides(get_bulk(), get_meta().universal_part(), {});

    stk::mesh::Selector ownedOrShared = get_meta().locally_owned_part() | get_meta().globally_shared_part();

    if(stk::parallel_machine_rank(get_comm()) == 0)
      EXPECT_EQ(6u, stk::mesh::count_selected_entities(ownedOrShared, get_bulk().buckets(stk::topology::FACE_RANK)));
    else
      EXPECT_EQ(1u, stk::mesh::count_selected_entities(ownedOrShared, get_bulk().buckets(stk::topology::FACE_RANK)));
  }
}

void expect_id_ord_perm(stk::mesh::BulkData &bulk,
                        stk::mesh::Entity face,
                        unsigned elemOffset,
                        const stk::mesh::EntityId elemId,
                        const stk::mesh::ConnectivityOrdinal ord,
                        const stk::mesh::Permutation perm)
{
  EXPECT_EQ(elemId, bulk.identifier(bulk.begin_elements(face)[elemOffset])) << "elemOffset=" << elemOffset;
  EXPECT_EQ(ord, bulk.begin_element_ordinals(face)[elemOffset]) << "elemOffset=" << elemOffset;
  EXPECT_EQ(perm, bulk.begin_element_permutations(face)[elemOffset]) << "elemOffset=" << elemOffset;
}

TEST_F(HexShellShell, SideConnections)
{
  if(stk::parallel_machine_size(get_comm()) == 1u)
  {
    setup_hex_shell_shell_on_procs({0, 0, 0});

    stk::mesh::ElemElemGraph elemElemGraph(get_bulk());
    stk::mesh::SideConnector& sideConnector = elemElemGraph.get_side_connector();

    get_bulk().modification_begin();
    stk::mesh::Entity shell2 = get_bulk().get_entity(stk::topology::ELEM_RANK, 2);
    int shellSide = 1;
    stk::mesh::Entity face = get_bulk().declare_element_side(shell2, shellSide, stk::mesh::ConstPartVector{});
    sideConnector.connect_side_to_all_elements(face, shell2, shellSide);
    get_bulk().modification_end();

    ASSERT_EQ(3u, get_bulk().num_elements(face));
    expect_id_ord_perm(get_bulk(), face, 0, 2, stk::mesh::ConnectivityOrdinal(1), stk::mesh::Permutation(0));
    expect_id_ord_perm(get_bulk(), face, 1, 3, stk::mesh::ConnectivityOrdinal(1), stk::mesh::Permutation(0));
    expect_id_ord_perm(get_bulk(), face, 2, 1, stk::mesh::ConnectivityOrdinal(5), stk::mesh::Permutation(4));
  }
}

void expect_correct_connected_element_via_side(stk::mesh::ElemElemGraph& elemElemGraph, stk::mesh::Entity elem, int k, stk::mesh::Entity otherElem, int viaSide)
{
  stk::mesh::impl::ElementViaSidePair elem_via_side = elemElemGraph.get_connected_element_and_via_side(elem, k);
  EXPECT_EQ(viaSide,      elem_via_side.side);
  EXPECT_EQ(otherElem, elem_via_side.element);
}


TEST( ElementGraph, HexAddShellAddShellSerial )
{
  //  ID.proc
  //
  //          3.0------------7.0
  //          /|             /|
  //         / |            / |
  //        /  |           /  |
  //      4.0------------8.0  |
  //       |   |          |   |
  //       |   |   1.0    |2.0|
  //       |   |          |3.0|
  //       |  2.0---------|--6.0
  //       |  /           |  /
  //       | /            | /
  //       |/             |/
  //      1.0------------5.0
  //                      ^
  //                      |
  //                       ---- Two stacked shell elements

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  int p_size = stk::parallel_machine_size(pm);

  if(p_size > 1)
  {
    return;
  }

  const unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, pm);
  stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();
  stk::mesh::BulkData& mesh = *bulkPtr;

  stk::mesh::Part * hexPart   = &meta.declare_part_with_topology("hex_part",   stk::topology::HEX_8);
  stk::mesh::Part * shellPart = &meta.declare_part_with_topology("shell_part", stk::topology::SHELL_QUAD_4);
  meta.commit();

  std::vector<stk::mesh::EntityIdVector> hexNodeIDs {
    { 1, 2, 3, 4, 5,  6,  7,  8 }
  };
  stk::mesh::EntityId hexElemIDs[] = { 1 };

  std::vector<stk::mesh::EntityIdVector> shellNodeIDs {
    { 5, 6, 7, 8 },
    { 5, 6, 7, 8 }
  };
  stk::mesh::EntityId shellElemIDs[] = { 2, 3 };

  mesh.modification_begin();
  for (size_t i = 0; i < hexNodeIDs.size(); ++i) {
    stk::mesh::declare_element(mesh, *hexPart, hexElemIDs[i], hexNodeIDs[i]);
  }
  mesh.modification_end();

  stk::mesh::ElemElemGraph elemElemGraph(mesh);

  mesh.modification_begin();
  stk::mesh::EntityVector added_shells;
  for (size_t i = 0; i < shellNodeIDs.size(); ++i) {
    added_shells.push_back( stk::mesh::declare_element(mesh, *shellPart, shellElemIDs[i], shellNodeIDs[i]) );
  }
  mesh.modification_end();

  elemElemGraph.add_elements(added_shells);

  const stk::mesh::Entity hex1   = mesh.get_entity(stk::topology::ELEM_RANK, 1);
  const stk::mesh::Entity shell2 = mesh.get_entity(stk::topology::ELEM_RANK, 2);
  const stk::mesh::Entity shell3 = mesh.get_entity(stk::topology::ELEM_RANK, 3);

  // Connectivity for Hex Element 1
  EXPECT_EQ(2u,     elemElemGraph.get_num_connected_elems(hex1));
  expect_correct_connected_element_via_side(elemElemGraph, hex1, 0, shell2, 5);
  expect_correct_connected_element_via_side(elemElemGraph, hex1, 1, shell3, 5);
  EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));
  EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 1));

  // Connectivity for Shell Element 2
  EXPECT_EQ(1u,   elemElemGraph.get_num_connected_elems(shell2));
  expect_correct_connected_element_via_side(elemElemGraph, shell2, 0, hex1, 1);
  EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell2, 0));

  // Connectivity for Shell Element 3
  EXPECT_EQ(1u,   elemElemGraph.get_num_connected_elems(shell3));
  expect_correct_connected_element_via_side(elemElemGraph, shell3, 0, hex1, 1);
  EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell3, 0));

  EXPECT_EQ(0u, elemElemGraph.num_parallel_edges());
}

TEST( ElementGraph, HexAddShellAddShellHexSerial )
{
  //  ID.proc
  //
  //          3.0------------7.0-----------11.0
  //          /|             /|             /|
  //         / |            / |            / |
  //        /  |           /  |           /  |
  //      4.0------------8.0-----------12.0  |
  //       |   |          |   |          |   |
  //       |   |   1.0    |3.0|   2.0    |   |
  //       |   |          |4.0|          |   |
  //       |  2.0---------|--6.0---------|-10.0
  //       |  /           |  /           |  /
  //       | /            | /            | /
  //       |/             |/             |/
  //      1.0------------5.0------------9.0
  //                      ^
  //                      |
  //                       ---- Added two stacked shell elements

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  int p_size = stk::parallel_machine_size(pm);

  if(p_size > 1)
  {
    return;
  }

  const unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, pm);
  stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();
  stk::mesh::BulkData& mesh = *bulkPtr;

  stk::mesh::Part * hexPart   = &meta.declare_part_with_topology("hex_part",   stk::topology::HEX_8);
  stk::mesh::Part * shellPart = &meta.declare_part_with_topology("shell_part", stk::topology::SHELL_QUAD_4);
  meta.commit();

  std::vector<stk::mesh::EntityIdVector> hexNodeIDs {
    { 1, 2, 3, 4, 5,  6,  7,  8 },
    { 5, 6, 7, 8, 9, 10, 11, 12 }
  };
  stk::mesh::EntityId hexElemIDs[] = { 1, 2 };

  std::vector<stk::mesh::EntityIdVector> shellNodeIDs {
    { 5, 6, 7, 8 },
    { 5, 6, 7, 8 }
  };
  stk::mesh::EntityId shellElemIDs[] = { 3, 4 };

  mesh.modification_begin();
  for (size_t i = 0; i < hexNodeIDs.size(); ++i) {
    stk::mesh::declare_element(mesh, *hexPart, hexElemIDs[i], hexNodeIDs[i]);
  }
  mesh.modification_end();

  stk::mesh::ElemElemGraph elemElemGraph(mesh);

  mesh.modification_begin();
  stk::mesh::EntityVector added_shells;
  for (size_t i = 0; i < shellNodeIDs.size(); ++i) {
    added_shells.push_back( stk::mesh::declare_element(mesh, *shellPart, shellElemIDs[i], shellNodeIDs[i]) );
  }
  mesh.modification_end();

  elemElemGraph.add_elements(added_shells);

  const stk::mesh::Entity hex1   = mesh.get_entity(stk::topology::ELEM_RANK, 1);
  const stk::mesh::Entity hex2   = mesh.get_entity(stk::topology::ELEM_RANK, 2);
  const stk::mesh::Entity shell3 = mesh.get_entity(stk::topology::ELEM_RANK, 3);
  const stk::mesh::Entity shell4 = mesh.get_entity(stk::topology::ELEM_RANK, 4);

  // Connectivity for Hex Element 1
  EXPECT_EQ(2u,     elemElemGraph.get_num_connected_elems(hex1));
  expect_correct_connected_element_via_side(elemElemGraph, hex1, 0, shell3, 5);
  expect_correct_connected_element_via_side(elemElemGraph, hex1, 1, shell4, 5);
  EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));
  EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 1));

  // Connectivity for Hex Element 2
  EXPECT_EQ(2u,     elemElemGraph.get_num_connected_elems(hex2));
  expect_correct_connected_element_via_side(elemElemGraph, hex2, 0, shell3, 4);
  expect_correct_connected_element_via_side(elemElemGraph, hex2, 1, shell4, 4);
  EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex2, 0));
  EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex2, 1));

  // Connectivity for Shell Element 3
  EXPECT_EQ(2u,   elemElemGraph.get_num_connected_elems(shell3));
  expect_correct_connected_element_via_side(elemElemGraph, shell3, 0, hex1, 1);
  expect_correct_connected_element_via_side(elemElemGraph, shell3, 1, hex2, 0);
  EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell3, 0));
  EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell3, 1));

  // Connectivity for Shell Element 4
  EXPECT_EQ(2u,   elemElemGraph.get_num_connected_elems(shell4));
  expect_correct_connected_element_via_side(elemElemGraph, shell4, 0, hex1, 1);
  expect_correct_connected_element_via_side(elemElemGraph, shell4, 1, hex2, 0);
  EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell4, 0));
  EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell4, 1));

  EXPECT_EQ(0u, elemElemGraph.num_parallel_edges());
}

namespace {
class ElemElemGraphTester : public stk::mesh::ElemElemGraph
{
public:
  ElemElemGraphTester(stk::mesh::BulkData& bulkData)
    :ElemElemGraph(bulkData) {}
  const stk::mesh::impl::SparseGraph& my_get_coincident_graph() {return m_coincidentGraph; }
};

class ShellMeshModification : public stk::unit_test_util::MeshFixture
{
protected:
  ShellMeshModification()
  {
  }

  ~ShellMeshModification()
  {
    delete elemElemGraph;
  }

  void initialize(stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
    setup_empty_mesh(auraOption);
    shellPart = &get_meta().declare_part_with_topology("shell_part", stk::topology::SHELL_QUAD_4);
  }

  void create_shell2_on_lowest_rank_proc_and_setup_node_sharing()
  {
    get_bulk().modification_begin();
    declare_shell2_on_lowest_rank_proc();
    if (get_bulk().parallel_size()>1)
      share_nodes();
    get_bulk().modification_end();
  }

  void create_shell3_on_highest_rank_proc()
  {
    get_bulk().modification_begin();
    declare_shell3_on_highest_rank_proc();
    get_bulk().modification_end();
  }

  void create_stacked_shells()
  {
    get_bulk().modification_begin();
    declare_shell2_on_lowest_rank_proc();
    declare_shell3_on_highest_rank_proc();
    if (get_bulk().parallel_size()>1)
      share_nodes();
    get_bulk().modification_end();
  }

  void delete_shell2()
  {
    get_bulk().modification_begin();
    destroy_shell2();
    get_bulk().modification_end();
  }

  void create_elem_elem_graph()
  {
    elemElemGraph = new ElemElemGraphTester(get_bulk());
    updater = std::make_shared<stk::mesh::ElemElemGraphUpdater>(get_bulk(), *elemElemGraph);
    get_bulk().register_observer(updater);
    coincident_graph = &elemElemGraph->my_get_coincident_graph();
  }

  void verify_graph_for_stacked_shells()
  {
    size_t expectedNumElementsThisProcessor = 2u;
    size_t expectedNumCoincidentElementsThisProcessor = 2u;
    size_t expectedNumParallelEdges = 0u;
    if (get_bulk().parallel_size()>1) {
      expectedNumElementsThisProcessor = 1u;
      expectedNumCoincidentElementsThisProcessor = 1u;
      expectedNumParallelEdges = 2u;
    }
    verify_graph_has_num_elements_and_num_edges(expectedNumElementsThisProcessor, 0u);
    verify_num_parallel_edges(expectedNumParallelEdges);
    verify_num_elements_and_num_edges_in_coincident_graph(expectedNumCoincidentElementsThisProcessor,2u);
  }

  void verify_graph_for_single_shell3_on_highest_rank_proc()
  {
    size_t numElementsThisProcessor = 1u;
    if (get_bulk().parallel_size()>1 && get_bulk().parallel_rank()==0)
      numElementsThisProcessor = 0u;
    const size_t numCoincidentElementsThisProcessor = 0u;
    verify_graph_has_num_elements_and_num_edges(numElementsThisProcessor, 0u);
    verify_num_parallel_edges(0u);
    verify_num_elements_and_num_edges_in_coincident_graph(numCoincidentElementsThisProcessor,2u);
  }

  void test_create_stacked_shells_then_delete_one(stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
    initialize(auraOption);
    create_stacked_shells();
    create_elem_elem_graph();
    //        write_graph("-------------- Graph after creating stacked shells --------------\n");
    verify_graph_for_stacked_shells();
    delete_shell2();
    //        write_graph("-------------- Graph after deleting one shell --------------\n");
    verify_graph_for_single_shell3_on_highest_rank_proc();
  }

  void test_create_shell_then_create_another_one(stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
    initialize(auraOption);
    create_shell3_on_highest_rank_proc();
    create_elem_elem_graph();
    //        write_graph("-------------- Graph after creating one shell --------------\n");
    verify_graph_for_single_shell3_on_highest_rank_proc();
    create_shell2_on_lowest_rank_proc_and_setup_node_sharing();
    //        write_graph("-------------- Graph after creating second (stacked) shell --------------\n");
    verify_graph_for_stacked_shells();
  }

  void create_and_delete_shell_in_same_mod_cycle()
  {
    get_bulk().modification_begin();
    declare_shell2_on_lowest_rank_proc();
    destroy_shell2();
    get_bulk().modification_end();
  }

  void verify_empty_graph()
  {
    verify_graph_has_num_elements_and_num_edges(0u,0u);
    verify_num_elements_and_num_edges_in_coincident_graph(0u,0u);
  }

private:
  void verify_local_and_parallel_edges_from_entity_to_remote_entity(stk::mesh::Entity entity, stk::mesh::EntityId remoteEntityId, size_t numParallelEdges)
  {
    if(get_bulk().is_valid(entity) && get_bulk().bucket(entity).owned())
    {
      EXPECT_EQ(0u, elemElemGraph->get_num_connected_elems(entity));
    }
  }

  void verify_num_parallel_edges(size_t numParallelEdges)
  {
    verify_local_and_parallel_edges_from_entity_to_remote_entity(shell2, 3u, numParallelEdges);
    verify_local_and_parallel_edges_from_entity_to_remote_entity(shell3, 2u, numParallelEdges);
  }

  void verify_num_elements_and_num_edges_in_coincident_graph(size_t numCoincidentElements, size_t numEdgesToCoincidentElements)
  {
    EXPECT_EQ(numCoincidentElements, coincident_graph->get_num_elements_in_graph());
    for(const stk::mesh::impl::SparseGraph::value_type& map_iter : *coincident_graph)
    {
      const std::vector<stk::mesh::GraphEdge>& coincident_edges = map_iter.second;
      EXPECT_EQ(numEdgesToCoincidentElements, coincident_edges.size());
    }
  }

  void verify_graph_has_num_elements_and_num_edges(size_t numElementsInGraph, size_t numEdgesInGraph)
  {
    EXPECT_EQ(numElementsInGraph, elemElemGraph->size());
    EXPECT_EQ(numEdgesInGraph, elemElemGraph->num_edges());
  }

  void share_nodes()
  {
    int other_proc = get_bulk().parallel_rank() == 0 ? 1 : 0;
    get_bulk().add_node_sharing(get_bulk().get_entity(stk::topology::NODE_RANK,shellNodeIDs[0][0]),other_proc);
    get_bulk().add_node_sharing(get_bulk().get_entity(stk::topology::NODE_RANK,shellNodeIDs[0][1]),other_proc);
    get_bulk().add_node_sharing(get_bulk().get_entity(stk::topology::NODE_RANK,shellNodeIDs[0][2]),other_proc);
    get_bulk().add_node_sharing(get_bulk().get_entity(stk::topology::NODE_RANK,shellNodeIDs[0][3]),other_proc);
  }

  void declare_shell2_on_lowest_rank_proc()
  {
    if (get_bulk().parallel_size()==1 || get_bulk().parallel_rank()==0 )
      shell2 = stk::mesh::declare_element(get_bulk(), *shellPart, shellElemIDs[0], shellNodeIDs[0]);
  }

  void declare_shell3_on_highest_rank_proc()
  {
    if (get_bulk().parallel_size()==1 || get_bulk().parallel_rank()==1 )
      shell3 = stk::mesh::declare_element(get_bulk(), *shellPart, shellElemIDs[1], shellNodeIDs[1]);
  }

  void destroy_shell2()
  {
    if (get_bulk().parallel_size()==1 || get_bulk().parallel_rank()==0 )
      get_bulk().destroy_entity(shell2);
  }


private:
  stk::mesh::Part* shellPart = nullptr;
  std::vector<stk::mesh::EntityIdVector> shellNodeIDs {
    { 5, 6, 7, 8 },
    { 5, 6, 7, 8 }
  };
  stk::mesh::EntityId shellElemIDs[2] = { 2, 3 };
  stk::mesh::Entity shell2;
  stk::mesh::Entity shell3;
  ElemElemGraphTester* elemElemGraph = nullptr;
  std::shared_ptr<stk::mesh::ElemElemGraphUpdater> updater;
  const stk::mesh::impl::SparseGraph* coincident_graph = nullptr;
};
} // namespace


TEST_F(ShellMeshModification, CreateStackedShellsThenTestDeleteOne)
{
  if (stk::parallel_machine_size(get_comm()) == 1)
  {
    test_create_stacked_shells_then_delete_one(stk::mesh::BulkData::AUTO_AURA);
  }
}

// disabled due to split coincident elements
TEST_F( ShellMeshModification, DISABLED_CreateShellThenTestCreateAnotherShellAutoAura)
{
  if (stk::parallel_machine_size(get_comm()) <= 2)
  {
    test_create_shell_then_create_another_one(stk::mesh::BulkData::AUTO_AURA);
  }
}

// disabled due to split coincident elements
TEST_F( ShellMeshModification, DISABLED_CreateShellThenTestCreateAnotherShellNoAura)
{
  if (stk::parallel_machine_size(get_comm()) == 2)
  {
    test_create_shell_then_create_another_one(stk::mesh::BulkData::NO_AUTO_AURA);
  }
}

TEST_F( ShellMeshModification, CreateAndDeleteShellInSameModCycle)
{
  if (stk::parallel_machine_size(get_comm()) == 1)
  {
    initialize(stk::mesh::BulkData::NO_AUTO_AURA);
    create_elem_elem_graph();
    EXPECT_NO_THROW(create_and_delete_shell_in_same_mod_cycle());
    verify_empty_graph();
  }
}


TEST( ElementGraph, HexShellShellSerial )
{
  //  ID.proc
  //
  //          3.0------------7.0
  //          /|             /|
  //         / |            / |
  //        /  |           /  |
  //      4.0------------8.0  |
  //       |   |          |   |
  //       |   |   1.0    |2.0|
  //       |   |          |3.0|
  //       |  2.0---------|--6.0
  //       |  /           |  /
  //       | /            | /
  //       |/             |/
  //      1.0------------5.0
  //                      ^
  //                      |
  //                       ---- Two stacked shell elements

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  int p_size = stk::parallel_machine_size(pm);

  if(p_size > 1)
  {
    return;
  }

  const unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, pm);
  stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();
  stk::mesh::BulkData& mesh = *bulkPtr;

  stk::mesh::Part * hexPart   = &meta.declare_part_with_topology("hex_part",   stk::topology::HEX_8);
  stk::mesh::Part * shellPart = &meta.declare_part_with_topology("shell_part", stk::topology::SHELL_QUAD_4);
  meta.commit();

  std::vector<stk::mesh::EntityIdVector> hexNodeIDs {
    { 1, 2, 3, 4, 5,  6,  7,  8 }
  };
  stk::mesh::EntityId hexElemIDs[] = { 1 };

  std::vector<stk::mesh::EntityIdVector> shellNodeIDs {
    { 5, 6, 7, 8 },
    { 5, 6, 7, 8 }
  };
  stk::mesh::EntityId shellElemIDs[] = { 2, 3 };

  mesh.modification_begin();
  for (size_t i = 0; i < hexNodeIDs.size(); ++i) {
    stk::mesh::declare_element(mesh, *hexPart, hexElemIDs[i], hexNodeIDs[i]);
  }
  for (size_t i = 0; i < shellNodeIDs.size(); ++i) {
    stk::mesh::declare_element(mesh, *shellPart, shellElemIDs[i], shellNodeIDs[i]);
  }
  mesh.modification_end();

  stk::mesh::ElemElemGraph elemElemGraph(mesh);

  const stk::mesh::Entity hex1   = mesh.get_entity(stk::topology::ELEM_RANK, 1);
  const stk::mesh::Entity shell2 = mesh.get_entity(stk::topology::ELEM_RANK, 2);
  const stk::mesh::Entity shell3 = mesh.get_entity(stk::topology::ELEM_RANK, 3);

  // Connectivity for Hex Element 1
  EXPECT_EQ(2u,     elemElemGraph.get_num_connected_elems(hex1));
  expect_correct_connected_element_via_side(elemElemGraph, hex1, 0, shell2, 5);
  expect_correct_connected_element_via_side(elemElemGraph, hex1, 1, shell3, 5);
  EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));
  EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 1));

  // Connectivity for Shell Element 2
  EXPECT_EQ(1u,   elemElemGraph.get_num_connected_elems(shell2));
  expect_correct_connected_element_via_side(elemElemGraph, shell2, 0, hex1, 1);
  EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell2, 0));

  // Connectivity for Shell Element 3
  EXPECT_EQ(1u,   elemElemGraph.get_num_connected_elems(shell3));
  expect_correct_connected_element_via_side(elemElemGraph, shell3, 0, hex1, 1);
  EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell3, 0));

  EXPECT_EQ(4u, elemElemGraph.num_edges());
  EXPECT_EQ(0u, elemElemGraph.num_parallel_edges());
}

TEST( ElementGraph, HexShellShellHexSerial )
{
  //  ID.proc
  //
  //          3.0------------7.0-----------11.0
  //          /|             /|             /|
  //         / |            / |            / |
  //        /  |           /  |           /  |
  //      4.0------------8.0-----------12.0  |
  //       |   |          |   |          |   |
  //       |   |   1.0    |3.0|   2.0    |   |
  //       |   |          |4.0|          |   |
  //       |  2.0---------|--6.0---------|-10.0
  //       |  /           |  /           |  /
  //       | /            | /            | /
  //       |/             |/             |/
  //      1.0------------5.0------------9.0
  //                      ^
  //                      |
  //                       ---- Two stacked shell elements

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  int p_size = stk::parallel_machine_size(pm);

  if(p_size > 1)
  {
    return;
  }

  const unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, pm);
  stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();
  stk::mesh::BulkData& mesh = *bulkPtr;

  stk::mesh::Part * hexPart   = &meta.declare_part_with_topology("hex_part",   stk::topology::HEX_8);
  stk::mesh::Part * shellPart = &meta.declare_part_with_topology("shell_part", stk::topology::SHELL_QUAD_4);
  meta.commit();

  std::vector<stk::mesh::EntityIdVector> hexNodeIDs {
    { 1, 2, 3, 4, 5,  6,  7,  8 },
    { 5, 6, 7, 8, 9, 10, 11, 12 }
  };
  stk::mesh::EntityId hexElemIDs[] = { 1, 2 };

  std::vector<stk::mesh::EntityIdVector> shellNodeIDs {
    { 5, 6, 7, 8 },
    { 5, 6, 7, 8 }
  };
  stk::mesh::EntityId shellElemIDs[] = { 3, 4 };

  mesh.modification_begin();
  for (size_t i = 0; i < hexNodeIDs.size(); ++i) {
    stk::mesh::declare_element(mesh, *hexPart, hexElemIDs[i], hexNodeIDs[i]);
  }
  for (size_t i = 0; i < shellNodeIDs.size(); ++i) {
    stk::mesh::declare_element(mesh, *shellPart, shellElemIDs[i], shellNodeIDs[i]);
  }
  mesh.modification_end();

  stk::mesh::ElemElemGraph elemElemGraph(mesh);

  const stk::mesh::Entity hex1   = mesh.get_entity(stk::topology::ELEM_RANK, 1);
  const stk::mesh::Entity hex2   = mesh.get_entity(stk::topology::ELEM_RANK, 2);
  const stk::mesh::Entity shell3 = mesh.get_entity(stk::topology::ELEM_RANK, 3);
  const stk::mesh::Entity shell4 = mesh.get_entity(stk::topology::ELEM_RANK, 4);

  // Connectivity for Hex Element 1
  EXPECT_EQ(2u,     elemElemGraph.get_num_connected_elems(hex1));
  expect_correct_connected_element_via_side(elemElemGraph, hex1, 0, shell3, 5);
  expect_correct_connected_element_via_side(elemElemGraph, hex1, 1, shell4, 5);
  EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));
  EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 1));

  // Connectivity for Hex Element 2
  EXPECT_EQ(2u,     elemElemGraph.get_num_connected_elems(hex2));
  expect_correct_connected_element_via_side(elemElemGraph, hex2, 0, shell3, 4);
  expect_correct_connected_element_via_side(elemElemGraph, hex2, 1, shell4, 4);
  EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex2, 0));
  EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex2, 1));

  // Connectivity for Shell Element 3
  EXPECT_EQ(2u,   elemElemGraph.get_num_connected_elems(shell3));
  expect_correct_connected_element_via_side(elemElemGraph, shell3, 0, hex1, 1);
  expect_correct_connected_element_via_side(elemElemGraph, shell3, 1, hex2, 0);
  EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell3, 0));
  EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell3, 1));

  // Connectivity for Shell Element 4
  EXPECT_EQ(2u,   elemElemGraph.get_num_connected_elems(shell4));
  expect_correct_connected_element_via_side(elemElemGraph, shell4, 0, hex1, 1);
  expect_correct_connected_element_via_side(elemElemGraph, shell4, 1, hex2, 0);
  EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell4, 0));
  EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell4, 1));

  EXPECT_EQ(8u, elemElemGraph.num_edges());
  EXPECT_EQ(0u, elemElemGraph.num_parallel_edges());
}

TEST( ElementGraph, HexShellReversedShellHexSerial )
{
  //  ID.proc
  //
  //          3.0------------7.0-----------11.0
  //          /|             /|             /|
  //         / |            / |            / |
  //        /  |           /  |           /  |
  //      4.0------------8.0-----------12.0  |
  //       |   |          |   |          |   |
  //       |   |   1.0    |3.0|   2.0    |   |
  //       |   |          |4.0|          |   |
  //       |  2.0---------|--6.0---------|-10.0
  //       |  /           |  /           |  /
  //       | /            | /            | /
  //       |/             |/             |/
  //      1.0------------5.0------------9.0
  //                      ^
  //                      |
  //                       ---- Two stacked shell elements

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  int p_size = stk::parallel_machine_size(pm);

  if(p_size > 1)
  {
    return;
  }

  const unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, pm);
  stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();
  stk::mesh::BulkData& mesh = *bulkPtr;

  stk::mesh::Part * hexPart   = &meta.declare_part_with_topology("hex_part",   stk::topology::HEX_8);
  stk::mesh::Part * shellPart = &meta.declare_part_with_topology("shell_part", stk::topology::SHELL_QUAD_4);
  meta.commit();

  std::vector<stk::mesh::EntityIdVector> hexNodeIDs {
    { 1, 2, 3, 4, 5,  6,  7,  8 },
    { 5, 6, 7, 8, 9, 10, 11, 12 }
  };
  stk::mesh::EntityId hexElemIDs[] = { 1, 2 };

  std::vector<stk::mesh::EntityIdVector> shellNodeIDs {
    { 5, 6, 7, 8 },
    { 5, 8, 7, 6 }
  };
  stk::mesh::EntityId shellElemIDs[] = { 3, 4 };

  mesh.modification_begin();
  for (size_t i = 0; i < hexNodeIDs.size(); ++i) {
    stk::mesh::declare_element(mesh, *hexPart, hexElemIDs[i], hexNodeIDs[i]);
  }
  for (size_t i = 0; i < shellNodeIDs.size(); ++i) {
    stk::mesh::declare_element(mesh, *shellPart, shellElemIDs[i], shellNodeIDs[i]);
  }
  mesh.modification_end();

  stk::mesh::ElemElemGraph elemElemGraph(mesh);

  const stk::mesh::Entity hex1   = mesh.get_entity(stk::topology::ELEM_RANK, 1);
  const stk::mesh::Entity hex2   = mesh.get_entity(stk::topology::ELEM_RANK, 2);
  const stk::mesh::Entity shell3 = mesh.get_entity(stk::topology::ELEM_RANK, 3);
  const stk::mesh::Entity shell4 = mesh.get_entity(stk::topology::ELEM_RANK, 4);

  // Connectivity for Hex Element 1
  EXPECT_EQ(2u,     elemElemGraph.get_num_connected_elems(hex1));
  expect_correct_connected_element_via_side(elemElemGraph, hex1, 0, shell3, 5);
  expect_correct_connected_element_via_side(elemElemGraph, hex1, 1, shell4, 5);
  EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));
  EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 1));

  // Connectivity for Hex Element 2
  EXPECT_EQ(2u,     elemElemGraph.get_num_connected_elems(hex2));
  expect_correct_connected_element_via_side(elemElemGraph, hex2, 0, shell3, 4);
  expect_correct_connected_element_via_side(elemElemGraph, hex2, 1, shell4, 4);
  EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex2, 0));
  EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex2, 1));

  // Connectivity for Shell Element 3
  EXPECT_EQ(2u,   elemElemGraph.get_num_connected_elems(shell3));
  expect_correct_connected_element_via_side(elemElemGraph, shell3, 0, hex1, 1);
  expect_correct_connected_element_via_side(elemElemGraph, shell3, 1, hex2, 0);
  EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell3, 0));
  EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell3, 1));

  // Connectivity for Shell Element 4
  EXPECT_EQ(2u,   elemElemGraph.get_num_connected_elems(shell4));
  expect_correct_connected_element_via_side(elemElemGraph, shell4, 0, hex1, 0);
  expect_correct_connected_element_via_side(elemElemGraph, shell4, 1, hex2, 1);
  EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell4, 0));
  EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell4, 1));

  EXPECT_EQ(8u, elemElemGraph.num_edges());
  EXPECT_EQ(0u, elemElemGraph.num_parallel_edges());
}

TEST( ElementGraph, Hex0Shell0Shell0Hex1Parallel )
{
  //  ID.proc
  //
  //          3.0------------7.0-----------11.1
  //          /|             /|             /|
  //         / |            / |            / |
  //        /  |           /  |           /  |
  //      4.0------------8.0-----------12.1  |
  //       |   |          |   |          |   |
  //       |   |   1.0    |3.0|   2.1    |   |
  //       |   |          |4.0|          |   |
  //       |  2.0---------|--6.0---------|-10.1
  //       |  /           |  /           |  /
  //       | /            | /            | /
  //       |/             |/             |/
  //      1.0------------5.0------------9.1
  //                      ^
  //                      |
  //                       ---- Two stacked shells

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  unsigned p_size = stk::parallel_machine_size(pm);
  unsigned p_rank = stk::parallel_machine_rank(pm);

  if(p_size != 2u)
  {
    return;
  }

  const unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, pm);
  stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();
  stk::mesh::BulkData& mesh = *bulkPtr;

  stk::mesh::Part * hexPart   = &meta.declare_part_with_topology("hex_part",   stk::topology::HEX_8);
  stk::mesh::Part * shellPart = &meta.declare_part_with_topology("shell_part", stk::topology::SHELL_QUAD_4);
  meta.commit();

  std::vector<stk::mesh::EntityIdVector> hexNodeIDs {
    { 1, 2, 3, 4, 5,  6,  7,  8 },
    { 5, 6, 7, 8, 9, 10, 11, 12 }
  };
  stk::mesh::EntityId hexElemIDs[] = { 1, 2 };
  stk::mesh::EntityId hexElemOwningProc[] = { 0, 1 };

  std::vector<stk::mesh::EntityIdVector> shellNodeIDs {
    { 5, 6, 7, 8 },
    { 5, 6, 7, 8 }
  };
  stk::mesh::EntityId shellElemIDs[] = { 3, 4 };
  stk::mesh::EntityId shellElemOwningProc[] = { 0, 0 };

  // list of triplets: (owner-proc, shared-nodeID, sharing-proc)
  std::vector< std::vector<unsigned> > shared_nodeIDs_and_procs
  {
    { 0, 5, 1 },  // proc 0
    { 0, 6, 1 },
    { 0, 7, 1 },
    { 0, 8, 1 },
    { 1, 5, 0 },  // proc 1
    { 1, 6, 0 },
    { 1, 7, 0 },
    { 1, 8, 0 }
  };

  mesh.modification_begin();
  for (size_t i = 0; i < hexNodeIDs.size(); ++i) {
    if (hexElemOwningProc[i] == p_rank) {
      stk::mesh::declare_element(mesh, *hexPart, hexElemIDs[i], hexNodeIDs[i]);
    }
  }
  for (size_t i = 0; i < shellNodeIDs.size(); ++i) {
    if (shellElemOwningProc[i] == p_rank) {
      stk::mesh::declare_element(mesh, *shellPart, shellElemIDs[i], shellNodeIDs[i]);
    }
  }
  setup_node_sharing(mesh, shared_nodeIDs_and_procs );
  mesh.modification_end();

  stk::mesh::ElemElemGraph elemElemGraph(mesh);

  const stk::mesh::Entity hex1   = mesh.get_entity(stk::topology::ELEM_RANK, 1);
  const stk::mesh::Entity hex2   = mesh.get_entity(stk::topology::ELEM_RANK, 2);
  const stk::mesh::Entity shell3 = mesh.get_entity(stk::topology::ELEM_RANK, 3);
  const stk::mesh::Entity shell4 = mesh.get_entity(stk::topology::ELEM_RANK, 4);

  if (p_rank == 0) {
    // Connectivity for Hex Element 1
    EXPECT_EQ(2u,     elemElemGraph.get_num_connected_elems(hex1));
    expect_correct_connected_element_via_side(elemElemGraph, hex1, 0, shell3, 5);
    expect_correct_connected_element_via_side(elemElemGraph, hex1, 1, shell4, 5);
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 1));

    // Connectivity for Shell Element 3
    EXPECT_EQ(2u,   elemElemGraph.get_num_connected_elems(shell3));
    expect_correct_connected_element_via_side(elemElemGraph, shell3, 0, hex1, 1);
    EXPECT_EQ(0,    elemElemGraph.get_connected_remote_id_and_via_side(shell3, 1).side);
    EXPECT_EQ(2u,   elemElemGraph.get_connected_remote_id_and_via_side(shell3, 1).id);
    EXPECT_TRUE (elemElemGraph.is_connected_elem_locally_owned(shell3, 0));
    EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell3, 1));

    // Connectivity for Shell Element 4
    EXPECT_EQ(2u,   elemElemGraph.get_num_connected_elems(shell4));
    expect_correct_connected_element_via_side(elemElemGraph, shell4, 0, hex1, 1);
    EXPECT_EQ(0,    elemElemGraph.get_connected_remote_id_and_via_side(shell4, 1).side);
    EXPECT_EQ(2u,   elemElemGraph.get_connected_remote_id_and_via_side(shell4, 1).id);
    EXPECT_TRUE (elemElemGraph.is_connected_elem_locally_owned(shell4, 0));
    EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell4, 1));

    EXPECT_EQ(6u, elemElemGraph.num_edges());
    EXPECT_EQ(2u, elemElemGraph.num_parallel_edges());
  }
  else if (p_rank == 1) {
    // Connectivity for Hex Element 2
    EXPECT_EQ(2u,     elemElemGraph.get_num_connected_elems(hex2));
    EXPECT_EQ(4,      elemElemGraph.get_connected_remote_id_and_via_side(hex2, 0).side);
    EXPECT_EQ(4,      elemElemGraph.get_connected_remote_id_and_via_side(hex2, 1).side);
    EXPECT_EQ(3u,   elemElemGraph.get_connected_remote_id_and_via_side(hex2, 0).id);
    EXPECT_EQ(4u,   elemElemGraph.get_connected_remote_id_and_via_side(hex2, 1).id);
    EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex2, 0));
    EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex2, 1));

    EXPECT_EQ(2u, elemElemGraph.num_edges());
    EXPECT_EQ(2u, elemElemGraph.num_parallel_edges());
  }
}

// disabled due to split coincident elements
TEST( ElementGraph, DISABLED_Hex0Shell0Shell1Hex1Parallel )
{
  //  ID.proc
  //
  //          3.0------------7.0-----------11.1
  //          /|             /|             /|
  //         / |            / |            / |
  //        /  |           /  |           /  |
  //      4.0------------8.0-----------12.1  |
  //       |   |          |   |          |   |
  //       |   |   1.0    |3.0|   2.1    |   |
  //       |   |          |4.1|          |   |
  //       |  2.0---------|--6.0---------|-10.1
  //       |  /           |  /           |  /
  //       | /            | /            | /
  //       |/             |/             |/
  //      1.0------------5.0------------9.1
  //                      ^
  //                      |
  //                       ---- Two stacked shells

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  unsigned p_size = stk::parallel_machine_size(pm);
  unsigned p_rank = stk::parallel_machine_rank(pm);

  if(p_size != 2u)
  {
    return;
  }

  const unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, pm);
  stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();
  stk::mesh::BulkData& mesh = *bulkPtr;

  stk::mesh::Part * hexPart   = &meta.declare_part_with_topology("hex_part",   stk::topology::HEX_8);
  stk::mesh::Part * shellPart = &meta.declare_part_with_topology("shell_part", stk::topology::SHELL_QUAD_4);
  meta.commit();

  std::vector<stk::mesh::EntityIdVector> hexNodeIDs {
    { 1, 2, 3, 4, 5,  6,  7,  8 },
    { 5, 6, 7, 8, 9, 10, 11, 12 }
  };
  stk::mesh::EntityId hexElemIDs[] = { 1, 2 };
  stk::mesh::EntityId hexElemOwningProc[] = { 0, 1 };

  std::vector<stk::mesh::EntityIdVector> shellNodeIDs {
    { 5, 6, 7, 8 },
    { 5, 6, 7, 8 }
  };
  stk::mesh::EntityId shellElemIDs[] = { 3, 4 };
  stk::mesh::EntityId shellElemOwningProc[] = { 0, 1 };

  // list of triplets: (owner-proc, shared-nodeID, sharing-proc)
  std::vector< std::vector<unsigned> > shared_nodeIDs_and_procs
  {
    { 0, 5, 1 },  // proc 0
    { 0, 6, 1 },
    { 0, 7, 1 },
    { 0, 8, 1 },
    { 1, 5, 0 },  // proc 1
    { 1, 6, 0 },
    { 1, 7, 0 },
    { 1, 8, 0 }
  };

  mesh.modification_begin();
  for (size_t i = 0; i < hexNodeIDs.size(); ++i) {
    if (hexElemOwningProc[i] == p_rank) {
      stk::mesh::declare_element(mesh, *hexPart, hexElemIDs[i], hexNodeIDs[i]);
    }
  }
  for (size_t i = 0; i < shellNodeIDs.size(); ++i) {
    if (shellElemOwningProc[i] == p_rank) {
      stk::mesh::declare_element(mesh, *shellPart, shellElemIDs[i], shellNodeIDs[i]);
    }
  }
  setup_node_sharing(mesh, shared_nodeIDs_and_procs );
  mesh.modification_end();

  stk::mesh::ElemElemGraph elemElemGraph(mesh);

  const stk::mesh::Entity hex1   = mesh.get_entity(stk::topology::ELEM_RANK, 1);
  const stk::mesh::Entity hex2   = mesh.get_entity(stk::topology::ELEM_RANK, 2);
  const stk::mesh::Entity shell3 = mesh.get_entity(stk::topology::ELEM_RANK, 3);
  const stk::mesh::Entity shell4 = mesh.get_entity(stk::topology::ELEM_RANK, 4);

  if (p_rank == 0) {
    // Connectivity for Hex Element 1
    EXPECT_EQ(2u,     elemElemGraph.get_num_connected_elems(hex1));
    expect_correct_connected_element_via_side(elemElemGraph, hex1, 0, shell3, 5);
    EXPECT_EQ(5,      elemElemGraph.get_connected_remote_id_and_via_side(hex1, 1).side);
    EXPECT_EQ(4u,     elemElemGraph.get_connected_remote_id_and_via_side(hex1, 1).id);
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));
    EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex1, 1));

    // Connectivity for Shell Element 3
    EXPECT_EQ(2u,   elemElemGraph.get_num_connected_elems(shell3));
    expect_correct_connected_element_via_side(elemElemGraph, shell3, 0, hex1, 1);
    EXPECT_EQ(0,    elemElemGraph.get_connected_remote_id_and_via_side(shell3, 1).side);
    EXPECT_EQ(2u,   elemElemGraph.get_connected_remote_id_and_via_side(shell3, 1).id);
    EXPECT_TRUE (elemElemGraph.is_connected_elem_locally_owned(shell3, 0));
    EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell3, 1));
  }
  else if (p_rank == 1) {
    // Connectivity for Shell Element 4
    EXPECT_EQ(2u,   elemElemGraph.get_num_connected_elems(shell4));
    expect_correct_connected_element_via_side(elemElemGraph, shell4, 0, hex2, 0);
    EXPECT_EQ(1,    elemElemGraph.get_connected_remote_id_and_via_side(shell4, 1).side);
    EXPECT_EQ(1u,   elemElemGraph.get_connected_remote_id_and_via_side(shell4, 1).id);
    EXPECT_TRUE (elemElemGraph.is_connected_elem_locally_owned(shell4, 0));
    EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell4, 1));

    // Connectivity for Hex Element 2
    EXPECT_EQ(2u,     elemElemGraph.get_num_connected_elems(hex2));
    expect_correct_connected_element_via_side(elemElemGraph, hex2, 0, shell4, 4);
    EXPECT_EQ(4,      elemElemGraph.get_connected_remote_id_and_via_side(hex2, 1).side);
    EXPECT_EQ(3u,     elemElemGraph.get_connected_remote_id_and_via_side(hex2, 1).id);
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex2, 0));
    EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex2, 1));
  }

  EXPECT_EQ(4u, elemElemGraph.num_edges());
  EXPECT_EQ(2u, elemElemGraph.num_parallel_edges());
}

TEST( ElementGraph, Hex0Shell0ReversedShell0Hex1Parallel )
{
  //  ID.proc
  //
  //          3.0------------7.0-----------11.1
  //          /|             /|             /|
  //         / |            / |            / |
  //        /  |           /  |           /  |
  //      4.0------------8.0-----------12.1  |
  //       |   |          |   |          |   |
  //       |   |   1.0    |3.0|   2.1    |   |
  //       |   |          |4.0|          |   |
  //       |  2.0---------|--6.0---------|-10.1
  //       |  /           |  /           |  /
  //       | /            | /            | /
  //       |/             |/             |/
  //      1.0------------5.0------------9.1
  //                      ^
  //                      |
  //                       ---- Two stacked shells, opposite orientation

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  unsigned p_size = stk::parallel_machine_size(pm);
  unsigned p_rank = stk::parallel_machine_rank(pm);

  if(p_size != 2u)
  {
    return;
  }

  const unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, pm);
  stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();
  stk::mesh::BulkData& mesh = *bulkPtr;

  stk::mesh::Part * hexPart   = &meta.declare_part_with_topology("hex_part",   stk::topology::HEX_8);
  stk::mesh::Part * shellPart = &meta.declare_part_with_topology("shell_part", stk::topology::SHELL_QUAD_4);
  meta.commit();

  std::vector<stk::mesh::EntityIdVector> hexNodeIDs {
    { 1, 2, 3, 4, 5,  6,  7,  8 },
    { 5, 6, 7, 8, 9, 10, 11, 12 }
  };
  stk::mesh::EntityId hexElemIDs[] = { 1, 2 };
  stk::mesh::EntityId hexElemOwningProc[] = { 0, 1 };

  std::vector<stk::mesh::EntityIdVector> shellNodeIDs {
    { 5, 6, 7, 8 },
    { 5, 8, 7, 6 }
  };
  stk::mesh::EntityId shellElemIDs[] = { 3, 4 };
  stk::mesh::EntityId shellElemOwningProc[] = { 0, 0 };

  // list of triplets: (owner-proc, shared-nodeID, sharing-proc)
  std::vector< std::vector<unsigned> > shared_nodeIDs_and_procs
  {
    { 0, 5, 1 },  // proc 0
    { 0, 6, 1 },
    { 0, 7, 1 },
    { 0, 8, 1 },
    { 1, 5, 0 },  // proc 1
    { 1, 6, 0 },
    { 1, 7, 0 },
    { 1, 8, 0 }
  };

  mesh.modification_begin();
  for (size_t i = 0; i < hexNodeIDs.size(); ++i) {
    if (hexElemOwningProc[i] == p_rank) {
      stk::mesh::declare_element(mesh, *hexPart, hexElemIDs[i], hexNodeIDs[i]);
    }
  }
  for (size_t i = 0; i < shellNodeIDs.size(); ++i) {
    if (shellElemOwningProc[i] == p_rank) {
      stk::mesh::declare_element(mesh, *shellPart, shellElemIDs[i], shellNodeIDs[i]);
    }
  }
  setup_node_sharing(mesh, shared_nodeIDs_and_procs );
  mesh.modification_end();

  stk::mesh::ElemElemGraph elemElemGraph(mesh);

  const stk::mesh::Entity hex1   = mesh.get_entity(stk::topology::ELEM_RANK, 1);
  const stk::mesh::Entity hex2   = mesh.get_entity(stk::topology::ELEM_RANK, 2);
  const stk::mesh::Entity shell3 = mesh.get_entity(stk::topology::ELEM_RANK, 3);
  const stk::mesh::Entity shell4 = mesh.get_entity(stk::topology::ELEM_RANK, 4);

  if (p_rank == 0) {
    // Connectivity for Hex Element 1
    EXPECT_EQ(2u,     elemElemGraph.get_num_connected_elems(hex1));
    expect_correct_connected_element_via_side(elemElemGraph, hex1, 0, shell3, 5);
    expect_correct_connected_element_via_side(elemElemGraph, hex1, 1, shell4, 5);
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 1));

    // Connectivity for Shell Element 3
    EXPECT_EQ(2u,   elemElemGraph.get_num_connected_elems(shell3));
    expect_correct_connected_element_via_side(elemElemGraph, shell3, 0, hex1, 1);
    EXPECT_EQ(0,    elemElemGraph.get_connected_remote_id_and_via_side(shell3, 1).side);
    EXPECT_EQ(2u,   elemElemGraph.get_connected_remote_id_and_via_side(shell3, 1).id);
    EXPECT_TRUE (elemElemGraph.is_connected_elem_locally_owned(shell3, 0));
    EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell3, 1));

    // Connectivity for Shell Element 4
    EXPECT_EQ(2u,   elemElemGraph.get_num_connected_elems(shell4));
    expect_correct_connected_element_via_side(elemElemGraph, shell4, 0, hex1, 0);
    EXPECT_EQ(1,    elemElemGraph.get_connected_remote_id_and_via_side(shell4, 1).side);
    EXPECT_EQ(2u,   elemElemGraph.get_connected_remote_id_and_via_side(shell4, 1).id);
    EXPECT_TRUE (elemElemGraph.is_connected_elem_locally_owned(shell4, 0));
    EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell4, 1));
    EXPECT_EQ(6u, elemElemGraph.num_edges());
    EXPECT_EQ(2u, elemElemGraph.num_parallel_edges());
  }
  else if (p_rank == 1) {
    // Connectivity for Hex Element 2
    EXPECT_EQ(2u, elemElemGraph.get_num_connected_elems(hex2));
    EXPECT_EQ(3u,  elemElemGraph.get_connected_remote_id_and_via_side(hex2, 0).id);
    EXPECT_EQ(4,  elemElemGraph.get_connected_remote_id_and_via_side(hex2, 0).side);
    EXPECT_EQ(4u,  elemElemGraph.get_connected_remote_id_and_via_side(hex2, 1).id);
    EXPECT_EQ(4,  elemElemGraph.get_connected_remote_id_and_via_side(hex2, 1).side);
    EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex2, 0));
    EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex2, 1));
    EXPECT_EQ(2u, elemElemGraph.num_edges());
    EXPECT_EQ(2u, elemElemGraph.num_parallel_edges());
  }
}

TEST( ElementGraph, Hex1Shell0Shell0Hex1Parallel )
{
  //  ID.proc
  //
  //          3.0------------7.1-----------11.1
  //          /|             /|             /|
  //         / |            / |            / |
  //        /  |           /  |           /  |
  //      4.0------------8.1-----------12.1  |
  //       |   |          |   |          |   |
  //       |   |   1.1    |3.0|   2.1    |   |
  //       |   |          |4.0|          |   |
  //       |  2.0---------|--6.1---------|-10.1
  //       |  /           |  /           |  /
  //       | /            | /            | /
  //       |/             |/             |/
  //      1.0------------5.1------------9.1
  //                      ^
  //                      |
  //                       ---- Two stacked shells

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  unsigned p_size = stk::parallel_machine_size(pm);
  unsigned p_rank = stk::parallel_machine_rank(pm);

  if(p_size != 2u)
  {
    return;
  }

  const unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, pm);
  stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();
  stk::mesh::BulkData& mesh = *bulkPtr;

  stk::mesh::Part * hexPart   = &meta.declare_part_with_topology("hex_part",   stk::topology::HEX_8);
  stk::mesh::Part * shellPart = &meta.declare_part_with_topology("shell_part", stk::topology::SHELL_QUAD_4);
  meta.commit();

  std::vector<stk::mesh::EntityIdVector> hexNodeIDs {
    { 1, 2, 3, 4, 5,  6,  7,  8 },
    { 5, 6, 7, 8, 9, 10, 11, 12 }
  };
  stk::mesh::EntityId hexElemIDs[] = { 1, 2 };
  stk::mesh::EntityId hexElemOwningProc[] = { 1, 1 };

  std::vector<stk::mesh::EntityIdVector> shellNodeIDs {
    { 5, 6, 7, 8 },
    { 5, 6, 7, 8 }
  };
  stk::mesh::EntityId shellElemIDs[] = { 3, 4 };
  stk::mesh::EntityId shellElemOwningProc[] = { 0, 0 };

  // list of triplets: (owner-proc, shared-nodeID, sharing-proc)
  std::vector< std::vector<unsigned> > shared_nodeIDs_and_procs
  {
    { 0, 5, 1 },  // proc 0
    { 0, 6, 1 },
    { 0, 7, 1 },
    { 0, 8, 1 },
    { 1, 5, 0 },  // proc 1
    { 1, 6, 0 },
    { 1, 7, 0 },
    { 1, 8, 0 }
  };

  mesh.modification_begin();
  for (size_t i = 0; i < hexNodeIDs.size(); ++i) {
    if (hexElemOwningProc[i] == p_rank) {
      stk::mesh::declare_element(mesh, *hexPart, hexElemIDs[i], hexNodeIDs[i]);
    }
  }
  for (size_t i = 0; i < shellNodeIDs.size(); ++i) {
    if (shellElemOwningProc[i] == p_rank) {
      stk::mesh::declare_element(mesh, *shellPart, shellElemIDs[i], shellNodeIDs[i]);
    }
  }
  setup_node_sharing(mesh, shared_nodeIDs_and_procs );
  mesh.modification_end();

  stk::mesh::ElemElemGraph elemElemGraph(mesh);

  const stk::mesh::Entity hex1   = mesh.get_entity(stk::topology::ELEM_RANK, 1);
  const stk::mesh::Entity hex2   = mesh.get_entity(stk::topology::ELEM_RANK, 2);
  const stk::mesh::Entity shell3 = mesh.get_entity(stk::topology::ELEM_RANK, 3);
  const stk::mesh::Entity shell4 = mesh.get_entity(stk::topology::ELEM_RANK, 4);

  if (p_rank == 0) {
    // Connectivity for Shell Element 3
    size_t hex1Index = 0;
    size_t hex2Index = 1;
    EXPECT_EQ(2u, elemElemGraph.get_num_connected_elems(shell3));
    EXPECT_EQ(0,  elemElemGraph.get_connected_remote_id_and_via_side(shell3, hex2Index).side);
    EXPECT_EQ(1,  elemElemGraph.get_connected_remote_id_and_via_side(shell3, hex1Index).side);
    EXPECT_EQ(2u, elemElemGraph.get_connected_remote_id_and_via_side(shell3, hex2Index).id);
    EXPECT_EQ(1u, elemElemGraph.get_connected_remote_id_and_via_side(shell3, hex1Index).id);
    EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell3, hex2Index));
    EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell3, hex1Index));

    // Connectivity for Shell Element 4
    EXPECT_EQ(2u, elemElemGraph.get_num_connected_elems(shell4));
    EXPECT_EQ(0,  elemElemGraph.get_connected_remote_id_and_via_side(shell4, hex2Index).side);
    EXPECT_EQ(1,  elemElemGraph.get_connected_remote_id_and_via_side(shell4, hex1Index).side);
    EXPECT_EQ(2u, elemElemGraph.get_connected_remote_id_and_via_side(shell4, hex2Index).id);
    EXPECT_EQ(1u, elemElemGraph.get_connected_remote_id_and_via_side(shell4, hex1Index).id);
    EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell4, hex2Index));
    EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell4, hex1Index));
  }
  else if (p_rank == 1) {
    size_t shell3Index = 0;
    size_t shell4Index = 1;
    // Connectivity for Hex Element 1
    EXPECT_EQ(2u, elemElemGraph.get_num_connected_elems(hex1));
    EXPECT_EQ(5,  elemElemGraph.get_connected_remote_id_and_via_side(hex1, shell3Index).side);
    EXPECT_EQ(5,  elemElemGraph.get_connected_remote_id_and_via_side(hex1, shell4Index).side);
    EXPECT_EQ(3u, elemElemGraph.get_connected_remote_id_and_via_side(hex1, shell3Index).id);
    EXPECT_EQ(4u, elemElemGraph.get_connected_remote_id_and_via_side(hex1, shell4Index).id);
    EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex1, shell3Index));
    EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex1, shell4Index));

    // Connectivity for Hex Element 2
    EXPECT_EQ(2u, elemElemGraph.get_num_connected_elems(hex2));
    EXPECT_EQ(4,  elemElemGraph.get_connected_remote_id_and_via_side(hex2, shell3Index).side);
    EXPECT_EQ(4,  elemElemGraph.get_connected_remote_id_and_via_side(hex2, shell4Index).side);
    EXPECT_EQ(3u, elemElemGraph.get_connected_remote_id_and_via_side(hex2, shell3Index).id);
    EXPECT_EQ(4u, elemElemGraph.get_connected_remote_id_and_via_side(hex2, shell4Index).id);
    EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex2, shell3Index));
    EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex2, shell4Index));
  }

  EXPECT_EQ(4u, elemElemGraph.num_edges());
  EXPECT_EQ(4u, elemElemGraph.num_parallel_edges());
}

}
