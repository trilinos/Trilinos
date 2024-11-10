#include <gtest/gtest.h>

#include <vector>
#include <algorithm>
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

#include <stk_mesh/base/SkinMesh.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_mesh/base/CreateFaces.hpp>
#include <stk_mesh/base/SkinMeshUtil.hpp>
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp>
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraphImpl.hpp>
#include <stk_mesh/baseImpl/elementGraph/ProcessKilledElements.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelVectorConcat.hpp>
#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/environment/memory_util.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include <stk_util/util/ReportHandler.hpp>

#include <stk_io/IossBridge.hpp>
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker

#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_unit_test_utils/getOption.h>

#include "stk_unit_test_utils/ElemGraphTestUtils.hpp"
#include "stk_unit_test_utils/unittestMeshUtils.hpp"
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/BuildMesh.hpp>

#include "SetupKeyholeMesh.hpp"
#include <stk_unit_test_utils/stk_mesh_fixtures/QuadFixture.hpp>  // for QuadFixture
#include <stk_unit_test_utils/stk_mesh_fixtures/degenerate_mesh.hpp>
#include <stk_unit_test_utils/stk_mesh_fixtures/heterogeneous_mesh.hpp>

#include "BulkDataElementGraphTester.hpp"
#include "ElementGraphTester.hpp"

namespace {

using namespace stk::mesh::impl;
using namespace stk::mesh;
using stk::unit_test_util::build_mesh;

// Not to be used with ElemElemGraph or ElemElemGraphTester class.
bool is_valid_graph_element(const impl::ElementGraph &elem_graph, stk::mesh::impl::LocalId elem_id);

// Not to be used with ElemElemGraph or ElemElemGraphTester class.
int check_connectivity(const impl::ElementGraph &elem_graph, const impl::SidesForElementGraph &via_sides,
                       stk::mesh::impl::LocalId element_id1, stk::mesh::impl::LocalId element_id2);

// Not to be used with ElemElemGraph or ElemElemGraphTester class.
int get_side_from_element1_to_element2(const impl::ElementGraph &elem_graph,
                                       const impl::SidesForElementGraph &via_sides,
                                       stk::mesh::impl::LocalId element1_local_id,
                                       stk::mesh::impl::LocalId other_element_id);

ElemElemGraphTester create_base_1x1x4_elem_graph(stk::ParallelMachine &comm, stk::mesh::BulkData &bulk);

void expect_correct_connected_element_via_side(ElemElemGraphTester& elemElemGraph, stk::mesh::Entity elem, int k, stk::mesh::Entity otherElem, int viaSide)
{
  impl::ElementViaSidePair elem_via_side = elemElemGraph.get_connected_element_and_via_side(elem, k);
  EXPECT_EQ(viaSide,      elem_via_side.side);
  EXPECT_EQ(otherElem, elem_via_side.element);
}

int get_side_from_element1_to_element2(const stk::mesh::BulkData &bulkData,
                                       stk::mesh::ElemElemGraph &elem_graph,
                                       stk::mesh::EntityId element1_local_id,
                                       stk::mesh::impl::LocalId other_element_id)
{
  int side = -1;
  stk::mesh::impl::LocalId elemLocal = elem_graph.get_local_element_id(bulkData.get_entity(stk::topology::ELEM_RANK, element1_local_id));
  for(const stk::mesh::GraphEdge &graphEdge : elem_graph.get_edges_for_element(elemLocal))
    if(graphEdge.elem2() == -other_element_id)
      side = graphEdge.side1();
  return side;
}

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

void test_add_elements_to_pre_existing_graph_and_mesh(stk::mesh::BulkData &bulkData)
{
  const int p_size  = bulkData.parallel_size();
  const int p_rank = bulkData.parallel_rank();

  STK_ThrowRequire(2 == p_size);

  stk::mesh::MetaData &meta = bulkData.mesh_meta_data();

  stk::mesh::Part * hexPart = &meta.declare_part_with_topology("hex_part", stk::topology::HEX_8);

  stk::io::fill_mesh("generated:1x1x2", bulkData);

  std::vector<size_t> counts;
  stk::mesh::count_entities(bulkData.mesh_meta_data().locally_owned_part(), bulkData, counts);

  stk::mesh::ElemElemGraph &elem_graph = bulkData.get_face_adjacent_element_graph();

  EXPECT_EQ(1u, elem_graph.size());
  EXPECT_EQ(1u, counts[stk::topology::ELEM_RANK]);

  const std::vector<size_t> numHex{1,1};
  stk::mesh::EntityIdVector hexNodeIDs[] {
    { 12, 11, 18, 17, 8, 7, 16, 15 },
    {  8,  7, 16, 15, 4, 3, 14, 13 }
  };
  stk::mesh::EntityIdVector hexElemIDs {3, 4};

  // list of triplets: (owner-proc, shared-nodeID, sharing-proc)
  std::vector< std::vector<unsigned> > shared_nodeIDs_and_procs
  {
    { 0,  3, 1 },  // proc 0
    { 0,  4, 1 },
    { 0, 11, 1 },
    { 0, 12, 1 },
    { 0, 15, 1 },
    { 0, 16, 1 },
    { 1,  3, 0 },  // proc 1
    { 1,  4, 0 },
    { 1, 11, 0 },
    { 1, 12, 0 },
    { 1, 15, 0 },
    { 1, 16, 0 },
  };

  bulkData.modification_begin();
  stk::mesh::declare_element(bulkData, *hexPart, hexElemIDs[p_rank], hexNodeIDs[p_rank]);
  setup_node_sharing(bulkData, shared_nodeIDs_and_procs );
  bulkData.modification_end();

  stk::mesh::EntityVector elements_to_add;

  if (0 == p_rank)
  {
    EXPECT_EQ(5, get_side_from_element1_to_element2(bulkData, elem_graph, 1, 2));
    EXPECT_EQ(2, get_side_from_element1_to_element2(bulkData, elem_graph, 1, 4));
    EXPECT_EQ(0, get_side_from_element1_to_element2(bulkData, elem_graph, 3, 2));
    EXPECT_EQ(5, get_side_from_element1_to_element2(bulkData, elem_graph, 3, 4));
  }
  else
  {
    EXPECT_EQ(4, get_side_from_element1_to_element2(bulkData, elem_graph, 2, 1));
    EXPECT_EQ(2, get_side_from_element1_to_element2(bulkData, elem_graph, 2, 3));
    EXPECT_EQ(0, get_side_from_element1_to_element2(bulkData, elem_graph, 4, 1));
    EXPECT_EQ(4, get_side_from_element1_to_element2(bulkData, elem_graph, 4, 3));

  }
  EXPECT_EQ(2u, elem_graph.size());
  EXPECT_EQ(4u, elem_graph.num_edges());
  EXPECT_EQ(4u, elem_graph.num_parallel_edges());
}

void test_delete_elements_from_graph(ElemElemGraphTester &elem_graph, std::vector<stk::mesh::EntityId> &ids_to_delete)
{
  std::set<stk::mesh::EntityId> currentElements = {1,2,3,4};
  stk::mesh::BulkData &bulkData = elem_graph.get_bulk_data();

  stk::mesh::impl::DeletedElementInfoVector elements_to_delete;
  for (size_t i = 0; i < ids_to_delete.size(); ++i)
  {
    stk::mesh::EntityId id = ids_to_delete[i];
    currentElements.erase(id);
    stk::mesh::Entity elem = bulkData.get_entity(stk::topology::ELEM_RANK, id);
    if (bulkData.is_valid(elem) && bulkData.bucket(elem).owned())
    {
      elements_to_delete.push_back({elem, id, bulkData.bucket(elem).topology().is_shell()});
    }
  }

  bulkData.modification_begin();
  for (stk::mesh::impl::DeletedElementInfo elem : elements_to_delete)
  {
    bulkData.destroy_entity(elem.entity);
    stk::mesh::Entity elemCheck = bulkData.get_entity(stk::topology::ELEM_RANK, elem.identifier);
    EXPECT_FALSE(bulkData.is_valid(elemCheck));
  }
  bulkData.modification_end();

  elem_graph.delete_elements(elements_to_delete);

  std::vector<size_t> counts;
  stk::mesh::count_entities(bulkData.mesh_meta_data().locally_owned_part(), bulkData, counts);
  unsigned numLocallyOwnedElems = counts[stk::topology::ELEM_RANK];

  EXPECT_EQ(numLocallyOwnedElems, elem_graph.size());

  unsigned numEdges = 0;
  unsigned numParallelEdges = 0;
  for (stk::mesh::EntityId elem_id : currentElements)
  {
    stk::mesh::Entity elem = bulkData.get_entity(stk::topology::ELEM_RANK, elem_id);

    if (!bulkData.is_valid(elem) || !bulkData.bucket(elem).owned())
    {
      continue;
    }

    stk::mesh::EntityId leftNeighbor = elem_id - 1;
    if (currentElements.find(leftNeighbor) != currentElements.end())
    {
      EXPECT_EQ(4, elem_graph.get_side_from_element1_to_element2(elem_id, leftNeighbor));
      ++numEdges;

      stk::mesh::Entity leftElem = bulkData.get_entity(stk::topology::ELEM_RANK, leftNeighbor);
      bool ownedLeftNeighbor = bulkData.is_valid(leftElem) && bulkData.bucket(leftElem).owned();
      if(!ownedLeftNeighbor)
      {
        numParallelEdges++;
      }
    }

    stk::mesh::EntityId rightNeighbor = elem_id + 1;
    if (currentElements.find(rightNeighbor) != currentElements.end())
    {
      EXPECT_EQ(5, elem_graph.get_side_from_element1_to_element2(elem_id, rightNeighbor));
      ++numEdges;

      stk::mesh::Entity rightElem = bulkData.get_entity(stk::topology::ELEM_RANK, rightNeighbor);
      bool ownedRightNeighbor = bulkData.is_valid(rightElem) && bulkData.bucket(rightElem).owned();
      if(!ownedRightNeighbor)
      {
        numParallelEdges++;
      }
    }
  }
  EXPECT_EQ(numEdges, elem_graph.num_edges());
  EXPECT_EQ(numParallelEdges, elem_graph.num_parallel_edges());
}

void test_element_graph_delete_elements_from_graph(stk::mesh::BulkData::AutomaticAuraOption auto_aura_option)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int pSize = stk::parallel_machine_size(comm);

  // Test designed for 1, 2, or 4 procs.
  if (pSize > 4)
  {
    return;
  }

  {
    std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(comm, auto_aura_option);

    ElemElemGraphTester elem_graph = create_base_1x1x4_elem_graph(comm, *bulkPtr);
    std::vector<stk::mesh::EntityId> ids_to_delete = {1, 4};
    test_delete_elements_from_graph(elem_graph, ids_to_delete);
  }

  {
    std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(comm, auto_aura_option);

    ElemElemGraphTester elem_graph = create_base_1x1x4_elem_graph(comm, *bulkPtr);
    std::vector<stk::mesh::EntityId> ids_to_delete = {2, 3};
    test_delete_elements_from_graph(elem_graph, ids_to_delete);
  }

  {
    std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(comm, auto_aura_option);

    ElemElemGraphTester elem_graph = create_base_1x1x4_elem_graph(comm, *bulkPtr);
    std::vector<stk::mesh::EntityId> ids_to_delete = {1, 3};
    test_delete_elements_from_graph(elem_graph, ids_to_delete);
  }
}

TEST(ElementGraph, delete_elements_from_graph_aura_on)
{
  test_element_graph_delete_elements_from_graph(stk::mesh::BulkData::AUTO_AURA);
}

TEST(ElementGraph, delete_elements_from_graph_aura_off)
{
  test_element_graph_delete_elements_from_graph(stk::mesh::BulkData::NO_AUTO_AURA);
}

TEST(ElementGraph, add_and_delete_elements_from_graph_serial)
{
  MPI_Comm comm = MPI_COMM_WORLD;

  if(stk::parallel_machine_size(comm) == 1)
  {
    std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(comm);
    stk::mesh::BulkData& bulkData = *bulkPtr;

    ElemElemGraphTester elem_graph(bulkData);

    EXPECT_EQ(0u, elem_graph.size());
    EXPECT_EQ(0u, elem_graph.num_edges());
    EXPECT_EQ(0u, elem_graph.num_parallel_edges());

    stk::io::fill_mesh("generated:1x1x4", bulkData);

    std::vector<size_t> counts;
    stk::mesh::count_entities(bulkData.mesh_meta_data().locally_owned_part(), bulkData, counts);
    int numLocallyOwnedElems = counts[stk::topology::ELEM_RANK];
    EXPECT_EQ(4, numLocallyOwnedElems);

    stk::mesh::Entity elem1 = bulkData.get_entity(stk::topology::ELEM_RANK, 1);
    stk::mesh::Entity elem2 = bulkData.get_entity(stk::topology::ELEM_RANK, 2);
    stk::mesh::Entity elem3 = bulkData.get_entity(stk::topology::ELEM_RANK, 3);
    stk::mesh::Entity elem4 = bulkData.get_entity(stk::topology::ELEM_RANK, 4);

    stk::mesh::EntityVector elements_to_add;
    elements_to_add.push_back(elem1);
    elements_to_add.push_back(elem2);
    elements_to_add.push_back(elem3);
    elements_to_add.push_back(elem4);

    for (unsigned i=0; i<elements_to_add.size(); i++)
    {
      EXPECT_TRUE(bulkData.is_valid(elements_to_add[i]));
      EXPECT_EQ(0, bulkData.parallel_owner_rank(elements_to_add[i]));
    }

    elem_graph.add_elements(elements_to_add);

    EXPECT_EQ(4u, elem_graph.size());
    EXPECT_EQ(6u, elem_graph.num_edges());
    EXPECT_EQ(0u, elem_graph.num_parallel_edges());

    EXPECT_EQ(5, elem_graph.check_local_connectivity(elem1, elem2));
    EXPECT_EQ(4, elem_graph.check_local_connectivity(elem2, elem1));
    EXPECT_EQ(5, elem_graph.check_local_connectivity(elem2, elem3));
    EXPECT_EQ(4, elem_graph.check_local_connectivity(elem3, elem2));
    EXPECT_EQ(5, elem_graph.check_local_connectivity(elem3, elem4));
    EXPECT_EQ(4, elem_graph.check_local_connectivity(elem4, elem3));
    EXPECT_EQ(-1, elem_graph.check_local_connectivity(elem1, elem3));
    EXPECT_EQ(-1, elem_graph.check_local_connectivity(elem2, elem4));

    stk::mesh::impl::DeletedElementInfoVector elems_to_delete;
    elems_to_delete.push_back({elem2, 2, false});
    elems_to_delete.push_back({elem3, 3, false});

    elem_graph.delete_elements(elems_to_delete);

    EXPECT_EQ(2u, elem_graph.size());
    EXPECT_EQ(0u, elem_graph.num_edges());
    EXPECT_EQ(0u, elem_graph.num_parallel_edges());

    EXPECT_EQ(-1, elem_graph.check_local_connectivity(elem1, elem2));
    EXPECT_EQ(-1, elem_graph.check_local_connectivity(elem2, elem1));
    EXPECT_EQ(-1, elem_graph.check_local_connectivity(elem2, elem3));
    EXPECT_EQ(-1, elem_graph.check_local_connectivity(elem3, elem2));
    EXPECT_EQ(-1, elem_graph.check_local_connectivity(elem3, elem4));
    EXPECT_EQ(-1, elem_graph.check_local_connectivity(elem4, elem3));
    EXPECT_EQ(-1, elem_graph.check_local_connectivity(elem1, elem3));
    EXPECT_EQ(-1, elem_graph.check_local_connectivity(elem2, elem4));

    elements_to_add.clear();
    elements_to_add.push_back(elem2);
    elements_to_add.push_back(elem3);

    elem_graph.add_elements(elements_to_add);

    EXPECT_EQ(4u, elem_graph.size());
    EXPECT_EQ(6u, elem_graph.num_edges());
    EXPECT_EQ(0u, elem_graph.num_parallel_edges());

    EXPECT_EQ(5, elem_graph.check_local_connectivity(elem1, elem2));
    EXPECT_EQ(4, elem_graph.check_local_connectivity(elem2, elem1));
    EXPECT_EQ(5, elem_graph.check_local_connectivity(elem2, elem3));
    EXPECT_EQ(4, elem_graph.check_local_connectivity(elem3, elem2));
    EXPECT_EQ(5, elem_graph.check_local_connectivity(elem3, elem4));
    EXPECT_EQ(4, elem_graph.check_local_connectivity(elem4, elem3));
    EXPECT_EQ(-1, elem_graph.check_local_connectivity(elem1, elem3));
    EXPECT_EQ(-1, elem_graph.check_local_connectivity(elem2, elem4));

    elems_to_delete.clear();
    elems_to_delete.push_back({elem4, 4, false});
    elems_to_delete.push_back({elem2, 2, false});
    elems_to_delete.push_back({elem1, 1, false});
    elems_to_delete.push_back({elem3, 3, false});

    elem_graph.delete_elements(elems_to_delete);

    EXPECT_EQ(0u, elem_graph.size());
    for(unsigned i=0; i<elem_graph.get_graph().get_num_elements_in_graph(); ++i)
    {
      EXPECT_EQ(0u, elem_graph.get_graph().get_num_edges_for_element(i));
    }
  }
}

TEST(ElementGraph, HexAddShellSerial)
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
  //       |   |          |   |
  //       |  2.0---------|--6.0
  //       |  /           |  /
  //       | /            | /
  //       |/             |/
  //      1.0------------5.0
  //                      ^
  //                      |
  //                       ---- Added single shell element

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  int p_size = stk::parallel_machine_size(pm);

  if(p_size != 1) {
    return;
  }

  const unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, pm);
  stk::mesh::BulkData& mesh = *bulkPtr;
  stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();

  stk::mesh::Part * hexPart   = &meta.declare_part_with_topology("hex_part",   stk::topology::HEX_8);
  stk::mesh::Part * shellPart = &meta.declare_part_with_topology("shell_part", stk::topology::SHELL_QUAD_4);
  meta.commit();

  std::vector<stk::mesh::EntityIdVector> hexNodeIDs {
    { 1, 2, 3, 4, 5, 6, 7, 8 }
  };
  stk::mesh::EntityId hexElemIDs[] = { 1 };

  std::vector<stk::mesh::EntityIdVector> shellNodeIDs {
    { 5, 6, 7, 8 }
  };
  stk::mesh::EntityId shellElemIDs[] = { 2 };

  // Build the base hex mesh
  mesh.modification_begin();
  for (size_t i = 0; i < hexNodeIDs.size(); ++i) {
    stk::mesh::declare_element(mesh, *hexPart, hexElemIDs[i], hexNodeIDs[i]);
  }
  mesh.modification_end();

  // Initialize the graph based on the existing hex mesh
  ElemElemGraphTester elem_graph(mesh);
  EXPECT_EQ(1u, elem_graph.size());

  // Add a shell
  mesh.modification_begin();
  stk::mesh::EntityVector added_shells;
  for (size_t i = 0; i < shellNodeIDs.size(); ++i) {
    added_shells.push_back( stk::mesh::declare_element(mesh, *shellPart, shellElemIDs[i], shellNodeIDs[i]) );
  }
  mesh.modification_end();

  elem_graph.add_elements(added_shells);

  EXPECT_EQ(2u, elem_graph.size());
  EXPECT_EQ(0u, elem_graph.num_parallel_edges());

  const Entity hex1   = mesh.get_entity(stk::topology::ELEM_RANK, 1);
  const Entity shell2 = mesh.get_entity(stk::topology::ELEM_RANK, 2);

  {
    // Connectivity for Hex Element 1
    EXPECT_EQ(1u,     elem_graph.get_num_connected_elems(hex1));
    impl::ElementViaSidePair elem_via_side = elem_graph.get_connected_element_and_via_side(hex1, 0);
    EXPECT_EQ(5,      elem_via_side.side);
    EXPECT_EQ(shell2, elem_via_side.element);
    EXPECT_TRUE(elem_graph.is_connected_elem_locally_owned(hex1, 0));
  }
  {
    // Connectivity for Shell Element 2
    EXPECT_EQ(1u,             elem_graph.get_num_connected_elems(shell2));
    impl::ElementViaSidePair elem_via_side = elem_graph.get_connected_element_and_via_side(shell2, 0);
    EXPECT_EQ(1,              elem_via_side.side);
    EXPECT_EQ(hex1, elem_via_side.element);
    EXPECT_TRUE(elem_graph.is_connected_elem_locally_owned(shell2, 0));
  }
  EXPECT_EQ(0u, elem_graph.num_parallel_edges());
}

TEST( ElementGraph, HexDelShellSerial )
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
  //       |   |          |   |
  //       |  2.0---------|--6.0
  //       |  /           |  /
  //       | /            | /
  //       |/             |/
  //      1.0------------5.0
  //                      ^
  //                      |
  //                       ---- Deleting shell element

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  int p_size = stk::parallel_machine_size(pm);

  if(p_size > 1)
  {
    return;
  }

  const unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, pm);
  stk::mesh::BulkData& mesh = *bulkPtr;
  stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();

  stk::mesh::Part * hexPart   = &meta.declare_part_with_topology("hex_part",   stk::topology::HEX_8);
  stk::mesh::Part * shellPart = &meta.declare_part_with_topology("shell_part", stk::topology::SHELL_QUAD_4);
  meta.commit();

  std::vector<stk::mesh::EntityIdVector> hexNodeIDs {
    { 1, 2, 3, 4, 5,  6,  7,  8 }
  };
  stk::mesh::EntityId hexElemIDs[] = { 1 };

  std::vector<stk::mesh::EntityIdVector> shellNodeIDs {
    { 5, 6, 7, 8 }
  };
  stk::mesh::EntityId shellElemIDs[] = { 2 };

  mesh.modification_begin();
  for (size_t i = 0; i < hexNodeIDs.size(); ++i) {
    stk::mesh::declare_element(mesh, *hexPart, hexElemIDs[i], hexNodeIDs[i]);
  }
  for (size_t i = 0; i < shellNodeIDs.size(); ++i) {
    stk::mesh::declare_element(mesh, *shellPart, shellElemIDs[i], shellNodeIDs[i]);
  }
  mesh.modification_end();

  ElemElemGraphTester elemElemGraph(mesh);

  const Entity hex1   = mesh.get_entity(stk::topology::ELEM_RANK, 1);
  const Entity shell2 = mesh.get_entity(stk::topology::ELEM_RANK, 2);

  stk::mesh::impl::DeletedElementInfoVector elements_to_delete;
  elements_to_delete.push_back({shell2, 2, mesh.bucket(shell2).topology().is_shell()});

  elemElemGraph.delete_elements(elements_to_delete);

  mesh.modification_begin();
  mesh.destroy_entity(shell2);
  mesh.modification_end();

  // Connectivity for Hex Element 1
  EXPECT_EQ(0u, elemElemGraph.get_num_connected_elems(hex1));
  EXPECT_EQ(0u, elemElemGraph.num_parallel_edges());

}

TEST( ElementGraph, HexAddShellHexSerial )
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
  //       |   |          |   |          |   |
  //       |  2.0---------|--6.0---------|-10.0
  //       |  /           |  /           |  /
  //       | /            | /            | /
  //       |/             |/             |/
  //      1.0------------5.0------------9.0
  //                      ^
  //                      |
  //                       ---- Added single shell element

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  int p_size = stk::parallel_machine_size(pm);

  if(p_size > 1)
  {
    return;
  }

  const unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, pm);
  stk::mesh::BulkData& mesh = *bulkPtr;
  stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();

  stk::mesh::Part * hexPart   = &meta.declare_part_with_topology("hex_part",   stk::topology::HEX_8);
  stk::mesh::Part * shellPart = &meta.declare_part_with_topology("shell_part", stk::topology::SHELL_QUAD_4);
  meta.commit();

  std::vector<stk::mesh::EntityIdVector> hexNodeIDs {
    { 1, 2, 3, 4, 5,  6,  7,  8 },
    { 5, 6, 7, 8, 9, 10, 11, 12 }
  };
  stk::mesh::EntityId hexElemIDs[] = { 1, 2 };

  std::vector<stk::mesh::EntityIdVector> shellNodeIDs {
    { 5, 6, 7, 8 }
  };
  stk::mesh::EntityId shellElemIDs[] = { 3 };

  mesh.modification_begin();
  for (size_t i = 0; i < hexNodeIDs.size(); ++i) {
    stk::mesh::declare_element(mesh, *hexPart, hexElemIDs[i], hexNodeIDs[i]);
  }
  mesh.modification_end();

  ElemElemGraphTester elemElemGraph(mesh);

  mesh.modification_begin();
  stk::mesh::EntityVector added_shells;
  for (size_t i = 0; i < shellNodeIDs.size(); ++i) {
    added_shells.push_back( stk::mesh::declare_element(mesh, *shellPart, shellElemIDs[i], shellNodeIDs[i]) );
  }
  mesh.modification_end();

  elemElemGraph.add_elements(added_shells);

  const Entity hex1   = mesh.get_entity(stk::topology::ELEM_RANK, 1);
  const Entity hex2   = mesh.get_entity(stk::topology::ELEM_RANK, 2);
  const Entity shell3 = mesh.get_entity(stk::topology::ELEM_RANK, 3);

  {
    // Connectivity for Hex Element 1
    EXPECT_EQ(1u,     elemElemGraph.get_num_connected_elems(hex1));
    impl::ElementViaSidePair elem_via_side = elemElemGraph.get_connected_element_and_via_side(hex1, 0);
    EXPECT_EQ(5,      elem_via_side.side);
    EXPECT_EQ(shell3, elem_via_side.element);
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));
  }
  {
    // Connectivity for Hex Element 2
    EXPECT_EQ(1u,     elemElemGraph.get_num_connected_elems(hex2));
    impl::ElementViaSidePair elem_via_side = elemElemGraph.get_connected_element_and_via_side(hex2, 0);
    EXPECT_EQ(4,      elem_via_side.side);
    EXPECT_EQ(shell3, elem_via_side.element);
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex2, 0));
  }
  // Connectivity for Shell Element 3
  EXPECT_EQ(2u,   elemElemGraph.get_num_connected_elems(shell3));
  impl::ElementViaSidePair elem0_via_side = elemElemGraph.get_connected_element_and_via_side(shell3, 0);
  impl::ElementViaSidePair elem1_via_side = elemElemGraph.get_connected_element_and_via_side(shell3, 1);
  EXPECT_EQ(1,    elem0_via_side.side);
  EXPECT_EQ(0,    elem1_via_side.side);
  EXPECT_EQ(hex1, elem0_via_side.element);
  EXPECT_EQ(hex2, elem1_via_side.element);
  EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell3, 0));
  EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell3, 1));

  EXPECT_EQ(0u, elemElemGraph.num_parallel_edges());
}

stk::mesh::Entity get_element_side(stk::mesh::BulkData& bulkData, stk::mesh::Entity element, stk::mesh::ConnectivityOrdinal side_ordinal)
{
  stk::mesh::Entity side = stk::mesh::Entity();
  stk::mesh::EntityRank side_rank = bulkData.mesh_meta_data().side_rank();

  unsigned elem_num_faces = bulkData.num_connectivity(element, side_rank);
  const stk::mesh::Entity * elem_sides = bulkData.begin(element, side_rank);
  const stk::mesh::ConnectivityOrdinal * elem_ord_it = bulkData.begin_ordinals(element, side_rank);
  for (unsigned i=0 ; i<elem_num_faces ; ++i)
  {
    if (elem_ord_it[i] == static_cast<unsigned>(side_ordinal))
    {
      side = elem_sides[i];
      break;
    }
  }

  return side;
}

void test_parallel_graph_info(const stk::mesh::Graph& graph, const ParallelInfoForGraphEdges& parallel_graph,
                              LocalId this_element, LocalId other_element, int other_proc, int permutation)
{
  bool did_find = false;

  for(size_t i=0;i<graph.get_num_elements_in_graph();++i)
  {
    size_t numConnected = graph.get_num_edges_for_element(i);
    for(size_t j=0;j<numConnected;++j)
    {
      const stk::mesh::GraphEdge& graphEdge = graph.get_edge_for_element(i, j);
      if(graphEdge.elem2()==-1*other_element && static_cast<LocalId>(i) == this_element)
      {
        const stk::mesh::impl::ParallelInfo& parallelInfo= parallel_graph.get_parallel_info_for_graph_edge(graphEdge);
        did_find = true;
        EXPECT_EQ(other_proc, parallelInfo.get_proc_rank_of_neighbor());
        EXPECT_EQ(permutation, parallelInfo.m_permutation);
      }
    }
  }
  ASSERT_TRUE(did_find);
}

template<class T>
void print_graph(const std::string &title, int proc_id, T& elem_graph)
{
  std::ostringstream os;

  os << title << " for processor " << proc_id << std::endl;
  for (size_t i=0;i<elem_graph.size();++i)
  {
    os << "Element " << i << "::\t";
    for (size_t j=0;j<elem_graph[i].size();++j)
    {
      os << elem_graph[i][j] << "\t";
    }
    os << std::endl;
  }
  std::cerr << os.str();
}

//BeginDocExample1

stk::mesh::EntityVector get_killed_elements(stk::mesh::BulkData& bulkData, const int killValue, const stk::mesh::Part& active)
{
  stk::mesh::EntityVector killedElements;
  const stk::mesh::BucketVector& buckets = bulkData.buckets(stk::topology::ELEMENT_RANK);
  for(size_t b = 0; b < buckets.size(); ++b)
  {
    const stk::mesh::Bucket &bucket = *buckets[b];
    if(bucket.owned() && bucket.member(active))
    {
      for(size_t e = 0; e < bucket.size(); ++e)
      {
        stk::mesh::Entity entity = bucket[e];
        bool should_element_be_killed = bulkData.identifier(entity) < static_cast<stk::mesh::EntityId>(killValue);
        if(bulkData.bucket(entity).member(active) && should_element_be_killed == true)
        {
          killedElements.push_back(bucket[e]);
        }
      }
    }
  }
  return killedElements;
}


TEST(ElementGraph, check_graph_connectivity)
{
  // element0 --> element1 --> element2
  ElementGraph elem_graph = {
    {1},
    {0,2},
    {1}
  };

  SidesForElementGraph via_side = {
    {4},
    {1,5},
    {3}
  };

  EXPECT_EQ(4, check_connectivity(elem_graph, via_side, 0, 1));
  EXPECT_EQ(1, check_connectivity(elem_graph, via_side, 1, 0));
  EXPECT_EQ(5, check_connectivity(elem_graph, via_side, 1, 2));
  EXPECT_EQ(3, check_connectivity(elem_graph, via_side, 2, 1));

  EXPECT_EQ(-1, check_connectivity(elem_graph, via_side, 0, 2));
  EXPECT_EQ(-1, check_connectivity(elem_graph, via_side, 2, 0));
  EXPECT_EQ(-1, check_connectivity(elem_graph, via_side, 3, 0));
  EXPECT_EQ(-1, check_connectivity(elem_graph, via_side, 0, 3));
  EXPECT_EQ(-1, check_connectivity(elem_graph, via_side, 0, 0));
}

TEST(ElementGraph, create_element_graph_serial)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  std::vector<double> wall_times;
  wall_times.reserve(10);
  std::vector<std::string> msgs;
  msgs.reserve(10);

  std::vector<size_t> mem_usage;

  wall_times.push_back(stk::wall_time());
  msgs.push_back("program-start");
  mem_usage.push_back(stk::get_memory_usage_now());

  if(stk::parallel_machine_size(comm) == 1)
  {
    std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(comm);
    stk::mesh::BulkData& bulkData = *bulkPtr;
    std::ostringstream os;

    stk::io::fill_mesh("generated:1x1x4", bulkData);

    wall_times.push_back(stk::wall_time());
    msgs.push_back("after mesh-read");
    mem_usage.push_back(stk::get_memory_usage_now());

    std::vector<size_t> counts;
    stk::mesh::count_entities(bulkData.mesh_meta_data().locally_owned_part(), bulkData, counts);
    int numElems = counts[stk::topology::ELEM_RANK];
    EXPECT_EQ(4, numElems);

    ElemElemGraphTester elemElemGraph(bulkData);

    size_t expectedNumElems = counts[stk::topology::ELEM_RANK];
    ASSERT_EQ(expectedNumElems, elemElemGraph.get_graph_size());

    wall_times.push_back(stk::wall_time());
    msgs.push_back("after ElemElemGraphTester constructor");
    mem_usage.push_back(stk::get_memory_usage_now());

    int left_side_id = 4;
    int right_side_id = 5;

    for(size_t i=0; i<elemElemGraph.get_graph().get_num_elements_in_graph(); ++i)
    {
      size_t numConnected = elemElemGraph.get_graph().get_num_edges_for_element(i);
      if (i == 0)
      {
        ASSERT_EQ(1u, numConnected);
        const stk::mesh::GraphEdge & graphEdge = elemElemGraph.get_graph().get_edge_for_element(i, 0);
        EXPECT_EQ(1, graphEdge.elem2());
        EXPECT_EQ(right_side_id, graphEdge.side1());
      }
      else if (i == elemElemGraph.get_graph().get_num_elements_in_graph() - 1)
      {
        LocalId second_to_last_element_index = elemElemGraph.get_graph().get_num_elements_in_graph() - 2;
        ASSERT_EQ(1u, numConnected);
        const stk::mesh::GraphEdge & graphEdge = elemElemGraph.get_graph().get_edge_for_element(i, 0);
        EXPECT_EQ(second_to_last_element_index, graphEdge.elem2());
        EXPECT_EQ(left_side_id, graphEdge.side1());
      }
      else
      {
        ASSERT_EQ(2u, numConnected);
        LocalId element_to_the_left = i-1;
        LocalId element_to_the_right = i+1;
        const stk::mesh::GraphEdge & graphEdge0 = elemElemGraph.get_graph().get_edge_for_element(i, 0);
        const stk::mesh::GraphEdge & graphEdge1 = elemElemGraph.get_graph().get_edge_for_element(i, 1);
        EXPECT_EQ(element_to_the_left, graphEdge0.elem2());
        EXPECT_EQ(element_to_the_right, graphEdge1.elem2());
        EXPECT_EQ(left_side_id, graphEdge0.side1());
        EXPECT_EQ(right_side_id, graphEdge1.side1());
      }
    }

    EXPECT_EQ(0u, elemElemGraph.num_parallel_edges());

    if (stk::parallel_machine_rank(comm) == 0)
    {
      for(size_t i=0;i<wall_times.size();++i)
      {
        std::cerr << "Wall time " << msgs[i] << ":\t" << wall_times[i] - wall_times[0] << std::endl;
      }

      for(size_t i=0;i<mem_usage.size();++i)
      {
        std::cerr << "Memory usage " << msgs[i] << ":\t" << mem_usage[i] - mem_usage[0] << std::endl;
      }
    }
  }
}

TEST(ElementGraph, create_element_graph_parallel)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  std::vector<double> wall_times;
  wall_times.reserve(10);
  std::vector<std::string> msgs;
  msgs.reserve(10);

  std::vector<size_t> mem_usage;

  wall_times.push_back(stk::wall_time());
  msgs.push_back("program-start");
  mem_usage.push_back(stk::get_memory_usage_now());

  if(stk::parallel_machine_size(comm) == 2)
  {
    std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(comm);
    stk::mesh::BulkData& bulkData = *bulkPtr;

    stk::io::fill_mesh("generated:1x1x4", bulkData);

    std::vector<size_t> counts;
    stk::mesh::count_entities(bulkData.mesh_meta_data().locally_owned_part(), bulkData, counts);
    int numLocallyOwnedElems = counts[stk::topology::ELEM_RANK];
    EXPECT_EQ(2, numLocallyOwnedElems);

    ElemElemGraphTester elemElemGraph(bulkData);

    size_t expectedNumElems = counts[stk::topology::ELEM_RANK];
    ASSERT_EQ(expectedNumElems, elemElemGraph.get_graph_size());

    wall_times.push_back(stk::wall_time());
    msgs.push_back("after ElemElemGraphTester constructor");
    mem_usage.push_back(stk::get_memory_usage_now());

    LocalId element_to_test_local_id = std::numeric_limits<LocalId>::max();
    int side_id = -1;
    int left_side_id = 4;
    int right_side_id = 5;

    EXPECT_EQ(2u, elemElemGraph.get_graph().get_num_elements_in_graph());
    EXPECT_EQ(1u, elemElemGraph.num_parallel_edges());

    if (stk::parallel_machine_rank(comm) == 0)
    {
      element_to_test_local_id = 1;
      side_id = right_side_id;
    }
    else
    {
      element_to_test_local_id = 0;
      side_id = left_side_id;
    }

    for(size_t i=0;i<elemElemGraph.get_graph().get_num_elements_in_graph();++i)
    {
      size_t numConnectedElements = elemElemGraph.get_graph().get_num_edges_for_element(i);
      if (static_cast<LocalId>(i) == element_to_test_local_id)
      {
        // Element on parallel boundary
        ASSERT_EQ(2u, numConnectedElements);
        const stk::mesh::GraphEdge & graphEdge = elemElemGraph.get_graph().get_edge_for_element(i, 1);
        ASSERT_GE(-1, graphEdge.elem2());
        ASSERT_EQ(side_id, graphEdge.side1());
      }
      else
      {
        EXPECT_EQ(1u, numConnectedElements);
      }
    }

    if (stk::parallel_machine_rank(comm) == 0)
    {
      for(size_t i=0;i<wall_times.size();++i)
      {
        std::cerr << "Wall time " << msgs[i] << ":\t" << wall_times[i] - wall_times[0] << std::endl;
      }

      for(size_t i=0;i<mem_usage.size();++i)
      {
        std::cerr << "Memory usage " << msgs[i] << ":\t" << mem_usage[i] - mem_usage[0] << std::endl;
      }
    }
  }
}

TEST(ElementSide, get_or_create_element_side_with_permutation)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  if(stk::parallel_machine_size(comm) == 1)
  {
    unsigned spatialDim = 3;
    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::Part& new_faces_part = meta.declare_part_with_topology("surface_5", stk::topology::QUAD_4);
    stk::mesh::PartVector face_parts {&new_faces_part};
    stk::io::put_io_part_attribute(new_faces_part);
    BulkDataElementGraphTester bulkData(meta, comm);

    stk::io::fill_mesh("generated:1x1x3", bulkData);

    //////////////// Make first side

    bulkData.modification_begin();

    //get_or_create_element_side_with_permutation(bulkData, element);
    stk::mesh::Entity element1 = bulkData.get_entity(stk::topology::ELEM_RANK, 1);
    stk::mesh::EntityId side_global_id = 11;
    stk::mesh::EntityRank side_rank = bulkData.mesh_meta_data().side_rank();
    stk::mesh::Permutation side_permutation = static_cast<stk::mesh::Permutation>(4);
    stk::mesh::ConnectivityOrdinal side_ordinal = static_cast<stk::mesh::ConnectivityOrdinal>(1);

    EXPECT_FALSE(impl::is_id_already_in_use_locally(bulkData, side_rank, side_global_id));
    EXPECT_FALSE(impl::does_side_exist_with_different_permutation(bulkData, element1, side_ordinal, side_permutation));
    EXPECT_FALSE(impl::does_element_side_exist(bulkData, element1, side_ordinal));

    impl::connect_side_to_element(bulkData, element1, side_global_id, side_ordinal, side_permutation, face_parts);

    bulkData.modification_end();

    stk::mesh::Entity side1 = get_element_side(bulkData, element1, side_ordinal);
    EXPECT_TRUE(bulkData.is_valid(side1));

    EXPECT_TRUE(impl::is_id_already_in_use_locally(bulkData, side_rank, side_global_id));

    stk::mesh::Permutation side_permutation1 = static_cast<stk::mesh::Permutation>(0);
    EXPECT_TRUE(impl::does_side_exist_with_different_permutation(bulkData, element1, side_ordinal, side_permutation1));

    size_t num_sides = stk::mesh::count_selected_entities(new_faces_part, bulkData.buckets(side_rank));
    EXPECT_EQ(1u, num_sides);
  }
}

TEST(ElementGraph, test_parallel_graph_info_data_structure)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1)
  {
    const int other_side_ord = 2;

    stk::mesh::Graph graph;
    graph.set_num_local_elements(2);
    stk::mesh::GraphEdge graphEdge(1, 5, -3, other_side_ord);
    std::vector<stk::mesh::GraphEdge> edges =
    { stk::mesh::GraphEdge(0, 4, 1, 1),
      stk::mesh::GraphEdge(1, 1, 0, 4),
      graphEdge };
    graph.add_sorted_edges(edges);

    stk::mesh::ParallelInfoForGraphEdges parallel_graph(stk::parallel_machine_rank(MPI_COMM_WORLD));
    int other_proc = 1;
    LocalId local_element = 1;
    LocalId other_element = 3;
    int permutation = 0;

    parallel_graph.insert_parallel_info_for_graph_edge(graphEdge, stk::mesh::impl::ParallelInfo(other_proc, permutation,
                                                                                                stk::topology::INVALID_TOPOLOGY));

    size_t num_elems_this_proc = graph.get_num_elements_in_graph();
    EXPECT_EQ(2u, num_elems_this_proc);

    test_parallel_graph_info(graph, parallel_graph, local_element, other_element, other_proc, permutation);
  }
}

TEST(ElementGraph, test_parallel_graph_info_with_parallel_element_graph)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;

  if(stk::parallel_machine_size(comm) == 2)
  {
    std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(comm);
    stk::mesh::BulkData& bulkData = *bulkPtr;

    stk::io::fill_mesh("generated:1x1x4", bulkData);

    std::vector<size_t> counts;
    stk::mesh::count_entities(bulkData.mesh_meta_data().locally_owned_part(), bulkData, counts);
    size_t numLocallyOwnedElems = counts[stk::topology::ELEM_RANK];

    ElemElemGraphTester elemElemGraph(bulkData);

    ASSERT_EQ(numLocallyOwnedElems, elemElemGraph.get_graph_size());

    ParallelInfoForGraphEdges & parallel_graph = elemElemGraph.get_parallel_graph();

    if(stk::parallel_machine_rank(comm)==0)
    {
      LocalId local_element = 1;
      LocalId other_element = 3;
      int other_proc = 1;
      int permutation = 4;

      test_parallel_graph_info(elemElemGraph.get_graph(), parallel_graph, local_element, other_element, other_proc, permutation);
    }
    else
    {
      LocalId local_element = 0;
      LocalId other_element = 2;
      int other_proc = 0;
      int permutation = 4;

      test_parallel_graph_info(elemElemGraph.get_graph(), parallel_graph, local_element, other_element, other_proc, permutation);
    }

    EXPECT_EQ(3u, elemElemGraph.num_edges());
    EXPECT_EQ(1u, elemElemGraph.num_parallel_edges());
  }
}

TEST(ElementGraph, create_faces_using_element_graph_parallel)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  std::vector<double> wall_times;
  wall_times.reserve(10);
  std::vector<std::string> msgs;
  msgs.reserve(10);

  std::vector<size_t> mem_usage;

  wall_times.push_back(stk::wall_time());
  msgs.push_back("program-start");
  mem_usage.push_back(stk::get_memory_usage_now());

  if(stk::parallel_machine_size(comm) <= 2)
  {
    unsigned spatialDim = 3;
    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::Part& new_faces_part = meta.declare_part_with_topology("surface_5", stk::topology::QUAD_4);
    stk::io::put_io_part_attribute(new_faces_part);
    BulkDataElementGraphTester bulkData(meta, comm);

    stk::io::fill_mesh("generated:1x1x4", bulkData);

    wall_times.push_back(stk::wall_time());
    msgs.push_back("after mesh-read");
    mem_usage.push_back(stk::get_memory_usage_now());

    stk::mesh::experimental::create_faces(bulkData);
    const bool connectFacesToEdges = false;
    stk::mesh::create_all_sides(bulkData, meta.universal_part(), {&new_faces_part}, connectFacesToEdges);

    const stk::mesh::BucketVector& sharedNodeBuckets = bulkData.get_buckets(stk::topology::NODE_RANK, meta.globally_shared_part());
    for(size_t bucket_index=0; bucket_index<sharedNodeBuckets.size(); ++bucket_index)
    {
      const stk::mesh::Bucket& bucket = *sharedNodeBuckets[bucket_index];
      EXPECT_TRUE(bucket.member(new_faces_part));
    }

    wall_times.push_back(stk::wall_time());
    msgs.push_back("after create-faces");
    mem_usage.push_back(stk::get_memory_usage_now());

    std::vector<size_t> entity_counts;
    stk::mesh::comm_mesh_counts(bulkData, entity_counts);
    stk::mesh::EntityRank side_rank = meta.side_rank();
    unsigned num_faces = entity_counts[side_rank];
    EXPECT_EQ(21u, num_faces);

    if (stk::parallel_machine_rank(comm) == 0)
    {
      for(size_t i=0;i<wall_times.size();++i)
      {
        std::cerr << "Wall time " << msgs[i] << ":\t" << wall_times[i] - wall_times[0] << std::endl;
      }

      for(size_t i=0;i<mem_usage.size();++i)
      {
        std::cerr << "Memory usage " << msgs[i] << ":\t" << mem_usage[i] - mem_usage[0] << std::endl;
      }
    }
  }
}

TEST(ElementGraph, create_faces_using_element_graph_parallel_block_membership)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  std::vector<double> wall_times;
  wall_times.reserve(10);
  std::vector<std::string> msgs;
  msgs.reserve(10);

  std::vector<size_t> mem_usage;

  wall_times.push_back(stk::wall_time());
  msgs.push_back("program-start");
  mem_usage.push_back(stk::get_memory_usage_now());

  if(stk::parallel_machine_size(comm) <= 2)
  {
    unsigned spatialDim = 3;
    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::Part& new_faces_part = meta.declare_part_with_topology("surface_5", stk::topology::QUAD_4);

    stk::mesh::Part& block_2 = meta.declare_part_with_topology("block_2", stk::topology::HEX_8);
    EXPECT_EQ(stk::topology::ELEM_RANK, block_2.primary_entity_rank());

    stk::io::put_io_part_attribute(new_faces_part);
    BulkDataElementGraphTester bulkData(meta, comm);

    stk::io::fill_mesh("generated:1x1x4", bulkData);

    wall_times.push_back(stk::wall_time());
    msgs.push_back("after mesh-read");
    mem_usage.push_back(stk::get_memory_usage_now());

    stk::mesh::Part& block_1 = *meta.get_part("block_1");

    bulkData.modification_begin();

    if (bulkData.parallel_rank() == 1)
    {
      //Move elem 3 from block_1 to block_2 so that the boundary between block_1 and block_2
      //will coincide with the proc boundary, and the shared face between elems 2 & 3
      //should have both block_1 and block_2.
      stk::mesh::Entity elem3 = bulkData.get_entity(stk::topology::ELEM_RANK, 3);
      ASSERT_TRUE(bulkData.is_valid(elem3));
      stk::mesh::PartVector add_parts(1, &block_2);
      stk::mesh::PartVector rem_parts(1, &block_1);
      bulkData.change_entity_parts(elem3, add_parts, rem_parts);
    }

    bulkData.modification_end();

    const bool connectFacesToEdges = false;
    stk::mesh::create_all_sides(bulkData, meta.universal_part(), {&new_faces_part}, connectFacesToEdges);

    wall_times.push_back(stk::wall_time());
    msgs.push_back("after create-faces");
    mem_usage.push_back(stk::get_memory_usage_now());

    std::vector<size_t> entity_counts;
    stk::mesh::comm_mesh_counts(bulkData, entity_counts);
    stk::mesh::EntityRank side_rank = bulkData.mesh_meta_data().side_rank();
    unsigned num_faces = entity_counts[side_rank];
    EXPECT_EQ(21u, num_faces);

    const stk::mesh::BucketVector& sharedFaceBuckets = bulkData.get_buckets(stk::topology::FACE_RANK, meta.globally_shared_part());
    if (bulkData.parallel_size() == 2)
    {
      ASSERT_EQ(1u, sharedFaceBuckets.size());

      const stk::mesh::Bucket& bucket = *sharedFaceBuckets[0];
      ASSERT_EQ(1u, bucket.size());
      EXPECT_TRUE(bucket.member(block_2));
      EXPECT_TRUE(bucket.member(block_1));
    }

    if (stk::parallel_machine_rank(comm) == 0)
    {
      for(size_t i=0;i<wall_times.size();++i)
      {
        std::cerr << "Wall time " << msgs[i] << ":\t" << wall_times[i] - wall_times[0] << std::endl;
      }

      for(size_t i=0;i<mem_usage.size();++i)
      {
        std::cerr << "Memory usage " << msgs[i] << ":\t" << mem_usage[i] - mem_usage[0] << std::endl;
      }
    }
  }
}

TEST(ElementGraph, compare_performance_create_faces)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;

  int xdim = stk::unit_test_util::get_command_line_option("--xdim", 3);
  int ydim = xdim;
  int zdim = xdim * stk::parallel_machine_size(comm);

  std::string filename = stk::unit_test_util::get_name_of_generated_mesh(xdim, ydim, zdim, "|nodeset:zZ");

  {
    unsigned spatialDim = 3;
    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::Part& faces_part = meta.declare_part_with_topology("surface_5", stk::topology::QUAD_4);
    stk::io::put_io_part_attribute(faces_part);
    BulkDataElementGraphTester bulkData(meta, comm);

    stk::io::fill_mesh(filename, bulkData);

    {
      double wall_time_start = stk::wall_time();

      stk::mesh::PartVector parts(1, &faces_part);
      stk::mesh::create_faces(bulkData);

      double elapsed_time = stk::wall_time() - wall_time_start;

      if (stk::parallel_machine_rank(comm) == 0)
      {
        std::cerr << "Create faces without graph time: " << elapsed_time << std::endl;
      }
    }
  }

  {
    unsigned spatialDim = 3;
    stk::mesh::MetaData meta(spatialDim);
    bool force_no_induce = true;
    stk::mesh::Part& faces_part = meta.declare_part_with_topology("surface_5", stk::topology::QUAD_4, force_no_induce);
    stk::io::put_io_part_attribute(faces_part);
    const bool connectFacesToEdges = false;
    BulkDataElementGraphTester bulkData(meta, comm);

    stk::io::fill_mesh(filename, bulkData);

    {
      double wall_time_start = stk::wall_time();

      stk::mesh::create_all_sides(bulkData, meta.universal_part(), {&faces_part}, connectFacesToEdges);

      double elapsed_time = stk::wall_time() - wall_time_start;

      std::vector<size_t> counts;
      stk::mesh::comm_mesh_counts(bulkData, counts);

      if (stk::parallel_machine_rank(comm) == 0)
      {
        std::cerr << "Total time for creating graph and faces: " << elapsed_time << std::endl;
        std::cerr << "Total # of elements: " << counts[stk::topology::ELEM_RANK] << std::endl;
      }
    }
  }
}

TEST(ElementGraph, make_items_inactive)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;

  unsigned nProc = stk::parallel_machine_size(comm);

  if(nProc <= 2)
  {
    unsigned spatialDim = 3;

    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::Part& faces_part = meta.declare_part_with_topology("surface_5", stk::topology::QUAD_4);
    stk::mesh::PartVector boundary_mesh_parts { &faces_part };
    stk::io::put_io_part_attribute(faces_part);
    ElemGraphTestUtils::ElementDeathBulkDataTester bulkData(meta, comm, stk::mesh::BulkData::AUTO_AURA);

    stk::mesh::Part& active = meta.declare_part("active"); // can't specify rank, because it gets checked against size of rank_names

    ASSERT_TRUE(active.primary_entity_rank() == stk::topology::INVALID_RANK);

    stk::io::fill_mesh("generated:1x1x4", bulkData);

    stk::unit_test_util::put_mesh_into_part(bulkData, active);

    ElemElemGraphTester graph(bulkData);

    size_t num_gold_edges =  6/bulkData.parallel_size();
    ASSERT_EQ(num_gold_edges, graph.num_edges());
    ASSERT_EQ((nProc-1), graph.num_parallel_edges());

    stk::mesh::EntityVector deactivated_elems;

    stk::mesh::Entity elem1 = bulkData.get_entity(stk::topology::ELEM_RANK, 1);
    if (bulkData.is_valid(elem1) && bulkData.bucket(elem1).owned() )
    {
      deactivated_elems.push_back(elem1);
    }
    stk::mesh::Entity elem3 = bulkData.get_entity(stk::topology::ELEM_RANK, 3);
    if (bulkData.is_valid(elem3) && bulkData.bucket(elem3).owned())
    {
      deactivated_elems.push_back(elem3);
    }

    ElemGraphTestUtils::deactivate_elements(deactivated_elems, bulkData, active);
    bulkData.modification_begin();
    bulkData.my_de_induce_unranked_part_from_nodes(deactivated_elems, active);

    for(size_t i=0; i<deactivated_elems.size(); ++i)
    {
      EXPECT_FALSE(bulkData.bucket(deactivated_elems[i]).member(active));
    }

    const stk::mesh::BucketVector& all_node_buckets = bulkData.buckets(stk::topology::NODE_RANK);
    for(size_t i=0; i<all_node_buckets.size(); ++i)
    {
      const stk::mesh::Bucket& bucket = *all_node_buckets[i];
      for(size_t node_index=0; node_index<bucket.size(); ++node_index)
      {
        stk::mesh::EntityId id = bulkData.identifier(bucket[node_index]);
        if (id >=1 && id <= 4)
        {
          EXPECT_FALSE(bucket.member(active));
        }
        else
        {
          EXPECT_TRUE(bucket.member(active)) << "for node id " << id << std::endl;
        }
      }
    }

    bulkData.modification_end();
    ASSERT_EQ((nProc-1), graph.num_parallel_edges());
  }
}

TEST(ElementGraph, test_element_death)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;

  if(stk::parallel_machine_size(comm) <= 2)
  {
    //IO error when this is <4.  Shared face being attached to the wrong element
    int xdim = stk::unit_test_util::get_command_line_option("--zdim", 4);
    int ydim = xdim;
    int zdim = xdim; //  * stk::parallel_machine_size(comm);

    std::string filename = stk::unit_test_util::get_name_of_generated_mesh(xdim, ydim, zdim, "|nodeset:zZ");

    {
      unsigned spatialDim = 3;
      stk::mesh::MetaData meta(spatialDim);
      stk::mesh::Part& faces_part = meta.declare_part_with_topology("surface_5", stk::topology::QUAD_4);
      stk::mesh::PartVector boundary_mesh_parts { &faces_part };
      stk::io::put_io_part_attribute(faces_part);
      BulkDataElementGraphTester bulkData(meta, comm);

      stk::mesh::Part& active = meta.declare_part("active", stk::topology::ELEMENT_RANK);
      stk::io::fill_mesh(filename, bulkData);

      double start_graph = stk::wall_time();

      ASSERT_TRUE(meta.get_part("block_1") != NULL);

      stk::mesh::Part& block_1 = *meta.get_part("block_1");

      stk::unit_test_util::put_mesh_into_part(bulkData, active);

      std::ostringstream os;
      os << "Proc id: " << bulkData.parallel_rank() << std::endl;

      stk::mesh::ElemElemGraph &elementGraph = bulkData.get_face_adjacent_element_graph();

      double elapsed_graph_time = stk::wall_time() - start_graph;
      os << "Time to create graph: " << elapsed_graph_time << std::endl;

      stk::mesh::EntityRank side_rank = meta.side_rank();

      int num_time_steps = xdim * ydim * zdim;
      double elapsed_death_time = 0;

      boundary_mesh_parts.push_back(&active);

      for(int i = 0; i < num_time_steps; ++i)
      {
        stk::mesh::EntityVector killedElements = get_killed_elements(bulkData, i, active);
        stk::unit_test_util::move_killed_elements_out_of_parts(bulkData, killedElements, {&block_1, &active});

        stk::mesh::impl::ParallelSelectedInfo remoteActiveSelector;
        stk::mesh::impl::populate_selected_value_for_remote_elements(bulkData, elementGraph, active, remoteActiveSelector);

        double start_time = stk::wall_time();
        process_killed_elements(bulkData, killedElements, active, remoteActiveSelector, boundary_mesh_parts, &boundary_mesh_parts);
        elapsed_death_time += (stk::wall_time() - start_time);
      }

      stk::mesh::Selector sel = block_1;
      std::vector<size_t> counts1;
      stk::mesh::comm_mesh_counts(bulkData, counts1, &sel);

      size_t num_active = counts1[stk::topology::ELEM_RANK];

      stk::mesh::Selector sel2 = faces_part;
      stk::mesh::comm_mesh_counts(bulkData, counts1, &sel2);

      size_t num_faces = counts1[side_rank];

      EXPECT_EQ(2u, num_active);
      EXPECT_EQ(5u, num_faces);

      if(stk::parallel_machine_rank(comm) == 0)
      {
        os << "Total time: " << elapsed_death_time << std::endl;
        os << "Total # of alive elements: " << num_active << std::endl;
        os << "Total # of faces: " << num_faces << std::endl;
      }

      std::cerr << os.str();
    }
  }
}


class ElementDeathRestartTest
{
public:
  typedef stk::mesh::Field<double> ActiveFieldType;

  ElementDeathRestartTest() {
  }

  virtual ~ElementDeathRestartTest() {
  }

  void load_without_restart()
  {
    initializeObjects();
    std::string filename = stk::unit_test_util::get_name_of_generated_mesh(1, 1, 2, "|sideset:zZ");
    stk::io::fill_mesh(filename, *bulkData);
    stk::unit_test_util::put_mesh_into_part(*bulkData, *activePart);
  }

  void read_restart_file(const std::string& filename)
  {
    stk::io::StkMeshIoBroker stkIo(bulkData->parallel());
    stkIo.set_bulk_data(*bulkData);
    stkIo.add_mesh_database(filename, stk::io::READ_RESTART);
    stkIo.create_input_mesh();
    stkIo.populate_bulk_data();
    stk::io::MeshField ioActiveField(*deathStatusField);
    stkIo.add_input_field(ioActiveField);
    stkIo.read_defined_input_fields(restartTime);
  }

  void load_with_restart(const std::string& filename)
  {
    initializeObjects();
    read_restart_file(filename);

    stk::mesh::EntityVector deactivatedElements = set_active_part_from_field();

    stk::mesh::impl::ParallelSelectedInfo remoteActiveSelector;
    stk::mesh::impl::populate_selected_value_for_remote_elements(*bulkData, bulkData->get_face_adjacent_element_graph(), *activePart, remoteActiveSelector);

    stk::mesh::process_killed_elements(*bulkData, deactivatedElements, *activePart, remoteActiveSelector, partsForCreatingSides, NULL);
  }

  void initializeObjects()
  {
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    stk::mesh::MeshBuilder builder(comm);
    builder.set_spatial_dimension(3);
    builder.set_aura_option(stk::mesh::BulkData::NO_AUTO_AURA);
    bulkData = builder.create();
    metaData = &(bulkData->mesh_meta_data());
    stk::io::put_io_part_attribute(metaData->universal_part());
    deathStatusField = &metaData->declare_field<double>(stk::topology::ELEM_RANK,deathStatusFieldName);
    stk::mesh::put_field_on_mesh(*deathStatusField,metaData->universal_part(), nullptr);

    activePart = &metaData->declare_part("active");
    deathPart = &metaData->declare_part("death_1", metaData->side_rank());
    const bool forceNoInduce = true;
    sidesCreatedDuringDeath = &metaData->declare_part("sides_created_during_death", metaData->side_rank(),forceNoInduce);

    partsForCreatingSides.push_back(activePart);
    partsForCreatingSides.push_back(sidesCreatedDuringDeath);
    partsForCreatingSides.push_back(deathPart);
    boundaryMeshParts.push_back(deathPart);
    activePartVector.push_back(activePart);

  }

  void kill_element(int globalId)
  {
    stk::mesh::Entity element = bulkData->get_entity(stk::topology::ELEM_RANK,globalId);
    stk::mesh::EntityVector elementsToKill;
    if (bulkData->is_valid(element))
    {
      elementsToKill.push_back(element);
    }
    stk::unit_test_util::move_killed_elements_out_of_parts(*bulkData, elementsToKill, activePartVector);

    stk::mesh::impl::ParallelSelectedInfo remoteActiveSelector;
    stk::mesh::impl::populate_selected_value_for_remote_elements(*bulkData, bulkData->get_face_adjacent_element_graph(), *activePart, remoteActiveSelector);

    process_killed_elements(*bulkData, elementsToKill, *activePart, remoteActiveSelector, partsForCreatingSides, &boundaryMeshParts);
  }

  void verify_mesh_before_death() const
  {
    stk::mesh::Entity element1 = bulkData->get_entity(stk::topology::ELEM_RANK,1);
    if (bulkData->is_valid(element1)) {
      EXPECT_EQ(bulkData->num_faces(element1), 1u);
      stk::mesh::Entity side4OfElement1 = bulkData->begin_faces(element1)[0];
      EXPECT_TRUE(bulkData->is_valid(side4OfElement1));
      EXPECT_EQ(bulkData->begin_face_ordinals(element1)[0], 4u);
      EXPECT_FALSE(bulkData->bucket(side4OfElement1).member(*sidesCreatedDuringDeath));
      EXPECT_TRUE(bulkData->bucket(side4OfElement1).member(*activePart));
    }
  }

  void verify_mesh_after_killing_element_1() const
  {
    stk::mesh::Entity element1 = bulkData->get_entity(stk::topology::ELEM_RANK,1);
    if (bulkData->is_valid(element1)) {
      EXPECT_FALSE( bulkData->bucket(element1).member(*activePart) );
      EXPECT_EQ(bulkData->num_faces(element1), 2u);
      const stk::mesh::Entity* sides = bulkData->begin_faces(element1);
      const stk::mesh::ConnectivityOrdinal* sideOrdinals = bulkData->begin_face_ordinals(element1);
      for (size_t sideI=0 ; sideI<2u ; ++sideI) {
        stk::mesh::Entity side = sides[sideI];
        int sideOrdinal = sideOrdinals[sideI];
        if (sideOrdinal == 5) {
          // Face between element1 and element2
          EXPECT_TRUE(bulkData->bucket(side).member(*activePart));
          EXPECT_TRUE(bulkData->bucket(side).member(*sidesCreatedDuringDeath));
        }
        else {
          // Side from generated mesh
          EXPECT_FALSE(bulkData->bucket(side).member(*activePart));
          EXPECT_FALSE(bulkData->bucket(side).member(*sidesCreatedDuringDeath));
        }
      }
    }
  }

  void verify_mesh_after_killing_element_2() const
  {
    stk::mesh::Entity element1 = bulkData->get_entity(stk::topology::ELEM_RANK,1);
    if (bulkData->is_valid(element1)) {
      EXPECT_EQ(bulkData->num_faces(element1), 1u);
      stk::mesh::Entity side = bulkData->begin_faces(element1)[0];
      stk::mesh::ConnectivityOrdinal sideOrdinalForSideFromGeneratedMesh = bulkData->begin_face_ordinals(element1)[0];
      EXPECT_EQ(sideOrdinalForSideFromGeneratedMesh, 4u);
      EXPECT_FALSE(bulkData->bucket(side).member(*activePart));
      EXPECT_FALSE(bulkData->bucket(side).member(*sidesCreatedDuringDeath));
    }

    stk::mesh::Entity element2 = bulkData->get_entity(stk::topology::ELEM_RANK,2);
    if (bulkData->is_valid(element2)) {
      EXPECT_FALSE( bulkData->bucket(element2).member(*activePart) );
      EXPECT_EQ(bulkData->num_faces(element2), 1u);
      stk::mesh::Entity side = bulkData->begin_faces(element2)[0];
      stk::mesh::ConnectivityOrdinal sideOrdinalForSideFromGeneratedMesh = bulkData->begin_face_ordinals(element2)[0];
      EXPECT_EQ(sideOrdinalForSideFromGeneratedMesh, 5u);
      EXPECT_FALSE(bulkData->bucket(side).member(*activePart));
      EXPECT_FALSE(bulkData->bucket(side).member(*sidesCreatedDuringDeath));
    }
  }

  stk::mesh::EntityVector set_active_part_from_field()
  {
    stk::mesh::EntityVector deactivatedElements;
    stk::mesh::EntityVector entities;
    std::vector<stk::mesh::PartVector> addParts;
    std::vector<stk::mesh::PartVector> removeParts;
    const stk::mesh::BucketVector& elements = bulkData->buckets(stk::topology::ELEM_RANK);
    const int deadElementStatus = 1;
    for (size_t bucketI=0 ; bucketI<elements.size() ; ++bucketI) {
      stk::mesh::Bucket& bucket = *(elements[bucketI]);
      for (size_t elementI=0 ; elementI<bucket.size() ; ++elementI) {
        stk::mesh::Entity element = bucket[elementI];
        double* deathStatus = stk::mesh::field_data(*deathStatusField,element);
        if (static_cast<int>(deathStatus[0]) == deadElementStatus) {
          entities.push_back(element);
          addParts.push_back({});
          removeParts.push_back({activePart});
          deactivatedElements.push_back(element);
        } else {
          entities.push_back(element);
          addParts.push_back({activePart});
          removeParts.push_back({});
        }
      }
    }
    bulkData->batch_change_entity_parts( entities, addParts, removeParts);
    return deactivatedElements;
  }

  void set_active_field_from_part()
  {
    const stk::mesh::BucketVector& elements = bulkData->buckets(stk::topology::ELEM_RANK);
    for (size_t bucketI=0 ; bucketI<elements.size() ; ++bucketI) {
      stk::mesh::Bucket& bucket = *(elements[bucketI]);
      for (size_t elementI=0 ; elementI<bucket.size() ; ++elementI) {
        stk::mesh::Entity element = bucket[elementI];
        double* deathStatus = stk::mesh::field_data(*deathStatusField,element);
        deathStatus[0] = !bulkData->bucket(element).member(*activePart);
      }
    }
  }

  void write_restart_file(std::string restartFileName)
  {
    stk::mesh::FieldBase* deathStatusFieldBase = deathStatusField;
    set_active_field_from_part();

    stk::io::StkMeshIoBroker stkIo(bulkData->parallel());
    stkIo.set_bulk_data(*bulkData);
    size_t fileHandle = stkIo.create_output_mesh(restartFileName, stk::io::WRITE_RESTART);
    stkIo.write_output_mesh(fileHandle);
    stkIo.add_field(fileHandle, *deathStatusFieldBase);
    stkIo.begin_output_step(fileHandle, restartTime);
    stkIo.write_defined_output_fields(fileHandle);
    stkIo.end_output_step(fileHandle);
  }
public:
  stk::mesh::MetaData* metaData = nullptr;
  std::shared_ptr<stk::mesh::BulkData> bulkData;
  stk::mesh::Part* activePart;
  stk::mesh::PartVector activePartVector;
  stk::mesh::PartVector partsForCreatingSides;
  stk::mesh::PartVector boundaryMeshParts;
  stk::mesh::Part* deathPart;
  stk::mesh::Part* sidesCreatedDuringDeath;
  ActiveFieldType* deathStatusField;
  std::string deathStatusFieldName = "death_status";
  double restartTime = 1.0;

};


TEST(ElementDeath, test_element_death_without_restart)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  if(stk::parallel_machine_size(comm) <= 2)
  {
    ElementDeathRestartTest elementDeathTest;
    elementDeathTest.load_without_restart();
    elementDeathTest.verify_mesh_before_death();
    elementDeathTest.kill_element(1);
    elementDeathTest.verify_mesh_after_killing_element_1();
    elementDeathTest.kill_element(2);
    elementDeathTest.verify_mesh_after_killing_element_2();
  }
}

TEST(ElementDeath, test_element_death_with_restart)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  if(stk::parallel_machine_size(comm) <= 2)
  {
    {
      ElementDeathRestartTest elementDeathTest;
      elementDeathTest.load_without_restart();
      elementDeathTest.verify_mesh_before_death();
      elementDeathTest.kill_element(1);
      elementDeathTest.verify_mesh_after_killing_element_1();
      elementDeathTest.write_restart_file("elemDeathRestartFile.exo");
    }
    {
      ElementDeathRestartTest elementDeathTest;
      elementDeathTest.load_with_restart("elemDeathRestartFile.exo");
      elementDeathTest.verify_mesh_after_killing_element_1();
      elementDeathTest.kill_element(2);
      elementDeathTest.verify_mesh_after_killing_element_2();
    }
    stk::unit_test_util::delete_mesh("elemDeathRestartFile.exo");
  }
}

//EndDocExample1

TEST( ElementGraph, HexHexHexSerial )
{
  //  ID.proc
  //
  //          3.0------------7.0-----------11.0-----------15.0
  //          /|             /|             /|             /|
  //         / |            / |            / |            / |
  //        /  |           /  |           /  |           /  |
  //      4.0------------8.0-----------12.0-----------16.0  |
  //       |   |          |   |          |   |          |   |
  //       |   |   1.0    |   |   2.0    |   |   3.0    |   |
  //       |   |          |   |          |   |          |   |
  //       |  2.0---------|--6.0---------|-10.0---------|-14.0
  //       |  /           |  /           |  /           |  /
  //       | /            | /            | /            | /
  //       |/             |/             |/             |/
  //      1.0------------5.0------------9.0-----------13.0
  //

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  unsigned p_size = stk::parallel_machine_size(pm);
  unsigned p_rank = stk::parallel_machine_rank(pm);

  if(p_size != 1u)
  {
    return;
  }

  const unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, pm);
  stk::mesh::BulkData& mesh = *bulkPtr;
  stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();

  stk::mesh::Part * hexPart = &meta.declare_part_with_topology("hex_part", stk::topology::HEX_8);
  meta.commit();

  std::vector<stk::mesh::EntityIdVector> hexNodeIDs {
    { 1,  2,  3,  4,  5,  6,  7,  8 },
    { 5,  6,  7,  8,  9, 10, 11, 12 },
    { 9, 10, 11, 12, 13, 14, 15, 16 }
  };
  stk::mesh::EntityId hexElemIDs[] = { 1, 2, 3 };
  stk::mesh::EntityId hexElemOwningProc[] = { 0, 0, 0 };

  mesh.modification_begin();
  for (size_t i = 0; i < hexNodeIDs.size(); ++i) {
    if (hexElemOwningProc[i] == p_rank) {
      stk::mesh::declare_element(mesh, *hexPart, hexElemIDs[i], hexNodeIDs[i]);
    }
  }
  mesh.modification_end();

  ElemElemGraphTester elemElemGraph(mesh);

  const Entity hex1 = mesh.get_entity(stk::topology::ELEM_RANK, 1);
  const Entity hex2 = mesh.get_entity(stk::topology::ELEM_RANK, 2);
  const Entity hex3 = mesh.get_entity(stk::topology::ELEM_RANK, 3);

  {
    // Connectivity for Hex Element 1
    EXPECT_EQ(1u,   elemElemGraph.get_num_connected_elems(hex1));
    impl::ElementViaSidePair elem_via_side = elemElemGraph.get_connected_element_and_via_side(hex1, 0);
    EXPECT_EQ(5,    elem_via_side.side);
    EXPECT_EQ(hex2, elem_via_side.element);
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));
  }

  // Connectivity for Hex Element 2
  EXPECT_EQ(2u,   elemElemGraph.get_num_connected_elems(hex2));
  impl::ElementViaSidePair elem0_via_side = elemElemGraph.get_connected_element_and_via_side(hex2, 0);
  impl::ElementViaSidePair elem1_via_side = elemElemGraph.get_connected_element_and_via_side(hex2, 1);
  EXPECT_EQ(4,    elem0_via_side.side);
  EXPECT_EQ(5,    elem1_via_side.side);
  EXPECT_EQ(hex1, elem0_via_side.element);
  EXPECT_EQ(hex3, elem1_via_side.element);
  EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex2, 0));
  EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex2, 1));

  {
    // Connectivity for Hex Element 3
    EXPECT_EQ(1u,   elemElemGraph.get_num_connected_elems(hex3));
    impl::ElementViaSidePair elem_via_side = elemElemGraph.get_connected_element_and_via_side(hex3, 0);
    EXPECT_EQ(4,    elem_via_side.side);
    EXPECT_EQ(hex2, elem_via_side.element);
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex3, 0));
  }
  EXPECT_EQ(4u, elemElemGraph.num_edges());
  EXPECT_EQ(0u, elemElemGraph.num_parallel_edges());
}

TEST( ElementGraph, HexShellSerial )
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
  //       |   |          |   |
  //       |  2.0---------|--6.0
  //       |  /           |  /
  //       | /            | /
  //       |/             |/
  //      1.0------------5.0
  //                      ^
  //                      |
  //                       ---- Single shell element

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  int p_size = stk::parallel_machine_size(pm);

  if(p_size > 1)
  {
    return;
  }

  const unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, pm);
  stk::mesh::BulkData& mesh = *bulkPtr;
  stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();

  stk::mesh::Part * hexPart   = &meta.declare_part_with_topology("hex_part",   stk::topology::HEX_8);
  stk::mesh::Part * shellPart = &meta.declare_part_with_topology("shell_part", stk::topology::SHELL_QUAD_4);
  meta.commit();

  std::vector<stk::mesh::EntityIdVector> hexNodeIDs {
    { 1, 2, 3, 4, 5,  6,  7,  8 }
  };
  stk::mesh::EntityId hexElemIDs[] = { 1 };

  std::vector<stk::mesh::EntityIdVector> shellNodeIDs {
    { 5, 6, 7, 8 }
  };
  stk::mesh::EntityId shellElemIDs[] = { 2 };

  mesh.modification_begin();
  for (size_t i = 0; i < hexNodeIDs.size(); ++i) {
    stk::mesh::declare_element(mesh, *hexPart, hexElemIDs[i], hexNodeIDs[i]);
  }
  for (size_t i = 0; i < shellNodeIDs.size(); ++i) {
    stk::mesh::declare_element(mesh, *shellPart, shellElemIDs[i], shellNodeIDs[i]);
  }
  mesh.modification_end();

  ElemElemGraphTester elemElemGraph(mesh);

  const Entity hex1   = mesh.get_entity(stk::topology::ELEM_RANK, 1);
  const Entity shell2 = mesh.get_entity(stk::topology::ELEM_RANK, 2);

  {
    // Connectivity for Hex Element 1
    EXPECT_EQ(1u,     elemElemGraph.get_num_connected_elems(hex1));
    impl::ElementViaSidePair elem_via_side = elemElemGraph.get_connected_element_and_via_side(hex1, 0);
    EXPECT_EQ(5,      elem_via_side.side);
    EXPECT_EQ(shell2, elem_via_side.element);
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));
  }
  {
    // Connectivity for Shell Element 2
    EXPECT_EQ(1u,             elemElemGraph.get_num_connected_elems(shell2));
    impl::ElementViaSidePair elem_via_side = elemElemGraph.get_connected_element_and_via_side(shell2, 0);
    EXPECT_EQ(1,              elem_via_side.side);
    EXPECT_EQ(hex1, elem_via_side.element);
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell2, 0));
  }
  EXPECT_EQ(2u, elemElemGraph.num_edges());
  EXPECT_EQ(0u, elemElemGraph.num_parallel_edges());
}

TEST( ElementGraph, AdjacentHexShellSerial )
{
  //  ID.proc
  //
  //         12.0-----------11.0
  //          /|             /|
  //         / |            / |
  //        /  |           /  |
  //      9.0-----------10.0  |
  //       |   |          |   |
  //       |   |   2.0    |4.0|
  //       |   |          |   |
  //       |  8.0---------|--7.0
  //       |  /|          |  /|
  //       | / |          | / |
  //       |/  |          |/  |
  //      5.0------------6.0  |
  //       |   |          |   |
  //       |   |   1.0    |3.0|
  //       |   |          |   |
  //       |  4.0---------|--3.0
  //       |  /           |  /
  //       | /            | /
  //       |/             |/
  //      1.0------------2.0
  //                      ^
  //                      |
  //                       ---- Single shell element

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  int p_size = stk::parallel_machine_size(pm);

  if(p_size > 1)
  {
    return;
  }

  const unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, pm);
  stk::mesh::BulkData& mesh = *bulkPtr;
  stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();

  stk::mesh::Part * hexPart   = &meta.declare_part_with_topology("hex_part",   stk::topology::HEX_8);
  stk::mesh::Part * shellPart = &meta.declare_part_with_topology("shell_part", stk::topology::SHELL_QUAD_4);
  meta.commit();

  std::vector<stk::mesh::EntityIdVector> hexNodeIDs {
    { 1, 2, 3, 4, 5,  6,  7,  8 },
    { 5, 6, 7, 8, 9, 10, 11, 12 }
  };
  stk::mesh::EntityId hexElemIDs[] = { 1, 2 };

  std::vector<stk::mesh::EntityIdVector> shellNodeIDs {
    { 2, 3,  7,  6 },
    { 6, 7, 11, 10 },
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

  ElemElemGraphTester elemElemGraph(mesh);

  const Entity hex1   = mesh.get_entity(stk::topology::ELEM_RANK, 1);
  const Entity hex2   = mesh.get_entity(stk::topology::ELEM_RANK, 2);
  const Entity shell3 = mesh.get_entity(stk::topology::ELEM_RANK, 3);
  const Entity shell4 = mesh.get_entity(stk::topology::ELEM_RANK, 4);

  // Connectivity for Hex Element 1
  EXPECT_EQ(2u,     elemElemGraph.get_num_connected_elems(hex1));
  expect_correct_connected_element_via_side(elemElemGraph, hex1, 0, hex2, 5);
  expect_correct_connected_element_via_side(elemElemGraph, hex1, 1, shell3, 1);
  EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));
  EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 1));

  // Connectivity for Hex Element 2
  EXPECT_EQ(2u,     elemElemGraph.get_num_connected_elems(hex2));
  expect_correct_connected_element_via_side(elemElemGraph, hex2, 0, hex1, 4);
  expect_correct_connected_element_via_side(elemElemGraph, hex2, 1, shell4, 1);
  EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex2, 0));
  EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex2, 1));

  // Connectivity for Shell Element 3
  EXPECT_EQ(1u,   elemElemGraph.get_num_connected_elems(shell3));
  expect_correct_connected_element_via_side(elemElemGraph, shell3, 0, hex1, 1);
  EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell3, 0));

  // Connectivity for Shell Element 4
  EXPECT_EQ(1u,   elemElemGraph.get_num_connected_elems(shell4));
  expect_correct_connected_element_via_side(elemElemGraph, shell4, 0, hex2, 1);
  EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell4, 0));

  EXPECT_EQ(6u, elemElemGraph.num_edges());
  EXPECT_EQ(0u, elemElemGraph.num_parallel_edges());
}

TEST( ElementGraph, HexShellHexSerial )
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
  //       |   |          |   |          |   |
  //       |  2.0---------|--6.0---------|-10.0
  //       |  /           |  /           |  /
  //       | /            | /            | /
  //       |/             |/             |/
  //      1.0------------5.0------------9.0
  //                      ^
  //                      |
  //                       ---- Single shell element

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  int p_size = stk::parallel_machine_size(pm);

  if(p_size > 1)
  {
    return;
  }

  const unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, pm);
  stk::mesh::BulkData& mesh = *bulkPtr;
  stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();

  stk::mesh::Part * hexPart   = &meta.declare_part_with_topology("hex_part",   stk::topology::HEX_8);
  stk::mesh::Part * shellPart = &meta.declare_part_with_topology("shell_part", stk::topology::SHELL_QUAD_4);
  meta.commit();

  std::vector<stk::mesh::EntityIdVector> hexNodeIDs {
    { 1, 2, 3, 4, 5,  6,  7,  8 },
    { 5, 6, 7, 8, 9, 10, 11, 12 }
  };
  stk::mesh::EntityId hexElemIDs[] = { 1, 2 };

  std::vector<stk::mesh::EntityIdVector> shellNodeIDs {
    { 5, 6, 7, 8 }
  };
  stk::mesh::EntityId shellElemIDs[] = { 3 };

  mesh.modification_begin();
  for (size_t i = 0; i < hexNodeIDs.size(); ++i) {
    stk::mesh::declare_element(mesh, *hexPart, hexElemIDs[i], hexNodeIDs[i]);
  }
  for (size_t i = 0; i < shellNodeIDs.size(); ++i) {
    stk::mesh::declare_element(mesh, *shellPart, shellElemIDs[i], shellNodeIDs[i]);
  }
  mesh.modification_end();

  ElemElemGraphTester elemElemGraph(mesh);

  const Entity hex1   = mesh.get_entity(stk::topology::ELEM_RANK, 1);
  const Entity hex2   = mesh.get_entity(stk::topology::ELEM_RANK, 2);
  const Entity shell3 = mesh.get_entity(stk::topology::ELEM_RANK, 3);

  {
    // Connectivity for Hex Element 1
    EXPECT_EQ(1u,     elemElemGraph.get_num_connected_elems(hex1));
    impl::ElementViaSidePair elem_via_side = elemElemGraph.get_connected_element_and_via_side(hex1, 0);
    EXPECT_EQ(5,      elem_via_side.side);
    EXPECT_EQ(shell3, elem_via_side.element);
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));
  }
  {
    // Connectivity for Hex Element 2
    EXPECT_EQ(1u,     elemElemGraph.get_num_connected_elems(hex2));
    impl::ElementViaSidePair elem_via_side = elemElemGraph.get_connected_element_and_via_side(hex2, 0);
    EXPECT_EQ(4,      elem_via_side.side);
    EXPECT_EQ(shell3, elem_via_side.element);
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex2, 0));
  }
  // Connectivity for Shell Element 3
  EXPECT_EQ(2u,   elemElemGraph.get_num_connected_elems(shell3));
  impl::ElementViaSidePair elem0_via_side = elemElemGraph.get_connected_element_and_via_side(shell3, 0);
  impl::ElementViaSidePair elem1_via_side = elemElemGraph.get_connected_element_and_via_side(shell3, 1);
  EXPECT_EQ(1,    elem0_via_side.side);
  EXPECT_EQ(0,    elem1_via_side.side);
  EXPECT_EQ(hex1, elem0_via_side.element);
  EXPECT_EQ(hex2, elem1_via_side.element);
  EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell3, 0));
  EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(shell3, 1));

  EXPECT_EQ(4u, elemElemGraph.num_edges());
  EXPECT_EQ(0u, elemElemGraph.num_parallel_edges());
}

TEST( ElementGraph, Hex0Hex0Hex1Parallel )
{
  //  ID.proc
  //
  //          3.0------------7.0-----------11.0-----------15.1
  //          /|             /|             /|             /|
  //         / |            / |            / |            / |
  //        /  |           /  |           /  |           /  |
  //      4.0------------8.0-----------12.0-----------16.1  |
  //       |   |          |   |          |   |          |   |
  //       |   |   1.0    |   |   2.0    |   |   3.1    |   |
  //       |   |          |   |          |   |          |   |
  //       |  2.0---------|--6.0---------|-10.0---------|-14.1
  //       |  /           |  /           |  /           |  /
  //       | /            | /            | /            | /
  //       |/             |/             |/             |/
  //      1.0------------5.0------------9.0-----------13.1
  //

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  unsigned p_size = stk::parallel_machine_size(pm);
  unsigned p_rank = stk::parallel_machine_rank(pm);

  if(p_size != 2u)
  {
    return;
  }

  const unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, pm);
  stk::mesh::BulkData& mesh = *bulkPtr;
  stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();

  stk::mesh::Part * hexPart = &meta.declare_part_with_topology("hex_part", stk::topology::HEX_8);
  meta.commit();

  std::vector<stk::mesh::EntityIdVector> hexNodeIDs {
    { 1,  2,  3,  4,  5,  6,  7,  8 },
    { 5,  6,  7,  8,  9, 10, 11, 12 },
    { 9, 10, 11, 12, 13, 14, 15, 16 }
  };
  stk::mesh::EntityId hexElemIDs[] = { 1, 2, 3 };
  stk::mesh::EntityId hexElemOwningProc[] = { 0, 0, 1 };

  // list of triplets: (owner-proc, shared-nodeID, sharing-proc)
  std::vector< std::vector<unsigned> > shared_nodeIDs_and_procs
  {
    { 0,  9, 1 },  // proc 0
    { 0, 10, 1 },
    { 0, 11, 1 },
    { 0, 12, 1 },
    { 1,  9, 0 },  // proc 1
    { 1, 10, 0 },
    { 1, 11, 0 },
    { 1, 12, 0 }
  };

  mesh.modification_begin();
  for (size_t i = 0; i < hexNodeIDs.size(); ++i) {
    if (hexElemOwningProc[i] == p_rank) {
      stk::mesh::declare_element(mesh, *hexPart, hexElemIDs[i], hexNodeIDs[i]);
    }
  }
  setup_node_sharing(mesh, shared_nodeIDs_and_procs );
  mesh.modification_end();

  ElemElemGraphTester elemElemGraph(mesh);

  const Entity hex1 = mesh.get_entity(stk::topology::ELEM_RANK, 1);
  const Entity hex2 = mesh.get_entity(stk::topology::ELEM_RANK, 2);
  const Entity hex3 = mesh.get_entity(stk::topology::ELEM_RANK, 3);

  if (p_rank == 0) {
    {
      // Connectivity for Hex Element 1
      EXPECT_EQ(1u,   elemElemGraph.get_num_connected_elems(hex1));
      impl::ElementViaSidePair elem_via_side = elemElemGraph.get_connected_element_and_via_side(hex1, 0);
      EXPECT_EQ(5,    elem_via_side.side);
      EXPECT_EQ(hex2, elem_via_side.element);
      EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));
    }
    {
      // Connectivity for Hex Element 2
      EXPECT_EQ(2u,   elemElemGraph.get_num_connected_elems(hex2));
      impl::ElementViaSidePair elem_via_side = elemElemGraph.get_connected_element_and_via_side(hex2, 0);
      EXPECT_EQ(4,    elem_via_side.side);
      EXPECT_EQ(hex1, elem_via_side.element);
      EXPECT_EQ(5,    elemElemGraph.get_connected_remote_id_and_via_side(hex2, 1).side);
      EXPECT_EQ(3u,   elemElemGraph.get_connected_remote_id_and_via_side(hex2, 1).id);
      EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex2, 0));
      EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex2, 1));
    }
    EXPECT_EQ(3u, elemElemGraph.num_edges());
    EXPECT_EQ(1u, elemElemGraph.num_parallel_edges());
  }
  else if (p_rank == 1) {
    // Connectivity for Hex Element 3
    EXPECT_EQ(1u,   elemElemGraph.get_num_connected_elems(hex3));
    impl::IdViaSidePair remote_elem_id_via_side = elemElemGraph.get_connected_remote_id_and_via_side(hex3, 0);
    EXPECT_EQ(2u, remote_elem_id_via_side.id);
    EXPECT_EQ(4, remote_elem_id_via_side.side);
    EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex3, 0));

    EXPECT_EQ(1u, elemElemGraph.num_edges());
    EXPECT_EQ(1u, elemElemGraph.num_parallel_edges());
  }
}

TEST( ElementGraph, OneHex )
{
  if (1 == stk::parallel_machine_size(MPI_COMM_WORLD))
  {
    const unsigned spatialDim = 3;
    std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, MPI_COMM_WORLD, stk::mesh::BulkData::NO_AUTO_AURA);
    stk::mesh::BulkData& bulk = *bulkPtr;
    stk::io::fill_mesh("generated:1x1x1", bulk);

    bulk.get_face_adjacent_element_graph().write_graph(std::cerr, "Before");

    bulk.modification_begin();
    stk::mesh::Entity elem1 = bulk.get_entity(stk::topology::ELEM_RANK,1);
    bulk.destroy_entity(elem1);
    bulk.modification_end();

    bulk.get_face_adjacent_element_graph().write_graph(std::cerr, "After");
    EXPECT_EQ(0u, bulk.get_face_adjacent_element_graph().size());
  }
}

TEST( ElementGraph, Hex0Hex1Hex0Parallel )
{
  //  ID.proc
  //
  //          3.0------------7.0-----------11.0-----------15.0
  //          /|             /|             /|             /|
  //         / |            / |            / |            / |
  //        /  |           /  |           /  |           /  |
  //      4.0------------8.0-----------12.0-----------16.0  |
  //       |   |          |   |          |   |          |   |
  //       |   |   1.0    |   |   2.1    |   |   3.0    |   |
  //       |   |          |   |          |   |          |   |
  //       |  2.0---------|--6.0---------|-10.0---------|-14.0
  //       |  /           |  /           |  /           |  /
  //       | /            | /            | /            | /
  //       |/             |/             |/             |/
  //      1.0------------5.0------------9.0-----------13.0
  //

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  unsigned p_size = stk::parallel_machine_size(pm);
  unsigned p_rank = stk::parallel_machine_rank(pm);

  if(p_size != 2u)
  {
    return;
  }

  const unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, pm);
  stk::mesh::BulkData& mesh = *bulkPtr;
  stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();

  stk::mesh::Part * hexPart = &meta.declare_part_with_topology("hex_part", stk::topology::HEX_8);
  meta.commit();

  std::vector<stk::mesh::EntityIdVector> hexNodeIDs {
    { 1,  2,  3,  4,  5,  6,  7,  8 },
    { 5,  6,  7,  8,  9, 10, 11, 12 },
    { 9, 10, 11, 12, 13, 14, 15, 16 }
  };
  stk::mesh::EntityId hexElemIDs[] = { 1, 2, 3 };
  stk::mesh::EntityId hexElemOwningProc[] = { 0, 1, 0 };

  // list of triplets: (owner-proc, shared-nodeID, sharing-proc)
  std::vector< std::vector<unsigned> > shared_nodeIDs_and_procs
  {
    { 0,  5, 1 },  // proc 0
    { 0,  6, 1 },
    { 0,  7, 1 },
    { 0,  8, 1 },
    { 0,  9, 1 },
    { 0, 10, 1 },
    { 0, 11, 1 },
    { 0, 12, 1 },
    { 1,  5, 0 },  // proc 1
    { 1,  6, 0 },
    { 1,  7, 0 },
    { 1,  8, 0 },
    { 1,  9, 0 },
    { 1, 10, 0 },
    { 1, 11, 0 },
    { 1, 12, 0 }
  };

  mesh.modification_begin();
  for (size_t i = 0; i < hexNodeIDs.size(); ++i) {
    if (hexElemOwningProc[i] == p_rank) {
      stk::mesh::declare_element(mesh, *hexPart, hexElemIDs[i], hexNodeIDs[i]);
    }
  }
  setup_node_sharing(mesh, shared_nodeIDs_and_procs );
  mesh.modification_end();

  ElemElemGraphTester elemElemGraph(mesh);

  const Entity hex1 = mesh.get_entity(stk::topology::ELEM_RANK, 1);
  const Entity hex2 = mesh.get_entity(stk::topology::ELEM_RANK, 2);
  const Entity hex3 = mesh.get_entity(stk::topology::ELEM_RANK, 3);

  if (p_rank == 0) {
    // Connectivity for Hex Element 1
    EXPECT_EQ(1u, elemElemGraph.get_num_connected_elems(hex1));
    EXPECT_EQ(5,  elemElemGraph.get_connected_remote_id_and_via_side(hex1, 0).side);
    EXPECT_EQ(2u, elemElemGraph.get_connected_remote_id_and_via_side(hex1, 0).id);
    EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));

    // Connectivity for Hex Element 3
    EXPECT_EQ(1u, elemElemGraph.get_num_connected_elems(hex3));
    EXPECT_EQ(4,  elemElemGraph.get_connected_remote_id_and_via_side(hex3, 0).side);
    EXPECT_EQ(2u, elemElemGraph.get_connected_remote_id_and_via_side(hex3, 0).id);
    EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex3, 0));
  }
  else if (p_rank == 1) {
    // Connectivity for Hex Element 2
    EXPECT_EQ(2u, elemElemGraph.get_num_connected_elems(hex2));
    EXPECT_EQ(4,  elemElemGraph.get_connected_remote_id_and_via_side(hex2, 0).side);
    EXPECT_EQ(5,  elemElemGraph.get_connected_remote_id_and_via_side(hex2, 1).side);
    EXPECT_EQ(1u, elemElemGraph.get_connected_remote_id_and_via_side(hex2, 0).id);
    EXPECT_EQ(3u, elemElemGraph.get_connected_remote_id_and_via_side(hex2, 1).id);
    EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex2, 0));
    EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex2, 1));
  }

  EXPECT_EQ(2u, elemElemGraph.num_edges());
  EXPECT_EQ(2u, elemElemGraph.num_parallel_edges());
}

TEST( ElementGraph, Hex0Hex1Hex2Parallel )
{
  //  ID.proc
  //
  //          3.0------------7.0-----------11.1-----------15.2
  //          /|             /|             /|             /|
  //         / |            / |            / |            / |
  //        /  |           /  |           /  |           /  |
  //      4.0------------8.0-----------12.1-----------16.2  |
  //       |   |          |   |          |   |          |   |
  //       |   |   1.0    |   |   2.1    |   |   3.2    |   |
  //       |   |          |   |          |   |          |   |
  //       |  2.0---------|--6.0---------|-10.1---------|-14.2
  //       |  /           |  /           |  /           |  /
  //       | /            | /            | /            | /
  //       |/             |/             |/             |/
  //      1.0------------5.0------------9.1-----------13.2
  //

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  unsigned p_size = stk::parallel_machine_size(pm);
  unsigned p_rank = stk::parallel_machine_rank(pm);

  if(p_size != 3u)
  {
    return;
  }

  const unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, pm);
  stk::mesh::BulkData& mesh = *bulkPtr;
  stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();

  stk::mesh::Part * hexPart = &meta.declare_part_with_topology("hex_part", stk::topology::HEX_8);
  meta.commit();

  std::vector<stk::mesh::EntityIdVector> hexNodeIDs {
    { 1,  2,  3,  4,  5,  6,  7,  8 },
    { 5,  6,  7,  8,  9, 10, 11, 12 },
    { 9, 10, 11, 12, 13, 14, 15, 16 }
  };
  stk::mesh::EntityId hexElemIDs[] = { 1, 2, 3 };
  stk::mesh::EntityId hexElemOwningProc[] = { 0, 1, 2 };

  // list of triplets: (owner-proc, shared-nodeID, sharing-proc)
  std::vector< std::vector<unsigned> > shared_nodeIDs_and_procs
  {
    { 0,  5, 1 },  // proc 0
    { 0,  6, 1 },
    { 0,  7, 1 },
    { 0,  8, 1 },
    { 1,  5, 0 },  // proc 1
    { 1,  6, 0 },
    { 1,  7, 0 },
    { 1,  8, 0 },
    { 1,  9, 2 },
    { 1, 10, 2 },
    { 1, 11, 2 },
    { 1, 12, 2 },
    { 2,  9, 1 },  // proc 2
    { 2, 10, 1 },
    { 2, 11, 1 },
    { 2, 12, 1 }
  };

  mesh.modification_begin();
  for (size_t i = 0; i < hexNodeIDs.size(); ++i) {
    if (hexElemOwningProc[i] == p_rank) {
      stk::mesh::declare_element(mesh, *hexPart, hexElemIDs[i], hexNodeIDs[i]);
    }
  }
  setup_node_sharing(mesh, shared_nodeIDs_and_procs );
  mesh.modification_end();

  ElemElemGraphTester elemElemGraph(mesh);

  const Entity hex1 = mesh.get_entity(stk::topology::ELEM_RANK, 1);
  const Entity hex2 = mesh.get_entity(stk::topology::ELEM_RANK, 2);
  const Entity hex3 = mesh.get_entity(stk::topology::ELEM_RANK, 3);

  if (p_rank == 0) {
    // Connectivity for Hex Element 1
    EXPECT_EQ(1u, elemElemGraph.get_num_connected_elems(hex1));
    EXPECT_EQ(5,  elemElemGraph.get_connected_remote_id_and_via_side(hex1, 0).side);
    EXPECT_EQ(2u, elemElemGraph.get_connected_remote_id_and_via_side(hex1, 0).id);
    EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));
    EXPECT_EQ(1u, elemElemGraph.num_edges());
    EXPECT_EQ(1u, elemElemGraph.num_parallel_edges());
  }
  else if (p_rank == 1) {
    // Connectivity for Hex Element 2
    EXPECT_EQ(2u, elemElemGraph.get_num_connected_elems(hex2));
    EXPECT_EQ(4,  elemElemGraph.get_connected_remote_id_and_via_side(hex2, 0).side);
    EXPECT_EQ(5,  elemElemGraph.get_connected_remote_id_and_via_side(hex2, 1).side);
    EXPECT_EQ(1u, elemElemGraph.get_connected_remote_id_and_via_side(hex2, 0).id);
    EXPECT_EQ(3u, elemElemGraph.get_connected_remote_id_and_via_side(hex2, 1).id);
    EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex2, 0));
    EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex2, 1));
    EXPECT_EQ(2u, elemElemGraph.num_edges());
    EXPECT_EQ(2u, elemElemGraph.num_parallel_edges());
  }
  else if (p_rank == 2) {
    // Connectivity for Hex Element 3
    EXPECT_EQ(1u, elemElemGraph.get_num_connected_elems(hex3));
    EXPECT_EQ(4,  elemElemGraph.get_connected_remote_id_and_via_side(hex3, 0).side);
    EXPECT_EQ(2u, elemElemGraph.get_connected_remote_id_and_via_side(hex3, 0).id);
    EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex3, 0));
    EXPECT_EQ(1u, elemElemGraph.num_edges());
    EXPECT_EQ(1u, elemElemGraph.num_parallel_edges());
  }
}

TEST( ElementGraph, Hex0Shell1Parallel )
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
  //       |   |          |   |
  //       |  2.0---------|--6.0
  //       |  /           |  /
  //       | /            | /
  //       |/             |/
  //      1.0------------5.0
  //                      ^
  //                      |
  //                       ---- Single shell element

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  unsigned p_size = stk::parallel_machine_size(pm);
  unsigned p_rank = stk::parallel_machine_rank(pm);

  if(p_size != 2u)
  {
    return;
  }

  const unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, pm);
  stk::mesh::BulkData& mesh = *bulkPtr;
  stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();

  stk::mesh::Part * hexPart   = &meta.declare_part_with_topology("hex_part",   stk::topology::HEX_8);
  stk::mesh::Part * shellPart = &meta.declare_part_with_topology("shell_part", stk::topology::SHELL_QUAD_4);
  meta.commit();

  std::vector<stk::mesh::EntityIdVector> hexNodeIDs {
    { 1, 2, 3, 4, 5,  6,  7,  8 }
  };
  stk::mesh::EntityId hexElemIDs[] = { 1 };
  stk::mesh::EntityId hexElemOwningProc[] = { 0 };

  std::vector<stk::mesh::EntityIdVector> shellNodeIDs {
    { 5, 6, 7, 8 }
  };
  stk::mesh::EntityId shellElemIDs[] = { 2 };
  stk::mesh::EntityId shellElemOwningProc[] = { 1 };

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

  ElemElemGraphTester elemElemGraph(mesh);

  const Entity hex1   = mesh.get_entity(stk::topology::ELEM_RANK, 1);
  const Entity shell2 = mesh.get_entity(stk::topology::ELEM_RANK, 2);

  if (p_rank == 0) {
    // Connectivity for Hex Element 1
    EXPECT_EQ(1u, elemElemGraph.num_edges());
    EXPECT_EQ(1u, elemElemGraph.num_parallel_edges());
    EXPECT_EQ(1u, elemElemGraph.get_num_connected_elems(hex1));
    EXPECT_EQ(5,  elemElemGraph.get_connected_remote_id_and_via_side(hex1, 0).side);
    EXPECT_EQ(2u, elemElemGraph.get_connected_remote_id_and_via_side(hex1, 0).id);
    EXPECT_EQ(1,  elemElemGraph.get_owning_proc_id_of_remote_element(hex1, 0));
    EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));
  }
  else if (p_rank == 1) {
    // Connectivity for Shell Element 2
    EXPECT_EQ(1u, elemElemGraph.num_edges());
    EXPECT_EQ(1u, elemElemGraph.num_parallel_edges());
    EXPECT_EQ(1u, elemElemGraph.get_num_connected_elems(shell2));
    EXPECT_EQ(1,  elemElemGraph.get_connected_remote_id_and_via_side(shell2, 0).side);
    EXPECT_EQ(1u, elemElemGraph.get_connected_remote_id_and_via_side(shell2, 0).id);
    EXPECT_EQ(0,  elemElemGraph.get_owning_proc_id_of_remote_element(shell2, 0));
    EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell2, 0));
  }

  EXPECT_EQ(1u, elemElemGraph.num_edges());
  EXPECT_EQ(1u, elemElemGraph.num_parallel_edges());
}

TEST( ElementGraph, Hex0DelShell1Parallel )
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
  //       |   |          |   |
  //       |  2.0---------|--6.0
  //       |  /           |  /
  //       | /            | /
  //       |/             |/
  //      1.0------------5.0
  //                      ^
  //                      |
  //                       ---- Deleted single shell element

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  unsigned p_size = stk::parallel_machine_size(pm);
  unsigned p_rank = stk::parallel_machine_rank(pm);

  if(p_size != 2u)
  {
    return;
  }

  const unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, pm);
  stk::mesh::BulkData& mesh = *bulkPtr;
  stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();

  stk::mesh::Part * hexPart   = &meta.declare_part_with_topology("hex_part",   stk::topology::HEX_8);
  stk::mesh::Part * shellPart = &meta.declare_part_with_topology("shell_part", stk::topology::SHELL_QUAD_4);
  meta.commit();

  std::vector<stk::mesh::EntityIdVector> hexNodeIDs {
    { 1, 2, 3, 4, 5,  6,  7,  8 }
  };
  stk::mesh::EntityId hexElemIDs[] = { 1 };
  stk::mesh::EntityId hexElemOwningProc[] = { 0 };

  std::vector<stk::mesh::EntityIdVector> shellNodeIDs {
    { 5, 6, 7, 8 }
  };
  stk::mesh::EntityId shellElemIDs[] = { 2 };
  stk::mesh::EntityId shellElemOwningProc[] = { 1 };

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

  ElemElemGraphTester elemElemGraph(mesh);

  const Entity hex1   = mesh.get_entity(stk::topology::ELEM_RANK, 1);
  const Entity shell2 = mesh.get_entity(stk::topology::ELEM_RANK, 2);

  stk::mesh::impl::DeletedElementInfoVector elements_to_delete;
  if (p_rank == 1) {
    elements_to_delete.push_back({shell2, 2, mesh.bucket(shell2).topology().is_shell()});
  }

  elemElemGraph.delete_elements( elements_to_delete );

  mesh.modification_begin();
  if (p_rank == 1) {
    mesh.destroy_entity(shell2);
  }
  mesh.modification_end();

  if (p_rank == 0) {
    // Connectivity for Hex Element 1
    EXPECT_EQ(0u, elemElemGraph.num_edges());
    EXPECT_EQ(0u, elemElemGraph.num_parallel_edges());
    EXPECT_EQ(0u, elemElemGraph.get_num_connected_elems(hex1));
  }
  else if (p_rank == 1) {
    EXPECT_EQ(0u, elemElemGraph.num_edges());
    EXPECT_EQ(0u, elemElemGraph.num_parallel_edges());
    EXPECT_EQ(0u, elemElemGraph.size());
  }

  EXPECT_EQ(0u, elemElemGraph.num_edges());
  EXPECT_EQ(0u, elemElemGraph.num_parallel_edges());
}

TEST( ElementGraph, Hex0AddShell1Parallel )
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
  //       |   |          |   |
  //       |  2.0---------|--6.0
  //       |  /           |  /
  //       | /            | /
  //       |/             |/
  //      1.0------------5.0
  //                      ^
  //                      |
  //                       ---- Add a single shell element

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  unsigned p_size = stk::parallel_machine_size(pm);
  unsigned p_rank = stk::parallel_machine_rank(pm);

  if(p_size != 2u)
  {
    return;
  }

  const unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, pm);
  stk::mesh::BulkData& mesh = *bulkPtr;
  stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();

  stk::mesh::Part * hexPart   = &meta.declare_part_with_topology("hex_part",   stk::topology::HEX_8);
  stk::mesh::Part * shellPart = &meta.declare_part_with_topology("shell_part", stk::topology::SHELL_QUAD_4);
  meta.commit();

  std::vector<stk::mesh::EntityIdVector> hexNodeIDs {
    { 1, 2, 3, 4, 5, 6, 7, 8 }
  };
  stk::mesh::EntityId hexElemIDs[] = { 1 };
  stk::mesh::EntityId hexElemOwningProc[] = { 0 };

  std::vector<stk::mesh::EntityIdVector> shellNodeIDs {
    { 5, 6, 7, 8 }
  };
  stk::mesh::EntityId shellElemIDs[] = { 2 };
  stk::mesh::EntityId shellElemOwningProc[] = { 1 };

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
  mesh.modification_end();

  ElemElemGraphTester elemElemGraph(mesh);

  EXPECT_EQ(0u, elemElemGraph.num_edges());
  EXPECT_EQ(0u, elemElemGraph.num_parallel_edges());

  mesh.modification_begin();
  for (size_t i = 0; i < shellNodeIDs.size(); ++i) {
    if (shellElemOwningProc[i] == p_rank) {
      stk::mesh::declare_element(mesh, *shellPart, shellElemIDs[i], shellNodeIDs[i]);
    }
  }
  setup_node_sharing(mesh, shared_nodeIDs_and_procs);
  mesh.modification_end();

  const Entity hex1   = mesh.get_entity(stk::topology::ELEM_RANK, 1);
  const Entity shell2 = mesh.get_entity(stk::topology::ELEM_RANK, 2);


  stk::mesh::EntityVector elements_to_add;
  if (p_rank == 1) {
    elements_to_add.push_back(shell2);
  }
  elemElemGraph.add_elements(elements_to_add);

  if (p_rank == 0) {
    // Connectivity for Hex Element 1
    EXPECT_EQ(1u, elemElemGraph.num_edges());
    EXPECT_EQ(1u, elemElemGraph.num_parallel_edges());
    EXPECT_EQ(1u, elemElemGraph.get_num_connected_elems(hex1));
    EXPECT_EQ(5,  elemElemGraph.get_connected_remote_id_and_via_side(hex1, 0).side);
    EXPECT_EQ(2u, elemElemGraph.get_connected_remote_id_and_via_side(hex1, 0).id);
    EXPECT_EQ(1,  elemElemGraph.get_owning_proc_id_of_remote_element(hex1, 0));
    EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));
  }
  else if (p_rank == 1) {
    // Connectivity for Shell Element 2
    EXPECT_EQ(1u, elemElemGraph.num_edges());
    EXPECT_EQ(1u, elemElemGraph.num_parallel_edges());
    EXPECT_EQ(1u, elemElemGraph.get_num_connected_elems(shell2));
    EXPECT_EQ(1,  elemElemGraph.get_connected_remote_id_and_via_side(shell2, 0).side);
    EXPECT_EQ(1u, elemElemGraph.get_connected_remote_id_and_via_side(shell2, 0).id);
    EXPECT_EQ(0,  elemElemGraph.get_owning_proc_id_of_remote_element(shell2, 0));
    EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell2, 0));
  }
}


TEST( ElementGraph, Hex0AddShell0Hex1Parallel )
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
  //       |   |          |   |          |   |
  //       |  2.0---------|--6.0---------|-10.1
  //       |  /           |  /           |  /
  //       | /            | /            | /
  //       |/             |/             |/
  //      1.0------------5.0------------9.1
  //                      ^
  //                      |
  //                       ---- Add a Single shell element here
  //
  //      side 0 of the shell points right to element 2
  //      side 1 of the shell points left to element 1

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  unsigned p_size = stk::parallel_machine_size(pm);
  unsigned p_rank = stk::parallel_machine_rank(pm);

  if(p_size != 2u)
  {
    return;
  }

  const unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, pm);
  stk::mesh::BulkData& mesh = *bulkPtr;
  stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();

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
    { 5, 6, 7, 8 }
  };
  stk::mesh::EntityId shellElemIDs[] = { 3 };
  stk::mesh::EntityId shellElemOwningProc[] = { 0 };

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
  setup_node_sharing(mesh, shared_nodeIDs_and_procs );
  mesh.modification_end();


  ElemElemGraphTester elemElemGraph(mesh);

  EXPECT_EQ(1u, elemElemGraph.num_edges());
  EXPECT_EQ(1u, elemElemGraph.num_parallel_edges());

  mesh.modification_begin();
  for (size_t i = 0; i < shellNodeIDs.size(); ++i) {
    if (shellElemOwningProc[i] == p_rank) {
      stk::mesh::declare_element(mesh, *shellPart, shellElemIDs[i], shellNodeIDs[i]);
    }
  }
  mesh.modification_end();


  const Entity hex1   = mesh.get_entity(stk::topology::ELEM_RANK, 1);
  const Entity hex2   = mesh.get_entity(stk::topology::ELEM_RANK, 2);
  const Entity shell3 = mesh.get_entity(stk::topology::ELEM_RANK, 3);

  stk::mesh::EntityVector addVector;
  if (0 == p_rank) {
    addVector.push_back(shell3);
  }

  elemElemGraph.add_elements(addVector);

  if (p_rank == 0) {
    // Connectivity for Hex Element 1
    EXPECT_EQ(1u,     elemElemGraph.get_num_connected_elems(hex1));
    expect_correct_connected_element_via_side(elemElemGraph, hex1, 0, shell3, 5);
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));

    // Connectivity for Shell Element 3
    EXPECT_EQ(2u,   elemElemGraph.get_num_connected_elems(shell3));
    expect_correct_connected_element_via_side(elemElemGraph, shell3, 0, hex1, 1);
    EXPECT_EQ(0,    elemElemGraph.get_connected_remote_id_and_via_side(shell3, 1).side);
    EXPECT_EQ(2u,   elemElemGraph.get_connected_remote_id_and_via_side(shell3, 1).id);
    EXPECT_TRUE (elemElemGraph.is_connected_elem_locally_owned(shell3, 0));
    EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell3, 1));
    EXPECT_EQ(3u, elemElemGraph.num_edges());
    EXPECT_EQ(1u, elemElemGraph.num_parallel_edges());
  }
  else if (p_rank == 1) {
    // Connectivity for Hex Element 2
    EXPECT_EQ(1u, elemElemGraph.get_num_connected_elems(hex2));
    EXPECT_EQ(4,  elemElemGraph.get_connected_remote_id_and_via_side(hex2, 0).side);
    EXPECT_EQ(3u, elemElemGraph.get_connected_remote_id_and_via_side(hex2, 0).id);
    EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex2, 0));
    EXPECT_EQ(1u, elemElemGraph.num_edges());
    EXPECT_EQ(1u, elemElemGraph.num_parallel_edges());
  }
}

TEST( ElementGraph, Hex0AddShell1Hex2Parallel )
{
  //  ID.proc
  //
  //          3.0------------7.0-----------11.2
  //          /|             /|             /|
  //         / |            / |            / |
  //        /  |           /  |           /  |
  //      4.0------------8.0-----------12.2  |
  //       |   |          |   |          |   |
  //       |   |   1.0    |3.1|   2.2    |   |
  //       |   |          |   |          |   |
  //       |  2.0---------|--6.0---------|-10.2
  //       |  /           |  /           |  /
  //       | /            | /            | /
  //       |/             |/             |/
  //      1.0------------5.0------------9.2
  //                      ^
  //                      |
  //                       ---- Add a Single shell element here

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  unsigned p_size = stk::parallel_machine_size(pm);
  unsigned p_rank = stk::parallel_machine_rank(pm);

  if(p_size != 3u)
  {
    return;
  }

  const unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, pm);
  stk::mesh::BulkData& mesh = *bulkPtr;
  stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();

  stk::mesh::Part * hexPart   = &meta.declare_part_with_topology("hex_part",   stk::topology::HEX_8);
  stk::mesh::Part * shellPart = &meta.declare_part_with_topology("shell_part", stk::topology::SHELL_QUAD_4);
  meta.commit();

  std::vector<stk::mesh::EntityIdVector> hexNodeIDs {
    { 1, 2, 3, 4, 5,  6,  7,  8 },
    { 5, 6, 7, 8, 9, 10, 11, 12 }
  };
  stk::mesh::EntityId hexElemIDs[] = { 1, 2 };
  stk::mesh::EntityId hexElemOwningProc[] = { 0, 2 };

  std::vector<stk::mesh::EntityIdVector> shellNodeIDs {
    { 5, 6, 7, 8 }
  };
  stk::mesh::EntityId shellElemIDs[] = { 3 };
  stk::mesh::EntityId shellElemOwningProc[] = { 1 };

  // list of triplets: (owner-proc, shared-nodeID, sharing-proc)
  std::vector< std::vector<unsigned> > shared_nodeIDs_and_procs
  {
    { 0, 5, 1 },  // proc 0
    { 0, 6, 1 },
    { 0, 7, 1 },
    { 0, 8, 1 },
    { 0, 5, 2 },  // proc 0
    { 0, 6, 2 },
    { 0, 7, 2 },
    { 0, 8, 2 },
    { 1, 5, 0 },  // proc 1
    { 1, 6, 0 },
    { 1, 7, 0 },
    { 1, 8, 0 },
    { 1, 5, 2 },
    { 1, 6, 2 },
    { 1, 7, 2 },
    { 1, 8, 2 },
    { 2, 5, 1 },  // proc 2
    { 2, 6, 1 },
    { 2, 7, 1 },
    { 2, 8, 1 },
    { 2, 5, 0 },  // proc 2
    { 2, 6, 0 },
    { 2, 7, 0 },
    { 2, 8, 0 }
  };
  mesh.modification_begin();
  for (size_t i = 0; i < hexNodeIDs.size(); ++i) {
    if (hexElemOwningProc[i] == p_rank) {
      stk::mesh::declare_element(mesh, *hexPart, hexElemIDs[i], hexNodeIDs[i]);
    }
  }
  if (1 == p_rank) {
    mesh.declare_node(5);
    mesh.declare_node(6);
    mesh.declare_node(7);
    mesh.declare_node(8);
  }
  setup_node_sharing(mesh, shared_nodeIDs_and_procs );
  mesh.modification_end();

  ElemElemGraphTester elemElemGraph(mesh);

  if (0 == p_rank || 2 == p_rank) {
    EXPECT_EQ(1u, elemElemGraph.num_edges());
    EXPECT_EQ(1u, elemElemGraph.num_parallel_edges());
  }
  else {
    EXPECT_EQ(0u, elemElemGraph.num_edges());
    EXPECT_EQ(0u, elemElemGraph.num_parallel_edges());
  }

  mesh.modification_begin();
  for (size_t i = 0; i < shellNodeIDs.size(); ++i) {
    if (shellElemOwningProc[i] == p_rank) {
      stk::mesh::declare_element(mesh, *shellPart, shellElemIDs[i], shellNodeIDs[i]);
    }
  }
  mesh.modification_end();

  const Entity hex1   = mesh.get_entity(stk::topology::ELEM_RANK, 1);
  const Entity hex2   = mesh.get_entity(stk::topology::ELEM_RANK, 2);
  const Entity shell3 = mesh.get_entity(stk::topology::ELEM_RANK, 3);

  stk::mesh::EntityVector addVector;
  if (1 == p_rank) {
    addVector.push_back(shell3);
  }

  elemElemGraph.add_elements(addVector);

  if (p_rank == 0) {
    // Connectivity for Hex Element 1
    EXPECT_EQ(1u, elemElemGraph.get_num_connected_elems(hex1));
    EXPECT_EQ(5,  elemElemGraph.get_connected_remote_id_and_via_side(hex1, 0).side);
    EXPECT_EQ(3u, elemElemGraph.get_connected_remote_id_and_via_side(hex1, 0).id);
    EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));
    EXPECT_EQ(1u, elemElemGraph.num_edges());
    EXPECT_EQ(1u, elemElemGraph.num_parallel_edges());
  }
  else if (p_rank == 1) {
    // Connectivity for Shell Element 3
    EXPECT_EQ(2u, elemElemGraph.get_num_connected_elems(shell3));
    EXPECT_EQ(1,  elemElemGraph.get_connected_remote_id_and_via_side(shell3, 0).side);
    EXPECT_EQ(0,  elemElemGraph.get_connected_remote_id_and_via_side(shell3, 1).side);
    EXPECT_EQ(1u, elemElemGraph.get_connected_remote_id_and_via_side(shell3, 0).id);
    EXPECT_EQ(2u, elemElemGraph.get_connected_remote_id_and_via_side(shell3, 1).id);
    EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell3, 0));
    EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell3, 1));
    EXPECT_EQ(2u, elemElemGraph.num_edges());
    EXPECT_EQ(2u, elemElemGraph.num_parallel_edges());
  }
  else if (p_rank == 2) {
    // Connectivity for Hex Element 2
    EXPECT_EQ(1u, elemElemGraph.get_num_connected_elems(hex2));
    EXPECT_EQ(4,  elemElemGraph.get_connected_remote_id_and_via_side(hex2, 0).side);
    EXPECT_EQ(3u, elemElemGraph.get_connected_remote_id_and_via_side(hex2, 0).id);
    EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex2, 0));
    EXPECT_EQ(1u, elemElemGraph.num_edges());
    EXPECT_EQ(1u, elemElemGraph.num_parallel_edges());
  }
}

TEST( ElementGraph, Hex0Shell1AddHex2Parallel )
{
  //  ID.proc
  //
  //          3.0------------7.0-----------11.2
  //          /|             /|             /|
  //         / |            / |            / |
  //        /  |           /  |           /  |
  //      4.0------------8.0-----------12.2  |
  //       |   |          |   |          |   |
  //       |   |   1.0    |3.1|   2.2    |   |  <-- Add this hex on the right
  //       |   |          |   |          |   |
  //       |  2.0---------|--6.0---------|-10.2
  //       |  /           |  /           |  /
  //       | /            | /            | /
  //       |/             |/             |/
  //      1.0------------5.0------------9.2
  //                      ^
  //                      |
  //                       ---- Single shell element here

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  unsigned p_size = stk::parallel_machine_size(pm);
  unsigned p_rank = stk::parallel_machine_rank(pm);

  if(p_size != 3u)
  {
    return;
  }

  const unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, pm);
  stk::mesh::BulkData& mesh = *bulkPtr;
  stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();

  stk::mesh::Part * hexPart   = &meta.declare_part_with_topology("hex_part",   stk::topology::HEX_8);
  stk::mesh::Part * shellPart = &meta.declare_part_with_topology("shell_part", stk::topology::SHELL_QUAD_4);
  meta.commit();

  std::vector<stk::mesh::EntityIdVector> hexNodeIDs {
    { 1, 2, 3, 4, 5,  6,  7,  8 },
    { 5, 6, 7, 8, 9, 10, 11, 12 }
  };
  stk::mesh::EntityId hexElemIDs[] = { 1, 2 };
  stk::mesh::EntityId hexElemOwningProc[] = { 0, 2 };

  std::vector<stk::mesh::EntityIdVector> shellNodeIDs {
    { 5, 6, 7, 8 }
  };
  stk::mesh::EntityId shellElemIDs[] = { 3 };
  stk::mesh::EntityId shellElemOwningProc[] = { 1 };

  // list of triplets: (owner-proc, shared-nodeID, sharing-proc)
  std::vector< std::vector<unsigned> > shared_nodeIDs_and_procs
  {
    { 0, 5, 1 },  // proc 0
    { 0, 6, 1 },
    { 0, 7, 1 },
    { 0, 8, 1 },
    { 0, 5, 2 },
    { 0, 6, 2 },
    { 0, 7, 2 },
    { 0, 8, 2 },
    { 1, 5, 0 },  // proc 1
    { 1, 6, 0 },
    { 1, 7, 0 },
    { 1, 8, 0 },
    { 1, 5, 2 },
    { 1, 6, 2 },
    { 1, 7, 2 },
    { 1, 8, 2 },
    { 2, 5, 1 },  // proc 2
    { 2, 6, 1 },
    { 2, 7, 1 },
    { 2, 8, 1 },
    { 2, 5, 0 },
    { 2, 6, 0 },
    { 2, 7, 0 },
    { 2, 8, 0 }
  };
  mesh.modification_begin();
  if (hexElemOwningProc[0] == p_rank) {
    stk::mesh::declare_element(mesh, *hexPart, hexElemIDs[0], hexNodeIDs[0]);
  }
  for (size_t i = 0; i < shellNodeIDs.size(); ++i) {
    if (shellElemOwningProc[i] == p_rank) {
      stk::mesh::declare_element(mesh, *shellPart, shellElemIDs[i], shellNodeIDs[i]);
    }
  }
  if (2 == p_rank) {
    mesh.declare_node(5);
    mesh.declare_node(6);
    mesh.declare_node(7);
    mesh.declare_node(8);
  }
  setup_node_sharing(mesh, shared_nodeIDs_and_procs );
  mesh.modification_end();

  ElemElemGraphTester elemElemGraph(mesh);

  if (0 == p_rank || 1 == p_rank) {
    EXPECT_EQ(1u, elemElemGraph.num_edges());
    EXPECT_EQ(1u, elemElemGraph.num_parallel_edges());
  }
  else {
    EXPECT_EQ(0u, elemElemGraph.num_edges());
    EXPECT_EQ(0u, elemElemGraph.num_parallel_edges());
  }

  mesh.modification_begin();
  if (hexElemOwningProc[1] == p_rank) {
    stk::mesh::declare_element(mesh, *hexPart, hexElemIDs[1], hexNodeIDs[1]);
  }
  mesh.modification_end();

  const Entity hex1   = mesh.get_entity(stk::topology::ELEM_RANK, 1);
  const Entity hex2   = mesh.get_entity(stk::topology::ELEM_RANK, 2);
  const Entity shell3 = mesh.get_entity(stk::topology::ELEM_RANK, 3);

  stk::mesh::EntityVector addVector;
  if (2 == p_rank) {
    addVector.push_back(hex2);
  }

  elemElemGraph.add_elements(addVector);

  if (p_rank == 0) {
    // Connectivity for Hex Element 1
    EXPECT_EQ(1u, elemElemGraph.get_num_connected_elems(hex1));
    EXPECT_EQ(5,  elemElemGraph.get_connected_remote_id_and_via_side(hex1, 0).side);
    EXPECT_EQ(3u, elemElemGraph.get_connected_remote_id_and_via_side(hex1, 0).id);
    EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));
    EXPECT_EQ(1u, elemElemGraph.num_edges());
    EXPECT_EQ(1u, elemElemGraph.num_parallel_edges());
  }
  else if (p_rank == 1) {
    // Connectivity for Shell Element 3
    EXPECT_EQ(2u, elemElemGraph.get_num_connected_elems(shell3));
    EXPECT_EQ(1,  elemElemGraph.get_connected_remote_id_and_via_side(shell3, 0).side);
    EXPECT_EQ(0,  elemElemGraph.get_connected_remote_id_and_via_side(shell3, 1).side);
    EXPECT_EQ(1u, elemElemGraph.get_connected_remote_id_and_via_side(shell3, 0).id);
    EXPECT_EQ(2u, elemElemGraph.get_connected_remote_id_and_via_side(shell3, 1).id);
    EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell3, 0));
    EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell3, 1));
    EXPECT_EQ(2u, elemElemGraph.num_edges());
    EXPECT_EQ(2u, elemElemGraph.num_parallel_edges());
  }
  else if (p_rank == 2) {
    // Connectivity for Hex Element 2
    EXPECT_EQ(1u, elemElemGraph.get_num_connected_elems(hex2));
    EXPECT_EQ(4,  elemElemGraph.get_connected_remote_id_and_via_side(hex2, 0).side);
    EXPECT_EQ(3u, elemElemGraph.get_connected_remote_id_and_via_side(hex2, 0).id);
    EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex2, 0));
    EXPECT_EQ(1u, elemElemGraph.num_edges());
    EXPECT_EQ(1u, elemElemGraph.num_parallel_edges());
  }
}

TEST( ElementGraph, AdjacentHex0Shell1Parallel )
{
  //  ID.proc
  //
  //         12.0-----------11.0
  //          /|             /|
  //         / |            / |
  //        /  |           /  |
  //      9.0-----------10.0  |
  //       |   |          |   |
  //       |   |   2.0    |4.1|
  //       |   |          |   |
  //       |  8.0---------|--7.0
  //       |  /|          |  /|
  //       | / |          | / |
  //       |/  |          |/  |
  //      5.0------------6.0  |
  //       |   |          |   |
  //       |   |   1.0    |3.1|
  //       |   |          |   |
  //       |  4.0---------|--3.0
  //       |  /           |  /
  //       | /            | /
  //       |/             |/
  //      1.0------------2.0
  //                      ^
  //                      |
  //                       ---- Two adjacent shell elements

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  unsigned p_size = stk::parallel_machine_size(pm);
  unsigned p_rank = stk::parallel_machine_rank(pm);

  if(p_size != 2u)
  {
    return;
  }

  const unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, pm);
  stk::mesh::BulkData& mesh = *bulkPtr;
  stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();

  stk::mesh::Part * hexPart   = &meta.declare_part_with_topology("hex_part",   stk::topology::HEX_8);
  stk::mesh::Part * shellPart = &meta.declare_part_with_topology("shell_part", stk::topology::SHELL_QUAD_4);
  meta.commit();

  std::vector<stk::mesh::EntityIdVector> hexNodeIDs {
    { 1, 2, 3, 4, 5,  6,  7,  8 },
    { 5, 6, 7, 8, 9, 10, 11, 12 }
  };
  stk::mesh::EntityId hexElemIDs[] = { 1, 2 };
  stk::mesh::EntityId hexElemOwningProc[] = { 0, 0 };

  std::vector<stk::mesh::EntityIdVector> shellNodeIDs {
    { 2, 3,  7,  6 },
    { 6, 7, 11, 10 }
  };
  stk::mesh::EntityId shellElemIDs[] = { 3, 4 };
  stk::mesh::EntityId shellElemOwningProc[] = { 1, 1 };

  // list of triplets: (owner-proc, shared-nodeID, sharing-proc)
  std::vector< std::vector<unsigned> > shared_nodeIDs_and_procs
  {
    { 0,  2, 1 },  // proc 0
    { 0,  3, 1 },
    { 0,  7, 1 },
    { 0,  6, 1 },
    { 0, 11, 1 },
    { 0, 10, 1 },
    { 1,  2, 0 },  // proc 1
    { 1,  3, 0 },
    { 1,  7, 0 },
    { 1,  6, 0 },
    { 1, 11, 0 },
    { 1, 10, 0 }
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

  ElemElemGraphTester elemElemGraph(mesh);

  const Entity hex1   = mesh.get_entity(stk::topology::ELEM_RANK, 1);
  const Entity hex2   = mesh.get_entity(stk::topology::ELEM_RANK, 2);
  const Entity shell3 = mesh.get_entity(stk::topology::ELEM_RANK, 3);
  const Entity shell4 = mesh.get_entity(stk::topology::ELEM_RANK, 4);

  if (p_rank == 0u) {
    // Connectivity for Hex Element 1
    EXPECT_EQ(2u,   elemElemGraph.get_num_connected_elems(hex1));
    expect_correct_connected_element_via_side(elemElemGraph, hex1, 0, hex2, 5);
    EXPECT_EQ(1,    elemElemGraph.get_connected_remote_id_and_via_side(hex1, 1).side);
    EXPECT_EQ(3u,   elemElemGraph.get_connected_remote_id_and_via_side(hex1, 1).id);
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));
    EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex1, 1));

    // Connectivity for Hex Element 2
    EXPECT_EQ(2u,   elemElemGraph.get_num_connected_elems(hex2));
    expect_correct_connected_element_via_side(elemElemGraph, hex2, 0, hex1, 4);
    EXPECT_EQ(1,    elemElemGraph.get_connected_remote_id_and_via_side(hex2, 1).side);
    EXPECT_EQ(4u,   elemElemGraph.get_connected_remote_id_and_via_side(hex2, 1).id);
    EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex2, 0));
    EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex2, 1));
    EXPECT_EQ(4u, elemElemGraph.num_edges());
    EXPECT_EQ(2u, elemElemGraph.num_parallel_edges());
  }
  else if (p_rank == 1u) {
    // Connectivity for Shell Element 3
    EXPECT_EQ(1u,   elemElemGraph.get_num_connected_elems(shell3));
    EXPECT_EQ(1,    elemElemGraph.get_connected_remote_id_and_via_side(shell3, 0).side);
    EXPECT_EQ(1u,   elemElemGraph.get_connected_remote_id_and_via_side(shell3, 0).id);
    EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell3, 0));

    // Connectivity for Shell Element 4
    EXPECT_EQ(1u,   elemElemGraph.get_num_connected_elems(shell4));
    EXPECT_EQ(1,    elemElemGraph.get_connected_remote_id_and_via_side(shell4, 0).side);
    EXPECT_EQ(2u,   elemElemGraph.get_connected_remote_id_and_via_side(shell4, 0).id);
    EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell4, 0));
    EXPECT_EQ(2u, elemElemGraph.num_edges());
    EXPECT_EQ(2u, elemElemGraph.num_parallel_edges());
  }
}


TEST( ElementGraph, Hex0Shell0Hex1Parallel )
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
  //       |   |          |   |          |   |
  //       |  2.0---------|--6.0---------|-10.1
  //       |  /           |  /           |  /
  //       | /            | /            | /
  //       |/             |/             |/
  //      1.0------------5.0------------9.1
  //                      ^
  //                      |
  //                       ---- Single shell element

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  unsigned p_size = stk::parallel_machine_size(pm);
  unsigned p_rank = stk::parallel_machine_rank(pm);

  if(p_size != 2u)
  {
    return;
  }

  const unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, pm);
  stk::mesh::BulkData& mesh = *bulkPtr;
  stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();

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
    { 5, 6, 7, 8 }
  };
  stk::mesh::EntityId shellElemIDs[] = { 3 };
  stk::mesh::EntityId shellElemOwningProc[] = { 0 };

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

  ElemElemGraphTester elemElemGraph(mesh);

  const Entity hex1   = mesh.get_entity(stk::topology::ELEM_RANK, 1);
  const Entity hex2   = mesh.get_entity(stk::topology::ELEM_RANK, 2);
  const Entity shell3 = mesh.get_entity(stk::topology::ELEM_RANK, 3);

  if (p_rank == 0) {
    {
      // Connectivity for Hex Element 1
      EXPECT_EQ(1u,     elemElemGraph.get_num_connected_elems(hex1));
      impl::ElementViaSidePair elem_via_side = elemElemGraph.get_connected_element_and_via_side(hex1, 0);
      EXPECT_EQ(5,      elem_via_side.side);
      EXPECT_EQ(shell3, elem_via_side.element);
      EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));
    }
    {
      // Connectivity for Shell Element 3
      EXPECT_EQ(2u,   elemElemGraph.get_num_connected_elems(shell3));
      impl::ElementViaSidePair elem_via_side = elemElemGraph.get_connected_element_and_via_side(shell3, 0);
      EXPECT_EQ(1,    elem_via_side.side);
      EXPECT_EQ(hex1, elem_via_side.element);
      EXPECT_EQ(0,    elemElemGraph.get_connected_remote_id_and_via_side(shell3, 1).side);
      EXPECT_EQ(2u,   elemElemGraph.get_connected_remote_id_and_via_side(shell3, 1).id);
      EXPECT_TRUE (elemElemGraph.is_connected_elem_locally_owned(shell3, 0));
      EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell3, 1));
      EXPECT_EQ(3u, elemElemGraph.num_edges());
      EXPECT_EQ(1u, elemElemGraph.num_parallel_edges());
    }
  }
  else if (p_rank == 1) {
    // Connectivity for Hex Element 2
    EXPECT_EQ(1u, elemElemGraph.get_num_connected_elems(hex2));
    EXPECT_EQ(4,  elemElemGraph.get_connected_remote_id_and_via_side(hex2, 0).side);
    EXPECT_EQ(3u, elemElemGraph.get_connected_remote_id_and_via_side(hex2, 0).id);
    EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex2, 0));
    EXPECT_EQ(1u, elemElemGraph.num_edges());
    EXPECT_EQ(1u, elemElemGraph.num_parallel_edges());
  }
}

class HexShellHexMesh {
public:
  virtual ~HexShellHexMesh(){}
  HexShellHexMesh(MPI_Comm comm, const std::vector<int>& hexOwningProcs, const std::vector<int>& shellOwningProcs,
                  const std::vector<std::vector<unsigned> >& sharedNodesAndProcs)
    : meshPtr(build_mesh(3, comm, stk::mesh::BulkData::NO_AUTO_AURA)), mesh(*meshPtr), meta(mesh.mesh_meta_data())
  {
    //  ID.proc
    //
    //          3.0------------7.0-----------11.1
    //          /|             /|             /|
    //         / |            / |            / |
    //        /  |           /  |           /  |
    //      4.0------------8.0-----------12.1  |
    //       |   |          |   |          |   |
    //       |   |   1      | 3 |   2      |   |
    //       |   |          |   |          |   |
    //       |  2.0---------|--6.0---------|-10.1
    //       |  /           |  /           |  /
    //       | /            | /            | /
    //       |/             |/             |/
    //      1.0------------5.0------------9.1

    stk::mesh::Part * hexPart   = &meta.declare_part_with_topology("hex_part",   stk::topology::HEX_8);
    stk::mesh::Part * shellPart = &meta.declare_part_with_topology("shell_part", stk::topology::SHELL_QUAD_4);
    meta.commit();

    std::vector<stk::mesh::EntityIdVector> hexNodeIDs {
      { 1, 2, 3, 4, 5,  6,  7,  8 },
      { 5, 6, 7, 8, 9, 10, 11, 12 }
    };
    stk::mesh::EntityId hexElemIDs[] = { 1, 2 };

    std::vector<stk::mesh::EntityIdVector> shellNodeIDs {
      { 5, 6, 7, 8 }
    };
    stk::mesh::EntityId shellElemIDs[] = { 3 };

    mesh.modification_begin();
    for (size_t i = 0; i < hexNodeIDs.size(); ++i) {
      if (hexOwningProcs[i] == mesh.parallel_rank()) {
        stk::mesh::declare_element(mesh, *hexPart, hexElemIDs[i], hexNodeIDs[i]);
      }
    }
    for (size_t i = 0; i < shellNodeIDs.size(); ++i) {
      if (shellOwningProcs[i] == mesh.parallel_rank()) {
        stk::mesh::declare_element(mesh, *shellPart, shellElemIDs[i], shellNodeIDs[i]);
      }
    }
    setup_node_sharing(mesh, sharedNodesAndProcs );
    mesh.modification_end();

    mesh.initialize_face_adjacent_element_graph();
  }

  stk::mesh::BulkData& get_bulk() { return mesh; }

private:
  std::shared_ptr<stk::mesh::BulkData> meshPtr;
  stk::mesh::BulkData& mesh;
  stk::mesh::MetaData& meta;
};

TEST( ElementGraph, HexDelShellHexSerial )
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  int p_size = stk::parallel_machine_size(pm);

  if(p_size > 1)
  {
    return;
  }

  std::vector<int> hexElemOwningProc = {0, 0};
  std::vector<int> shellElemOwningProc = {0};
  std::vector<std::vector<unsigned> > sharedNodeIDsAndProcs = {};

  HexShellHexMesh mesh(pm, hexElemOwningProc, shellElemOwningProc, sharedNodeIDsAndProcs);

  const Entity hex1   = mesh.get_bulk().get_entity(stk::topology::ELEM_RANK, 1);
  const Entity hex2   = mesh.get_bulk().get_entity(stk::topology::ELEM_RANK, 2);
  const Entity shell3 = mesh.get_bulk().get_entity(stk::topology::ELEM_RANK, 3);

  mesh.get_bulk().modification_begin();
  mesh.get_bulk().destroy_entity(shell3);
  mesh.get_bulk().destroy_entity(hex2);
  EXPECT_NO_THROW(mesh.get_bulk().modification_end());

  ElemElemGraph& elemElemGraph = mesh.get_bulk().get_face_adjacent_element_graph();
  // Connectivity for Hex Element 1
  EXPECT_EQ(0u,   elemElemGraph.get_num_connected_elems(hex1));

  EXPECT_EQ(0u, elemElemGraph.num_parallel_edges());
}

TEST( ElementGraph, Hex0DelShell0DelHex1Parallel )
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  unsigned p_size = stk::parallel_machine_size(pm);
  unsigned p_rank = stk::parallel_machine_rank(pm);

  if(p_size != 2u)
  {
    return;
  }

  std::vector<int> hexElemOwningProc = {0, 1};
  std::vector<int> shellElemOwningProc = {0};
  std::vector<std::vector<unsigned> > sharedNodeIDsAndProcs =
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


  HexShellHexMesh mesh(pm, hexElemOwningProc, shellElemOwningProc, sharedNodeIDsAndProcs);

  const Entity hex1   = mesh.get_bulk().get_entity(stk::topology::ELEM_RANK, 1);
  const Entity hex2   = mesh.get_bulk().get_entity(stk::topology::ELEM_RANK, 2);
  const Entity shell3 = mesh.get_bulk().get_entity(stk::topology::ELEM_RANK, 3);

  mesh.get_bulk().modification_begin();
  if (p_rank == 0) {
    mesh.get_bulk().destroy_entity(shell3);
  }
  else if (p_rank == 1) {
    mesh.get_bulk().destroy_entity(hex2);
  }
  EXPECT_NO_THROW(mesh.get_bulk().modification_end());

  ElemElemGraph& elemElemGraph = mesh.get_bulk().get_face_adjacent_element_graph();
  if (p_rank == 0) {
    // Connectivity for Hex Element 1
    EXPECT_EQ(0u, elemElemGraph.get_num_connected_elems(hex1));
    EXPECT_EQ(1u, elemElemGraph.size());
  }
  else if (p_rank == 1) {
    EXPECT_EQ(0u, elemElemGraph.size());
  }
}

TEST( ElementGraph, Hex0DelShell1Hex2Parallel )
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  unsigned p_size = stk::parallel_machine_size(pm);
  unsigned p_rank = stk::parallel_machine_rank(pm);

  if(p_size != 3u)
  {
    return;
  }

  std::vector<int> hexElemOwningProc = {0, 2};
  std::vector<int> shellElemOwningProc = {1};
  std::vector<std::vector<unsigned> > sharedNodeIDsAndProcs =
  {
    { 0, 5, 1 },  // proc 0
    { 0, 6, 1 },
    { 0, 7, 1 },
    { 0, 8, 1 },
    { 0, 5, 2 },
    { 0, 6, 2 },
    { 0, 7, 2 },
    { 0, 8, 2 },
    { 1, 5, 0 },  // proc 1
    { 1, 6, 0 },
    { 1, 7, 0 },
    { 1, 8, 0 },
    { 1, 5, 2 },
    { 1, 6, 2 },
    { 1, 7, 2 },
    { 1, 8, 2 },
    { 2, 5, 0 },  // proc 2
    { 2, 6, 0 },
    { 2, 7, 0 },
    { 2, 8, 0 },
    { 2, 5, 1 },
    { 2, 6, 1 },
    { 2, 7, 1 },
    { 2, 8, 1 }
  };


  HexShellHexMesh mesh(pm, hexElemOwningProc, shellElemOwningProc, sharedNodeIDsAndProcs);

  ElemElemGraphTester elemElemGraph(mesh.get_bulk());

  const Entity hex1   = mesh.get_bulk().get_entity(stk::topology::ELEM_RANK, 1);
  const Entity hex2   = mesh.get_bulk().get_entity(stk::topology::ELEM_RANK, 2);
  const Entity shell3 = mesh.get_bulk().get_entity(stk::topology::ELEM_RANK, 3);

  stk::mesh::impl::DeletedElementInfoVector elements_to_delete;
  if (p_rank == 1) {
    elements_to_delete.push_back({shell3, 3, mesh.get_bulk().bucket(shell3).topology().is_shell()});
  }

  elemElemGraph.delete_elements( elements_to_delete );

  mesh.get_bulk().modification_begin();
  if (p_rank == 1) {
    mesh.get_bulk().destroy_entity(shell3);
  }
  mesh.get_bulk().modification_end();

  if (p_rank == 0) {
    // Connectivity for Hex Element 1
    EXPECT_EQ(1u, elemElemGraph.num_edges());
    EXPECT_EQ(1u, elemElemGraph.num_parallel_edges());
    EXPECT_EQ(1u, elemElemGraph.get_num_connected_elems(hex1));
    EXPECT_EQ(5,  elemElemGraph.get_connected_remote_id_and_via_side(hex1, 0).side);
    EXPECT_EQ(2u, elemElemGraph.get_connected_remote_id_and_via_side(hex1, 0).id);
    EXPECT_EQ(2,  elemElemGraph.get_owning_proc_id_of_remote_element(hex1, 0));
    EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));
    EXPECT_EQ(1u, elemElemGraph.num_edges());
    EXPECT_EQ(1u, elemElemGraph.num_parallel_edges());
  }
  else if (p_rank == 1) {
    // Connectivity for Shell Element 3
    EXPECT_EQ(0u, elemElemGraph.num_edges());
    EXPECT_EQ(0u, elemElemGraph.num_parallel_edges());
    EXPECT_EQ(0u, elemElemGraph.size());
    EXPECT_EQ(0u, elemElemGraph.num_edges());
    EXPECT_EQ(0u, elemElemGraph.num_parallel_edges());
  }
  else if (p_rank == 2) {
    // Connectivity for Hex Element 2
    EXPECT_EQ(1u, elemElemGraph.num_edges());
    EXPECT_EQ(1u, elemElemGraph.num_parallel_edges());
    EXPECT_EQ(1u, elemElemGraph.get_num_connected_elems(hex2));
    EXPECT_EQ(4,  elemElemGraph.get_connected_remote_id_and_via_side(hex2, 0).side);
    EXPECT_EQ(1u, elemElemGraph.get_connected_remote_id_and_via_side(hex2, 0).id);
    EXPECT_EQ(0,  elemElemGraph.get_owning_proc_id_of_remote_element(hex2, 0));
    EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex2, 0));
    EXPECT_EQ(1u, elemElemGraph.num_edges());
    EXPECT_EQ(1u, elemElemGraph.num_parallel_edges());
  }
  elemElemGraph.write_graph(std::cerr);
}

TEST( ElementGraph, Hex0Shell1Hex2Parallel )
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  unsigned p_size = stk::parallel_machine_size(pm);
  unsigned p_rank = stk::parallel_machine_rank(pm);

  if(p_size != 3u)
  {
    return;
  }

  std::vector<int> hexElemOwningProc = {0, 2};
  std::vector<int> shellElemOwningProc = {1};
  std::vector<std::vector<unsigned> > sharedNodeIDsAndProcs =
  {
    { 0, 5, 1 },  // proc 0
    { 0, 6, 1 },
    { 0, 7, 1 },
    { 0, 8, 1 },
    { 0, 5, 2 },
    { 0, 6, 2 },
    { 0, 7, 2 },
    { 0, 8, 2 },
    { 1, 5, 0 },  // proc 1
    { 1, 6, 0 },
    { 1, 7, 0 },
    { 1, 8, 0 },
    { 1, 5, 2 },
    { 1, 6, 2 },
    { 1, 7, 2 },
    { 1, 8, 2 },
    { 2, 5, 0 },  // proc 2
    { 2, 6, 0 },
    { 2, 7, 0 },
    { 2, 8, 0 },
    { 2, 5, 1 },
    { 2, 6, 1 },
    { 2, 7, 1 },
    { 2, 8, 1 }
  };


  HexShellHexMesh mesh(pm, hexElemOwningProc, shellElemOwningProc, sharedNodeIDsAndProcs);

  ElemElemGraphTester elemElemGraph(mesh.get_bulk());

  const Entity hex1   = mesh.get_bulk().get_entity(stk::topology::ELEM_RANK, 1);
  const Entity hex2   = mesh.get_bulk().get_entity(stk::topology::ELEM_RANK, 2);
  const Entity shell3 = mesh.get_bulk().get_entity(stk::topology::ELEM_RANK, 3);

  if (p_rank == 0) {
    // Connectivity for Hex Element 1
    EXPECT_EQ(1u, elemElemGraph.get_num_connected_elems(hex1));
    EXPECT_EQ(5,  elemElemGraph.get_connected_remote_id_and_via_side(hex1, 0).side);
    EXPECT_EQ(3u, elemElemGraph.get_connected_remote_id_and_via_side(hex1, 0).id);
    EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));
    EXPECT_EQ(1u, elemElemGraph.num_edges());
    EXPECT_EQ(1u, elemElemGraph.num_parallel_edges());
  }
  else if (p_rank == 1) {
    // Connectivity for Shell Element 3
    unsigned hex1Index = 0;
    unsigned hex2Index = 1;
    EXPECT_EQ(2u, elemElemGraph.get_num_connected_elems(shell3));
    EXPECT_EQ(0,  elemElemGraph.get_connected_remote_id_and_via_side(shell3, hex2Index).side);
    EXPECT_EQ(1,  elemElemGraph.get_connected_remote_id_and_via_side(shell3, hex1Index).side);
    EXPECT_EQ(2u, elemElemGraph.get_connected_remote_id_and_via_side(shell3, hex2Index).id);
    EXPECT_EQ(1u, elemElemGraph.get_connected_remote_id_and_via_side(shell3, hex1Index).id);
    EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell3, hex2Index));
    EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell3, hex1Index));
    EXPECT_EQ(2u, elemElemGraph.num_edges());
    EXPECT_EQ(2u, elemElemGraph.num_parallel_edges());
  }
  else if (p_rank == 2) {
    // Connectivity for Hex Element 2
    EXPECT_EQ(1u, elemElemGraph.get_num_connected_elems(hex2));
    EXPECT_EQ(4,  elemElemGraph.get_connected_remote_id_and_via_side(hex2, 0).side);
    EXPECT_EQ(3u, elemElemGraph.get_connected_remote_id_and_via_side(hex2, 0).id);
    EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex2, 0));
    EXPECT_EQ(1u, elemElemGraph.num_edges());
    EXPECT_EQ(1u, elemElemGraph.num_parallel_edges());
  }
}

TEST( ElementGraph, Hex0Shell1Hex0Parallel )
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  unsigned p_size = stk::parallel_machine_size(pm);
  unsigned p_rank = stk::parallel_machine_rank(pm);

  if(p_size != 2u)
  {
    return;
  }

  std::vector<int> hexElemOwningProc = {0, 0};
  std::vector<int> shellElemOwningProc = {1};
  std::vector<std::vector<unsigned> > sharedNodeIDsAndProcs =
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


  HexShellHexMesh mesh(pm, hexElemOwningProc, shellElemOwningProc, sharedNodeIDsAndProcs);

  ElemElemGraphTester elemElemGraph(mesh.get_bulk());

  const Entity hex1   = mesh.get_bulk().get_entity(stk::topology::ELEM_RANK, 1);
  const Entity hex2   = mesh.get_bulk().get_entity(stk::topology::ELEM_RANK, 2);
  const Entity shell3 = mesh.get_bulk().get_entity(stk::topology::ELEM_RANK, 3);

  if (p_rank == 0) {
    // Connectivity for Hex Element 1
    EXPECT_EQ(1u, elemElemGraph.get_num_connected_elems(hex1));
    EXPECT_EQ(5,  elemElemGraph.get_connected_remote_id_and_via_side(hex1, 0).side);
    EXPECT_EQ(3u, elemElemGraph.get_connected_remote_id_and_via_side(hex1, 0).id);
    EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));

    // Connectivity for Hex Element 2
    EXPECT_EQ(1u, elemElemGraph.get_num_connected_elems(hex2));
    EXPECT_EQ(4,  elemElemGraph.get_connected_remote_id_and_via_side(hex2, 0).side);
    EXPECT_EQ(3u, elemElemGraph.get_connected_remote_id_and_via_side(hex2, 0).id);
    EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex2, 0));
  }
  else if (p_rank == 1) {
    // Connectivity for Shell Element 3
    unsigned hex1Index = 0;
    unsigned hex2Index = 1;
    EXPECT_EQ(2u,   elemElemGraph.get_num_connected_elems(shell3));
    EXPECT_EQ(0,    elemElemGraph.get_connected_remote_id_and_via_side(shell3, hex2Index).side);
    EXPECT_EQ(1,    elemElemGraph.get_connected_remote_id_and_via_side(shell3, hex1Index).side);
    EXPECT_EQ(2u,   elemElemGraph.get_connected_remote_id_and_via_side(shell3, hex2Index).id);
    EXPECT_EQ(1u,   elemElemGraph.get_connected_remote_id_and_via_side(shell3, hex1Index).id);
    EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell3, hex2Index));
    EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell3, hex1Index));
  }

  EXPECT_EQ(2u, elemElemGraph.num_edges());
  EXPECT_EQ(2u, elemElemGraph.num_parallel_edges());
}

TEST( ElementGraph, Hex0AddShell1Hex0Parallel )
{
  //  ID.proc
  //
  //          3.0------------7.0-----------11.0
  //          /|             /|             /|
  //         / |            / |            / |
  //        /  |           /  |           /  |
  //      4.0------------8.0-----------12.0  |
  //       |   |          |   |          |   |
  //       |   |   1.0    |3.1|   2.0    |   |
  //       |   |          |   |          |   |
  //       |  2.0---------|--6.0---------|-10.0
  //       |  /           |  /           |  /
  //       | /            | /            | /
  //       |/             |/             |/
  //      1.0------------5.0------------9.0
  //                      ^
  //                      |
  //                       ---- Single shell element

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  unsigned p_size = stk::parallel_machine_size(pm);
  unsigned p_rank = stk::parallel_machine_rank(pm);

  if(p_size != 2u)
  {
    return;
  }

  const unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, pm);
  stk::mesh::BulkData& mesh = *bulkPtr;
  stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();

  stk::mesh::Part * hexPart   = &meta.declare_part_with_topology("hex_part",   stk::topology::HEX_8);
  stk::mesh::Part * shellPart = &meta.declare_part_with_topology("shell_part", stk::topology::SHELL_QUAD_4);
  meta.commit();

  std::vector<stk::mesh::EntityIdVector> hexNodeIDs {
    { 1, 2, 3, 4, 5,  6,  7,  8 },
    { 5, 6, 7, 8, 9, 10, 11, 12 }
  };
  stk::mesh::EntityId hexElemIDs[] = { 1, 2 };
  stk::mesh::EntityId hexElemOwningProc[] = { 0, 0 };

  std::vector<stk::mesh::EntityIdVector> shellNodeIDs {
    { 5, 6, 7, 8 }
  };
  stk::mesh::EntityId shellElemIDs[] = { 3 };
  stk::mesh::EntityId shellElemOwningProc[] = { 1 };

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
  if (1 == p_rank) {
    mesh.declare_node(5);
    mesh.declare_node(6);
    mesh.declare_node(7);
    mesh.declare_node(8);
  }
  setup_node_sharing(mesh, shared_nodeIDs_and_procs );
  mesh.modification_end();

  ElemElemGraphTester elemElemGraph(mesh);

  mesh.modification_begin();
  for (size_t i = 0; i < shellNodeIDs.size(); ++i) {
    if (shellElemOwningProc[i] == p_rank) {
      stk::mesh::declare_element(mesh, *shellPart, shellElemIDs[i], shellNodeIDs[i]);
    }
  }
  mesh.modification_end();

  const Entity hex1   = mesh.get_entity(stk::topology::ELEM_RANK, 1);
  const Entity hex2   = mesh.get_entity(stk::topology::ELEM_RANK, 2);
  const Entity shell3 = mesh.get_entity(stk::topology::ELEM_RANK, 3);

  stk::mesh::EntityVector addVector;
  if (1 == p_rank) {
    addVector.push_back(shell3);
  }
  elemElemGraph.add_elements(addVector);

  if (p_rank == 0) {
    // Connectivity for Hex Element 1
    EXPECT_EQ(1u, elemElemGraph.get_num_connected_elems(hex1));
    EXPECT_EQ(5,  elemElemGraph.get_connected_remote_id_and_via_side(hex1, 0).side);
    EXPECT_EQ(3u, elemElemGraph.get_connected_remote_id_and_via_side(hex1, 0).id);
    EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));

    // Connectivity for Hex Element 2
    EXPECT_EQ(1u, elemElemGraph.get_num_connected_elems(hex2));
    EXPECT_EQ(4,  elemElemGraph.get_connected_remote_id_and_via_side(hex2, 0).side);
    EXPECT_EQ(3u, elemElemGraph.get_connected_remote_id_and_via_side(hex2, 0).id);
    EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(hex2, 0));
  }
  else if (p_rank == 1) {
    // Connectivity for Shell Element 3
    unsigned hex1Index = 0;
    unsigned hex2Index = 1;
    EXPECT_EQ(2u,   elemElemGraph.get_num_connected_elems(shell3));
    EXPECT_EQ(1,    elemElemGraph.get_connected_remote_id_and_via_side(shell3, hex1Index).side);
    EXPECT_EQ(0,    elemElemGraph.get_connected_remote_id_and_via_side(shell3, hex2Index).side);
    EXPECT_EQ(1u,   elemElemGraph.get_connected_remote_id_and_via_side(shell3, hex1Index).id);
    EXPECT_EQ(2u,   elemElemGraph.get_connected_remote_id_and_via_side(shell3, hex2Index).id);
    EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell3, hex1Index));
    EXPECT_FALSE(elemElemGraph.is_connected_elem_locally_owned(shell3, hex2Index));
  }

  EXPECT_EQ(2u, elemElemGraph.num_edges());
  EXPECT_EQ(2u, elemElemGraph.num_parallel_edges());
}

TEST( ElementGraph, Hex0DelShell1Hex0Parallel )
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  unsigned p_size = stk::parallel_machine_size(pm);
  unsigned p_rank = stk::parallel_machine_rank(pm);

  if(p_size != 2u)
  {
    return;
  }

  std::vector<int> hexElemOwningProc = { 0, 0 };
  std::vector<int> shellElemOwningProc = { 1 };

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

  HexShellHexMesh mesh(pm, hexElemOwningProc, shellElemOwningProc, shared_nodeIDs_and_procs);

  ElemElemGraphTester elemElemGraph(mesh.get_bulk());

  EXPECT_EQ(2u, elemElemGraph.num_edges());
  EXPECT_EQ(2u, elemElemGraph.num_parallel_edges());

  const Entity hex1   = mesh.get_bulk().get_entity(stk::topology::ELEM_RANK, 1);
  const Entity hex2   = mesh.get_bulk().get_entity(stk::topology::ELEM_RANK, 2);
  const Entity shell3 = mesh.get_bulk().get_entity(stk::topology::ELEM_RANK, 3);

  stk::mesh::impl::DeletedElementInfoVector elements_to_delete;
  if (p_rank == 1) {
    elements_to_delete.push_back({shell3, 3, mesh.get_bulk().bucket(shell3).topology().is_shell()});
  }

  elemElemGraph.delete_elements( elements_to_delete );

  mesh.get_bulk().modification_begin();
  if (p_rank == 1) {
    mesh.get_bulk().destroy_entity(shell3);
  }
  mesh.get_bulk().modification_end();

  if (p_rank == 0) {
    EXPECT_EQ(2u, elemElemGraph.num_edges());
    EXPECT_EQ(0u, elemElemGraph.num_parallel_edges());

    {
      // Connectivity for Hex Element 1
      EXPECT_EQ(1u,   elemElemGraph.get_num_connected_elems(hex1));
      impl::ElementViaSidePair elem_via_side = elemElemGraph.get_connected_element_and_via_side(hex1, 0);
      EXPECT_EQ(5,    elem_via_side.side);
      EXPECT_EQ(hex2, elem_via_side.element);
      EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex1, 0));
    }
    {
      // Connectivity for Hex Element 2
      EXPECT_EQ(1u, elemElemGraph.get_num_connected_elems(hex2));
      impl::ElementViaSidePair elem_via_side = elemElemGraph.get_connected_element_and_via_side(hex2, 0);
      EXPECT_EQ(4,    elem_via_side.side);
      EXPECT_EQ(hex1, elem_via_side.element);
      EXPECT_TRUE(elemElemGraph.is_connected_elem_locally_owned(hex2, 0));
    }
  }
  else if (p_rank == 1) {
    EXPECT_EQ(0u, elemElemGraph.size());
    EXPECT_EQ(0u, elemElemGraph.num_edges());
    EXPECT_EQ(0u, elemElemGraph.num_parallel_edges());
  }
}

void test_add_element_to_graph_with_element_death(stk::mesh::BulkData::AutomaticAuraOption autoAuraOption)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  int pSize = stk::parallel_machine_size(comm);
  int pRank = stk::parallel_machine_rank(comm);

  if(2 == pSize)
  {
    const unsigned spatialDim = 3;
    std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, comm);
    stk::mesh::BulkData& bulkData = *bulkPtr;
    stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();

    stk::mesh::Part& faces_part = meta.declare_part_with_topology("surface_5", stk::topology::QUAD_4);
    stk::mesh::PartVector boundary_mesh_parts {&faces_part};
    stk::io::put_io_part_attribute(faces_part);

    stk::mesh::Part& active = meta.declare_part("active"); // can't specify rank, because it gets checked against size of rank_names

    test_add_elements_to_pre_existing_graph_and_mesh(bulkData);
    stk::mesh::ElemElemGraph &graph = bulkData.get_face_adjacent_element_graph();

    stk::unit_test_util::put_mesh_into_part(bulkData, active);

    stk::mesh::EntityVector deactivated_elems;

    if (0 == pRank)
      deactivated_elems.push_back(bulkData.get_entity(stk::topology::ELEM_RANK, 1));
    else
      deactivated_elems.push_back(bulkData.get_entity(stk::topology::ELEM_RANK, 4));

    std::vector<size_t> entity_counts;
    stk::mesh::comm_mesh_counts(bulkData, entity_counts);
    ASSERT_TRUE(entity_counts[stk::topology::FACE_RANK] == 0);

    boundary_mesh_parts.push_back(&active);

    ElemGraphTestUtils::deactivate_elements(deactivated_elems, bulkData,  active);

    stk::mesh::impl::ParallelSelectedInfo remoteActiveSelector;
    stk::mesh::impl::populate_selected_value_for_remote_elements(bulkData, graph, active, remoteActiveSelector);

    stk::mesh::process_killed_elements(bulkData, deactivated_elems, active, remoteActiveSelector, boundary_mesh_parts, &boundary_mesh_parts);

    if (0 == pRank)
    {
      stk::mesh::EntityId elem1Id = 1;
      stk::mesh::EntityId elem2Id = 2;
      stk::mesh::Entity face_between_elem1_and_elem2 = ElemGraphTestUtils::get_face_between_element_ids(graph, bulkData, elem1Id, elem2Id);

      ASSERT_TRUE(bulkData.is_valid(face_between_elem1_and_elem2));
      EXPECT_TRUE(bulkData.bucket(face_between_elem1_and_elem2).member(active));

      stk::mesh::Entity elem1 = bulkData.get_entity(stk::topology::ELEM_RANK, 1);
      stk::mesh::Entity elem3 = bulkData.get_entity(stk::topology::ELEM_RANK, 3);
      EXPECT_FALSE(bulkData.bucket(elem1).member(active));
      EXPECT_TRUE(bulkData.bucket(elem3).member(active));
    }
    else
    {
      stk::mesh::EntityId elem3Id = 3;
      stk::mesh::EntityId elem4Id = 4;
      stk::mesh::Entity face_between_elem3_and_elem4 = ElemGraphTestUtils::get_face_between_element_ids(graph, bulkData, elem3Id, elem4Id);

      ASSERT_TRUE(bulkData.is_valid(face_between_elem3_and_elem4));
      EXPECT_TRUE(bulkData.bucket(face_between_elem3_and_elem4).member(active));

      stk::mesh::Entity elem2 = bulkData.get_entity(stk::topology::ELEM_RANK, 2);
      stk::mesh::Entity elem4 = bulkData.get_entity(stk::topology::ELEM_RANK, 4);
      EXPECT_TRUE(bulkData.bucket(elem2).member(active));
      EXPECT_FALSE(bulkData.bucket(elem4).member(active));
    }

    // bulkData.dump_all_mesh_info(std::cout, true);
    stk::mesh::comm_mesh_counts(bulkData, entity_counts);
    ASSERT_TRUE(entity_counts[stk::topology::FACE_RANK] == 2);

    EXPECT_EQ(4u, graph.num_edges());
    EXPECT_EQ(4u, graph.num_parallel_edges());
  }
}

TEST(ElementGraph, add_element_to_graph_with_element_death_aura_on)
{
  test_add_element_to_graph_with_element_death(stk::mesh::BulkData::AUTO_AURA);
}

TEST(ElementGraph, add_element_to_graph_with_element_death_aura_off)
{
  test_add_element_to_graph_with_element_death(stk::mesh::BulkData::NO_AUTO_AURA);
}

void test_delete_element_from_graph_with_element_death(stk::mesh::BulkData::AutomaticAuraOption autoAuraOption)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  int pSize = stk::parallel_machine_size(comm);
  int pRank = stk::parallel_machine_rank(comm);

  if(3 == pSize)
  {
    const int procRank = stk::parallel_machine_rank(comm);
    unsigned spatialDim = 3;
    std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, comm);
    stk::mesh::BulkData& bulkData = *bulkPtr;
    stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();

    stk::mesh::Part& faces_part = meta.declare_part_with_topology("surface_5", stk::topology::QUAD_4);
    stk::mesh::PartVector boundary_mesh_parts {&faces_part};
    stk::io::put_io_part_attribute(faces_part);

    stk::mesh::Part& active = meta.declare_part("active"); // can't specify rank, because it gets checked against size of rank_names

    stk::io::fill_mesh("generated:3x1x3", bulkData);

    bulkData.modification_begin();
    if (procRank == 0) {
      stk::mesh::Entity element1 = bulkData.get_entity(stk::topology::ELEM_RANK,1);
      stk::mesh::Entity element3 = bulkData.get_entity(stk::topology::ELEM_RANK,3);
      bulkData.destroy_entity(element1);
      bulkData.destroy_entity(element3);
    }
    if (procRank == 2) {
      stk::mesh::Entity element7 = bulkData.get_entity(stk::topology::ELEM_RANK,7);
      stk::mesh::Entity element9 = bulkData.get_entity(stk::topology::ELEM_RANK,9);
      bulkData.destroy_entity(element7);
      bulkData.destroy_entity(element9);
    }
    bulkData.modification_end();

    stk::mesh::ElemElemGraph &graph = bulkData.get_face_adjacent_element_graph();

    stk::unit_test_util::put_mesh_into_part(bulkData, active);

    stk::mesh::EntityVector deactivated_elems;

    if (1 == pRank)
    {
      deactivated_elems.push_back(bulkData.get_entity(stk::topology::ELEM_RANK, 5));
    }

    std::vector<size_t> entity_counts;
    stk::mesh::comm_mesh_counts(bulkData, entity_counts);
    ASSERT_TRUE(entity_counts[stk::topology::FACE_RANK] == 0);

    boundary_mesh_parts.push_back(&active);

    ElemGraphTestUtils::deactivate_elements(deactivated_elems, bulkData,  active);

    stk::mesh::impl::ParallelSelectedInfo remoteActiveSelector;
    stk::mesh::impl::populate_selected_value_for_remote_elements(bulkData, graph, active, remoteActiveSelector);

    stk::mesh::process_killed_elements(bulkData, deactivated_elems, active, remoteActiveSelector, boundary_mesh_parts, &boundary_mesh_parts);

    stk::mesh::comm_mesh_counts(bulkData, entity_counts);

    ASSERT_EQ(4u, entity_counts[stk::topology::FACE_RANK]);

    if (0 == pRank)
    {
      stk::mesh::EntityId elem2Id = 2;
      stk::mesh::EntityId elem5Id = 5;
      stk::mesh::Entity face_between_elem2_and_elem5 = ElemGraphTestUtils::get_face_between_element_ids(graph, bulkData, elem2Id, elem5Id);

      ASSERT_TRUE(bulkData.is_valid(face_between_elem2_and_elem5));
      EXPECT_TRUE(bulkData.bucket(face_between_elem2_and_elem5).member(active));

      stk::mesh::Entity elem2 = bulkData.get_entity(stk::topology::ELEM_RANK, 2);
      EXPECT_TRUE(bulkData.bucket(elem2).member(active));

      EXPECT_EQ(1u, graph.num_edges());
      EXPECT_EQ(1u, graph.num_parallel_edges());
    }
    else if(1 == pRank)
    {
      stk::mesh::EntityId elem4Id = 4;
      stk::mesh::EntityId elem5Id = 5;
      stk::mesh::EntityId elem6Id = 6;
      stk::mesh::Entity face_between_elem4_and_elem5 = ElemGraphTestUtils::get_face_between_element_ids(graph, bulkData, elem4Id, elem5Id);
      stk::mesh::Entity face_between_elem6_and_elem5 = ElemGraphTestUtils::get_face_between_element_ids(graph, bulkData, elem6Id, elem5Id);

      ASSERT_TRUE(bulkData.is_valid(face_between_elem4_and_elem5));
      EXPECT_TRUE(bulkData.bucket(face_between_elem4_and_elem5).member(active));

      ASSERT_TRUE(bulkData.is_valid(face_between_elem6_and_elem5));
      EXPECT_TRUE(bulkData.bucket(face_between_elem6_and_elem5).member(active));

      stk::mesh::Entity elem4 = bulkData.get_entity(stk::topology::ELEM_RANK, 4);
      stk::mesh::Entity elem5 = bulkData.get_entity(stk::topology::ELEM_RANK, 5);
      stk::mesh::Entity elem6 = bulkData.get_entity(stk::topology::ELEM_RANK, 6);

      EXPECT_TRUE(bulkData.bucket(elem4).member(active));
      EXPECT_TRUE(bulkData.bucket(elem6).member(active));
      EXPECT_FALSE(bulkData.bucket(elem5).member(active));

      EXPECT_EQ(6u, graph.num_edges());
      EXPECT_EQ(2u, graph.num_parallel_edges());
    }
    else if (2 == pRank)
    {
      stk::mesh::EntityId elem8Id = 8;
      stk::mesh::EntityId elem5Id = 5;
      stk::mesh::Entity face_between_elem8_and_elem5 = ElemGraphTestUtils::get_face_between_element_ids(graph, bulkData, elem8Id, elem5Id);

      ASSERT_TRUE(bulkData.is_valid(face_between_elem8_and_elem5));
      EXPECT_TRUE(bulkData.bucket(face_between_elem8_and_elem5).member(active));

      stk::mesh::Entity elem8 = bulkData.get_entity(stk::topology::ELEM_RANK, 8);
      EXPECT_TRUE(bulkData.bucket(elem8).member(active));

      EXPECT_EQ(1u, graph.num_edges());
      EXPECT_EQ(1u, graph.num_parallel_edges());
    }
  }
}


TEST(ElementGraph, delete_element_from_graph_with_element_death_aura_on)
{
  test_delete_element_from_graph_with_element_death(stk::mesh::BulkData::AUTO_AURA);
}

TEST(ElementGraph, delete_element_from_graph_with_element_death_aura_off)
{
  test_delete_element_from_graph_with_element_death(stk::mesh::BulkData::NO_AUTO_AURA);
}

ElemElemGraphTester create_base_1x1x4_elem_graph(stk::ParallelMachine &comm, stk::mesh::BulkData &bulkData)
{
  unsigned nProc  = stk::parallel_machine_size(comm);
  unsigned myProc = stk::parallel_machine_rank(comm);

  STK_ThrowRequire(nProc <= 4);

  stk::io::fill_mesh("generated:1x1x4", bulkData);

  std::vector<size_t> counts;
  stk::mesh::count_entities(bulkData.mesh_meta_data().locally_owned_part(), bulkData, counts);
  unsigned numLocallyOwnedElems = counts[stk::topology::ELEM_RANK];

  if(myProc != 0)
  {
    EXPECT_EQ(4u/nProc, numLocallyOwnedElems);
  }
  else
  {
    EXPECT_EQ((4u/nProc + 4u%nProc), numLocallyOwnedElems);
  }

  ElemElemGraphTester elem_graph(bulkData);

  stk::mesh::Entity elem1 = bulkData.get_entity(stk::topology::ELEM_RANK, 1);
  stk::mesh::Entity elem2 = bulkData.get_entity(stk::topology::ELEM_RANK, 2);
  stk::mesh::Entity elem3 = bulkData.get_entity(stk::topology::ELEM_RANK, 3);
  stk::mesh::Entity elem4 = bulkData.get_entity(stk::topology::ELEM_RANK, 4);

  // We know that element 1 has 1 connection
  //                      2 has 2 connections
  //                      3 has 2 connections
  //                      4 has 1 connection

  bool ownedElem1 = false;
  bool ownedElem2 = false;
  bool ownedElem3 = false;
  bool ownedElem4 = false;

  unsigned numEdges = 0;

  if(bulkData.is_valid(elem1) && bulkData.bucket(elem1).owned())
  {
    numEdges += 1;
    ownedElem1 = true;
  }
  if(bulkData.is_valid(elem2) && bulkData.bucket(elem2).owned())
  {
    numEdges += 2;
    ownedElem2 = true;
  }

  if(bulkData.is_valid(elem3) && bulkData.bucket(elem3).owned())
  {
    numEdges += 2;
    ownedElem3 = true;
  }
  if(bulkData.is_valid(elem4) && bulkData.bucket(elem4).owned())
  {
    numEdges += 1;
    ownedElem4 = true;
  }

  unsigned numParallelEdges = 0;
  for(unsigned elem_id=1; elem_id<=4; ++elem_id)
  {
    stk::mesh::Entity elem = bulkData.get_entity(stk::topology::ELEM_RANK, elem_id);

    if (!bulkData.is_valid(elem) || !bulkData.bucket(elem).owned())
    {
      continue;
    }

    stk::mesh::EntityId leftNeighbor = elem_id - 1;
    stk::mesh::Entity leftElem = bulkData.get_entity(stk::topology::ELEM_RANK, leftNeighbor);
    bool ownedLeftNeighbor = bulkData.is_valid(leftElem) && bulkData.bucket(leftElem).owned();
    if(!ownedLeftNeighbor && (elem_id > 1))
    {
      numParallelEdges++;
    }

    stk::mesh::EntityId rightNeighbor = elem_id + 1;
    stk::mesh::Entity rightElem = bulkData.get_entity(stk::topology::ELEM_RANK, rightNeighbor);
    bool ownedRightNeighbor = bulkData.is_valid(rightElem) && bulkData.bucket(rightElem).owned();
    if(!ownedRightNeighbor && (elem_id < 4))
    {
      numParallelEdges++;
    }
  }

  EXPECT_EQ(numLocallyOwnedElems, elem_graph.size());
  EXPECT_EQ(numEdges, elem_graph.num_edges());
  EXPECT_EQ(numParallelEdges, elem_graph.num_parallel_edges());

  if (ownedElem1)
  {
    EXPECT_EQ(5, elem_graph.get_side_from_element1_to_element2(1, 2));
    EXPECT_EQ(-1, elem_graph.get_side_from_element1_to_element2(1, 3));
    EXPECT_EQ(-1, elem_graph.get_side_from_element1_to_element2(1, 4));

  }
  if (ownedElem2)
  {
    EXPECT_EQ(4, elem_graph.get_side_from_element1_to_element2(2, 1));
    EXPECT_EQ(5, elem_graph.get_side_from_element1_to_element2(2, 3));
    EXPECT_EQ(-1, elem_graph.get_side_from_element1_to_element2(2, 4));
  }
  if (ownedElem3)
  {
    EXPECT_EQ(4, elem_graph.get_side_from_element1_to_element2(3, 2));
    EXPECT_EQ(5, elem_graph.get_side_from_element1_to_element2(3, 4));
    EXPECT_EQ(-1, elem_graph.get_side_from_element1_to_element2(3, 1));
  }
  if (ownedElem4)
  {
    EXPECT_EQ(4, elem_graph.get_side_from_element1_to_element2(4, 3));
    EXPECT_EQ(-1, elem_graph.get_side_from_element1_to_element2(4, 1));
    EXPECT_EQ(-1, elem_graph.get_side_from_element1_to_element2(4, 2));
  }

  return elem_graph;
}

// element ids / proc_id:
// |-------|-------|-------|
// |       |       |       |
// |  1/0  |  4/2  |  7/2  |
// |       |       |       |
// |-------|-------|-------|
// |       |       |       |
// |  2/0  |  5/1  |  8/2  |
// |       |       |       |
// |-------|-------|-------|
// |       |       |       |
// |  3/0  |  6/2  |  9/2  |
// |       |       |       |
// |-------|-------|-------|


TEST(ElementGraph, TestKeyHoleSimilarProblemAInParallel)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;

  if(stk::parallel_machine_size(comm) == 3)
  {
    const int procRank = stk::parallel_machine_rank(comm);
    unsigned spatialDim = 3;

    std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, comm, stk::mesh::BulkData::NO_AUTO_AURA);
    stk::mesh::BulkData& bulkData = *bulkPtr;
    stk::io::fill_mesh("generated:3x1x3", bulkData);

    stk::mesh::EntityProcVec elementProcChanges;
    if (procRank == 1) {
      elementProcChanges.push_back(stk::mesh::EntityProc(bulkData.get_entity(stk::topology::ELEM_RANK,4),2));
      elementProcChanges.push_back(stk::mesh::EntityProc(bulkData.get_entity(stk::topology::ELEM_RANK,6),2));
    }
    bulkData.change_entity_owner(elementProcChanges);

    ElemElemGraphTester graph(bulkData);
    if (procRank == 0) {
      stk::mesh::Entity local_element = bulkData.get_entity(stk::topology::ELEM_RANK,2);
      ASSERT_TRUE(bulkData.bucket(local_element).owned());
      ASSERT_EQ(3u, graph.get_num_connected_elems(local_element));

      EXPECT_EQ( 3, graph.get_connected_element_and_via_side(local_element,0).side);
      EXPECT_EQ( 1, graph.get_connected_element_and_via_side(local_element,1).side);
      EXPECT_EQ( 5, graph.get_connected_remote_id_and_via_side(local_element,2).side);

      EXPECT_TRUE(graph.is_connected_elem_locally_owned(local_element, 0));
      EXPECT_TRUE(graph.is_connected_elem_locally_owned(local_element, 1));
      EXPECT_FALSE(graph.is_connected_elem_locally_owned(local_element, 2));

      EXPECT_EQ( 1u, bulkData.identifier(graph.get_connected_element_and_via_side(local_element, 0).element));
      EXPECT_EQ( 3u, bulkData.identifier(graph.get_connected_element_and_via_side(local_element, 1).element));
      EXPECT_EQ( 5u, graph.get_connected_remote_id_and_via_side(local_element, 2).id);

      EXPECT_EQ( 1, graph.get_owning_proc_id_of_remote_element(local_element, 2));

      EXPECT_EQ(7u, graph.num_edges());
      EXPECT_EQ(3u, graph.num_parallel_edges());
    }
    if (procRank == 1) {
      stk::mesh::Entity local_element = bulkData.get_entity(stk::topology::ELEM_RANK,5);
      ASSERT_TRUE(bulkData.bucket(local_element).owned());
      size_t numConnectedElems = graph.get_num_connected_elems(local_element);
      ASSERT_EQ(4u, numConnectedElems);

      EXPECT_EQ( 4, graph.get_connected_remote_id_and_via_side(local_element,0).side);
      EXPECT_EQ( 3, graph.get_connected_remote_id_and_via_side(local_element,1).side);
      EXPECT_EQ( 1, graph.get_connected_remote_id_and_via_side(local_element,2).side);
      EXPECT_EQ( 5, graph.get_connected_remote_id_and_via_side(local_element,3).side);

      EXPECT_FALSE(graph.is_connected_elem_locally_owned(local_element, 0));
      EXPECT_FALSE(graph.is_connected_elem_locally_owned(local_element, 1));
      EXPECT_FALSE(graph.is_connected_elem_locally_owned(local_element, 2));
      EXPECT_FALSE(graph.is_connected_elem_locally_owned(local_element, 3));

      EXPECT_EQ( 2u, graph.get_connected_remote_id_and_via_side(local_element, 0).id);
      EXPECT_EQ( 4u, graph.get_connected_remote_id_and_via_side(local_element, 1).id);
      EXPECT_EQ( 6u, graph.get_connected_remote_id_and_via_side(local_element, 2).id);
      EXPECT_EQ( 8u, graph.get_connected_remote_id_and_via_side(local_element, 3).id);

      EXPECT_EQ( 0, graph.get_owning_proc_id_of_remote_element(local_element, 0));
      EXPECT_EQ( 2, graph.get_owning_proc_id_of_remote_element(local_element, 1));
      EXPECT_EQ( 2, graph.get_owning_proc_id_of_remote_element(local_element, 2));
      EXPECT_EQ( 2, graph.get_owning_proc_id_of_remote_element(local_element, 3));

      EXPECT_EQ(4u, graph.num_edges());
      EXPECT_EQ(4u, graph.num_parallel_edges());
    }
    if (procRank == 2) {
      EXPECT_EQ(13u, graph.num_edges());
      EXPECT_EQ( 5u, graph.num_parallel_edges());
    }
  }
}

TEST(TestAuraHex, TestKeyHoleSimilarProblemBInParallel)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;

  if(stk::parallel_machine_size(comm) != 3) { GTEST_SKIP(); }

  const int procRank = stk::parallel_machine_rank(comm);
  unsigned spatialDim = 3;

  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, comm);
  stk::mesh::BulkData& bulkData = *bulkPtr;
  stk::io::fill_mesh("generated:3x1x3", bulkData);

  {
    stk::mesh::EntityProcVec elementProcChanges;
    if (procRank == 1) {
      elementProcChanges.push_back(stk::mesh::EntityProc(bulkData.get_entity(stk::topology::ELEM_RANK,4),2));
      elementProcChanges.push_back(stk::mesh::EntityProc(bulkData.get_entity(stk::topology::ELEM_RANK,6),2));
    }
    bulkData.change_entity_owner(elementProcChanges);
  }

{
stk::parallel_machine_barrier(bulkData.parallel());
std::ostringstream os;
os<<"P"<<bulkData.parallel_rank()<<" ****************** about to do test mod"<<std::endl;
stk::mesh::EntityVector elems;
stk::mesh::get_entities(bulkData, stk::topology::ELEM_RANK, bulkData.mesh_meta_data().locally_owned_part(), elems);
for(stk::mesh::Entity elem : elems) {
os<<bulkData.entity_key(elem)<<" nodes: ";
const stk::mesh::Entity* nodes = bulkData.begin_nodes(elem);
for(int i=0; i<8; ++i) os<<bulkData.identifier(nodes[i])<<" o="<<bulkData.parallel_owner_rank(nodes[i])<<",s="<<bulkData.bucket(nodes[i]).shared()<<",a="<<bulkData.bucket(nodes[i]).in_aura()<<"; ";
os<<std::endl;
}
std::cerr<<os.str();
stk::parallel_machine_barrier(bulkData.parallel());
}
  bulkData.modification_begin();
  if (procRank == 1) {
    stk::mesh::Entity elem5 = bulkData.get_entity(stk::topology::ELEM_RANK,5);
    ASSERT_TRUE(bulkData.parallel_owner_rank(elem5)==1);
{
std::ostringstream os;
os<<"P"<<bulkData.parallel_rank()<<" elem5 nodes: ";
const stk::mesh::Entity* nodes = bulkData.begin_nodes(elem5);
for(int i=0; i<8; ++i) {
os<<bulkData.identifier(nodes[i])<<" o="<<bulkData.parallel_owner_rank(nodes[i])<<",s="<<bulkData.bucket(nodes[i]).shared()<<",a="<<bulkData.bucket(nodes[i]).in_aura()<<"; ";
}
os<<std::endl;
std::cerr<<os.str();
}
    EXPECT_TRUE(bulkData.destroy_entity(elem5));
  }
  bulkData.modification_end();
}

// element ids / proc_id:
// |-------|-------|-------|
// |       |       |       |
// |  1/0  |  4/2  |  7/2  |
// |       |       |       |
// |-------|-------|-------|
// |       |       |       |
// |  2/0  |  n/a  |  8/2  |
// |       |       |       |
// |-------|-------|-------|
// |       |       |       |
// |  3/0  |  6/2  |  9/2  |
// |       |       |       |
// |-------|-------|-------|
// The element in the middle has been deleted

TEST(ElementGraph, TestKeyHoleSimilarProblemBInParallel)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;

  if(stk::parallel_machine_size(comm) != 3) { GTEST_SKIP(); }

  const int procRank = stk::parallel_machine_rank(comm);
  unsigned spatialDim = 3;

  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, comm);
  stk::mesh::BulkData& bulkData = *bulkPtr;
  stk::io::fill_mesh("generated:3x1x3", bulkData);

  stk::mesh::EntityProcVec elementProcChanges;
  if (procRank == 1) {
    elementProcChanges.push_back(stk::mesh::EntityProc(bulkData.get_entity(stk::topology::ELEM_RANK,4),2));
    elementProcChanges.push_back(stk::mesh::EntityProc(bulkData.get_entity(stk::topology::ELEM_RANK,6),2));
  }

  bulkData.change_entity_owner(elementProcChanges);
{
stk::parallel_machine_barrier(bulkData.parallel());
std::ostringstream os;
os<<"P"<<bulkData.parallel_rank()<<" about to do test mod"<<std::endl;
std::cerr<<os.str();
stk::parallel_machine_barrier(bulkData.parallel());
}
  bulkData.modification_begin();
  if (procRank == 1) {
    stk::mesh::Entity local_element5 = bulkData.get_entity(stk::topology::ELEM_RANK,5);
    bulkData.destroy_entity(local_element5);
  }
  bulkData.modification_end();

  ElemElemGraphTester graph(bulkData);
  if (procRank == 0) {
    EXPECT_EQ(2u, graph.get_num_connected_elems(bulkData.get_entity(stk::topology::ELEM_RANK,1)));
    EXPECT_EQ(2u, graph.get_num_connected_elems(bulkData.get_entity(stk::topology::ELEM_RANK,2)));
    EXPECT_EQ(2u, graph.get_num_connected_elems(bulkData.get_entity(stk::topology::ELEM_RANK,3)));
    EXPECT_EQ(6u, graph.num_edges());
    EXPECT_EQ(2u, graph.num_parallel_edges());
  }
  if (procRank == 1) {
    EXPECT_EQ(0u, graph.size());
    EXPECT_EQ(0u, graph.num_edges());
    EXPECT_EQ(0u, graph.num_parallel_edges());
  }
  if (procRank == 2) {
    EXPECT_EQ(2u, graph.get_num_connected_elems(bulkData.get_entity(stk::topology::ELEM_RANK,4)));
    EXPECT_EQ(2u, graph.get_num_connected_elems(bulkData.get_entity(stk::topology::ELEM_RANK,6)));
    EXPECT_EQ(2u, graph.get_num_connected_elems(bulkData.get_entity(stk::topology::ELEM_RANK,7)));
    EXPECT_EQ(2u, graph.get_num_connected_elems(bulkData.get_entity(stk::topology::ELEM_RANK,8)));
    EXPECT_EQ(2u, graph.get_num_connected_elems(bulkData.get_entity(stk::topology::ELEM_RANK,9)));
    EXPECT_EQ(10u, graph.num_edges());
    EXPECT_EQ( 2u, graph.num_parallel_edges());
  }
}

void test_parallel_uniqueness(const std::vector<stk::mesh::EntityId> &ids_in_use, const std::vector<stk::mesh::EntityId>& requested_ids, stk::ParallelMachine comm)
{
  std::vector<stk::mesh::EntityId> global_ids_in_use;
  stk::parallel_vector_concat(comm, ids_in_use, global_ids_in_use);
  std::sort(global_ids_in_use.begin(), global_ids_in_use.end());

  std::vector<stk::mesh::EntityId> global_requested_ids;
  stk::parallel_vector_concat(comm, requested_ids, global_requested_ids);
  std::sort(global_requested_ids.begin(), global_requested_ids.end());

  std::vector<stk::mesh::EntityId> intersection;

  std::set_intersection(global_ids_in_use.begin(), global_ids_in_use.end(),
                        global_requested_ids.begin(), global_requested_ids.end(),
                        std::back_inserter(intersection));

  EXPECT_TRUE(intersection.empty());
}

TEST(ElemGraph, test_id_reservation)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  if(stk::parallel_machine_size(comm)==2)
  {
    unsigned spatialDim = 3;

    stk::mesh::MetaData meta(spatialDim);
    BulkDataElementGraphTester bulkData(meta, comm);
    stk::io::fill_mesh("generated:1x1x4", bulkData);

    std::vector<stk::mesh::EntityId> ids_in_use = bulkData.my_internal_get_ids_in_use_this_proc_for_locally_owned(stk::topology::ELEM_RANK);

    size_t num_ids_requested_per_proc = 10;
    std::vector<stk::mesh::EntityId> requested_ids;
    bulkData.generate_new_ids(stk::topology::ELEM_RANK, num_ids_requested_per_proc, requested_ids);

    test_parallel_uniqueness(ids_in_use, requested_ids, comm);

    std::vector<stk::mesh::EntityId> requested_ids_again;
    bulkData.generate_new_ids_given_reserved_ids(stk::topology::ELEM_RANK, num_ids_requested_per_proc, requested_ids, requested_ids_again);

    test_parallel_uniqueness(ids_in_use, requested_ids_again, comm);
    test_parallel_uniqueness(requested_ids, requested_ids_again, comm);
  }
}

bool is_valid_graph_element(const impl::ElementGraph &elem_graph, stk::mesh::impl::LocalId elem_id)
{
  stk::mesh::impl::LocalId max_elem_id = static_cast<stk::mesh::impl::LocalId>(elem_graph.size());
  return (elem_id >= 0 && elem_id < max_elem_id);
}

int check_connectivity(const impl::ElementGraph &elem_graph, const impl::SidesForElementGraph &via_sides,
                       stk::mesh::impl::LocalId element_id1, stk::mesh::impl::LocalId element_id2)
{
  int side=-1;
  if (is_valid_graph_element(elem_graph, element_id1) && is_valid_graph_element(elem_graph, element_id2)) {
    side = get_side_from_element1_to_element2(elem_graph, via_sides, element_id1, element_id2);
  }
  return side;
}

int get_side_from_element1_to_element2(const impl::ElementGraph &elem_graph,
                                       const impl::SidesForElementGraph &via_sides,
                                       stk::mesh::impl::LocalId element1_local_id,
                                       stk::mesh::impl::LocalId other_element_id)
{
  int side = -1;
  const std::vector<stk::mesh::impl::LocalId>& conn_elements = elem_graph[element1_local_id];

  std::vector<stk::mesh::impl::LocalId>::const_iterator iter = std::find(conn_elements.begin(), conn_elements.end(), other_element_id);
  if ( iter != conn_elements.end() )
  {
    int64_t index = iter - conn_elements.begin();
    side = via_sides[element1_local_id][index];
  }
  return side;
}
//EndDocExample1

void add_elem3_on_proc_1(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part& block1 = *bulk.mesh_meta_data().get_part("block_1");
  stk::mesh::PartVector parts = {&block1};

  bulk.modification_begin();

  if (bulk.parallel_rank() == 1)
  {
    stk::mesh::EntityId elemId = 3;

    stk::mesh::EntityIdVector nodeIds = {6, 13, 14, 8, 10, 15, 16, 12};
    stk::mesh::declare_element(bulk, parts, elemId, nodeIds);
  }

  bulk.modification_end();
}

TEST(ElemGraph, get_all_sides_sideset_including_ghosts)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  if(stk::parallel_machine_size(comm)==2)
  {
    {
      std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(3, comm, stk::mesh::BulkData::AUTO_AURA);
      stk::mesh::BulkData& bulk = *bulkPtr;
      stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();
      stk::io::fill_mesh("generated:1x1x2", bulk);
      add_elem3_on_proc_1(bulk);

      bool includeAuraElementSides = true;
      std::vector<stk::mesh::SideSetEntry> sides = stk::mesh::SkinMeshUtil::get_all_sides_sideset(bulk, meta.universal_part(), includeAuraElementSides);
      EXPECT_EQ(16u, sides.size());
    }
    {
      std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(3, comm, stk::mesh::BulkData::AUTO_AURA);
      stk::mesh::BulkData& bulk = *bulkPtr;
      stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();
      stk::io::fill_mesh("generated:1x2x2", bulk);
      bool includeAuraElementSides = true;
      std::vector<stk::mesh::SideSetEntry> sides = stk::mesh::SkinMeshUtil::get_all_sides_sideset(bulk, meta.universal_part(), includeAuraElementSides);
      EXPECT_EQ(20u, sides.size());
    }
    {
      std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(3, comm, stk::mesh::BulkData::AUTO_AURA);
      stk::mesh::BulkData& bulk = *bulkPtr;
      stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();
      stk::io::fill_mesh("generated:2x2x2", bulk);
      bool includeAuraElementSides = true;
      std::vector<stk::mesh::SideSetEntry> sides = stk::mesh::SkinMeshUtil::get_all_sides_sideset(bulk, meta.universal_part(), includeAuraElementSides);
      EXPECT_EQ(36u, sides.size());
    }
  }
}

TEST(ElemGraph, test_initial_graph_creation_with_deactivated_elements)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  if(stk::parallel_machine_size(comm)==2)
  {
    stk::mesh::MetaData meta(3);
    stk::mesh::Part &activePart = meta.declare_part("active");
    BulkDataElementGraphTester bulkData(meta, comm);
    stk::io::fill_mesh("generated:1x1x4", bulkData);
    stk::unit_test_util::put_mesh_into_part(bulkData, activePart);

    bulkData.modification_begin();
    if(bulkData.parallel_rank() == 0)
    {
      stk::mesh::Entity elem2 = bulkData.get_entity(stk::topology::ELEMENT_RANK, 2);
      bulkData.change_entity_parts(elem2, ConstPartVector{}, ConstPartVector{&activePart});
    }
    bulkData.modification_end();

    stk::mesh::ElemElemGraph graph(bulkData);
    stk::mesh::impl::ParallelSelectedInfo remoteSkinSelector;
    stk::mesh::impl::populate_selected_value_for_remote_elements(bulkData, graph, activePart, remoteSkinSelector);
    if (bulkData.parallel_rank() == 0) {
      stk::mesh::Entity elem1 = bulkData.get_entity(stk::topology::ELEM_RANK,1);
      stk::mesh::Entity elem2 = bulkData.get_entity(stk::topology::ELEM_RANK,2);


      ASSERT_EQ(1u, graph.get_num_connected_elems(elem1));

      EXPECT_TRUE(graph.is_connected_elem_locally_owned(elem1, 0));
      EXPECT_EQ(elem2, graph.get_connected_element_and_via_side(elem1, 0).element);


      ASSERT_EQ(2u, graph.get_num_connected_elems(elem2));

      EXPECT_TRUE(graph.is_connected_elem_locally_owned(elem2, 0));
      EXPECT_EQ(elem1, graph.get_connected_element_and_via_side(elem2, 0).element);

      EXPECT_TRUE(!graph.is_connected_elem_locally_owned(elem2, 1));
      stk::mesh::EntityId remoteElemId = 3;
      EXPECT_EQ(remoteElemId, graph.get_connected_remote_id_and_via_side(elem2, 1).id);
      stk::mesh::impl::ParallelInfo &parallelInfo = graph.get_parallel_edge_info(elem2, 5, remoteElemId, 4);
      EXPECT_EQ(1, parallelInfo.get_proc_rank_of_neighbor());
      EXPECT_TRUE(remoteSkinSelector[-remoteElemId]);
    }
    else if (bulkData.parallel_rank() == 1) {
      stk::mesh::Entity elem3 = bulkData.get_entity(stk::topology::ELEM_RANK,3);
      stk::mesh::Entity elem4 = bulkData.get_entity(stk::topology::ELEM_RANK,4);


      ASSERT_EQ(2u, graph.get_num_connected_elems(elem3));

      EXPECT_TRUE(graph.is_connected_elem_locally_owned(elem3, 0));
      EXPECT_EQ(elem4, graph.get_connected_element_and_via_side(elem3, 0).element);

      EXPECT_TRUE(!graph.is_connected_elem_locally_owned(elem3, 1));
      stk::mesh::EntityId remoteElemId = 2;
      EXPECT_EQ(remoteElemId, graph.get_connected_remote_id_and_via_side(elem3, 1).id);
      stk::mesh::impl::ParallelInfo &parallelInfo = graph.get_parallel_edge_info(elem3, 4, remoteElemId, 5);
      EXPECT_EQ(0, parallelInfo.get_proc_rank_of_neighbor());
      EXPECT_TRUE(!remoteSkinSelector[-remoteElemId]);


      ASSERT_EQ(1u, graph.get_num_connected_elems(elem4));

      EXPECT_TRUE(graph.is_connected_elem_locally_owned(elem4, 0));
      EXPECT_EQ(elem3, graph.get_connected_element_and_via_side(elem4, 0).element);
    }
  }
}

} // namespace
