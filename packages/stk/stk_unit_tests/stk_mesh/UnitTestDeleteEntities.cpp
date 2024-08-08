#include <gtest/gtest.h>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include "stk_mesh/base/FEMHelpers.hpp"
#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_mesh/base/MeshUtils.hpp"
#include "stk_unit_test_utils/BulkDataTester.hpp"

namespace
{

void expect_num_elements_and_nodes(stk::mesh::BulkData &bulk, size_t goldNumElems, size_t goldNumNodes)
{
  std::vector<size_t> entityCounts;
  stk::mesh::comm_mesh_counts(bulk, entityCounts);
  EXPECT_EQ(goldNumElems, entityCounts[stk::topology::ELEM_RANK]);
  EXPECT_EQ(goldNumNodes, entityCounts[stk::topology::NODE_RANK]);
}

void expect_num_elements_and_faces_and_nodes(stk::mesh::BulkData &bulk, size_t goldNumElems, size_t goldNumFaces, size_t goldNumNodes)
{
  std::vector<size_t> entityCounts;
  stk::mesh::comm_mesh_counts(bulk, entityCounts);
  EXPECT_EQ(goldNumElems, entityCounts[stk::topology::ELEM_RANK]);
  EXPECT_EQ(goldNumFaces, entityCounts[stk::topology::FACE_RANK]);
  EXPECT_EQ(goldNumNodes, entityCounts[stk::topology::NODE_RANK]);
}

class HexShellHexMesh : public stk::unit_test_util::MeshFixture
{
protected:
  HexShellHexMesh()
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
    std::string meshDesc;
    if (get_parallel_size() == 1) {
      meshDesc =
          "0,1,HEX_8,1,2,3,4,5,6,7,8\n\
          0,2,HEX_8,5,6,7,8,9,10,11,12\n\
          0,3,SHELL_QUAD_4,5,6,7,8";
    }
    else {
      meshDesc =
          "0,1,HEX_8,1,2,3,4,5,6,7,8\n\
          1,2,HEX_8,5,6,7,8,9,10,11,12\n\
          0,3,SHELL_QUAD_4,5,6,7,8";
    }
    stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);
  }
};

TEST_F(HexShellHexMesh, DeleteShell_OnlyHexesRemain)
{
  expect_num_elements_and_nodes(get_bulk(), 3u, 12u);
  get_bulk().destroy_elements_of_topology(stk::topology::SHELL_QUAD_4);
  expect_num_elements_and_nodes(get_bulk(), 2u, 12u);
}

class HexHexShellMesh : public stk::unit_test_util::MeshFixture
{
protected:
  HexHexShellMesh()
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
    std::string meshDesc;
    if (get_parallel_size() == 1) {
      meshDesc =
          "0,1,HEX_8,1,2,3,4,5,6,7,8\n\
          0,2,HEX_8,5,6,7,8,9,10,11,12\n\
          0,3,SHELL_QUAD_4,9,10,11,12";
    }
    else {
      meshDesc =
          "0,1,HEX_8,1,2,3,4,5,6,7,8\n\
          1,2,HEX_8,5,6,7,8,9,10,11,12\n\
          0,3,SHELL_QUAD_4,9,10,11,12";
    }
    stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);
  }
};

TEST_F(HexHexShellMesh, DeleteAllHexes_OnlyShellRemains)
{
  expect_num_elements_and_nodes(get_bulk(), 3u, 12u);
  get_bulk().destroy_elements_of_topology(stk::topology::HEX_8);
  expect_num_elements_and_nodes(get_bulk(), 1u, 4u);

  stk::mesh::EntityVector nodes;
  stk::mesh::get_entities(get_bulk(), stk::topology::NODE_RANK, nodes);
  for(stk::mesh::Entity node : nodes)
    EXPECT_TRUE(!get_bulk().bucket(node).shared());
}

class HexWedgeHexMesh : public stk::unit_test_util::MeshFixture
{
protected:
  HexWedgeHexMesh()
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
    std::string meshDesc;
    if (get_parallel_size() == 1) {
      meshDesc =
          "0,1,HEX_8,1,2,3,4,5,6,7,8\n\
          0,2,WEDGE_6,5,9,8,6,10,7\n\
          0,3,HEX_8,11,12,13,14,5,9,10,6";
    }
    else {
      meshDesc =
          "0,1,HEX_8,1,2,3,4,5,6,7,8\n\
          1,2,WEDGE_6,5,9,8,6,10,7\n\
          1,3,HEX_8,11,12,13,14,5,9,10,6";
    }
    stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);
  }
};

TEST_F(HexWedgeHexMesh, CreateFacesThenDeleteAllHexes_OnlyWedgeRemains)
{
  expect_num_elements_and_faces_and_nodes(get_bulk(), 3u, 0u, 14u);
  stk::mesh::create_all_sides(get_bulk(), get_meta().universal_part(), {}, false);
  expect_num_elements_and_faces_and_nodes(get_bulk(), 3u, 15u, 14u);
  get_bulk().destroy_elements_of_topology(stk::topology::HEX_8);
  expect_num_elements_and_faces_and_nodes(get_bulk(), 1u, 5u, 6u);
}

TEST_F(HexWedgeHexMesh, CreateFacesThenDeleteWedge_TwoHexesRemain)
{
  expect_num_elements_and_faces_and_nodes(get_bulk(), 3u, 0u, 14u);
  stk::mesh::create_all_sides(get_bulk(), get_meta().universal_part(), {}, false);
  expect_num_elements_and_faces_and_nodes(get_bulk(), 3u, 15u, 14u);
  get_bulk().destroy_elements_of_topology(stk::topology::WEDGE_6);
  expect_num_elements_and_faces_and_nodes(get_bulk(), 2u, 12u, 14u);
}

void add_hexes_back(stk::mesh::BulkData &bulk)
{
  bulk.modification_begin();
  stk::mesh::PartVector topologyParts = {&bulk.mesh_meta_data().get_topology_root_part(stk::topology::HEX_8)};
  if(bulk.parallel_rank() == 0 || bulk.parallel_size() == 1)
    stk::mesh::declare_element(bulk, topologyParts, 1, {1,2,3,4,5,6,7,8});
  if(bulk.parallel_rank() == 1 || bulk.parallel_size() == 1)
    stk::mesh::declare_element(bulk, topologyParts, 3, {11,12,13,14,5,9,10,6});
  if(bulk.parallel_size() > 1 && bulk.parallel_rank() < 2)
  {
    int otherProcRank = bulk.parallel_rank() == 0 ? 1 : 0;
    bulk.add_node_sharing(bulk.get_entity(stk::topology::NODE_RANK, 5), otherProcRank);
    bulk.add_node_sharing(bulk.get_entity(stk::topology::NODE_RANK, 6), otherProcRank);
    bulk.add_node_sharing(bulk.get_entity(stk::topology::NODE_RANK, 7), otherProcRank);
    bulk.add_node_sharing(bulk.get_entity(stk::topology::NODE_RANK, 8), otherProcRank);
  }
  bulk.modification_end();
}

TEST_F(HexWedgeHexMesh, CreateFacesThenDeleteAllHexesThenCreateFaces_OnlyWedgeRemains)
{
  expect_num_elements_and_faces_and_nodes(get_bulk(), 3u, 0u, 14u);

  stk::mesh::create_all_sides(get_bulk(), get_meta().universal_part(), {}, false);
  expect_num_elements_and_faces_and_nodes(get_bulk(), 3u, 15u, 14u);

  get_bulk().destroy_elements_of_topology(stk::topology::HEX_8);
  expect_num_elements_and_faces_and_nodes(get_bulk(), 1u, 5u, 6u);

  stk::mesh::create_all_sides(get_bulk(), get_meta().universal_part(), {}, false);
  expect_num_elements_and_faces_and_nodes(get_bulk(), 1u, 5u, 6u);

  add_hexes_back(get_bulk());
  expect_num_elements_and_faces_and_nodes(get_bulk(), 3u, 5u, 14u);

  //    stk::mesh::create_all_sides(get_bulk(), get_meta().universal_part(), {}, false);
  //    expect_num_elements_and_faces_and_nodes(get_bulk(), 3u, 15u, 14u);
}

TEST_F(HexWedgeHexMesh, CreateFacesThenDeleteWedgeThenCreateFaces_TwoHexesRemain)
{
  expect_num_elements_and_faces_and_nodes(get_bulk(), 3u, 0u, 14u);
  stk::mesh::create_all_sides(get_bulk(), get_meta().universal_part(), {}, false);
  expect_num_elements_and_faces_and_nodes(get_bulk(), 3u, 15u, 14u);
  get_bulk().destroy_elements_of_topology(stk::topology::WEDGE_6);
  expect_num_elements_and_faces_and_nodes(get_bulk(), 2u, 12u, 14u);
  stk::mesh::create_all_sides(get_bulk(), get_meta().universal_part(), {}, false);
  expect_num_elements_and_faces_and_nodes(get_bulk(), 2u, 12u, 14u);
}

class SingleHexMesh : public stk::unit_test_util::MeshFixture
{
protected:
  const stk::mesh::EntityId firstHexId = 1;
  const stk::mesh::EntityId secondHexId = 2;
  void create_hex_on_proc_zero()
  {
    get_bulk().modification_begin();
    if(get_bulk().parallel_rank() == 0)
      stk::mesh::declare_element(get_bulk(), get_elem_part(), firstHexId, {1,2,3,4,5,6,7,8});
    get_bulk().modification_end();
  }
  void create_face_on_proc_zero()
  {
    get_bulk().modification_begin();
    if(get_bulk().parallel_rank() == 0)
      get_bulk().declare_element_side(get_bulk().get_entity(stk::topology::ELEM_RANK, firstHexId), 5, stk::mesh::ConstPartVector{});
    get_bulk().modification_end();
  }
  void create_adjacent_hex_on_last_proc()
  {
    get_bulk().modification_begin();
    if(get_bulk().parallel_rank() == 1 || get_bulk().parallel_size() == 1)
      stk::mesh::declare_element(get_bulk(), get_elem_part(), secondHexId, {5,6,7,8,9,10,11,12});

    int otherProcRank = get_bulk().parallel_rank() == 0 ? 1 : 0;
    for(stk::mesh::EntityId nodeId : {5, 6, 7, 8})
      get_bulk().add_node_sharing(get_bulk().get_entity(stk::topology::NODE_RANK, nodeId), otherProcRank);
    get_bulk().modification_end();
  }
  void expect_one_face_connected_to_two_elements()
  {
    stk::mesh::EntityVector faces;
    stk::mesh::get_entities(get_bulk(), stk::topology::FACE_RANK, faces);
    ASSERT_EQ(1u, faces.size());
    unsigned numElems = get_bulk().num_elements(faces[0]);
    ASSERT_EQ(2u, numElems);
    const stk::mesh::Entity * elems = get_bulk().begin_elements(faces[0]);
    EXPECT_EQ(firstHexId, get_bulk().identifier(elems[0]));
    EXPECT_EQ(secondHexId, get_bulk().identifier(elems[1]));
  }
  stk::mesh::Part &get_elem_part()
  {
    return get_meta().get_topology_root_part(stk::topology::HEX_8);
  }
};
// pre-existing face is not being attached to newly created element
TEST_F(SingleHexMesh, DISABLED_CreateFacesThenCreateAnotherElement_ConnectivityIsWrong)
{
  if(stk::parallel_machine_size(get_comm()) <= 2)
  {
    setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
    create_hex_on_proc_zero();
    create_face_on_proc_zero();
    create_adjacent_hex_on_last_proc();
    expect_one_face_connected_to_two_elements();
  }
}

std::string get_many_block_mesh_desc(unsigned numBlocks, unsigned nProcs, bool allBlocksOnProc0 = true)
{
  std::ostringstream oss;
  unsigned proc = 0;
  for(unsigned i = 0; i < numBlocks; ++i) {
    unsigned elemId = i + 1;
    unsigned firstNodeId = i * 4 + 1;
    oss << proc << "," << elemId << ",HEX_8,";
    for(unsigned node = firstNodeId; node < firstNodeId + 8; ++node) {
      oss << node << ",";
    }
    unsigned blockId = i + 1;
    oss << "block_" << blockId;

    if(i < numBlocks - 1) {
      oss << '\n';
    }

    if(!allBlocksOnProc0){
      proc++;
      proc = proc%nProcs;
    }
  }

  return oss.str();
}

void create_mesh_from_mesh_read(stk::mesh::BulkData& bulk,
                                const std::string& sidesetSpec,
                                const std::vector<std::string>& omittedBlocks,
                                const std::string& fileName)
{
  size_t nProc = bulk.parallel_size();

  {
    stk::mesh::MeshBuilder builder(bulk.parallel());
    builder.set_spatial_dimension(bulk.mesh_meta_data().spatial_dimension());

    std::string meshDesc = get_many_block_mesh_desc(nProc, nProc, false);
    meshDesc += sidesetSpec;

    std::shared_ptr<stk::mesh::BulkData> outBulk = builder.create();
    stk::unit_test_util::setup_text_mesh(*outBulk, meshDesc);
    stk::io::write_mesh(fileName, *outBulk);
  }

  stk::io::StkMeshIoBroker broker;

  broker.set_bulk_data(bulk);
  size_t inputIndex = broker.add_mesh_database(fileName, stk::io::READ_MESH);
  Ioss::DatabaseIO* dbIo = broker.get_input_database(inputIndex);
  dbIo->set_block_omissions(omittedBlocks);
  std::shared_ptr<Ioss::Region> inputRegion = broker.get_input_ioss_region();
  inputRegion->property_add(Ioss::Property(stk::io::s_processAllInputNodes, true));
  inputRegion->property_add(Ioss::Property(stk::io::s_ignoreDisconnectedNodes, false));
  broker.create_input_mesh();
  broker.populate_bulk_data();

  unlink(dbIo->decoded_filename().c_str());
}

TEST(CleanupOrphans, deleteOnlyNodes_withBlockOmit)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  size_t nProc = stk::parallel_machine_size(communicator);

  const unsigned spatialDim = 3;

  stk::mesh::MeshBuilder builder(communicator);
  builder.set_spatial_dimension(spatialDim);
  std::shared_ptr<stk::mesh::BulkData> bulk = builder.create();

  const std::vector<std::string> omittedBlocks{"block_1"};
  const std::string fileName = "mesh.g";
  const std::string sidesetSpec = "";
  create_mesh_from_mesh_read(*bulk, sidesetSpec, omittedBlocks, fileName);

  std::vector<size_t> counts;
  stk::mesh::comm_mesh_counts(*bulk, counts);

  size_t baseNumNodes = (nProc+1)*4;
  size_t validNumberOfNodes = baseNumNodes;
  size_t validNumberOfEdges = 0u;
  size_t validNumberOfFaces = 0u;
  EXPECT_EQ(validNumberOfNodes, counts[stk::topology::NODE_RANK]);
  EXPECT_EQ(validNumberOfEdges, counts[stk::topology::EDGE_RANK]);
  EXPECT_EQ(validNumberOfFaces, counts[stk::topology::FACE_RANK]);
  EXPECT_EQ(nProc-1, counts[stk::topology::ELEM_RANK]);

  stk::mesh::cleanup_orphan_nodes(*bulk);

  validNumberOfNodes = (nProc == 1) ? 0 : nProc*4;
  stk::mesh::comm_mesh_counts(*bulk, counts);
  EXPECT_EQ(validNumberOfNodes, counts[stk::topology::NODE_RANK]);
  EXPECT_EQ(validNumberOfEdges, counts[stk::topology::EDGE_RANK]);
  EXPECT_EQ(validNumberOfFaces, counts[stk::topology::FACE_RANK]);
}


TEST(CleanupOrphans, deleteOnlyFaces_withoutBlockOmit)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  size_t nProc = stk::parallel_machine_size(communicator);

  const unsigned spatialDim = 3;

  stk::mesh::MeshBuilder builder(communicator);
  builder.set_spatial_dimension(spatialDim);
  std::shared_ptr<stk::mesh::BulkData> bulk = builder.create();

  const std::vector<std::string> omittedBlocks{};
  const std::string fileName = "mesh.g";
  const std::string sidesetSpec = "|sideset:data=1,1, 1,2, 1,3, 1,4";
  create_mesh_from_mesh_read(*bulk, sidesetSpec, omittedBlocks, fileName);

  // Destroy element 1 leaving the orphaned faces
  stk::mesh::Entity elem1 = bulk->get_entity(stk::topology::ELEM_RANK, 1u);
  bulk->modification_begin();
  if(bulk->is_valid(elem1)) {
    bulk->destroy_entity(elem1);
  }
  bulk->modification_end();

  std::vector<size_t> counts;
  stk::mesh::comm_mesh_counts(*bulk, counts);

  size_t baseNumNodes = (nProc+1)*4;
  size_t validNumberOfNodes = baseNumNodes;
  size_t validNumberOfEdges = 0u;
  size_t validNumberOfFaces = 4u;
  EXPECT_EQ(validNumberOfNodes, counts[stk::topology::NODE_RANK]);
  EXPECT_EQ(validNumberOfEdges, counts[stk::topology::EDGE_RANK]);
  EXPECT_EQ(validNumberOfFaces, counts[stk::topology::FACE_RANK]);
  EXPECT_EQ(nProc-1, counts[stk::topology::ELEM_RANK]);

  stk::mesh::cleanup_orphan_faces(*bulk);

  validNumberOfFaces = 0u;
  stk::mesh::comm_mesh_counts(*bulk, counts);
  EXPECT_EQ(validNumberOfNodes, counts[stk::topology::NODE_RANK]);
  EXPECT_EQ(validNumberOfEdges, counts[stk::topology::EDGE_RANK]);
  EXPECT_EQ(validNumberOfFaces, counts[stk::topology::FACE_RANK]);
}

TEST(CleanupOrphans, deleteAllOrphans_withoutBlockOmit)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  size_t nProc = stk::parallel_machine_size(communicator);

  const unsigned spatialDim = 3;

  stk::mesh::MeshBuilder builder(communicator);
  builder.set_spatial_dimension(spatialDim);
  std::shared_ptr<stk::mesh::BulkData> bulk = builder.create();

  const std::vector<std::string> omittedBlocks{};
  const std::string fileName = "mesh.g";
  const std::string sidesetSpec = "|sideset:data=1,1, 1,2, 1,3, 1,4";
  create_mesh_from_mesh_read(*bulk, sidesetSpec, omittedBlocks, fileName);

  // Destroy element 1 leaving the orphaned faces
  stk::mesh::Entity elem1 = bulk->get_entity(stk::topology::ELEM_RANK, 1u);
  bulk->modification_begin();
  if(bulk->is_valid(elem1)) {
    bulk->destroy_entity(elem1);
  }
  bulk->modification_end();

  std::vector<size_t> counts;
  stk::mesh::comm_mesh_counts(*bulk, counts);

  size_t baseNumNodes = (nProc+1)*4;
  size_t validNumberOfNodes = baseNumNodes;
  size_t validNumberOfEdges = 0u;
  size_t validNumberOfFaces = 4u;
  EXPECT_EQ(validNumberOfNodes, counts[stk::topology::NODE_RANK]);
  EXPECT_EQ(validNumberOfEdges, counts[stk::topology::EDGE_RANK]);
  EXPECT_EQ(validNumberOfFaces, counts[stk::topology::FACE_RANK]);
  EXPECT_EQ(nProc-1, counts[stk::topology::ELEM_RANK]);

  stk::mesh::cleanup_all_downward_orphan_entities(*bulk, stk::topology::FACE_RANK);

  validNumberOfNodes = (nProc == 1) ? 0 : nProc*4;
  validNumberOfFaces = 0u;
  stk::mesh::comm_mesh_counts(*bulk, counts);
  EXPECT_EQ(validNumberOfNodes, counts[stk::topology::NODE_RANK]);
  EXPECT_EQ(validNumberOfEdges, counts[stk::topology::EDGE_RANK]);
  EXPECT_EQ(validNumberOfFaces, counts[stk::topology::FACE_RANK]);
}

}
