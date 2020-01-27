#include "gtest/gtest.h"
#include "mpi.h"
#include <string>
#include <vector>
#include <stk_mesh/base/Field.hpp>      // for Field
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, put_field, etc
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_io/IossBridge.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>

namespace {

class TextMeshFixture : public stk::unit_test_util::MeshFixture
{
protected:
  TextMeshFixture(unsigned spatialDim) : stk::unit_test_util::MeshFixture(spatialDim)
  { }

  void verify_shared_nodes(const stk::mesh::EntityIdVector& nodeIds, int sharingProc)
  {
    std::vector<size_t> counts;
    stk::mesh::count_entities(get_meta().globally_shared_part(), get_bulk(), counts);
    EXPECT_EQ(nodeIds.size(), counts[stk::topology::NODE_RANK]);

    for (stk::mesh::EntityId nodeId : nodeIds) {
      EXPECT_TRUE(get_bulk().in_shared(stk::mesh::EntityKey(stk::topology::NODE_RANK, nodeId), sharingProc));
    }
  }

  void verify_num_elements(size_t goldCount)
  {
    std::vector<size_t> counts;
    stk::mesh::count_entities(get_meta().universal_part(), get_bulk(), counts);
    EXPECT_EQ(goldCount, counts[stk::topology::ELEM_RANK]);
  }

  void verify_single_element(stk::mesh::EntityId elemId,
                             stk::topology topology,
                             const stk::mesh::EntityIdVector& nodeIds)
  {
    stk::mesh::Entity element = get_bulk().get_entity(stk::topology::ELEM_RANK, elemId);
    EXPECT_TRUE(get_bulk().is_valid(element));
    EXPECT_EQ(topology, get_bulk().bucket(element).topology());
    verify_nodes_on_element(element, nodeIds);
  }

  struct PartInfo
  {
    std::string blockName;
    std::set<stk::mesh::EntityId> ids;
  };

  void verify_part_membership(const std::vector<PartInfo> golds)
  {
    for (const PartInfo& gold : golds) {
      stk::mesh::Part* blockPart = get_meta().get_part(gold.blockName);

      verify_part(blockPart);
      verify_elements_on_part(blockPart, gold.ids);
    }
  }

  void verify_coordinates(const stk::mesh::EntityIdVector& goldNodeIds, const std::vector<double>& goldCoordinates)
  {
    CoordinateVerifier cv(get_bulk(), goldNodeIds, goldCoordinates);
    cv.verify();
  }

private:
  void verify_nodes_on_element(stk::mesh::Entity element, const stk::mesh::EntityIdVector& goldNodeIds)
  {
    stk::mesh::EntityVector nodes(get_bulk().begin_nodes(element), get_bulk().end_nodes(element));
    EXPECT_EQ(goldNodeIds, get_node_ids(nodes));
  }

  stk::mesh::EntityIdVector get_node_ids(const stk::mesh::EntityVector& nodes)
  {
    stk::mesh::EntityIdVector nodeIds;
    for(const stk::mesh::Entity& node : nodes) {
      nodeIds.emplace_back(get_bulk().identifier(node));
    }
    return nodeIds;
  }

  void verify_part(stk::mesh::Part* blockPart)
  {
    ASSERT_TRUE(blockPart != nullptr);
    EXPECT_TRUE(stk::io::is_part_io_part(*blockPart));
  }

  void verify_elements_on_part(stk::mesh::Part* blockPart, const std::set<stk::mesh::EntityId>& goldIds)
  {
    stk::mesh::EntityVector elems;
    stk::mesh::get_selected_entities(*blockPart, get_bulk().buckets(stk::topology::ELEM_RANK), elems);

    ASSERT_EQ(goldIds.size(), elems.size());
    for (const stk::mesh::Entity& elem : elems) {
      stk::mesh::EntityId elemId = get_bulk().identifier(elem);
      EXPECT_EQ(1u, goldIds.count(elemId));
    }
  }

  class CoordinateVerifier
  {
  public:
    CoordinateVerifier(const stk::mesh::BulkData& b,
                       const stk::mesh::EntityIdVector& ids,
                       const std::vector<double>& coords)
      : bulk(b), meta(bulk.mesh_meta_data()),
        spatialDim(meta.spatial_dimension()),
        goldNodeIds(ids), goldCoordinates(coords)
    { }

    void verify()
    {
      verify_num_nodes();

      for(size_t nodeIndex=0; nodeIndex<goldNodeIds.size(); nodeIndex++) {
        stk::mesh::EntityId nodeId = goldNodeIds[nodeIndex];
        EXPECT_TRUE(bulk.is_valid(get_node(nodeId)));

        const double* nodalCoords = get_nodal_coordinates(nodeId);
        const double* goldCoords = &goldCoordinates[nodeIndex*spatialDim];

        verify_nodal_coordinates(nodeId, goldCoords, nodalCoords);
      }
    }

  private:
    void verify_num_nodes()
    {
      stk::mesh::EntityVector nodes;
      bulk.get_entities(stk::topology::NODE_RANK, meta.universal_part(), nodes);
      EXPECT_EQ(goldNodeIds.size(), nodes.size());
    }

    const double* get_nodal_coordinates(const stk::mesh::EntityId& nodeId)
    {
      const stk::mesh::CoordinatesField& coordsField =
          static_cast<const stk::mesh::CoordinatesField&>(*meta.coordinate_field());
      return stk::mesh::field_data(coordsField, get_node(nodeId));
    }

    const stk::mesh::Entity get_node(const stk::mesh::EntityId& nodeId)
    {
      return bulk.get_entity(stk::topology::NODE_RANK, nodeId);
    }

    void verify_nodal_coordinates(const stk::mesh::EntityId& nodeId,
                                  const double* goldCoords,
                                  const double* nodalCoords)
    {
      for (unsigned i=0; i<spatialDim; i++) {
        EXPECT_NEAR(goldCoords[i], nodalCoords[i], 1.0e-9) << error_message(nodeId, i);
      }
    }

    std::string error_message(const stk::mesh::EntityId& nodeId, unsigned coordIndex)
    {
      std::stringstream message;
      message << "Proc " << bulk.parallel_rank() << ", Node ID " << nodeId << ", coord index " << coordIndex;
      return message.str();
    }

    const stk::mesh::BulkData& bulk;
    const stk::mesh::MetaData& meta;

    const unsigned spatialDim;

    const stk::mesh::EntityIdVector& goldNodeIds;
    const std::vector<double>& goldCoordinates;
  };
};

class TestTextMesh : public TextMeshFixture
{
protected:
  TestTextMesh() : TextMeshFixture(3)
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  }
};

TEST_F(TestTextMesh, singleHex)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8";
  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

  verify_num_elements(1);
  verify_single_element(1u, stk::topology::HEX_8, stk::mesh::EntityIdVector{1,2,3,4,5,6,7,8});
}

TEST_F(TestTextMesh, singleHexWithCoordinates)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8";
  stk::mesh::EntityIdVector nodeIds = stk::mesh::EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8};
  std::vector<double> coordinates = { 0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1 };

  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc, coordinates);

  verify_num_elements(1);
  verify_single_element(1u, stk::topology::HEX_8, nodeIds);
  verify_coordinates(nodeIds, coordinates);
}

TEST_F(TestTextMesh, singleHexWithCoordinates_separatedNodeIds)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1,2,3,4,11,12,13,14";
  stk::mesh::EntityIdVector nodeIds = stk::mesh::EntityIdVector{1, 2, 3, 4, 11, 12, 13, 14};
  std::vector<double> coordinates = { 0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1 };

  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc, coordinates);

  verify_num_elements(1);
  verify_single_element(1u, stk::topology::HEX_8, nodeIds);
  verify_coordinates(nodeIds, coordinates);
}

TEST_F(TestTextMesh, twoHexesSerial)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                         "0,2,HEX_8,5,6,7,8,9,10,11,12";
  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

  verify_num_elements(2);
  verify_single_element(1u, stk::topology::HEX_8, stk::mesh::EntityIdVector{1,2,3,4,5,6,7,8});
  verify_single_element(2u, stk::topology::HEX_8, stk::mesh::EntityIdVector{5,6,7,8,9,10,11,12});
}

TEST_F(TestTextMesh, twoTet10Serial)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,TET_10,1,2,3,4,5,6,7,8,9,10\n"
                         "0,2,TET_10,2,11,3,4,12,13,6,9,14,10";
//                                       1       2      3        4          5          6         7           8
  std::vector<double> coordinates = { 0,0,0, 1,0,0, 0.5,1,0, 0.5,0.5,1, 0.5,0,0, 0.75,0.5,0, 0.25,0.5,0, 0.25,0.25,0.5,
//                                       9              10            11         12           13        14
                                      0.75,0.25,0.5, 0.5,0.75,0.5, 1.5,0.5,0, 1.25,0.25,0, 1,0.75,0, 1,0.5,0.5 };

  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc, coordinates);

  verify_num_elements(2);
  verify_single_element(1u, stk::topology::TET_10, stk::mesh::EntityIdVector{1,2,3,4,5,6,7,8,9,10});
  verify_single_element(2u, stk::topology::TET_10, stk::mesh::EntityIdVector{2,11,3,4,12,13,6,9,14,10});
}

TEST_F(TestTextMesh, twoHexDisconnectedWithCoordinatesAndParts)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,block_1\n"
                         "0,2,HEX_8,9,10,11,12,13,14,15,16,block_2";
  std::vector<double> coordinates = { 0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1,
                                      0,0,1, 1,0,1, 1,1,1, 0,1,1, 0,0,2, 1,0,2, 1,1,2, 0,1,2 };

  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc, coordinates);

  verify_num_elements(2);
  verify_single_element(1u, stk::topology::HEX_8, stk::mesh::EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8});
  verify_single_element(2u, stk::topology::HEX_8, stk::mesh::EntityIdVector{9, 10, 11, 12, 13, 14, 15, 16});
  verify_coordinates(stk::mesh::EntityIdVector{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16}, coordinates);
  verify_part_membership({{"block_1", {1u}}, {"block_2", {2u}}});
}

TEST_F(TestTextMesh, twoHexDisconnectedWithDefaultParts)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8\n"
                         "0,2,HEX_8,9,10,11,12,13,14,15,16";
  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

  verify_num_elements(2);
  verify_single_element(1u, stk::topology::HEX_8, stk::mesh::EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8});
  verify_single_element(2u, stk::topology::HEX_8, stk::mesh::EntityIdVector{9, 10, 11, 12, 13, 14, 15, 16});
  verify_part_membership({{"block_HEXAHEDRON_8", {1u, 2u}}});
}

TEST_F(TestTextMesh, threeTriShellsWithCoordinatesAndParts)
{
  //      4-----5             //
  //      |\  2 |\            //
  //      |  \  |  \          //
  //      | 1  \| 3  \        //
  //      1-----2-----3       //

  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,SHELL_TRI_3,1,2,4,block_1\n"
                         "0,2,SHELL_TRI_3,2,5,4,block_2\n"
                         "0,3,SHELL_TRI_3,2,3,5,block_2";
  std::vector<double> coordinates = {
    0,0,0, 1,0,0, 2,0,0,
    0,1,0, 1,1,0
  };

  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc, coordinates);

  verify_num_elements(3);
  verify_single_element(1u, stk::topology::SHELL_TRI_3, stk::mesh::EntityIdVector{1, 2, 4});
  verify_single_element(2u, stk::topology::SHELL_TRI_3, stk::mesh::EntityIdVector{2, 5, 4});
  verify_single_element(3u, stk::topology::SHELL_TRI_3, stk::mesh::EntityIdVector{2, 3, 5});
  verify_coordinates(stk::mesh::EntityIdVector{1, 2, 3, 4, 5}, coordinates);
  verify_part_membership({{"block_1", {1u}}, {"block_2", {2u, 3u}}});
}

TEST_F(TestTextMesh, threeTriShellsWithDefaultParts)
{
  //      4-----5             //
  //      |\  2 |\            //
  //      |  \  |  \          //
  //      | 1  \| 3  \        //
  //      1-----2-----3       //

  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,SHELL_TRI_3,1,2,4\n"
                         "0,2,SHELL_TRI_3,2,5,4\n"
                         "0,3,SHELL_TRI_3,2,3,5";
  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

  verify_num_elements(3);
  verify_single_element(1u, stk::topology::SHELL_TRI_3, stk::mesh::EntityIdVector{1, 2, 4});
  verify_single_element(2u, stk::topology::SHELL_TRI_3, stk::mesh::EntityIdVector{2, 5, 4});
  verify_single_element(3u, stk::topology::SHELL_TRI_3, stk::mesh::EntityIdVector{2, 3, 5});
  verify_part_membership({{"block_SHELL_TRIANGLE_3", {1u, 2u, 3u}}});
}

TEST_F(TestTextMesh, twoHexesParallel)
{
  if (get_parallel_size() != 2) return;
  int rank = get_parallel_rank();

  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                         "1,2,HEX_8,5,6,7,8,9,10,11,12";
  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

  if (rank == 0) {
    verify_num_elements(1);
    verify_single_element(1u, stk::topology::HEX_8, stk::mesh::EntityIdVector{1,2,3,4,5,6,7,8});
    verify_shared_nodes(stk::mesh::EntityIdVector{5,6,7,8}, 1);
  }
  else if (rank == 1) {
    verify_num_elements(1);
    verify_single_element(2u, stk::topology::HEX_8, stk::mesh::EntityIdVector{5,6,7,8,9,10,11,12});
    verify_shared_nodes(stk::mesh::EntityIdVector{5,6,7,8}, 0);
  }
}

TEST_F(TestTextMesh, twoQuadShellsWithCoordinatesParallel)
{
  //      4-----5-----6
  //      |     |     |
  //      |  1  |  2  |
  //      |     |     |
  //      1-----2-----3

  if (get_parallel_size() != 2) return;
  int rank = get_parallel_rank();

  std::string meshDesc = "0,1,SHELL_QUAD_4,1,2,5,4\n"
                         "1,2,SHELL_QUAD_4,2,3,6,5";
  std::vector<double> coordinates = {
    0,0,0, 1,0,0, 2,0,0,
    0,1,0, 1,1,0, 2,1,0
  };

  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc, coordinates);

  if (rank == 0) {
    verify_num_elements(1);
    verify_single_element(1u, stk::topology::SHELL_QUAD_4, stk::mesh::EntityIdVector{1, 2, 5, 4});
    verify_coordinates(stk::mesh::EntityIdVector{1, 2, 4, 5}, {0,0,0, 1,0,0, 0,1,0, 1,1,0});
  }
  else if (rank == 1) {
    verify_num_elements(1);
    verify_single_element(2u, stk::topology::SHELL_QUAD_4, stk::mesh::EntityIdVector{2, 3, 6, 5});
    verify_coordinates(stk::mesh::EntityIdVector{2, 3, 5, 6}, {1,0,0, 2,0,0, 1,1,0, 2,1,0});
  }
}

TEST_F(TestTextMesh, threeTriShellsWithCoordinatesParallel)
{
  //      4-----5             //
  //      |\  2 |\            //
  //      |  \  |  \          //
  //      | 1  \| 3  \        //
  //      1-----2-----3       //

  if (get_parallel_size() != 2) return;
  int rank = get_parallel_rank();

  std::string meshDesc = "0,1,SHELL_TRI_3,1,2,4\n"
                         "0,2,SHELL_TRI_3,2,5,4\n"
                         "1,3,SHELL_TRI_3,2,3,5";
  std::vector<double> coordinates = {
    0,0,0, 1,0,0, 2,0,0,
    0,1,0, 1,1,0
  };

  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc, coordinates);

  if (rank == 0) {
    verify_num_elements(2);
    verify_single_element(1u, stk::topology::SHELL_TRI_3, stk::mesh::EntityIdVector{1, 2, 4});
    verify_single_element(2u, stk::topology::SHELL_TRI_3, stk::mesh::EntityIdVector{2, 5, 4});
    verify_coordinates(stk::mesh::EntityIdVector{1, 2, 4, 5}, {0,0,0, 1,0,0, 0,1,0, 1,1,0});
  }
  else if (rank == 1) {
    verify_num_elements(1);
    verify_single_element(3u, stk::topology::SHELL_TRI_3, stk::mesh::EntityIdVector{2, 3, 5});
    verify_coordinates(stk::mesh::EntityIdVector{2, 3, 5}, {1,0,0, 2,0,0, 1,1,0});
  }
}

TEST_F(TestTextMesh, singleHexWithSpaces)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0, 1, HEX_8, 1, 2, 3, 4, 5, 6, 7, 8";
  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

  verify_num_elements(1);
  verify_single_element(1u, stk::topology::HEX_8, stk::mesh::EntityIdVector{1,2,3,4,5,6,7,8});
}

TEST_F(TestTextMesh, singleHexWithCoordinatesAndSpaces)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0, 1, HEX_8, 1, 2, 3, 4, 5, 6, 7, 8";
  stk::mesh::EntityIdVector nodeIds = stk::mesh::EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8};
  std::vector<double> coordinates = { 0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1 };

  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc, coordinates);

  verify_num_elements(1);
  verify_single_element(1u, stk::topology::HEX_8, nodeIds);
  verify_coordinates(nodeIds, coordinates);
}

TEST_F(TestTextMesh, singleHexWithLowerCase)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,Hex_8,1,2,3,4,5,6,7,8";
  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

  verify_num_elements(1);
  verify_single_element(1u, stk::topology::HEX_8, stk::mesh::EntityIdVector{1,2,3,4,5,6,7,8});
}

TEST_F(TestTextMesh, singleHexWithCoordinatesAndLowerCase)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,Hex_8,1,2,3,4,5,6,7,8";
  stk::mesh::EntityIdVector nodeIds = stk::mesh::EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8};
  std::vector<double> coordinates = { 0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1 };

  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc, coordinates);

  verify_num_elements(1);
  verify_single_element(1u, stk::topology::HEX_8, nodeIds);
  verify_coordinates(nodeIds, coordinates);
}

TEST_F(TestTextMesh, tooFewNodes)
{
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7";
  EXPECT_THROW(stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc), std::logic_error);
}

TEST_F(TestTextMesh, tooFewNodesWithCoordinates)
{
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7";
  std::vector<double> coordinates = { 0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1 };
  EXPECT_THROW(stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc, coordinates), std::logic_error);
}

TEST_F(TestTextMesh, tooFewCoordinates)
{
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8";
  std::vector<double> coordinates = { 0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1 };
  EXPECT_THROW(stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc, coordinates), std::logic_error);
}

TEST_F(TestTextMesh, tooManyNodes)
{
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,9,10";
  EXPECT_THROW(stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc), std::logic_error);
}

TEST_F(TestTextMesh, tooManyNodesWithCoordinates)
{
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,9,10";
  std::vector<double> coordinates = { 0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1, 2,0,0, 2,1,0 };
  EXPECT_THROW(stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc, coordinates), std::logic_error);
}

TEST_F(TestTextMesh, tooManyCoordinates)
{
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8";
  std::vector<double> coordinates = { 0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1, 52 };
  EXPECT_THROW(stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc, coordinates), std::logic_error);
}

TEST_F(TestTextMesh, tooLittleData)
{
  std::string meshDesc = "0,1,";
  EXPECT_THROW(stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc), std::logic_error);
}

TEST_F(TestTextMesh, tooLittleDataWithCoordinates)
{
  std::string meshDesc = "0,1,";
  std::vector<double> coordinates;
  EXPECT_THROW(stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc, coordinates), std::logic_error);
}

TEST_F(TestTextMesh, invalidTopology)
{
  std::string meshDesc = "0,1,invalid,1";
  EXPECT_THROW(stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc), std::logic_error);
}

TEST_F(TestTextMesh, invalidTopologyWithCoordinates)
{
  std::string meshDesc = "0,1,invalid,1";
  std::vector<double> coordinates = { 0,0,0 };
  EXPECT_THROW(stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc, coordinates), std::logic_error);
}

TEST_F(TestTextMesh, mixedSpatialDim)
{
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                         "0,2,QUAD_4_2D,5,6,7,8";
  EXPECT_THROW(stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc), std::logic_error);
}

TEST_F(TestTextMesh, mixedSpatialDimWithCoordinates)
{
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                         "0,2,QUAD_4_2D,5,6,7,8";
  std::vector<double> coordinates = { 0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1 };
  EXPECT_THROW(stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc, coordinates), std::logic_error);
}

TEST_F(TestTextMesh, spatialDimInconsistentWithMetaData)
{
  std::string meshDesc = "0,1,QUAD_4_2D,1,2,3,4";
  EXPECT_THROW(stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc), std::logic_error);
}

TEST_F(TestTextMesh, spatialDimInconsistentWithMetaDataWithCoordinates)
{
  std::string meshDesc = "0,1,QUAD_4_2D,1,2,3,4";
  std::vector<double> coordinates = { 0,0, 1,0, 1,1, 0,1 };
  EXPECT_THROW(stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc, coordinates), std::logic_error);
}

TEST_F(TestTextMesh, particleHex)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,PARTICLE,1\n"
                         "0,2,HEX_8,2,3,4,5,6,7,8,9";
  EXPECT_NO_THROW(stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc));

  verify_num_elements(2);
  verify_single_element(1u, stk::topology::PARTICLE, stk::mesh::EntityIdVector{1});
  verify_single_element(2u, stk::topology::HEX_8, stk::mesh::EntityIdVector{2,3,4,5,6,7,8,9});
}

TEST_F(TestTextMesh, particleHexWithCoordinates)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,PARTICLE,1\n"
                         "0,2,HEX_8,2,3,4,5,6,7,8,9";
  std::vector<double> coordinates = { 2,0,0, 0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1 };
  EXPECT_NO_THROW(stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc, coordinates));

  verify_num_elements(2);
  verify_single_element(1u, stk::topology::PARTICLE, stk::mesh::EntityIdVector{1});
  verify_single_element(2u, stk::topology::HEX_8, stk::mesh::EntityIdVector{2,3,4,5,6,7,8,9});
  verify_coordinates(stk::mesh::EntityIdVector{1,2,3,4,5,6,7,8,9}, coordinates);
}

class TestTextMeshAura : public TextMeshFixture
{
protected:
  TestTextMeshAura() : TextMeshFixture(3)
  {
    setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  }
};

TEST_F(TestTextMeshAura, twoQuadShellWithCoordinates)
{
  //      4-----5-----6
  //      |     |     |
  //      |  1  |  2  |
  //      |     |     |
  //      1-----2-----3

  if (get_parallel_size() != 2) return;
  int rank = get_parallel_rank();

  std::string meshDesc = "0,1,SHELL_QUAD_4,1,2,5,4\n"
                         "1,2,SHELL_QUAD_4,2,3,6,5";
  std::vector<double> coordinates = {
    0,0,0, 1,0,0, 2,0,0,
    0,1,0, 1,1,0, 2,1,0
  };

  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc, coordinates);

  verify_num_elements(2);
  verify_single_element(1u, stk::topology::SHELL_QUAD_4, stk::mesh::EntityIdVector{1, 2, 5, 4});
  verify_single_element(2u, stk::topology::SHELL_QUAD_4, stk::mesh::EntityIdVector{2, 3, 6, 5});
  verify_coordinates(stk::mesh::EntityIdVector{1, 2, 3, 4, 5, 6}, coordinates);

  if (rank == 0) verify_shared_nodes(stk::mesh::EntityIdVector{2,5}, 1);
  if (rank == 1) verify_shared_nodes(stk::mesh::EntityIdVector{2,5}, 0);
}

class TestTextMesh2d : public TextMeshFixture
{
protected:
  TestTextMesh2d() : TextMeshFixture(2)
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  }
};

TEST_F(TestTextMesh2d, singleQuad)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,QUAD_4_2D,1,2,3,4";
  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

  verify_num_elements(1);
  verify_single_element(1u, stk::topology::QUAD_4_2D, stk::mesh::EntityIdVector{1,2,3,4});
}

TEST_F(TestTextMesh2d, twoSprings)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,SPRING_2,1,2\n"
                         "0,2,SPRING_2,2,3";
  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

  verify_num_elements(2);
  verify_single_element(1u, stk::topology::SPRING_2, stk::mesh::EntityIdVector{1,2});
  verify_single_element(2u, stk::topology::SPRING_2, stk::mesh::EntityIdVector{2,3});
}

TEST_F(TestTextMesh2d, threeQuadsWithCoordinates)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,QUAD_4_2D,1,2,3,4\n"
                         "0,2,QUAD_4_2D,2,3,5,6\n"
                         "0,3,QUAD_4_2D,5,7,8,6";
  std::vector<double> coordinates = { 0,0, 1,0, 1,1, 0,1, 2,0, 2,1, 3,0, 3,1 };

  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc, coordinates);

  verify_num_elements(3);
  verify_single_element(1u, stk::topology::QUAD_4_2D, stk::mesh::EntityIdVector{1,2,3,4});
  verify_single_element(2u, stk::topology::QUAD_4_2D, stk::mesh::EntityIdVector{2,3,5,6});
  verify_single_element(3u, stk::topology::QUAD_4_2D, stk::mesh::EntityIdVector{5,7,8,6});
  verify_coordinates(stk::mesh::EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8}, coordinates);
}

TEST_F(TestTextMesh2d, twoQuadsWithCoordinatesParallel)
{
  //      4-----5-----6
  //      |     |     |
  //      |  1  |  2  |
  //      |     |     |
  //      1-----2-----3

  if (get_parallel_size() != 2) return;
  int rank = get_parallel_rank();

  std::string meshDesc = "0,1,QUAD_4_2D,1,2,5,4\n"
                         "1,2,QUAD_4_2D,2,3,6,5";
  std::vector<double> coordinates = {
    0,0, 1,0, 2,0,
    0,1, 1,1, 2,1
  };

  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc, coordinates);

  if (rank == 0) {
    verify_num_elements(1);
    verify_single_element(1u, stk::topology::QUAD_4_2D, stk::mesh::EntityIdVector{1, 2, 5, 4});
    verify_shared_nodes(stk::mesh::EntityIdVector{2,5}, 1);
    verify_coordinates(stk::mesh::EntityIdVector{1, 2, 4, 5}, {0,0, 1,0, 0,1, 1,1});
  }
  else if (rank == 1) {
    verify_num_elements(1);
    verify_single_element(2u, stk::topology::QUAD_4_2D, stk::mesh::EntityIdVector{2, 3, 6, 5});
    verify_shared_nodes(stk::mesh::EntityIdVector{2,5}, 0);
    verify_coordinates(stk::mesh::EntityIdVector{2, 3, 5, 6}, {1,0, 2,0, 1,1, 2,1});
  }
}

TEST_F(TestTextMesh2d, twoQuadOneShellParallel)
{
  if (get_parallel_size() != 3) return;
  int rank = get_parallel_rank();

  std::string meshDesc = "0,1,QUAD_4_2D,1,2,3,4\n"
                         "1,2,QUAD_4_2D,3,4,5,6\n"
                         "2,3,SHELL_LINE_2,3,4";
  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

  if (rank == 0) {
    verify_num_elements(1);
    verify_single_element(1u, stk::topology::QUAD_4_2D, stk::mesh::EntityIdVector{1,2,3,4});
    verify_shared_nodes(stk::mesh::EntityIdVector{3,4}, 1);
    verify_shared_nodes(stk::mesh::EntityIdVector{3,4}, 2);
  }
  else if (rank == 1) {
    verify_num_elements(1);
    verify_single_element(2u, stk::topology::QUAD_4_2D, stk::mesh::EntityIdVector{3,4,5,6});
    verify_shared_nodes(stk::mesh::EntityIdVector{3,4}, 0);
    verify_shared_nodes(stk::mesh::EntityIdVector{3,4}, 2);
  }
  else if (rank == 2) {
    verify_num_elements(1);
    verify_single_element(3u, stk::topology::SHELL_LINE_2, stk::mesh::EntityIdVector{3,4});
    verify_shared_nodes(stk::mesh::EntityIdVector{3,4}, 0);
    verify_shared_nodes(stk::mesh::EntityIdVector{3,4}, 1);
  }
}

class TestTextMesh1d : public TextMeshFixture
{
protected:
  TestTextMesh1d() : TextMeshFixture(1)
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  }
};

TEST_F(TestTextMesh1d, oneDimensionNotSupported)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,LINE_2_1D,1,2";
  EXPECT_THROW(stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc), std::logic_error);
}

} // namespace
