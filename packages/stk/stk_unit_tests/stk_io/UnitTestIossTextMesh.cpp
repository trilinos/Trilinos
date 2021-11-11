#include <init/Ionit_Initializer.h>  // for Initializer

#include <sstream>
#include <stk_io/IossBridge.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>  // for Field
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>  // for MetaData, put_field, etc
#include <stk_unit_test_utils/GetMeshSpec.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <stk_util/environment/FileUtils.hpp>
#include <string>
#include <vector>

#include "Ioss_DBUsage.h"         // for DatabaseUsage, etc
#include "Ioss_DatabaseIO.h"      // for DatabaseIO
#include "Ioss_EntityType.h"      // for EntityType, etc
#include "Ioss_Field.h"           // for Field, etc
#include "Ioss_GroupingEntity.h"  // for GroupingEntity
#include "Ioss_IOFactory.h"       // for IOFactory
#include "Ioss_MeshType.h"        // for MeshType, etc
#include "Ioss_NodeBlock.h"       // for NodeBlock
#include "Ioss_NodeSet.h"         // for NodeSet
#include "Ioss_Property.h"        // for Property
#include "Ioss_Region.h"          // for Region, etc
#include "Ioss_SideBlock.h"       // for NodeSet
#include "Ioss_SideSet.h"         // for NodeSet
#include "gtest/gtest.h"

#include "mpi.h"
#include "stk_unit_test_utils/Iotm_DatabaseIO.hpp"
#include "stk_unit_test_utils/TextMeshFixture.hpp"

namespace
{
class IossTextMeshFixture : public stk::unit_test_util::TextMeshFixture
{
 protected:
  IossTextMeshFixture(unsigned spatialDim) : stk::unit_test_util::TextMeshFixture(spatialDim)
  {
    Ioss::Init::Initializer io;
    Iotm::IOFactory::factory();
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  }

  void create_ioss_region()
  {
    if (m_region == nullptr) {
      EXPECT_TRUE(m_database != nullptr);
      m_region = new Ioss::Region(m_database, "input_model");
      EXPECT_TRUE(m_region != nullptr);
    }
  }

  void create_factory(const std::string& fileName, const std::string& meshType)
  {
    Ioss::DatabaseUsage db_usage = Ioss::READ_MODEL;

    std::string meshFileName = fileName;
    stk::util::filename_substitution(meshFileName);
    m_database = Ioss::IOFactory::create(meshType, meshFileName, db_usage, get_comm(), m_propertyManager);
    EXPECT_TRUE(m_database != nullptr);
    EXPECT_TRUE(m_database->ok(true));
    EXPECT_EQ(m_database->get_format(), "TextMesh");
  }

  std::pair<std::string, std::string> get_database_type_and_filename(const std::string& meshDesc)
  {
    std::string type;
    std::string filename;

    size_t colon = meshDesc.find(':');
    if (colon != std::string::npos && colon > 0) {
      type = meshDesc.substr(0, colon);
      filename = meshDesc.substr(colon + 1);
    } else {
      type = "textmesh";
      filename = meshDesc;
    }

    return std::make_pair(type, filename);
  }

  void create_factory(const std::string& meshDesc)
  {
    std::pair<std::string, std::string> result = get_database_type_and_filename(meshDesc);

    std::string type = result.first;
    std::string filename = result.second;

    EXPECT_EQ("textmesh", type);
    create_factory(filename, type);
  }

  std::string get_mesh_desc(const std::string& textMeshDesc)
  {
    std::string header = "textmesh:";
    std::string meshDesc = header + textMeshDesc;
    return meshDesc;
  }

  std::string get_mesh_desc(const std::string& textMeshDesc, std::vector<double>& coordVec)
  {
    std::stringstream coords;
    coords << "|coordinates:";

    for (double coord : coordVec) {
      coords << coord << ",";
    }

    std::string meshDesc = get_mesh_desc(textMeshDesc) + coords.str();
    return meshDesc;
  }

  std::string get_mesh_desc(const std::string& textMeshDesc, unsigned dimension)
  {
    std::stringstream dim;
    dim << "|dimension:" << dimension;

    std::string meshDesc = get_mesh_desc(textMeshDesc) + dim.str();
    return meshDesc;
  }

  std::string get_mesh_desc(const std::string& textMeshDesc, std::vector<double>& coordVec, unsigned dimension)
  {
    std::stringstream dim;
    dim << "|dimension:" << dimension;

    std::string meshDesc = get_mesh_desc(textMeshDesc, coordVec) + dim.str();
    return meshDesc;
  }

  void fill_mesh(const std::string& meshDesc)
  {
    stk::io::fill_mesh(meshDesc, get_bulk(), m_broker);
    Teuchos::RCP<Ioss::Region> iossRegion = m_broker.get_input_io_region();
    m_region = iossRegion.get();
    m_database = m_region->get_database();
    EXPECT_EQ(m_database->get_format(), "TextMesh");
  }

  Ioss::PropertyManager m_propertyManager;
  Ioss::DatabaseIO* m_database = nullptr;
  Ioss::Region* m_region = nullptr;
  stk::io::StkMeshIoBroker m_broker;
};

class TestIossTextMesh : public IossTextMeshFixture
{
 protected:
  TestIossTextMesh() : IossTextMeshFixture(3) {}
};

TEST_F(TestIossTextMesh, canCreateFactory)
{
  if (get_parallel_size() != 1) return;

  std::string textMeshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8";
  std::string meshDesc = get_mesh_desc(textMeshDesc);

  create_factory(meshDesc);
}

TEST_F(TestIossTextMesh, canCreateIossRegion_withoutCoordinates)
{
  if (get_parallel_size() != 1) return;

  std::string textMeshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8";
  std::string meshDesc = get_mesh_desc(textMeshDesc);

  create_factory(meshDesc);
  create_ioss_region();
}

TEST_F(TestIossTextMesh, canCreateIossRegion_withCoordinates)
{
  if (get_parallel_size() != 1) return;

  std::string textMeshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8";
  std::vector<double> coords = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1};
  std::string meshDesc = get_mesh_desc(textMeshDesc, coords);

  create_factory(meshDesc);
  create_ioss_region();
}

TEST_F(TestIossTextMesh, canFillSerialMesh_withCoordinates)
{
  if (get_parallel_size() != 1) return;

  std::string textMeshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8";
  std::vector<double> coords = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1};
  std::string meshDesc = get_mesh_desc(textMeshDesc, coords);

  fill_mesh(meshDesc);
}

TEST_F(TestIossTextMesh, canFillParallelMesh_withCoordinates)
{
  if (get_parallel_size() != 2) return;

  std::string textMeshDesc =
      "0,1,HEX_8,1,2,3,4,5, 6, 7, 8,block_1\n"
      "1,2,HEX_8,5,6,7,8,9,10,11,12,block_1";
  std::vector<double> coords = {
      0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 2, 1, 0, 2, 1, 1, 2, 0, 1, 2};
  std::string meshDesc = get_mesh_desc(textMeshDesc, coords);

  fill_mesh(meshDesc);
}

TEST_F(TestIossTextMesh, singleHex)
{
  if (get_parallel_size() != 1) return;

  std::string textMeshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8";
  stk::mesh::EntityIdVector nodeIds = stk::mesh::EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8};

  fill_mesh(get_mesh_desc(textMeshDesc));

  verify_num_elements(1);
  verify_single_element(1u, stk::topology::HEX_8, nodeIds);
}

TEST_F(TestIossTextMesh, singleHexWithCoordinates)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8";
  stk::mesh::EntityIdVector nodeIds = stk::mesh::EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8};
  std::vector<double> coordinates = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1};

  fill_mesh(get_mesh_desc(meshDesc, coordinates));

  verify_num_elements(1);
  verify_single_element(1u, stk::topology::HEX_8, nodeIds);
  verify_coordinates(nodeIds, coordinates);
}

TEST_F(TestIossTextMesh, singleHexWithCoordinates_separatedNodeIds)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1,2,3,4,11,12,13,14";
  stk::mesh::EntityIdVector nodeIds = stk::mesh::EntityIdVector{1, 2, 3, 4, 11, 12, 13, 14};
  std::vector<double> coordinates = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1};

  fill_mesh(get_mesh_desc(meshDesc, coordinates));

  verify_num_elements(1);
  verify_single_element(1u, stk::topology::HEX_8, nodeIds);
  verify_coordinates(nodeIds, coordinates);
}

TEST_F(TestIossTextMesh, twoHexesSerial)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc =
      "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
      "0,2,HEX_8,5,6,7,8,9,10,11,12";
  fill_mesh(get_mesh_desc(meshDesc));

  verify_num_elements(2);
  verify_single_element(1u, stk::topology::HEX_8, stk::mesh::EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8});
  verify_single_element(2u, stk::topology::HEX_8, stk::mesh::EntityIdVector{5, 6, 7, 8, 9, 10, 11, 12});
}

TEST_F(TestIossTextMesh, twoTet10Serial)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc =
      "0,1,TET_10,1,2,3,4,5,6,7,8,9,10\n"
      "0,2,TET_10,2,11,3,4,12,13,6,9,14,10";
  //                                       1       2      3        4          5          6         7           8
  std::vector<double> coordinates = {0, 0, 0, 1, 0, 0, 0.5, 1, 0, 0.5, 0.5, 1, 0.5, 0, 0, 0.75, 0.5, 0, 0.25, 0.5, 0,
      0.25, 0.25, 0.5,
      //                                       9              10            11         12           13        14
      0.75, 0.25, 0.5, 0.5, 0.75, 0.5, 1.5, 0.5, 0, 1.25, 0.25, 0, 1, 0.75, 0, 1, 0.5, 0.5};

  fill_mesh(get_mesh_desc(meshDesc, coordinates));

  verify_num_elements(2);
  verify_single_element(1u, stk::topology::TET_10, stk::mesh::EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8, 9, 10});
  verify_single_element(2u, stk::topology::TET_10, stk::mesh::EntityIdVector{2, 11, 3, 4, 12, 13, 6, 9, 14, 10});
}

TEST_F(TestIossTextMesh, twoHexDisconnectedWithCoordinatesAndParts)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc =
      "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,block_1\n"
      "0,2,HEX_8,9,10,11,12,13,14,15,16,block_2";
  std::vector<double> coordinates = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1,
      0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 2, 1, 0, 2, 1, 1, 2, 0, 1, 2};

  fill_mesh(get_mesh_desc(meshDesc, coordinates));

  verify_num_elements(2);
  verify_single_element(1u, stk::topology::HEX_8, stk::mesh::EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8});
  verify_single_element(2u, stk::topology::HEX_8, stk::mesh::EntityIdVector{9, 10, 11, 12, 13, 14, 15, 16});
  verify_coordinates(stk::mesh::EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16}, coordinates);
  verify_part_membership({{"block_1", {1u}}, {"block_2", {2u}}});
}

TEST_F(TestIossTextMesh, twoHexDisconnectedWithCoordinatesAndParts_reversedBlockmembership)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc =
      "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,block_2\n"
      "0,2,HEX_8,9,10,11,12,13,14,15,16,block_1";
  std::vector<double> coordinates = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1,
      0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 2, 1, 0, 2, 1, 1, 2, 0, 1, 2};

  fill_mesh(get_mesh_desc(meshDesc, coordinates));

  verify_num_elements(2);
  verify_single_element(1u, stk::topology::HEX_8, stk::mesh::EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8});
  verify_single_element(2u, stk::topology::HEX_8, stk::mesh::EntityIdVector{9, 10, 11, 12, 13, 14, 15, 16});
  verify_coordinates(stk::mesh::EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16}, coordinates);
  verify_part_membership({{"block_2", {1u}}, {"block_1", {2u}}});
}

TEST_F(TestIossTextMesh, twoHexDisconnectedWithCoordinatesAndSameParts_reversedElementListing)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc =
      "0,2,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,block_1\n"
      "0,1,HEX_8,9,10,11,12,13,14,15,16,block_1";
  std::vector<double> coordinates = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1,
      0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 2, 1, 0, 2, 1, 1, 2, 0, 1, 2};

  fill_mesh(get_mesh_desc(meshDesc, coordinates));

  verify_num_elements(2);
  verify_single_element(2u, stk::topology::HEX_8, stk::mesh::EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8});
  verify_single_element(1u, stk::topology::HEX_8, stk::mesh::EntityIdVector{9, 10, 11, 12, 13, 14, 15, 16});
  verify_coordinates(stk::mesh::EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16}, coordinates);
  verify_part_membership({{"block_1", {1u, 2u}}});
}

TEST_F(TestIossTextMesh, twoHexDisconnectedWithDefaultParts)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc =
      "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8\n"
      "0,2,HEX_8,9,10,11,12,13,14,15,16";
  fill_mesh(get_mesh_desc(meshDesc));

  verify_num_elements(2);
  verify_single_element(1u, stk::topology::HEX_8, stk::mesh::EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8});
  verify_single_element(2u, stk::topology::HEX_8, stk::mesh::EntityIdVector{9, 10, 11, 12, 13, 14, 15, 16});
  verify_part_membership({{"block_HEXAHEDRON_8", {1u, 2u}}});
}

TEST_F(TestIossTextMesh, threeTriShellsWithCoordinatesAndParts)
{
  //      4-----5             //
  //      |\  2 |\            //
  //      |  \  |  \          //
  //      | 1  \| 3  \        //
  //      1-----2-----3       //

  if (get_parallel_size() != 1) return;

  std::string meshDesc =
      "0,1,SHELL_TRI_3,1,2,4,block_1\n"
      "0,2,SHELL_TRI_3,2,5,4,block_2\n"
      "0,3,SHELL_TRI_3,2,3,5,block_2";
  std::vector<double> coordinates = {0, 0, 0, 1, 0, 0, 2, 0, 0, 0, 1, 0, 1, 1, 0};

  fill_mesh(get_mesh_desc(meshDesc, coordinates));

  verify_num_elements(3);
  verify_single_element(1u, stk::topology::SHELL_TRI_3, stk::mesh::EntityIdVector{1, 2, 4});
  verify_single_element(2u, stk::topology::SHELL_TRI_3, stk::mesh::EntityIdVector{2, 5, 4});
  verify_single_element(3u, stk::topology::SHELL_TRI_3, stk::mesh::EntityIdVector{2, 3, 5});
  verify_coordinates(stk::mesh::EntityIdVector{1, 2, 3, 4, 5}, coordinates);
  verify_part_membership({{"block_1", {1u}}, {"block_2", {2u, 3u}}});
}

TEST_F(TestIossTextMesh, threeTriShellsWithDefaultParts)
{
  //      4-----5             //
  //      |\  2 |\            //
  //      |  \  |  \          //
  //      | 1  \| 3  \        //
  //      1-----2-----3       //

  if (get_parallel_size() != 1) return;

  std::string meshDesc =
      "0,1,SHELL_TRI_3,1,2,4\n"
      "0,2,SHELL_TRI_3,2,5,4\n"
      "0,3,SHELL_TRI_3,2,3,5";
  fill_mesh(get_mesh_desc(meshDesc));

  verify_num_elements(3);
  verify_single_element(1u, stk::topology::SHELL_TRI_3, stk::mesh::EntityIdVector{1, 2, 4});
  verify_single_element(2u, stk::topology::SHELL_TRI_3, stk::mesh::EntityIdVector{2, 5, 4});
  verify_single_element(3u, stk::topology::SHELL_TRI_3, stk::mesh::EntityIdVector{2, 3, 5});
  verify_part_membership({{"block_SHELL_TRIANGLE_3", {1u, 2u, 3u}}});
}

TEST_F(TestIossTextMesh, partIds_oneDefaultPartOneElem)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8";
  fill_mesh(get_mesh_desc(meshDesc));

  verify_part_ids({{"block_HEXAHEDRON_8", 1u}});
}

TEST_F(TestIossTextMesh, partIds_oneDefaultPartTwoElems)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc =
      "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8\n"
      "0,2,HEX_8,9,10,11,12,13,14,15,16";
  fill_mesh(get_mesh_desc(meshDesc));

  verify_part_ids({{"block_HEXAHEDRON_8", 1u}});
}

TEST_F(TestIossTextMesh, partIds_twoDefaultParts)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc =
      "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8\n"
      "0,2,HEX_8,9,10,11,12,13,14,15,16\n"
      "0,3,SHELL_QUAD_4,17,18,19,20";
  fill_mesh(get_mesh_desc(meshDesc));

  verify_part_ids({{"block_HEXAHEDRON_8", 1u}, {"block_SHELL_QUADRILATERAL_4", 2u}});
}

TEST_F(TestIossTextMesh, partIds_oneDefaultPartOneUserSpecifiedPart)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc =
      "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8\n"
      "0,2,HEX_8,9,10,11,12,13,14,15,16,my_cool_part";
  fill_mesh(get_mesh_desc(meshDesc));

  verify_part_ids({{"block_HEXAHEDRON_8", 1u}, {"my_cool_part", 2u}});
}

TEST_F(TestIossTextMesh, partIds_samePartTwice)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc =
      "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,my_cool_part\n"
      "0,2,HEX_8,9,10,11,12,13,14,15,16,my_cool_part";
  fill_mesh(get_mesh_desc(meshDesc));

  verify_part_ids({{"my_cool_part", 1u}});
}

TEST_F(TestIossTextMesh, partIds_orderingIsByLine)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc =
      "0,2,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,partOne\n"
      "0,1,HEX_8,9,10,11,12,13,14,15,16,partTwo";
  fill_mesh(get_mesh_desc(meshDesc));

  verify_part_ids({{"partOne", 1u}, {"partTwo", 2u}});
}

TEST_F(TestIossTextMesh, partIds_respectExodusNamingConvention_onePart)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_101";
  fill_mesh(get_mesh_desc(meshDesc));

  verify_part_ids({{"block_101", 101u}});
}

TEST_F(TestIossTextMesh, partIds_respectExodusNamingConvention_twoParts)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc =
      "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,block_101\n"
      "0,2,HEX_8,9,10,11,12,13,14,15,16,block_201";
  fill_mesh(get_mesh_desc(meshDesc));

  verify_part_ids({{"block_101", 101u}, {"block_201", 201u}});
}

TEST_F(TestIossTextMesh, partIds_respectExodusNamingConvention_withDefaultPartFirst)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc =
      "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8\n"
      "0,2,HEX_8,9,10,11,12,13,14,15,16,block_101";
  fill_mesh(get_mesh_desc(meshDesc));

  verify_part_ids({{"block_HEXAHEDRON_8", 1u}, {"block_101", 101u}});
}

TEST_F(TestIossTextMesh, partIds_respectExodusNamingConvention_withDefaultPartSecond)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc =
      "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,block_101\n"
      "0,2,HEX_8,9,10,11,12,13,14,15,16";
  fill_mesh(get_mesh_desc(meshDesc));

  verify_part_ids({{"block_101", 101u}, {"block_HEXAHEDRON_8", 1u}});
}

TEST_F(TestIossTextMesh, partIds_respectExodusNamingConvention_withDefaultPartIdCollision)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc =
      "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8\n"
      "0,2,HEX_8,9,10,11,12,13,14,15,16,block_1";
  fill_mesh(get_mesh_desc(meshDesc));

  verify_part_ids({{"block_HEXAHEDRON_8", 2u}, {"block_1", 1u}});
}

TEST_F(TestIossTextMesh, partIds_userSpecifiedPartId_onePart)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,partThree,3";
  fill_mesh(get_mesh_desc(meshDesc));

  verify_part_ids({{"partThree", 3u}});
  verify_num_elements(1);
  verify_single_element(1u, stk::topology::HEX_8, stk::mesh::EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8});
  verify_part_membership({{"partThree", {1u}}});
}

TEST_F(TestIossTextMesh, partIds_userSpecifiedPartId_twoParts)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc =
      "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,partThree,3\n"
      "0,2,HEX_8,9,10,11,12,13,14,15,16,partFive,5";
  fill_mesh(get_mesh_desc(meshDesc));

  verify_part_ids({{"partThree", 3u}, {"partFive", 5u}});
}

TEST_F(TestIossTextMesh, partIds_userSpecifiedPartId_twoPartsSameId)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc =
      "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,partFour,4\n"
      "0,2,HEX_8,9,10,11,12,13,14,15,16,partFour,4";
  fill_mesh(get_mesh_desc(meshDesc));

  verify_part_ids({{"partFour", 4u}});
}

TEST_F(TestIossTextMesh, partIds_userSpecifiedPartId_samePartDifferentIds)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc =
      "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,partFour,4\n"
      "0,2,HEX_8,9,10,11,12,13,14,15,16,partFour,5";
  EXPECT_THROW(fill_mesh(get_mesh_desc(meshDesc)), std::logic_error);
}

TEST_F(TestIossTextMesh, partIds_userSpecifiedPartId_withDefaultPartIdCollision)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc =
      "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8\n"
      "0,2,HEX_8,9,10,11,12,13,14,15,16,partOne,1";
  fill_mesh(get_mesh_desc(meshDesc));

  verify_part_ids({{"block_HEXAHEDRON_8", 2u}, {"partOne", 1u}});
}

TEST_F(TestIossTextMesh, partIds_userSpecifiedPartId_withDefaultPart)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc =
      "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,partThree,3\n"
      "0,2,HEX_8,9,10,11,12,13,14,15,16";
  fill_mesh(get_mesh_desc(meshDesc));

  verify_part_ids({{"partThree", 3u}, {"block_HEXAHEDRON_8", 1u}});
}

TEST_F(TestIossTextMesh, partIds_userSpecifiedPartId_withExodusPart)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc =
      "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,partThree,3\n"
      "0,2,HEX_8,9,10,11,12,13,14,15,16,block_4";
  fill_mesh(get_mesh_desc(meshDesc));

  verify_part_ids({{"partThree", 3u}, {"block_4", 4u}});
}

TEST_F(TestIossTextMesh, partIds_userSpecifiedPartId_collidesWithExodusPart)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc =
      "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,partThree,3\n"
      "0,2,HEX_8,9,10,11,12,13,14,15,16,block_3";
  EXPECT_THROW(fill_mesh(get_mesh_desc(meshDesc)), std::logic_error);
}

TEST_F(TestIossTextMesh, partIds_userSpecifiedPartId_collidesWithPreviousSpec)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc =
      "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,partThreeA,3\n"
      "0,2,HEX_8,9,10,11,12,13,14,15,16,partThreeB,3";
  EXPECT_THROW(fill_mesh(get_mesh_desc(meshDesc)), std::logic_error);
}

TEST_F(TestIossTextMesh, partIds_userSpecifiedPartId_forExodusPart)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,block_2,3\n";
  EXPECT_THROW(fill_mesh(get_mesh_desc(meshDesc)), std::logic_error);
}

TEST_F(TestIossTextMesh, partIds_shortPartNamesAreValid)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc =
      "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,a\n"
      "0,2,HEX_8,9,10,11,12,13,14,15,16,b,3";
  fill_mesh(get_mesh_desc(meshDesc));

  verify_part_ids({{"a", 1u}, {"b", 3u}});
}

TEST_F(TestIossTextMesh, partIds_integerPartNamesAreInvalid)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,9";
  EXPECT_THROW(fill_mesh(get_mesh_desc(meshDesc)), std::logic_error);
}

TEST_F(TestIossTextMesh, twoHexesParallel)
{
  if (get_parallel_size() != 2) return;
  int rank = get_parallel_rank();

  std::string meshDesc =
      "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
      "1,2,HEX_8,5,6,7,8,9,10,11,12";
  fill_mesh(get_mesh_desc(meshDesc));

  if (rank == 0) {
    verify_num_elements(1);
    verify_single_element(1u, stk::topology::HEX_8, stk::mesh::EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8});
    verify_shared_nodes(stk::mesh::EntityIdVector{5, 6, 7, 8}, 1);
  } else if (rank == 1) {
    verify_num_elements(1);
    verify_single_element(2u, stk::topology::HEX_8, stk::mesh::EntityIdVector{5, 6, 7, 8, 9, 10, 11, 12});
    verify_shared_nodes(stk::mesh::EntityIdVector{5, 6, 7, 8}, 0);
  }
}

TEST_F(TestIossTextMesh, twoQuadShellsWithCoordinatesParallel)
{
  //      4-----5-----6
  //      |     |     |
  //      |  1  |  2  |
  //      |     |     |
  //      1-----2-----3

  if (get_parallel_size() != 2) return;
  int rank = get_parallel_rank();

  std::string meshDesc =
      "0,1,SHELL_QUAD_4,1,2,5,4\n"
      "1,2,SHELL_QUAD_4,2,3,6,5";
  std::vector<double> coordinates = {0, 0, 0, 1, 0, 0, 2, 0, 0, 0, 1, 0, 1, 1, 0, 2, 1, 0};

  fill_mesh(get_mesh_desc(meshDesc, coordinates));

  if (rank == 0) {
    verify_num_elements(1);
    verify_single_element(1u, stk::topology::SHELL_QUAD_4, stk::mesh::EntityIdVector{1, 2, 5, 4});
    verify_coordinates(stk::mesh::EntityIdVector{1, 2, 4, 5}, {0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0});
  } else if (rank == 1) {
    verify_num_elements(1);
    verify_single_element(2u, stk::topology::SHELL_QUAD_4, stk::mesh::EntityIdVector{2, 3, 6, 5});
    verify_coordinates(stk::mesh::EntityIdVector{2, 3, 5, 6}, {1, 0, 0, 2, 0, 0, 1, 1, 0, 2, 1, 0});
  }
}

TEST_F(TestIossTextMesh, threeTriShellsWithCoordinatesParallel)
{
  //      4-----5             //
  //      |\  2 |\            //
  //      |  \  |  \          //
  //      | 1  \| 3  \        //
  //      1-----2-----3       //

  if (get_parallel_size() != 2) return;
  int rank = get_parallel_rank();

  std::string meshDesc =
      "0,1,SHELL_TRI_3,1,2,4\n"
      "0,2,SHELL_TRI_3,2,5,4\n"
      "1,3,SHELL_TRI_3,2,3,5";
  std::vector<double> coordinates = {0, 0, 0, 1, 0, 0, 2, 0, 0, 0, 1, 0, 1, 1, 0};

  fill_mesh(get_mesh_desc(meshDesc, coordinates));

  if (rank == 0) {
    verify_num_elements(2);
    verify_single_element(1u, stk::topology::SHELL_TRI_3, stk::mesh::EntityIdVector{1, 2, 4});
    verify_single_element(2u, stk::topology::SHELL_TRI_3, stk::mesh::EntityIdVector{2, 5, 4});
    verify_coordinates(stk::mesh::EntityIdVector{1, 2, 4, 5}, {0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0});
  } else if (rank == 1) {
    verify_num_elements(1);
    verify_single_element(3u, stk::topology::SHELL_TRI_3, stk::mesh::EntityIdVector{2, 3, 5});
    verify_coordinates(stk::mesh::EntityIdVector{2, 3, 5}, {1, 0, 0, 2, 0, 0, 1, 1, 0});
  }
}

TEST_F(TestIossTextMesh, singleHexWithSpaces)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0, 1, HEX_8, 1, 2, 3, 4, 5, 6, 7, 8";
  fill_mesh(get_mesh_desc(meshDesc));

  verify_num_elements(1);
  verify_single_element(1u, stk::topology::HEX_8, stk::mesh::EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8});
}

TEST_F(TestIossTextMesh, singleHexWithCoordinatesAndSpaces)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0, 1, HEX_8, 1, 2, 3, 4, 5, 6, 7, 8";
  stk::mesh::EntityIdVector nodeIds = stk::mesh::EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8};
  std::vector<double> coordinates = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1};

  fill_mesh(get_mesh_desc(meshDesc, coordinates));

  verify_num_elements(1);
  verify_single_element(1u, stk::topology::HEX_8, nodeIds);
  verify_coordinates(nodeIds, coordinates);
}

TEST_F(TestIossTextMesh, singleHexWithLowerCase)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,Hex_8,1,2,3,4,5,6,7,8";
  fill_mesh(get_mesh_desc(meshDesc));

  verify_num_elements(1);
  verify_single_element(1u, stk::topology::HEX_8, stk::mesh::EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8});
}

TEST_F(TestIossTextMesh, singleHexWithCoordinatesAndLowerCase)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,Hex_8,1,2,3,4,5,6,7,8";
  stk::mesh::EntityIdVector nodeIds = stk::mesh::EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8};
  std::vector<double> coordinates = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1};

  fill_mesh(get_mesh_desc(meshDesc, coordinates));

  verify_num_elements(1);
  verify_single_element(1u, stk::topology::HEX_8, nodeIds);
  verify_coordinates(nodeIds, coordinates);
}

TEST_F(TestIossTextMesh, tooFewNodes)
{
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7";
  EXPECT_THROW(fill_mesh(get_mesh_desc(meshDesc)), std::logic_error);
}

TEST_F(TestIossTextMesh, tooFewNodesWithCoordinates)
{
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7";
  std::vector<double> coordinates = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1};
  EXPECT_THROW(fill_mesh(get_mesh_desc(meshDesc, coordinates)), std::logic_error);
}

TEST_F(TestIossTextMesh, tooFewCoordinates)
{
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8";
  std::vector<double> coordinates = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1};
  EXPECT_THROW(fill_mesh(get_mesh_desc(meshDesc, coordinates)), std::logic_error);
}

TEST_F(TestIossTextMesh, tooManyNodes)
{
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,9,10";
  EXPECT_THROW(fill_mesh(get_mesh_desc(meshDesc)), std::logic_error);
}

TEST_F(TestIossTextMesh, tooManyNodesWithCoordinates)
{
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,9,10";
  std::vector<double> coordinates = {
      0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 2, 0, 0, 2, 1, 0};
  EXPECT_THROW(fill_mesh(get_mesh_desc(meshDesc, coordinates)), std::logic_error);
}

TEST_F(TestIossTextMesh, tooManyCoordinates)
{
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8";
  std::vector<double> coordinates = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 52};
  EXPECT_THROW(fill_mesh(get_mesh_desc(meshDesc, coordinates)), std::logic_error);
}

TEST_F(TestIossTextMesh, tooLittleData_empty)
{
  std::string meshDesc = "";
  EXPECT_NO_THROW(fill_mesh(get_mesh_desc(meshDesc)));
}

TEST_F(TestIossTextMesh, tooLittleData_startsWithString)
{
  std::string meshDesc = "hi";
  EXPECT_THROW(fill_mesh(get_mesh_desc(meshDesc)), std::logic_error);
}

TEST_F(TestIossTextMesh, tooLittleData_noGlobalId)
{
  std::string meshDesc = "0 ";
  EXPECT_THROW(fill_mesh(get_mesh_desc(meshDesc)), std::logic_error);
}

TEST_F(TestIossTextMesh, tooLittleData_noTopology)
{
  std::string meshDesc = "0,1,";
  EXPECT_THROW(fill_mesh(get_mesh_desc(meshDesc)), std::logic_error);
}

TEST_F(TestIossTextMesh, tooLittleData_noNodes)
{
  std::string meshDesc = "0,1,HEX_8";
  EXPECT_THROW(fill_mesh(get_mesh_desc(meshDesc)), std::logic_error);
}

TEST_F(TestIossTextMesh, tooLittleDataWithCoordinates)
{
  std::string meshDesc = "0,1,";
  std::vector<double> coordinates;
  EXPECT_THROW(fill_mesh(get_mesh_desc(meshDesc, coordinates)), std::logic_error);
}

TEST_F(TestIossTextMesh, invalidTopology)
{
  std::string meshDesc = "0,1,invalid,1";
  EXPECT_THROW(fill_mesh(get_mesh_desc(meshDesc)), std::logic_error);
}

TEST_F(TestIossTextMesh, invalidTopologyWithCoordinates)
{
  std::string meshDesc = "0,1,invalid,1";
  std::vector<double> coordinates = {0, 0, 0};
  EXPECT_THROW(fill_mesh(get_mesh_desc(meshDesc, coordinates)), std::logic_error);
}

TEST_F(TestIossTextMesh, mixedSpatialDim)
{
  std::string meshDesc =
      "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
      "0,2,QUAD_4_2D,5,6,7,8";
  EXPECT_THROW(fill_mesh(get_mesh_desc(meshDesc)), std::logic_error);
}

TEST_F(TestIossTextMesh, mixedSpatialDimWithCoordinates)
{
  std::string meshDesc =
      "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
      "0,2,QUAD_4_2D,5,6,7,8";
  std::vector<double> coordinates = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1};
  EXPECT_THROW(fill_mesh(get_mesh_desc(meshDesc, coordinates)), std::logic_error);
}

TEST_F(TestIossTextMesh, spatialDimInconsistentWithMetaData)
{
  std::string meshDesc = "0,1,QUAD_4_2D,1,2,3,4";
  EXPECT_THROW(fill_mesh(get_mesh_desc(meshDesc)), std::logic_error);
}

TEST_F(TestIossTextMesh, spatialDimInconsistentWithMetaDataWithCoordinates)
{
  std::string meshDesc = "0,1,QUAD_4_2D,1,2,3,4";
  std::vector<double> coordinates = {0, 0, 1, 0, 1, 1, 0, 1};
  EXPECT_THROW(fill_mesh(get_mesh_desc(meshDesc, coordinates)), std::logic_error);
}

TEST_F(TestIossTextMesh, endingWithNewlineIsOk)
{
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n";
  EXPECT_NO_THROW(fill_mesh(get_mesh_desc(meshDesc)));
}

TEST_F(TestIossTextMesh, stringAfterPartNameIsError)
{
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1,bogus\n";
  EXPECT_THROW(fill_mesh(get_mesh_desc(meshDesc)), std::logic_error);
}

TEST_F(TestIossTextMesh, particleHex)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc =
      "0,1,PARTICLE,1\n"
      "0,2,HEX_8,2,3,4,5,6,7,8,9";
  EXPECT_NO_THROW(fill_mesh(get_mesh_desc(meshDesc)));

  verify_num_elements(2);
  verify_single_element(1u, stk::topology::PARTICLE, stk::mesh::EntityIdVector{1});
  verify_single_element(2u, stk::topology::HEX_8, stk::mesh::EntityIdVector{2, 3, 4, 5, 6, 7, 8, 9});
}

TEST_F(TestIossTextMesh, particleHexWithCoordinates)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc =
      "0,1,PARTICLE,1\n"
      "0,2,HEX_8,2,3,4,5,6,7,8,9";
  std::vector<double> coordinates = {2, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1};
  EXPECT_NO_THROW(fill_mesh(get_mesh_desc(meshDesc, coordinates)));

  verify_num_elements(2);
  verify_single_element(1u, stk::topology::PARTICLE, stk::mesh::EntityIdVector{1});
  verify_single_element(2u, stk::topology::HEX_8, stk::mesh::EntityIdVector{2, 3, 4, 5, 6, 7, 8, 9});
  verify_coordinates(stk::mesh::EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8, 9}, coordinates);
}

///////////////////////////////////////////////////////////////
class TestIossTextMesh2d : public IossTextMeshFixture
{
 protected:
  TestIossTextMesh2d() : IossTextMeshFixture(2) {}
};

TEST_F(TestIossTextMesh2d, singleQuad)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,QUAD_4_2D,1,2,3,4";
  fill_mesh(get_mesh_desc(meshDesc, 2u));

  verify_num_elements(1);
  verify_single_element(1u, stk::topology::QUAD_4_2D, stk::mesh::EntityIdVector{1, 2, 3, 4});
}

TEST_F(TestIossTextMesh2d, twoSprings)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc =
      "0,1,SPRING_2,1,2\n"
      "0,2,SPRING_2,2,3";
  fill_mesh(get_mesh_desc(meshDesc, 2u));

  verify_num_elements(2);
  verify_single_element(1u, stk::topology::SPRING_2, stk::mesh::EntityIdVector{1, 2});
  verify_single_element(2u, stk::topology::SPRING_2, stk::mesh::EntityIdVector{2, 3});
}

TEST_F(TestIossTextMesh2d, threeQuadsWithCoordinates)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc =
      "0,1,QUAD_4_2D,1,2,3,4\n"
      "0,2,QUAD_4_2D,2,3,5,6\n"
      "0,3,QUAD_4_2D,5,7,8,6";
  std::vector<double> coordinates = {0, 0, 1, 0, 1, 1, 0, 1, 2, 0, 2, 1, 3, 0, 3, 1};

  fill_mesh(get_mesh_desc(meshDesc, coordinates, 2u));

  verify_num_elements(3);
  verify_single_element(1u, stk::topology::QUAD_4_2D, stk::mesh::EntityIdVector{1, 2, 3, 4});
  verify_single_element(2u, stk::topology::QUAD_4_2D, stk::mesh::EntityIdVector{2, 3, 5, 6});
  verify_single_element(3u, stk::topology::QUAD_4_2D, stk::mesh::EntityIdVector{5, 7, 8, 6});
  verify_coordinates(stk::mesh::EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8}, coordinates);
}

TEST_F(TestIossTextMesh2d, twoQuadsWithCoordinatesParallel)
{
  //      4-----5-----6
  //      |     |     |
  //      |  1  |  2  |
  //      |     |     |
  //      1-----2-----3

  if (get_parallel_size() != 2) return;
  int rank = get_parallel_rank();

  std::string meshDesc =
      "0,1,QUAD_4_2D,1,2,5,4\n"
      "1,2,QUAD_4_2D,2,3,6,5";
  std::vector<double> coordinates = {0, 0, 1, 0, 2, 0, 0, 1, 1, 1, 2, 1};

  fill_mesh(get_mesh_desc(meshDesc, coordinates, 2u));

  if (rank == 0) {
    verify_num_elements(1);
    verify_single_element(1u, stk::topology::QUAD_4_2D, stk::mesh::EntityIdVector{1, 2, 5, 4});
    verify_shared_nodes(stk::mesh::EntityIdVector{2, 5}, 1);
    verify_coordinates(stk::mesh::EntityIdVector{1, 2, 4, 5}, {0, 0, 1, 0, 0, 1, 1, 1});
  } else if (rank == 1) {
    verify_num_elements(1);
    verify_single_element(2u, stk::topology::QUAD_4_2D, stk::mesh::EntityIdVector{2, 3, 6, 5});
    verify_shared_nodes(stk::mesh::EntityIdVector{2, 5}, 0);
    verify_coordinates(stk::mesh::EntityIdVector{2, 3, 5, 6}, {1, 0, 2, 0, 1, 1, 2, 1});
  }
}

TEST_F(TestIossTextMesh2d, twoQuadOneShellParallel)
{
  if (get_parallel_size() != 3) return;
  int rank = get_parallel_rank();

  std::string meshDesc =
      "0,1,QUAD_4_2D,1,2,3,4\n"
      "1,2,QUAD_4_2D,3,4,5,6\n"
      "2,3,SHELL_LINE_2,3,4";
  fill_mesh(get_mesh_desc(meshDesc, 2u));

  if (rank == 0) {
    verify_num_elements(1);
    verify_single_element(1u, stk::topology::QUAD_4_2D, stk::mesh::EntityIdVector{1, 2, 3, 4});
    verify_shared_nodes(stk::mesh::EntityIdVector{3, 4}, 1);
    verify_shared_nodes(stk::mesh::EntityIdVector{3, 4}, 2);
  } else if (rank == 1) {
    verify_num_elements(1);
    verify_single_element(2u, stk::topology::QUAD_4_2D, stk::mesh::EntityIdVector{3, 4, 5, 6});
    verify_shared_nodes(stk::mesh::EntityIdVector{3, 4}, 0);
    verify_shared_nodes(stk::mesh::EntityIdVector{3, 4}, 2);
  } else if (rank == 2) {
    verify_num_elements(1);
    verify_single_element(3u, stk::topology::SHELL_LINE_2, stk::mesh::EntityIdVector{3, 4});
    verify_shared_nodes(stk::mesh::EntityIdVector{3, 4}, 0);
    verify_shared_nodes(stk::mesh::EntityIdVector{3, 4}, 1);
  }
}

class TestIossTextMesh1d : public IossTextMeshFixture
{
 protected:
  TestIossTextMesh1d() : IossTextMeshFixture(1) {}
};

TEST_F(TestIossTextMesh1d, oneDimensionNotSupported)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,LINE_2_1D,1,2";
  EXPECT_THROW(fill_mesh(get_mesh_desc(meshDesc, 1u)), std::logic_error);
}

}  // namespace
