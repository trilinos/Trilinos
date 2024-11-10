#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stddef.h>                     // for size_t, nullptr
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_io/StkIoUtils.hpp>
#include <string>                       // for string
#include "mpi.h"                        // for MPI_COMM_WORLD
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/SideSetUtil.hpp>
#include "stk_unit_test_utils/FaceTestingUtils.hpp"

namespace
{

typedef std::pair<uint64_t, int> ElementIdToExodusSide;
typedef std::vector<ElementIdToExodusSide> ElementSidePairs;
typedef std::map<int,ElementSidePairs> ExodusSideSet;

void verify_element_side_pairs(stk::mesh::BulkData& bulkData, const ExodusSideSet& goldSideset)
{
  stk::mesh::create_bulkdata_sidesets(bulkData);

  std::map<int,ElementSidePairs>::const_iterator iter = goldSideset.begin();
  for(;iter!=goldSideset.end();++iter)
  {
    int id = iter->first;
    stk::mesh::Part *part = stk::unit_test_util::get_surface_part_with_id(bulkData.mesh_meta_data(), id);
    stk::mesh::SideSet &sset = bulkData.get_sideset(*part);
    ElementSidePairs goldSet = iter->second;
    ASSERT_EQ(goldSet.size(), sset.size());
    for(size_t j=0;j<goldSet.size();++j)
    {
      stk::mesh::Entity goldEntity = bulkData.get_entity(stk::topology::ELEMENT_RANK, goldSet[j].first);
      EXPECT_TRUE(sset.contains(goldEntity, goldSet[j].second));
    }
  }
}

struct TestData
{
  std::string filename = "";
  int sidesetId = -1;
  std::string surfaceName = "surface_" + std::to_string(sidesetId);
  std::vector<std::string> blockNames={};
  ExodusSideSet goldSideSetData;
  ExodusSideSet goldSideSetDataSide1;
  ExodusSideSet goldSideSetDataSide2;
  int sideOrdinal = -1;
};

struct ABTestData_2D : public TestData
{
  ABTestData_2D() : TestData()
  {
    filename = "2D_AB.e";
    sidesetId = 1;
    surfaceName = "surface_" + std::to_string(sidesetId);
    blockNames = {"block_1", "block_2" };
    goldSideSetData  =  { {1, { {1u, 0}, {2u, 2}} } };
    goldSideSetDataSide1  = { {1, {{1u, 0}} } };
    goldSideSetDataSide2  = { {1, {{2u, 2}} } };
    sideOrdinal = 0;
  }
};

struct ABTestData : public TestData
{
  ABTestData() : TestData()
  {
    filename = "AB.e";
    sidesetId = 1;
    surfaceName = "surface_" + std::to_string(sidesetId);
    blockNames = {"block_1", "block_2" };
    goldSideSetData  =  { {1, { {2u, 4}, {1u, 5}} } };
    goldSideSetDataSide1  = { {1, {{1u, 5}} } };
    goldSideSetDataSide2  = { {1, {{2u, 4}} } };
    sideOrdinal = 5;
  }
};

void testSidesetCreation(TestData &testData)
{
  std::shared_ptr<stk::mesh::BulkData> bulkData = stk::mesh::MeshBuilder(MPI_COMM_WORLD).create();
  stk::mesh::MetaData& meta = bulkData->mesh_meta_data();
  stk::io::StkMeshIoBroker stkIo;
  stkIo.property_add(Ioss::Property("DECOMPOSITION_METHOD", "RCB"));
  stkIo.set_bulk_data(bulkData);
  stkIo.add_mesh_database(testData.filename, stk::io::READ_MESH);
  stkIo.create_input_mesh();
  stkIo.add_all_mesh_fields_as_input_fields();

  stk::mesh::Part& sideSetPart = meta.declare_part(testData.surfaceName, meta.side_rank());
  meta.set_part_id(sideSetPart, testData.sidesetId);
  stk::io::put_io_part_attribute(sideSetPart);

  stkIo.populate_bulk_data();

  int sideOrdinal = testData.sideOrdinal;
  stk::mesh::Entity elem1 = bulkData->get_entity(stk::topology::ELEM_RANK, 1);
  bulkData->modification_begin();
  bulkData->declare_element_side(elem1, sideOrdinal, stk::mesh::PartVector{&sideSetPart});
  bulkData->modification_end();

  //BEGINDOC1
  std::vector<const stk::mesh::Part*> blocks;
  for(const std::string& blockName : testData.blockNames)
  {
    stk::mesh::Part *block = meta.get_part(blockName);
    blocks.push_back(block);
  }

  meta.set_surface_to_block_mapping(&sideSetPart, blocks);
  //ENDDOC1

  verify_element_side_pairs(*bulkData, testData.goldSideSetData);
}



TEST(StkIo, readABFromIoss_verifySideSetCreation_3D)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  if(stk::parallel_machine_size(comm) == 1)
  {
    ABTestData testDataADB;
    testDataADB.blockNames = {"block_1", "block_2" };
    testSidesetCreation(testDataADB);

    ABTestData testDataALB;
    testDataALB.blockNames = {"block_1"};
    testDataALB.goldSideSetData  = testDataALB.goldSideSetDataSide1;
    testSidesetCreation(testDataALB);

    ABTestData testDataARB;
    testDataARB.blockNames = {"block_2"};
    testDataARB.goldSideSetData  = testDataALB.goldSideSetDataSide2;
    testSidesetCreation(testDataARB);
  }
}

TEST(StkIo, readABFromIoss_verifySideSetCreation_2D)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  if(stk::parallel_machine_size(comm) == 1)
  {
    ABTestData_2D testDataADB;
    testDataADB.blockNames = {"block_1", "block_2" };
    testSidesetCreation(testDataADB);

    ABTestData_2D testDataALB;
    testDataALB.blockNames = {"block_1"};
    testDataALB.goldSideSetData  = testDataALB.goldSideSetDataSide1;
    testSidesetCreation(testDataALB);

    ABTestData_2D testDataARB;
    testDataARB.blockNames = {"block_2"};
    testDataARB.goldSideSetData  = testDataALB.goldSideSetDataSide2;
    testSidesetCreation(testDataARB);
  }
}

}
