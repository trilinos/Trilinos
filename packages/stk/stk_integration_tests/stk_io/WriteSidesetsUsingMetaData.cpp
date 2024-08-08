#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stddef.h>                     // for size_t, nullptr
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_io/StkIoUtils.hpp>

#include <string>                       // for string
#include "mpi.h"                        // for MPI_COMM_WORLD
#include <stk_io/DatabasePurpose.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_io/FillMesh.hpp>
#include <Ioss_SideSet.h>
#include <Ioss_SideBlock.h>
#include <stk_mesh/base/SideSetEntry.hpp>
#include <stk_mesh/base/SideSetUtil.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/ExodusTranslator.hpp>
#include "stk_io/WriteMesh.hpp"
#include "stk_unit_test_utils/FaceTestingUtils.hpp"
#include "stk_unit_test_utils/BuildMesh.hpp"

using stk::unit_test_util::build_mesh_no_simple_fields;

namespace
{

void verify_block_membership(const stk::io::IossBlockMembership &goldBlockMemberships, const stk::io::IossBlockMembership &blockMemberships)
{
    for(stk::io::IossBlockMembership::const_iterator goldIter = goldBlockMemberships.begin(); goldIter != goldBlockMemberships.end(); goldIter++)
    {
        std::string goldSideSetName = goldIter->first;
        stk::io::IossBlockMembership::const_iterator iossSideSetIter = blockMemberships.find(goldSideSetName);
        ASSERT_TRUE(iossSideSetIter != blockMemberships.end());
        ASSERT_EQ(goldIter->second.size(), iossSideSetIter->second.size());
        for(const std::string & goldBlockName : goldIter->second)
        {
            auto result = std::find(iossSideSetIter->second.begin(), iossSideSetIter->second.end(), goldBlockName);
            EXPECT_TRUE(result != iossSideSetIter->second.end()) << "failed to block, " << goldBlockName << ", in sideset, " << goldSideSetName << "." ;
        }
    }
}

void read_and_test_block_membership(const std::string& filename, const stk::io::IossBlockMembership &goldBlockMemberships)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    if(stk::parallel_machine_size(comm) <= 2)
    {
        stk::io::StkMeshIoBroker stkIo;
        stkIo.property_add(Ioss::Property("DECOMPOSITION_METHOD", "RCB"));
        stkIo.add_mesh_database(filename, stk::io::READ_MESH);
        stkIo.create_input_mesh();
        stkIo.add_all_mesh_fields_as_input_fields();
        stk::io::IossBlockMembership blockMemberships = stk::io::get_block_memberships(stkIo);
        verify_block_membership(goldBlockMemberships, blockMemberships);
    }
}

TEST(StkIo, readALBFromIoss_verifyBlockMembership)
{
    stk::io::IossBlockMembership goldBlockMem = { {"surface_1", {"block_1"}} };
    read_and_test_block_membership("ALB.e", goldBlockMem);
}

TEST(StkIo, readARBFromIoss_verifyBlockMembership)
{
    stk::io::IossBlockMembership goldBlockMem = { {"surface_1", {"block_2"}} };
    read_and_test_block_membership("ARB.e", goldBlockMem);
}

TEST(StkIo, readADBFromIoss_verifyBlockMembership)
{
    stk::io::IossBlockMembership goldBlockMem = { {"surface_1", {"block_2", "block_1"}} };
    read_and_test_block_membership("ADB.e", goldBlockMem);
}

TEST(StkIo, readALRBFromIoss_verifyBlockMembership)
{
    stk::io::IossBlockMembership goldBlockMem = { {"surface_1", {"block_1"}}, {"surface_2", {"block_2"}} };
    read_and_test_block_membership("ALRB.e", goldBlockMem);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

typedef std::pair<uint64_t, int> ElementIdToExodusSide;
typedef std::vector<ElementIdToExodusSide> ElementSidePairs;
typedef std::map<int,ElementSidePairs> ExodusSideSet;

ExodusSideSet set_gold_data(ExodusSideSet &goldSideSetData, ExodusSideSet &goldSideSetData0, ExodusSideSet &goldSideSetData1)
{
    ExodusSideSet goldData;
    if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1)
        goldData = goldSideSetData;
    else if(stk::parallel_machine_rank(MPI_COMM_WORLD) == 0)
        goldData = goldSideSetData0;
    else
        goldData = goldSideSetData1;
    return goldData;
}

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
        bool found_matching_sset_entry = false;
        for (auto && goldEntry : goldSet)
        {
          for (auto && sset_entry : sset)
          {
            if (goldEntry.first == bulkData.identifier(sset_entry.element) &&
                goldEntry.second == static_cast<int>(sset_entry.side))
            {
              found_matching_sset_entry = true;
              break;
            }
          }
          EXPECT_TRUE(found_matching_sset_entry);
        }
    }
}

void fill_mesh(stk::mesh::BulkData& bulkData, const std::string& filename)
{
    stk::io::StkMeshIoBroker stkIo;
    stkIo.property_add(Ioss::Property("DECOMPOSITION_METHOD", "RCB"));
    stkIo.set_bulk_data(bulkData);
    stkIo.add_mesh_database(filename, stk::io::READ_MESH);
    stkIo.create_input_mesh();
    stkIo.add_all_mesh_fields_as_input_fields();
    stkIo.populate_bulk_data();
}

std::string get_output_file_name(const std::string &input_file_name)
{
    size_t original_length = input_file_name.length();
    std::string base_file_name = input_file_name.substr(0, original_length-2); // assuming .e files here
    std::string output_file_name = base_file_name + "_new.e";
    return output_file_name;
}

void write_mesh(const std::string& filename, stk::mesh::BulkData &bulkData)
{
    std::string output_file_name = get_output_file_name(filename);
    stk::mesh::create_bulkdata_sidesets(bulkData);
    stk::io::write_mesh(output_file_name, bulkData);
}

void fill_sideset_data_structure_and_test(const std::string& filename, const ExodusSideSet &goldSideset)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    if(stk::parallel_machine_size(comm) <= 2)
    {
        std::shared_ptr<stk::mesh::BulkData> bulkData = build_mesh_no_simple_fields(comm);
        fill_mesh(*bulkData, filename);

        verify_element_side_pairs(*bulkData, goldSideset);
        write_mesh(filename, *bulkData);
    }
}

TEST(StkIo, readALBFromIoss_verifySideSetCreation)
{
    ExodusSideSet goldSideSetData  = { {1, {{1u, 5}} } };
    ExodusSideSet goldSideSetData0 = { {1, {{1u, 5}} } };
    ExodusSideSet goldSideSetData1 = { };
    ExodusSideSet goldData = set_gold_data(goldSideSetData, goldSideSetData0, goldSideSetData1);
    fill_sideset_data_structure_and_test("ALB.e", goldData);
}

TEST(StkIo, readARBFromIoss_verifySideSetCreation)
{
    ExodusSideSet goldSideSetData  = { {1, {{2u, 4}} } };
    ExodusSideSet goldSideSetData0 = { };
    ExodusSideSet goldSideSetData1 = { {1, {{2u, 4}} } };
    ExodusSideSet goldData = set_gold_data(goldSideSetData, goldSideSetData0, goldSideSetData1);
    fill_sideset_data_structure_and_test("ARB.e", goldData);
}

TEST(StkIo, readADBFromIoss_verifySideSetCreation)
{
    ExodusSideSet  goldSideSetData  =  { {1, {{2u, 4}, {1u, 5}} } };
    ExodusSideSet  goldSideSetData0 = { {1, {{1u, 5}} } };
    ExodusSideSet  goldSideSetData1 = { {1, {{2u, 4}} } };
    ExodusSideSet  goldData = set_gold_data(goldSideSetData, goldSideSetData0, goldSideSetData1);
    fill_sideset_data_structure_and_test("ADB.e", goldData);
}

TEST(StkIo, readALRBFromIoss_verifySideSetCreation)
{
    ExodusSideSet  goldSideSetData  = { {1, {{1u, 5}} }, {2, {{2u, 4}} } };
    ExodusSideSet  goldSideSetData0 = { {1, {{1u, 5}} } };
    ExodusSideSet  goldSideSetData1 = { {2, {{2u, 4}} } };
    ExodusSideSet  goldData = set_gold_data(goldSideSetData, goldSideSetData0, goldSideSetData1);
    fill_sideset_data_structure_and_test("ALRB.e", goldData);
}

void fill_sideset_data_structure_and_test_output_on_block_1(const std::string& filename, const ExodusSideSet &goldSideset)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    if(stk::parallel_machine_size(comm) <= 2)
    {
        std::string output_file_name = get_output_file_name(filename);
        {
            std::shared_ptr<stk::mesh::BulkData> bulkData = build_mesh_no_simple_fields(comm);
            stk::mesh::MetaData& meta = bulkData->mesh_meta_data();
            fill_mesh(*bulkData, filename);

            stk::mesh::Part* part = meta.get_part("block_1");
            ASSERT_TRUE(part!=nullptr);
            stk::mesh::Selector subset = *part;
            stk::io::write_mesh_subset(output_file_name, *bulkData, subset);
        }

        {
            std::shared_ptr<stk::mesh::BulkData> bulkData = build_mesh_no_simple_fields(comm);
            stk::io::fill_mesh(output_file_name, *bulkData);
            verify_element_side_pairs(*bulkData, goldSideset);
        }
    }
}

TEST(StkIo, readADBFromIoss_writeBlock1Subset_verifySideSetCreation)
{
    ExodusSideSet  goldSideSetData  = { {1, {{1u, 5}} } };
    ExodusSideSet  goldSideSetData0 = { {1, {{1u, 5}} } };
    ExodusSideSet  goldSideSetData1 = {};
    ExodusSideSet  goldData = set_gold_data(goldSideSetData, goldSideSetData0, goldSideSetData1);
    fill_sideset_data_structure_and_test_output_on_block_1("ADB.e", goldData);
}

}
