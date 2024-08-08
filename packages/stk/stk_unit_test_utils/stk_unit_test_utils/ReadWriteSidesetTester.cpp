// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include "stk_unit_test_utils/ReadWriteSidesetTester.hpp"
#include <gtest/gtest.h>                   // for AssertHelper, EXPECT_EQ, etc
#include <unistd.h>                        // for unlink
#include <ostream>                         // for operator<<
#include <string>                          // for string, operator<<
#include <vector>                          // for vector
#include "gtest/gtest-message.h"           // for Message
#include "stk_io/DatabasePurpose.hpp"
#include "stk_io/StkMeshIoBroker.hpp"      // for StkMeshIoBroker
#include "stk_mesh/base/BulkData.hpp"      // for BulkData
#include "stk_mesh/base/MetaData.hpp"      // for MetaData
#include "stk_mesh/base/SideSetEntry.hpp"  // for SideSet, SideSetEntry
#include "stk_mesh/base/Types.hpp"         // for EntityId, operator<<
#include "stk_util/parallel/Parallel.hpp"  // for ParallelMachine
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk{ namespace unit_test_util{ namespace sideset{

stk::mesh::SideSet* get_stk_side_set(stk::mesh::BulkData &bulk, const ElemIdSideVector &ss)
{
    stk::mesh::SideSet* sideSet = new stk::mesh::SideSet(bulk);
    for(size_t i=0; i<ss.size(); i++)
        sideSet->add(stk::mesh::SideSetEntry(bulk.get_entity(stk::topology::ELEM_RANK, ss[i].elem_id), ss[i].side_ordinal));

    return sideSet;
}

stk::unit_test_util::sideset::SideSetData get_stk_side_set_data(stk::mesh::BulkData &bulk, const SideSetIdAndElemIdSidesVector &ssData)
{
    stk::unit_test_util::sideset::SideSetData sideSetData(ssData.size());

    for(size_t i=0; i<ssData.size(); i++)
    {
        sideSetData[i].id = ssData[i].id;
        sideSetData[i].sideSet = get_stk_side_set(bulk, ssData[i].sideSet);
    }
    return sideSetData;
}

void compare_sidesets(const std::string& inputFileName,
                      stk::mesh::BulkData &bulk,
                      const SideSetIdAndElemIdSidesVector &expected)
{
    ASSERT_EQ(expected.size(), bulk.get_number_of_sidesets()) << "for file: " << inputFileName;
    for(size_t ss=0; ss<bulk.get_number_of_sidesets(); ++ss)
    {
        stk::mesh::Part *surface_part = get_surface_part_with_id(bulk.mesh_meta_data(), expected[ss].id);
        STK_ThrowRequire(surface_part != nullptr);
        const stk::mesh::SideSet& sideSet = bulk.get_sideset(*surface_part);
        const std::vector<ElemIdSide>& expectedSideSet = expected[ss].sideSet;
        EXPECT_EQ(expected[ss].id, surface_part->id());
        ASSERT_EQ(expectedSideSet.size(), sideSet.size()) << "for file: " << inputFileName;

        for(size_t i=0;i<sideSet.size();++i)
        {
            EXPECT_EQ(static_cast<stk::mesh::EntityId>(expectedSideSet[i].elem_id), bulk.identifier(sideSet[i].element)) << "for file: " << inputFileName;
            EXPECT_EQ(expectedSideSet[i].side_ordinal, sideSet[i].side) << "for file: " << inputFileName;
        }
    }
}

void setup_io_broker_for_read(stk::io::StkMeshIoBroker &stkIo, stk::mesh::BulkData &bulkData, std::string filename, ReadMode read_mode)
{
    if(read_mode == READ_SERIAL_AND_DECOMPOSE)
        stkIo.property_add(Ioss::Property("DECOMPOSITION_METHOD", "RCB"));
    stkIo.set_bulk_data(bulkData);
    stkIo.add_mesh_database(filename, stk::io::READ_MESH);
    stkIo.create_input_mesh();
    stkIo.add_all_mesh_fields_as_input_fields();
}

void load_mesh_and_fill_sideset_data(StkMeshIoBrokerTester &stkIo)
{
    stkIo.populate_bulk_data();
}


void read_exo_file( stk::mesh::BulkData &bulkData, std::string filename, ReadMode read_mode)
{
    StkMeshIoBrokerTester stkIo;
    setup_io_broker_for_read(stkIo, bulkData, filename, read_mode);
    stkIo.populate_bulk_data();
}

void write_exo_file(BulkDataTester &bulkData, const std::string &filename)
{
    StkMeshIoBrokerTester io_broker;
    io_broker.set_bulk_data(bulkData);
    size_t resultFileIndex = io_broker.create_output_mesh(filename, stk::io::WRITE_RESULTS);
    io_broker.write_output_mesh(resultFileIndex);
}

void compare_sidesets(const std::string& input_file_name,
                      BulkDataTester &bulk1,
                      BulkDataTester &bulk2)
{
    stk::mesh::SideSetVector ss1 = bulk1.get_sidesets();
    stk::mesh::SideSetVector ss2 = bulk2.get_sidesets();

    ASSERT_EQ(ss1.size(), ss2.size()) << "for file: " << input_file_name;
    for(size_t ss=0; ss<ss1.size(); ++ss)
    {
        const stk::mesh::SideSet& sideSet1 = *(ss1[ss]);
        const stk::mesh::SideSet& sideSet2 = *(ss2[ss]);

        ASSERT_EQ(sideSet1.size(), sideSet2.size()) << "for file: " << input_file_name;

        for(size_t i=0;i<sideSet1.size();++i)
        {
            EXPECT_EQ(bulk1.identifier(sideSet1[i].element), bulk2.identifier(sideSet2[i].element)) << "for file: " << input_file_name;
            EXPECT_EQ(sideSet1[i].side, sideSet2[i].side) << "for file: " << input_file_name;
        }
    }
}

void test_reading_writing_sideset_from_file(stk::ParallelMachine comm, const std::string& input_file_name, const std::string& output_file_name)
{
    stk::mesh::MetaData meta1;
    BulkDataTester bulk1(meta1, comm);
    read_exo_file( bulk1, input_file_name, READ_SERIAL_AND_DECOMPOSE);
    write_exo_file( bulk1, output_file_name);

    stk::mesh::MetaData meta2;
    BulkDataTester bulk2(meta2, comm);
    read_exo_file( bulk2, output_file_name, READ_ALREADY_DECOMPOSED);
    compare_sidesets(input_file_name, bulk1, bulk2);
    unlink(output_file_name.c_str());
}

namespace simple_fields {

stk::mesh::SideSet* get_stk_side_set(stk::mesh::BulkData &bulk, const ElemIdSideVector &ss)
{
    stk::mesh::SideSet* sideSet = new stk::mesh::SideSet(bulk);
    for(size_t i=0; i<ss.size(); i++)
        sideSet->add(stk::mesh::SideSetEntry(bulk.get_entity(stk::topology::ELEM_RANK, ss[i].elem_id), ss[i].side_ordinal));

    return sideSet;
}

stk::unit_test_util::sideset::SideSetData get_stk_side_set_data(stk::mesh::BulkData &bulk, const SideSetIdAndElemIdSidesVector &ssData)
{
    stk::unit_test_util::sideset::SideSetData sideSetData(ssData.size());

    for(size_t i=0; i<ssData.size(); i++)
    {
        sideSetData[i].id = ssData[i].id;
        sideSetData[i].sideSet = stk::unit_test_util::sideset::get_stk_side_set(bulk, ssData[i].sideSet);
    }
    return sideSetData;
}

void compare_sidesets(const std::string& inputFileName,
                      stk::mesh::BulkData &bulk,
                      const SideSetIdAndElemIdSidesVector &expected)
{
    ASSERT_EQ(expected.size(), bulk.get_number_of_sidesets()) << "for file: " << inputFileName;
    for(size_t ss=0; ss<bulk.get_number_of_sidesets(); ++ss)
    {
        stk::mesh::Part *surface_part = get_surface_part_with_id(bulk.mesh_meta_data(), expected[ss].id);
        STK_ThrowRequire(surface_part != nullptr);
        const stk::mesh::SideSet& sideSet = bulk.get_sideset(*surface_part);
        const std::vector<ElemIdSide>& expectedSideSet = expected[ss].sideSet;
        EXPECT_EQ(expected[ss].id, surface_part->id());
        ASSERT_EQ(expectedSideSet.size(), sideSet.size()) << "for file: " << inputFileName;

        for(size_t i=0;i<sideSet.size();++i)
        {
            EXPECT_EQ(static_cast<stk::mesh::EntityId>(expectedSideSet[i].elem_id), bulk.identifier(sideSet[i].element)) << "for file: " << inputFileName;
            EXPECT_EQ(expectedSideSet[i].side_ordinal, sideSet[i].side) << "for file: " << inputFileName;
        }
    }
}

void setup_io_broker_for_read(stk::io::StkMeshIoBroker &stkIo, stk::mesh::BulkData &bulkData, std::string filename, ReadMode read_mode)
{
    if(read_mode == READ_SERIAL_AND_DECOMPOSE)
        stkIo.property_add(Ioss::Property("DECOMPOSITION_METHOD", "RCB"));
    stkIo.set_bulk_data(bulkData);
    stkIo.add_mesh_database(filename, stk::io::READ_MESH);
    stkIo.create_input_mesh();
    stkIo.add_all_mesh_fields_as_input_fields();
}

void load_mesh_and_fill_sideset_data(StkMeshIoBrokerTester &stkIo)
{
    stkIo.populate_bulk_data();
}


void read_exo_file( stk::mesh::BulkData &bulkData, std::string filename, ReadMode read_mode)
{
    StkMeshIoBrokerTester stkIo;
    stk::unit_test_util::sideset::setup_io_broker_for_read(stkIo, bulkData, filename, read_mode);
    stkIo.populate_bulk_data();
}

void write_exo_file(BulkDataTester &bulkData, const std::string &filename)
{
    StkMeshIoBrokerTester io_broker;
    io_broker.set_bulk_data(bulkData);
    size_t resultFileIndex = io_broker.create_output_mesh(filename, stk::io::WRITE_RESULTS);
    io_broker.write_output_mesh(resultFileIndex);
}

void compare_sidesets(const std::string& input_file_name,
                      BulkDataTester &bulk1,
                      BulkDataTester &bulk2)
{
    stk::mesh::SideSetVector ss1 = bulk1.get_sidesets();
    stk::mesh::SideSetVector ss2 = bulk2.get_sidesets();

    ASSERT_EQ(ss1.size(), ss2.size()) << "for file: " << input_file_name;
    for(size_t ss=0; ss<ss1.size(); ++ss)
    {
        const stk::mesh::SideSet& sideSet1 = *(ss1[ss]);
        const stk::mesh::SideSet& sideSet2 = *(ss2[ss]);

        ASSERT_EQ(sideSet1.size(), sideSet2.size()) << "for file: " << input_file_name;

        for(size_t i=0;i<sideSet1.size();++i)
        {
            EXPECT_EQ(bulk1.identifier(sideSet1[i].element), bulk2.identifier(sideSet2[i].element)) << "for file: " << input_file_name;
            EXPECT_EQ(sideSet1[i].side, sideSet2[i].side) << "for file: " << input_file_name;
        }
    }
}

void test_reading_writing_sideset_from_file(stk::ParallelMachine comm, const std::string& input_file_name, const std::string& output_file_name)
{
    stk::mesh::MetaData meta1;
    stk::unit_test_util::sideset::BulkDataTester bulk1(meta1, comm);
    stk::unit_test_util::sideset::read_exo_file( bulk1, input_file_name, READ_SERIAL_AND_DECOMPOSE);
    stk::unit_test_util::sideset::write_exo_file( bulk1, output_file_name);

    stk::mesh::MetaData meta2;
    stk::unit_test_util::sideset::BulkDataTester bulk2(meta2, comm);
    stk::unit_test_util::sideset::read_exo_file( bulk2, output_file_name, READ_ALREADY_DECOMPOSED);
    stk::unit_test_util::sideset::compare_sidesets(input_file_name, bulk1, bulk2);
    unlink(output_file_name.c_str());
}

} // namespace simple_fields

}
}
}
