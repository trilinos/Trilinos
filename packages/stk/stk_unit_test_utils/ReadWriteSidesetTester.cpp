#include <gtest/gtest.h>
#include "stk_unit_test_utils/ReadWriteSidesetTester.hpp"

namespace stk{ namespace unit_test_util{ namespace sideset{

stk::mesh::SideSet get_stk_side_set(stk::mesh::BulkData &bulk, const ElemIdSideVector &ss)
{
    stk::mesh::SideSet sideSet(ss.size());
    for(size_t i=0; i<ss.size(); i++)
        sideSet[i] = stk::mesh::SideSetEntry(bulk.get_entity(stk::topology::ELEM_RANK, ss[i].elem_id), ss[i].side_ordinal);

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
    std::vector<int> ids = bulk.get_sideset_ids();
    ASSERT_EQ(expected.size(), ids.size()) << "for file: " << inputFileName;
    for(size_t ss=0; ss<ids.size(); ++ss)
    {
        const stk::mesh::SideSet& sideSet = bulk.get_sideset(ids[ss]);
        const std::vector<ElemIdSide>& expectedSideSet = expected[ss].sideSet;
        EXPECT_EQ(expected[ss].id, ids[ss]);
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
    std::vector<int> ids1 = bulk1.get_sideset_ids();
    std::vector<int> ids2 = bulk2.get_sideset_ids();

    ASSERT_EQ(ids1.size(), ids2.size()) << "for file: " << input_file_name;
    for(size_t ss=0; ss<ids1.size(); ++ss)
    {
        const stk::mesh::SideSet& sideSet1 = bulk1.get_sideset(ids1[ss]);
        const stk::mesh::SideSet& sideSet2 = bulk2.get_sideset(ids2[ss]);
        EXPECT_EQ(ids1[ss], ids2[ss]);
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

}
}
}
