#include <gtest/gtest.h>
#include "stk_unit_test_utils/ReadWriteSidesetTester.hpp"

namespace stk{ namespace unit_test_util{ namespace sideset{

void read_exo_file( stk::mesh::BulkData &bulkData, std::string filename, SideSetData &sideset_data, ReadMode read_mode)
{
    StkMeshIoBrokerTester stkIo;
    if(read_mode == READ_SERIAL_AND_DECOMPOSE)
        stkIo.property_add(Ioss::Property("DECOMPOSITION_METHOD", "RCB"));
    stkIo.set_bulk_data(bulkData);
    stkIo.add_mesh_database(filename, stk::io::READ_MESH);
    stkIo.create_input_mesh();
    stkIo.add_all_mesh_fields_as_input_fields();
    stkIo.populate_bulk_data();

    stkIo.fill_sideset_data(sideset_data);
}

void write_exo_file(BulkDataTester &bulkData, const std::string &filename, const SideSetData& sideset_data)
{
    for(const IdAndSideSet& sset : sideset_data)
        bulkData.tester_save_sideset_data(sset.id, sset.sideSet);

    StkMeshIoBrokerTester io_broker;
    io_broker.set_bulk_data(bulkData);
    size_t resultFileIndex = io_broker.create_output_mesh(filename, stk::io::WRITE_RESULTS);
    io_broker.write_output_mesh(resultFileIndex);
}

void fill_sideset_data_from_serial_input_file_and_write_decomposed_file(BulkDataTester &bulk, const std::string& input_filename, const std::string& output_file_name, SideSetData& sideset_data)
{
    read_exo_file( bulk, input_filename, sideset_data, READ_SERIAL_AND_DECOMPOSE);
    write_exo_file( bulk, output_file_name, sideset_data);
}

void fill_sideset_data_from_decomposed_input_file(BulkDataTester &bulk, const std::string& input_filename, SideSetData& sideset_data)
{
    read_exo_file( bulk, input_filename, sideset_data, READ_ALREADY_DECOMPOSED);
}

void compare_sidesets(const std::string& input_file_name,
                      BulkDataTester &bulk1,
                      const SideSetData &sideset_data1,
                      BulkDataTester &bulk2,
                      const SideSetData &sideset_data2)
{
    ASSERT_EQ(sideset_data1.size(), sideset_data2.size()) << "for file: " << input_file_name;
    for(size_t ss=0; ss<sideset_data1.size(); ++ss)
    {
        const stk::mesh::SideSet& sideSet1 = sideset_data1[ss].sideSet;
        const stk::mesh::SideSet& sideSet2 = sideset_data2[ss].sideSet;
        EXPECT_EQ(sideset_data1[ss].id, sideset_data2[ss].id);
        ASSERT_EQ(sideSet1.size(), sideSet2.size()) << "for file: " << input_file_name;

        for(size_t i=0;i<sideSet1.size();++i)
        {
            EXPECT_EQ(bulk1.identifier(sideSet1[i].element), bulk2.identifier(sideSet2[i].element)) << "for file: " << input_file_name;
            EXPECT_EQ(sideSet1[i].side, sideSet2[i].side) << "for file: " << input_file_name;
        }
    }
}

SideSetData get_sideset_data_from_file_read_and_write_new_file(BulkDataTester &bulk, const std::string& input_file_name, const std::string& output_file_name)
{
    SideSetData sideset_data1;
    fill_sideset_data_from_serial_input_file_and_write_decomposed_file(bulk, input_file_name, output_file_name, sideset_data1);
    return sideset_data1;
}

SideSetData get_sideset_data_from_written_file(BulkDataTester &bulk, const std::string& output_file_name)
{
    SideSetData sideset_data2;
    fill_sideset_data_from_decomposed_input_file(bulk, output_file_name, sideset_data2);
    return sideset_data2;
}

void test_reading_writing_sideset_from_file(stk::ParallelMachine comm, const std::string& input_file_name, const std::string& output_file_name)
{
    stk::mesh::MetaData meta1;
    BulkDataTester bulk1(meta1, comm);
    SideSetData sideset_data1 = get_sideset_data_from_file_read_and_write_new_file(bulk1, input_file_name, output_file_name);
    stk::mesh::MetaData meta2;
    BulkDataTester bulk2(meta2, comm);
    SideSetData sideset_data2 = get_sideset_data_from_written_file(bulk2, output_file_name);
    compare_sidesets(input_file_name, bulk1, sideset_data1, bulk2, sideset_data2);
    unlink(output_file_name.c_str());
}

}
}
}
