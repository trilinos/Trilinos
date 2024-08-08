#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_io/FillMesh.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <fstream>
#include <stk_unit_test_utils/BuildMesh.hpp>

namespace
{

using stk::unit_test_util::build_mesh;
using stk::unit_test_util::build_mesh_no_simple_fields;

class AttributesInFile : public stk::unit_test_util::MeshFixture
{
protected:
    const std::string filename = "fileWithAttr.e";

    void expect_ordered_attributes(const stk::io::StkMeshIoBroker &stkIo)
    {
        const stk::mesh::Part *blockPart = get_meta().get_part("block_10");
        ASSERT_NE(nullptr, blockPart);
        stk::mesh::FieldVector attrFields = stkIo.get_ordered_attribute_fields(blockPart);
        std::vector<std::string> expectedNames = {"area", "i1", "i2", "j", "orient", "offset"};
        ASSERT_EQ(expectedNames.size(), attrFields.size());
        for(size_t i=0; i<expectedNames.size(); i++)
            EXPECT_EQ(expectedNames[i], attrFields[i]->name());
    }

    void read_and_write(const std::string &outputFilename)
    {
        setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
        stk::io::StkMeshIoBroker stkIo(get_comm());
        stk::io::fill_mesh_preexisting(stkIo, filename, get_bulk());

        stk::io::StkMeshIoBroker outputStkIo(get_comm());
        outputStkIo.set_bulk_data(get_bulk());
        size_t outputFileIndex = outputStkIo.create_output_mesh(outputFilename, stk::io::WRITE_RESULTS);
        outputStkIo.set_attribute_field_ordering_stored_by_part_ordinal(stkIo.get_attribute_field_ordering_stored_by_part_ordinal());
        outputStkIo.write_output_mesh(outputFileIndex);
    }
};

TEST_F(AttributesInFile, reading_fieldsReturnedInOrder)
{
    if (!std::ifstream(filename)) {
      std::cout<<"file "<<filename<<" doesn't exist, skipping test."<<std::endl;
      return;
    }

    if(stk::parallel_machine_size(get_comm()) == 1)
    {
        setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
        stk::io::StkMeshIoBroker stkIo(get_comm());
        stk::io::fill_mesh_preexisting(stkIo, filename, get_bulk());

        expect_ordered_attributes(stkIo);
    }
}

TEST_F(AttributesInFile, readWriteRead_fieldsReturnedInOrder)
{
    if (!std::ifstream(filename)) {
      std::cout<<"file "<<filename<<" doesn't exist, skipping test."<<std::endl;
      return;
    }

    if(stk::parallel_machine_size(get_comm()) == 1)
    {
        const std::string outputFilename = "fileWithAttrOut.e";
        read_and_write(outputFilename);

        std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(get_comm());

        stk::io::StkMeshIoBroker stkIo(get_comm());
        stk::io::fill_mesh_preexisting(stkIo, outputFilename, *bulk);

        expect_ordered_attributes(stkIo);
    }
}

}
