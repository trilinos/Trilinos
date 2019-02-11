#include <unistd.h>
#include <gtest/gtest.h>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_unit_test_utils/getOption.h>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_io/WriteMesh.hpp>
#include <stk_unit_test_utils/ioUtils.hpp>

extern int gl_argc;
extern char** gl_argv;

namespace
{

TEST(StkIoHowTo, WriteMesh)
{
    std::string filename = "output.exo";
    {
        stk::mesh::MetaData meta;
        stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
        stk::io::fill_mesh("generated:1x1x4", bulk);

        stk::io::StkMeshIoBroker stkIo;
        stkIo.set_bulk_data(bulk);
        size_t outputFileIndex = stkIo.create_output_mesh(filename, stk::io::WRITE_RESULTS);
        stkIo.write_output_mesh(outputFileIndex);
        stkIo.write_defined_output_fields(outputFileIndex);
    }

    {
        stk::mesh::MetaData meta;
        stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
        stk::io::fill_mesh(filename, bulk);

        std::vector<size_t> entityCounts;
        stk::mesh::comm_mesh_counts(bulk, entityCounts);
        EXPECT_EQ(4u, entityCounts[stk::topology::ELEM_RANK]);
    }

    unlink(filename.c_str());
}

TEST(StkIoHowTo, generateHugeMesh)
{
    std::string meshSpec = stk::unit_test_util::get_option("-i", "1x1x4");
    std::string fullMeshSpec = "generated:"+meshSpec;

    std::string filename = "output.exo";
    stk::io::StkMeshIoBroker inputBroker;
    inputBroker.property_add(Ioss::Property("INTEGER_SIZE_API" , 8));
    inputBroker.property_add(Ioss::Property("INTEGER_SIZE_DB" , 8));
    stk::mesh::MetaData meta;
    stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
    stk::io::fill_mesh_preexisting(inputBroker, fullMeshSpec, bulk);

    stk::io::write_mesh_with_large_ids_and_fields(filename, bulk);

    if (gl_argc == 0)
    {
        unlink(filename.c_str());
    }
}

}
