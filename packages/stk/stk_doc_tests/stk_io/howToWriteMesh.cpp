#include <unistd.h>
#include <gtest/gtest.h>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_unit_test_utils/ioUtils.hpp>
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

}
