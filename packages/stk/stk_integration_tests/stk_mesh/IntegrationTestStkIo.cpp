#include <gtest/gtest.h>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include "stk_util/parallel/Parallel.hpp"
#include "stk_unit_test_utils/ioUtils.hpp"

namespace
{

TEST(StkIo, readingParallelFilesMissingParallelCommInfo)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    if(stk::parallel_machine_size(communicator) == 2)
    {
        stk::mesh::MetaData meta;
        stk::mesh::BulkData bulk(meta, communicator);
        // file from Salinas output is missing parallel info
        EXPECT_THROW(stk::unit_test_util::fill_mesh_using_stk_io("twoHexMissingParallelInfo.e", bulk, communicator), std::runtime_error);
    }
}

}
