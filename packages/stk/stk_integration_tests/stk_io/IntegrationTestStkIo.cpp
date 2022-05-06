#include <gtest/gtest.h>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include "stk_util/parallel/Parallel.hpp"
#include "stk_unit_test_utils/ioUtils.hpp"

namespace
{

TEST(StkIo, readingParallelFilesMissingParallelCommInfo)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    if(stk::parallel_machine_size(communicator) == 2)
    {
        std::shared_ptr<stk::mesh::BulkData> bulk = stk::mesh::MeshBuilder(communicator).create();
        // file from Salinas output is missing parallel info
        EXPECT_THROW(stk::io::fill_mesh("twoHexMissingParallelInfo.e", *bulk), std::runtime_error);
    }
}

}
