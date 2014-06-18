#include <gtest/gtest.h>
#include <sstream>
#include <mpi.h>

// macro for short term solution to not break builds for Trilinos.
// Dependency on unitTestUtils here is the main worry.

#if defined(STK_BUILT_IN_SIERRA)

#include <exampleMeshes/StkMeshFromGeneratedMesh.h>
#include <stk_mesh/base/Comm.hpp>

namespace
{

std::string getGeneratedMeshString(const int xdim, const int ydim, const int zdim)
{
    std::ostringstream oss;
    oss << "generated: " << xdim << "x" << ydim << "x" << zdim;
    return oss.str();
}

//DocTest1
TEST( CommMeshCounts, Serial )
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    int numprocs = -1;
    MPI_Comm_size(communicator, &numprocs);
    if ( numprocs == 1 )
    {
        const std::string generatedMeshSpec = getGeneratedMeshString(10,20,2);
        unitTestUtils::exampleMeshes::StkMeshCreator stkMesh(generatedMeshSpec, communicator);

        std::vector<size_t> comm_mesh_counts;
        stk::mesh::comm_mesh_counts(*stkMesh.getBulkData(), comm_mesh_counts);

        size_t goldNumElements = 10*20*2;
        EXPECT_EQ(goldNumElements, comm_mesh_counts[stk::topology::ELEMENT_RANK]);
    }
}
//DocTest2
TEST( CommMeshCounts, Parallel )
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    int numprocs = 1;
    MPI_Comm_size(communicator, &numprocs);

    const std::string generatedMeshSpec = getGeneratedMeshString(10,20,2*numprocs);
    unitTestUtils::exampleMeshes::StkMeshCreator stkMesh(generatedMeshSpec, communicator);

    std::vector<size_t> comm_mesh_counts;
    stk::mesh::comm_mesh_counts(*stkMesh.getBulkData(), comm_mesh_counts);

    size_t goldNumElements = 10*20*2*numprocs;
    EXPECT_EQ(goldNumElements, comm_mesh_counts[stk::topology::ELEMENT_RANK]);
}
//DocTest3
TEST( CommMeshCountsWithStats, Parallel )
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    int numprocs = 1;
    MPI_Comm_size(communicator, &numprocs);

    const std::string generatedMeshSpec = getGeneratedMeshString(10,20,2*numprocs);
    unitTestUtils::exampleMeshes::StkMeshCreator stkMesh(generatedMeshSpec, communicator);

    std::vector<size_t> comm_mesh_counts;
    std::vector<size_t> min_counts;
    std::vector<size_t> max_counts;

    stk::mesh::comm_mesh_counts(*stkMesh.getBulkData(), comm_mesh_counts, min_counts, max_counts);

    size_t goldNumElements = 10*20*2*numprocs;
    EXPECT_EQ(goldNumElements, comm_mesh_counts[stk::topology::ELEMENT_RANK]);

    size_t goldMinNumElements = 10*20*2;
    EXPECT_EQ(goldMinNumElements, min_counts[stk::topology::ELEMENT_RANK]);

    size_t goldMaxNumElements = goldMinNumElements;
    EXPECT_EQ(goldMaxNumElements, max_counts[stk::topology::ELEMENT_RANK]);
}
//EndDocTest
}
#endif
