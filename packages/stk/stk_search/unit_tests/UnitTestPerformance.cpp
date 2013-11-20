#include <stk_search/CoarseSearch.hpp>
#include <stk_search/BoundingBox.hpp>
#include <stk_search/IdentProc.hpp>

#include <mpi.h>
#include <gtest/gtest.h>
#include <vector>

#include <unit_tests/UnitTestUtils.hpp>

namespace
{

TEST(Performance, Spheres)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int proc = 0;
    MPI_Comm_rank(comm, &proc);
    int numProcs = 1;
    MPI_Comm_size(comm, &numProcs);

    size_t numColumnsPerProcessor = 1000;
    double boxSize = 1.0;

    std::vector<My3dAxisAlignedBoundingBox> smallBoxVector;
    std::vector<My3dAxisAlignedBoundingBox> bigBoxVector;

    std::vector<std::pair<MyBoxId, MyBoxId> > boxIdPairResults;

    if(proc % 2 == 0)
    {
        double startX = numColumnsPerProcessor * boxSize * proc / 2;

        for(size_t x = 0; x < numColumnsPerProcessor; x++)
        {
            for(size_t y = 0; y < numColumnsPerProcessor; y++)
            {
                double radius = boxSize / 2;
                double centerX = startX + x * boxSize + radius;
                double centerY = y * boxSize + radius;
                double centerZ = radius;

                int id = x * numColumnsPerProcessor + y;

                smallBoxVector.push_back(generate3dBoundingBox<My3dAxisAlignedBoundingBox>(centerX,
                                                                                           centerY,
                                                                                           centerZ,
                                                                                           radius,
                                                                                           id,
                                                                                           proc));
            }
        }
    }
    else
    {
        double radius = numColumnsPerProcessor * boxSize / 2;
        double startX = numColumnsPerProcessor * boxSize * (proc - 1) / 2;
        double centerX = startX + radius;
        double centerY = radius;
        double centerZ = boxSize / 2;
        int id = 1;
        bigBoxVector.push_back(generate3dBoundingBox<My3dAxisAlignedBoundingBox>(centerX,
                                                                                 centerY,
                                                                                 centerZ,
                                                                                 radius - boxSize / 2,
                                                                                 id,
                                                                                 proc));
    }

    stk::search::coarse_search(smallBoxVector, bigBoxVector, stk::search::OCTREE, comm, boxIdPairResults);

    size_t numExpectedResults = numColumnsPerProcessor * numColumnsPerProcessor;
    bool lastProcessorWithOddNumberOfProcs = numProcs % 2 != 0 && proc == numProcs - 1;
    if(lastProcessorWithOddNumberOfProcs)
    {
        numExpectedResults = 0;
    }

    EXPECT_EQ(numExpectedResults, boxIdPairResults.size());
}

}
