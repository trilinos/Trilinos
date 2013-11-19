#include <stk_search/CoarseSearch.hpp>
#include <stk_search/BoundingBox.hpp>
#include <stk_search/IdentProc.hpp>

#include <mpi.h>
#include <gtest/gtest.h>
#include <vector>

namespace
{

const int spatialDim = 3;
typedef stk::search::ident::IdentProc<unsigned, unsigned> MyBoxId;
typedef stk::search::box::SphereBoundingBox<MyBoxId, double, spatialDim> My3dSphereBoundingBox;

My3dSphereBoundingBox generate3dSphereBoundingBox(const double centerX,
                                                  const double centerY,
                                                  const double centerZ,
                                                  const double radius,
                                                  const int entityId,
                                                  const int procId)
{
    MyBoxId id(entityId, procId);
    double center[spatialDim];
    center[0] = centerX;
    center[1] = centerY;
    center[2] = centerZ;
    return My3dSphereBoundingBox(center, radius, id);
}

void runTwoSpheresTest(const double distanceBetweenSphereCenters, const double radius, std::vector< std::pair<MyBoxId, MyBoxId> > &boxIdPairResults)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int procId = -1;
    MPI_Comm_rank(comm, &procId);

    std::vector<My3dSphereBoundingBox> boxVector1;
    boxVector1.push_back(generate3dSphereBoundingBox(0, 0, 0, radius, 1, procId));

    std::vector<My3dSphereBoundingBox> boxVector2;
    boxVector2.push_back(generate3dSphereBoundingBox(distanceBetweenSphereCenters, 0, 0, radius, 2, procId));

    stk::search::coarse_search(boxVector1, boxVector2, stk::search::BOOST_RTREE, comm, boxIdPairResults);
}

const double radiusOfOneHalf = 0.5;

TEST(Verification, OverlappingSpheres)
{
    double distanceBetweenSphereCenters = 0.5;
    std::vector< std::pair<MyBoxId, MyBoxId> > boxIdPairResults;
    runTwoSpheresTest(distanceBetweenSphereCenters, radiusOfOneHalf, boxIdPairResults);

    ASSERT_EQ(1u, boxIdPairResults.size());
}

TEST(Verification, NonOverlappingSpheres)
{
    double distanceBetweenSphereCenters = 2.0;
    std::vector< std::pair<MyBoxId, MyBoxId> > boxIdPairResults;
    runTwoSpheresTest(distanceBetweenSphereCenters, radiusOfOneHalf, boxIdPairResults);

    ASSERT_EQ(0u, boxIdPairResults.size());
}

TEST(Verification, JustEdgeOverlappingSpheres)
{
    double distanceBetweenSphereCenters = 0.999999999;
    std::vector< std::pair<MyBoxId, MyBoxId> > boxIdPairResults;
    runTwoSpheresTest(distanceBetweenSphereCenters, radiusOfOneHalf, boxIdPairResults);

    ASSERT_EQ(1u, boxIdPairResults.size());
}

TEST(Verification, NotQuiteEdgeOverlappingSpheres)
{
    double distanceBetweenSphereCenters = 1.0000000001;
    std::vector< std::pair<MyBoxId, MyBoxId> > boxIdPairResults;
    runTwoSpheresTest(distanceBetweenSphereCenters, radiusOfOneHalf, boxIdPairResults);

    ASSERT_EQ(0u, boxIdPairResults.size());
}

void runSphereOverlappingEightSurroundingSpheres(const double radius, const unsigned numExpectedResults)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int procId = -1;
    MPI_Comm_rank(comm, &procId);

    std::vector<My3dSphereBoundingBox> boxVector1;
    boxVector1.push_back(generate3dSphereBoundingBox(0, 0, 0, radius, 1, procId));
    boxVector1.push_back(generate3dSphereBoundingBox(1, 0, 0, radius, 2, procId));
    boxVector1.push_back(generate3dSphereBoundingBox(2, 0, 0, radius, 3, procId));
    boxVector1.push_back(generate3dSphereBoundingBox(0, 1, 0, radius, 4, procId));
    //skip middle one
    boxVector1.push_back(generate3dSphereBoundingBox(2, 1, 0, radius, 6, procId));
    boxVector1.push_back(generate3dSphereBoundingBox(0, 2, 0, radius, 7, procId));
    boxVector1.push_back(generate3dSphereBoundingBox(1, 2, 0, radius, 8, procId));
    boxVector1.push_back(generate3dSphereBoundingBox(2, 2, 0, radius, 9, procId));

    std::vector<My3dSphereBoundingBox> boxVector2;
    boxVector2.push_back(generate3dSphereBoundingBox(1, 1, 0, radius, 5, procId));

    std::vector< std::pair<MyBoxId, MyBoxId> > boxIdPairResults;
    stk::search::coarse_search(boxVector1, boxVector2, stk::search::BOOST_RTREE, comm, boxIdPairResults);

    if(numExpectedResults != boxIdPairResults.size())
    {
        for(size_t i=0; i<boxIdPairResults.size(); i++)
        {
            std::cerr << boxIdPairResults[i].first << ", " << boxIdPairResults[i].second << std::endl;
        }
    }
    ASSERT_EQ(numExpectedResults, boxIdPairResults.size());
}

TEST(Verification, SphereOverlappingEightSurroundingSpheres)
{
    const double radius = 0.708;
    const unsigned numExpectedResults = 8;
    runSphereOverlappingEightSurroundingSpheres(radius, numExpectedResults);
}

TEST(Verification, SphereOverlappingFourOfEightSurroundingSpheres)
{
    const double radius = 0.706;
    const unsigned numExpectedResults = 4;
    runSphereOverlappingEightSurroundingSpheres(radius, numExpectedResults);

}

TEST(Verification, ParallelLineOfSpheres)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int procId = -1;
    MPI_Comm_rank(comm, &procId);

    const double radius = 0.708;
    const double distanceBetweenCenters = 1.0;

    const double xCoord = procId * distanceBetweenCenters;
    My3dSphereBoundingBox localBox = generate3dSphereBoundingBox(xCoord, 0, 0, radius, 1, procId);

    std::vector<My3dSphereBoundingBox> boxVector1;
    std::vector<My3dSphereBoundingBox> boxVector2;
    if(procId % 2 == 0)
    {
        boxVector1.push_back(localBox);
    }
    else
    {
        boxVector2.push_back(localBox);
    }

    std::vector< std::pair<MyBoxId, MyBoxId> > boxIdPairResults;
    stk::search::coarse_search(boxVector1, boxVector2, stk::search::BOOST_RTREE, comm, boxIdPairResults);

    int numProcs = -1;
    MPI_Comm_size(comm, &numProcs);

    unsigned numExpectedResults = 2;
    bool doOwnFirstOrLastSphereInLine = procId == 0 || procId == numProcs-1;
    if(doOwnFirstOrLastSphereInLine)
    {
        numExpectedResults = 1;
    }
    if(numProcs == 1)
    {
        numExpectedResults = 0;
    }
    ASSERT_EQ(numExpectedResults, boxIdPairResults.size()) << "on proc id " << procId;
}

}
