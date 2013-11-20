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
typedef stk::search::box::PointBoundingBox<MyBoxId, double, spatialDim> My3dPointBoundingBox;
typedef stk::search::box::AxisAlignedBoundingBox<MyBoxId, double, spatialDim> My3dAxisAlignedBoundingBox;

template<class BoundingBoxType>
BoundingBoxType generate3dBoundingBox(const double centerX,
                                                  const double centerY,
                                                  const double centerZ,
                                                  const double radius,
                                                  const int entityId,
                                                  const int procId);

template<>
My3dPointBoundingBox generate3dBoundingBox<My3dPointBoundingBox>(const double centerX,
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
    return My3dPointBoundingBox(center, id);
}

template<>
My3dSphereBoundingBox generate3dBoundingBox<My3dSphereBoundingBox>(const double centerX,
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

//       ------------
//      |            |
//      |      radius|
//      |      ------|
//      |            |
//      |            |
//       ------------
// width = 2*radius
template<>
My3dAxisAlignedBoundingBox generate3dBoundingBox<My3dAxisAlignedBoundingBox>(const double centerX,
                                                  const double centerY,
                                                  const double centerZ,
                                                  const double radius,
                                                  const int entityId,
                                                  const int procId)
{
    MyBoxId id(entityId, procId);
    double corners[2*spatialDim];
    corners[0] = centerX-radius;
    corners[1] = centerY-radius;
    corners[2] = centerZ-radius;
    corners[3] = centerX+radius;
    corners[4] = centerY+radius;
    corners[5] = centerZ+radius;
    return My3dAxisAlignedBoundingBox(corners, id);
}

void runTwoSpheresTest(stk::search::SearchMethod searchMethod, const double distanceBetweenSphereCenters, const double radius, std::vector< std::pair<MyBoxId, MyBoxId> > &boxIdPairResults)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int procId = -1;
    MPI_Comm_rank(comm, &procId);

    std::vector<My3dSphereBoundingBox> boxVector1;
    boxVector1.push_back(generate3dBoundingBox<My3dSphereBoundingBox>(0, 0, 0, radius, 1, procId));

    std::vector<My3dSphereBoundingBox> boxVector2;
    boxVector2.push_back(generate3dBoundingBox<My3dSphereBoundingBox>(distanceBetweenSphereCenters, 0, 0, radius, 2, procId));

    stk::search::coarse_search(boxVector1, boxVector2, searchMethod, comm, boxIdPairResults);
}

const double radiusOfOneHalf = 0.5;

int numProcessors()
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int numProcs = 1;
    MPI_Comm_size(comm, &numProcs);
    return numProcs;
}

TEST(Verification, OverlappingSpheres_BOOST_RTREE)
{
    if (numProcessors() > 1) {
        return;
    }

    double distanceBetweenSphereCenters = 0.5;
    std::vector< std::pair<MyBoxId, MyBoxId> > boxIdPairResults;
    runTwoSpheresTest(stk::search::BOOST_RTREE, distanceBetweenSphereCenters, radiusOfOneHalf, boxIdPairResults);

    ASSERT_EQ(1u, boxIdPairResults.size());
}

TEST(Verification, OverlappingSpheres_OCTREE)
{
    if (numProcessors() > 1) {
        return;
    }

    double distanceBetweenSphereCenters = 0.5;
    std::vector< std::pair<MyBoxId, MyBoxId> > boxIdPairResults;
    runTwoSpheresTest(stk::search::OCTREE, distanceBetweenSphereCenters, radiusOfOneHalf, boxIdPairResults);

    ASSERT_EQ(1u, boxIdPairResults.size());
}

TEST(Verification, NonOverlappingSpheres_BOOST_RTREE)
{
    double distanceBetweenSphereCenters = 2.0;
    std::vector< std::pair<MyBoxId, MyBoxId> > boxIdPairResults;
    runTwoSpheresTest(stk::search::BOOST_RTREE, distanceBetweenSphereCenters, radiusOfOneHalf, boxIdPairResults);

    ASSERT_EQ(0u, boxIdPairResults.size());
}

TEST(Verification, NonOverlappingSpheres_OCTREE)
{
    double distanceBetweenSphereCenters = 2.0;
    std::vector< std::pair<MyBoxId, MyBoxId> > boxIdPairResults;
    runTwoSpheresTest(stk::search::OCTREE, distanceBetweenSphereCenters, radiusOfOneHalf, boxIdPairResults);

    ASSERT_EQ(0u, boxIdPairResults.size());
}

TEST(Verification, JustEdgeOverlappingSpheres_BOOST_RTREE)
{
    if (numProcessors() > 1) {
        return;
    }

    double distanceBetweenSphereCenters = 0.999999999;
    std::vector< std::pair<MyBoxId, MyBoxId> > boxIdPairResults;
    runTwoSpheresTest(stk::search::BOOST_RTREE, distanceBetweenSphereCenters, radiusOfOneHalf, boxIdPairResults);

    ASSERT_EQ(1u, boxIdPairResults.size());
}

TEST(Verification, JustEdgeOverlappingSpheres_OCTREE)
{
    if (numProcessors() > 1) {
        return;
    }

    double distanceBetweenSphereCenters = 0.999999999;
    std::vector< std::pair<MyBoxId, MyBoxId> > boxIdPairResults;
    runTwoSpheresTest(stk::search::OCTREE, distanceBetweenSphereCenters, radiusOfOneHalf, boxIdPairResults);

    ASSERT_EQ(1u, boxIdPairResults.size());
}

TEST(Verification, NotQuiteEdgeOverlappingSpheres_BOOST_RTREE)
{
    double distanceBetweenSphereCenters = 1.0000000001;
    std::vector< std::pair<MyBoxId, MyBoxId> > boxIdPairResults;
    runTwoSpheresTest(stk::search::BOOST_RTREE, distanceBetweenSphereCenters, radiusOfOneHalf, boxIdPairResults);

    ASSERT_EQ(0u, boxIdPairResults.size());
}

TEST(Verification, NotQuiteEdgeOverlappingSpheres_OCTREE)
{
    double distanceBetweenSphereCenters = 1.0000000001;
    std::vector< std::pair<MyBoxId, MyBoxId> > boxIdPairResults;
    runTwoSpheresTest(stk::search::OCTREE, distanceBetweenSphereCenters, radiusOfOneHalf, boxIdPairResults);

    ASSERT_EQ(0u, boxIdPairResults.size());
}

template<class InnerBoundingBoxType, class OuterBoundingBoxType>
void runBoxOverlappingEightSurroundingBoxes(stk::search::SearchMethod searchMethod, const double radius, const unsigned numExpectedResults)
{
    int numProc = numProcessors();
    MPI_Comm comm = MPI_COMM_WORLD;
    int procId = -1;
    MPI_Comm_rank(comm, &procId);

    std::vector<OuterBoundingBoxType> boxVector1;
    if(procId == 0)
    {
        boxVector1.push_back(generate3dBoundingBox<OuterBoundingBoxType>(0, 0, 0, radius, 1, procId));
        boxVector1.push_back(generate3dBoundingBox<OuterBoundingBoxType>(1, 0, 0, radius, 2, procId));
        boxVector1.push_back(generate3dBoundingBox<OuterBoundingBoxType>(2, 0, 0, radius, 3, procId));
        boxVector1.push_back(generate3dBoundingBox<OuterBoundingBoxType>(0, 1, 0, radius, 4, procId));
        //skip middle one
        boxVector1.push_back(generate3dBoundingBox<OuterBoundingBoxType>(2, 1, 0, radius, 6, procId));
        boxVector1.push_back(generate3dBoundingBox<OuterBoundingBoxType>(0, 2, 0, radius, 7, procId));
        boxVector1.push_back(generate3dBoundingBox<OuterBoundingBoxType>(1, 2, 0, radius, 8, procId));
        boxVector1.push_back(generate3dBoundingBox<OuterBoundingBoxType>(2, 2, 0, radius, 9, procId));
    }

    std::vector<InnerBoundingBoxType> boxVector2;
    if(procId == numProc-1)
    {
        boxVector2.push_back(generate3dBoundingBox<InnerBoundingBoxType>(1, 1, 0, radius, 5, procId));
    }

    std::vector< std::pair<MyBoxId, MyBoxId> > boxIdPairResults;
    stk::search::coarse_search(boxVector1, boxVector2, searchMethod, comm, boxIdPairResults);

    if(!boxVector1.empty() || !boxVector2.empty())
    {
        if(numExpectedResults != boxIdPairResults.size())
        {
            for(size_t i=0; i<boxIdPairResults.size(); i++)
            {
                std::cerr << boxIdPairResults[i].first << ", " << boxIdPairResults[i].second << std::endl;
            }
        }
        ASSERT_EQ(numExpectedResults, boxIdPairResults.size());
    }
}

TEST(Verification, SphereOverlappingEightSurroundingSpheres_BOOST_RTREE)
{
    const double radius = 0.708;
    const unsigned numExpectedResults = 8;
    runBoxOverlappingEightSurroundingBoxes<My3dSphereBoundingBox,My3dSphereBoundingBox>(stk::search::BOOST_RTREE, radius, numExpectedResults);
}

TEST(Verification, SphereOverlappingEightSurroundingSpheres_OCTREE)
{
    const double radius = 0.708;
    const unsigned numExpectedResults = 8;
    runBoxOverlappingEightSurroundingBoxes<My3dSphereBoundingBox,My3dSphereBoundingBox>(stk::search::OCTREE, radius, numExpectedResults);
}

TEST(Verification, SphereOverlappingNoSurroundingPoints_BOOST_RTREE)
{
    const double radius = 0.99;
    const unsigned numExpectedResults = 0;
    runBoxOverlappingEightSurroundingBoxes<My3dSphereBoundingBox,My3dPointBoundingBox>(stk::search::BOOST_RTREE, radius, numExpectedResults);
}

TEST(Verification, SphereOverlappingFourSurroundingPoints_BOOST_RTREE)
{
    const double radius = 1.41;
    const unsigned numExpectedResults = 4;
    runBoxOverlappingEightSurroundingBoxes<My3dSphereBoundingBox,My3dPointBoundingBox>(stk::search::BOOST_RTREE, radius, numExpectedResults);
}

TEST(Verification, SphereOverlappingEightSurroundingPoints_BOOST_RTREE)
{
    const double radius = 1.42;
    const unsigned numExpectedResults = 8;
    runBoxOverlappingEightSurroundingBoxes<My3dSphereBoundingBox,My3dPointBoundingBox>(stk::search::BOOST_RTREE, radius, numExpectedResults);
}

TEST(Verification, SphereOverlappingNoSurroundingPoints_OCTREE)
{
    const double radius = 0.99;
    const unsigned numExpectedResults = 0;
    runBoxOverlappingEightSurroundingBoxes<My3dSphereBoundingBox,My3dPointBoundingBox>(stk::search::OCTREE, radius, numExpectedResults);
}

TEST(Verification, SphereOverlappingFourSurroundingPoints_OCTREE)
{
    const double radius = 1.41;
    const unsigned numExpectedResults = 4;
    runBoxOverlappingEightSurroundingBoxes<My3dSphereBoundingBox,My3dPointBoundingBox>(stk::search::OCTREE, radius, numExpectedResults);
}

TEST(Verification, SphereOverlappingEightSurroundingPoints_OCTREE)
{
    const double radius = 1.42;
    const unsigned numExpectedResults = 8;
    runBoxOverlappingEightSurroundingBoxes<My3dSphereBoundingBox,My3dPointBoundingBox>(stk::search::OCTREE, radius, numExpectedResults);
}

TEST(Verification, SphereOverlappingFourOfEightSurroundingSpheres_BOOST_RTREE)
{
    const double radius = 0.706;
    const unsigned numExpectedResults = 4;
    runBoxOverlappingEightSurroundingBoxes<My3dSphereBoundingBox,My3dSphereBoundingBox>(stk::search::BOOST_RTREE, radius, numExpectedResults);
}

TEST(Verification, SphereOverlappingFourOfEightSurroundingSpheres_OCTREE)
{
    const double radius = 0.706;
    const unsigned numExpectedResults = 4;
    runBoxOverlappingEightSurroundingBoxes<My3dSphereBoundingBox,My3dSphereBoundingBox>(stk::search::OCTREE, radius, numExpectedResults);
}

TEST(Verification, BoxOverlappingNoSurroundingPoints_BOOST_RTREE)
{
    const double radius = 0.99;
    const unsigned numExpectedResults = 0;
    runBoxOverlappingEightSurroundingBoxes<My3dAxisAlignedBoundingBox,My3dPointBoundingBox>(stk::search::BOOST_RTREE, radius, numExpectedResults);
}

TEST(Verification, BoxOverlappingNoSurroundingPoints_OCTREE)
{
    const double radius = 0.99;
    const unsigned numExpectedResults = 0;
    runBoxOverlappingEightSurroundingBoxes<My3dAxisAlignedBoundingBox,My3dPointBoundingBox>(stk::search::OCTREE, radius, numExpectedResults);
}

TEST(Verification, BoxOverlappingEightSurroundingPoints_BOOST_RTREE)
{
    const double radius = 1.01;
    const unsigned numExpectedResults = 8;
    runBoxOverlappingEightSurroundingBoxes<My3dAxisAlignedBoundingBox,My3dPointBoundingBox>(stk::search::BOOST_RTREE, radius, numExpectedResults);
}

TEST(Verification, BoxOverlappingEightSurroundingPoints_OCTREE)
{
    const double radius = 1.01;
    const unsigned numExpectedResults = 8;
    runBoxOverlappingEightSurroundingBoxes<My3dAxisAlignedBoundingBox,My3dPointBoundingBox>(stk::search::OCTREE, radius, numExpectedResults);
}

TEST(Verification, PointOverlappingNoSurroundingBoxes_BOOST_RTREE)
{
    const double radius = 0.99;
    const unsigned numExpectedResults = 0;
    runBoxOverlappingEightSurroundingBoxes<My3dPointBoundingBox,My3dAxisAlignedBoundingBox>(stk::search::BOOST_RTREE, radius, numExpectedResults);
}

TEST(Verification, PointOverlappingEightSurroundingBoxes_BOOST_RTREE)
{
    const double radius = 1.01;
    const unsigned numExpectedResults = 8;
    runBoxOverlappingEightSurroundingBoxes<My3dPointBoundingBox,My3dAxisAlignedBoundingBox>(stk::search::BOOST_RTREE, radius, numExpectedResults);
}

TEST(Verification, PointOverlappingNoSurroundingBoxes_OCTREE)
{
    const double radius = 0.99;
    const unsigned numExpectedResults = 0;
    runBoxOverlappingEightSurroundingBoxes<My3dPointBoundingBox,My3dAxisAlignedBoundingBox>(stk::search::OCTREE, radius, numExpectedResults);
}

TEST(Verification, PointOverlappingEightSurroundingBoxes_OCTREE)
{
    const double radius = 1.01;
    const unsigned numExpectedResults = 8;
    runBoxOverlappingEightSurroundingBoxes<My3dPointBoundingBox,My3dAxisAlignedBoundingBox>(stk::search::OCTREE, radius, numExpectedResults);
}

enum Axis
{
    xDim, yDim, zDim
};
template<class BoundingBoxType>
void runLineOfBoundingBoxes(stk::search::SearchMethod searchMethod, enum Axis axis)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int procId = -1;
    MPI_Comm_rank(comm, &procId);

    const double radius = 0.708;
    const double distanceBetweenCenters = 1.0;

    const double paramCoord = procId * distanceBetweenCenters;

    std::vector<BoundingBoxType> boxVector1;
    std::vector<BoundingBoxType> boxVector2;
    if(procId % 2 == 0)
    {
        switch(axis)
        {
            case xDim:
                boxVector1.push_back(generate3dBoundingBox<BoundingBoxType>(paramCoord, 0, 0, radius, 1, procId));
                break;
            case yDim:
                boxVector1.push_back(generate3dBoundingBox<BoundingBoxType>(0, paramCoord, 0, radius, 1, procId));
                break;
            case zDim:
                boxVector1.push_back(generate3dBoundingBox<BoundingBoxType>(0, 0, paramCoord, radius, 1, procId));
                break;
        }
    }
    else
    {
        switch(axis)
        {
            case xDim:
                boxVector2.push_back(generate3dBoundingBox<BoundingBoxType>(paramCoord, 0, 0, radius, 1, procId));
                break;
            case yDim:
                boxVector2.push_back(generate3dBoundingBox<BoundingBoxType>(0, paramCoord, 0, radius, 1, procId));
                break;
            case zDim:
                boxVector2.push_back(generate3dBoundingBox<BoundingBoxType>(0, 0, paramCoord, radius, 1, procId));
                break;
        }
    }

    std::vector< std::pair<MyBoxId, MyBoxId> > boxIdPairResults;
    stk::search::coarse_search(boxVector1, boxVector2, searchMethod, comm, boxIdPairResults);

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

TEST(Verification, LineOfSpheres_BOOST_RTREE)
{
    runLineOfBoundingBoxes<My3dSphereBoundingBox>(stk::search::BOOST_RTREE, xDim);
}

TEST(Verification, LineOfSpheres_OCTREE)
{
    runLineOfBoundingBoxes<My3dSphereBoundingBox>(stk::search::OCTREE, xDim);
}

TEST(Verification, LineOfBoxes_BOOST_RTREE)
{
    runLineOfBoundingBoxes<My3dAxisAlignedBoundingBox>(stk::search::BOOST_RTREE, yDim);
}

TEST(Verification, LineOfBoxes_OCTREE)
{
    runLineOfBoundingBoxes<My3dAxisAlignedBoundingBox>(stk::search::OCTREE, yDim);
}

TEST(Verification, LineOfSpheresZDimension_BOOST_RTREE)
{
    runLineOfBoundingBoxes<My3dSphereBoundingBox>(stk::search::BOOST_RTREE, zDim);
}

TEST(Verification, LineOfSpheresZDimension_OCTREE)
{
    runLineOfBoundingBoxes<My3dSphereBoundingBox>(stk::search::OCTREE, zDim);
}

}
