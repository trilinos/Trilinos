// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#include <stk_search/BoundingBox.hpp>
#include <stk_search/IdentProc.hpp>
#include <unit_tests/UnitTestUtils.hpp>
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine, etc

#include <gtest/gtest.h>
#include <vector>

namespace
{


void runTwoSpheresTest(stk::search::SearchMethod searchMethod, const double distanceBetweenSphereCenters, const double radius, std::vector< std::pair<Ident, Ident> > &boxIdPairResults)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int procId = stk::parallel_machine_rank(comm);

    std::vector< std::pair<Sphere,Ident> > boxVector1;
    boxVector1.push_back(generateBoundingVolume<Sphere>(0, 0, 0, radius, 1, procId));

    std::vector< std::pair<Sphere,Ident> > boxVector2;
    boxVector2.push_back(generateBoundingVolume<Sphere>(distanceBetweenSphereCenters, 0, 0, radius, 2, procId));

    stk::search::coarse_search(boxVector1, boxVector2, searchMethod, comm, boxIdPairResults);
}

const double radiusOfOneHalf = 0.5;

int numProcessors()
{
    MPI_Comm comm = MPI_COMM_WORLD;
    return stk::parallel_machine_size(comm);
}

TEST(Verification, OverlappingSpheres_BOOST_RTREE)
{
    if (numProcessors() > 1) {
        return;
    }

    double distanceBetweenSphereCenters = 0.5;
    std::vector< std::pair<Ident, Ident> > boxIdPairResults;
    runTwoSpheresTest(stk::search::BOOST_RTREE, distanceBetweenSphereCenters, radiusOfOneHalf, boxIdPairResults);

    ASSERT_EQ(1u, boxIdPairResults.size());
}

TEST(Verification, OverlappingSpheres_OCTREE)
{
    if (numProcessors() > 1) {
        return;
    }

    double distanceBetweenSphereCenters = 0.5;
    std::vector< std::pair<Ident, Ident> > boxIdPairResults;
    runTwoSpheresTest(stk::search::OCTREE, distanceBetweenSphereCenters, radiusOfOneHalf, boxIdPairResults);

    ASSERT_EQ(1u, boxIdPairResults.size());
}

TEST(Verification, NonOverlappingSpheres_BOOST_RTREE)
{
    double distanceBetweenSphereCenters = 2.0;
    std::vector< std::pair<Ident, Ident> > boxIdPairResults;
    runTwoSpheresTest(stk::search::BOOST_RTREE, distanceBetweenSphereCenters, radiusOfOneHalf, boxIdPairResults);

    ASSERT_EQ(0u, boxIdPairResults.size());
}

TEST(Verification, NonOverlappingSpheres_OCTREE)
{
    double distanceBetweenSphereCenters = 2.0;
    std::vector< std::pair<Ident, Ident> > boxIdPairResults;
    runTwoSpheresTest(stk::search::OCTREE, distanceBetweenSphereCenters, radiusOfOneHalf, boxIdPairResults);

    ASSERT_EQ(0u, boxIdPairResults.size());
}

TEST(Verification, JustEdgeOverlappingSpheres_BOOST_RTREE)
{
    if (numProcessors() > 1) {
        return;
    }

    double distanceBetweenSphereCenters = 0.999999999;
    std::vector< std::pair<Ident, Ident> > boxIdPairResults;
    runTwoSpheresTest(stk::search::BOOST_RTREE, distanceBetweenSphereCenters, radiusOfOneHalf, boxIdPairResults);

    ASSERT_EQ(1u, boxIdPairResults.size());
}

TEST(Verification, JustEdgeOverlappingSpheres_OCTREE)
{
    if (numProcessors() > 1) {
        return;
    }

    double distanceBetweenSphereCenters = 0.999999999;
    std::vector< std::pair<Ident, Ident> > boxIdPairResults;
    runTwoSpheresTest(stk::search::OCTREE, distanceBetweenSphereCenters, radiusOfOneHalf, boxIdPairResults);

    ASSERT_EQ(1u, boxIdPairResults.size());
}

TEST(Verification, NotQuiteEdgeOverlappingSpheres_BOOST_RTREE)
{
    double distanceBetweenSphereCenters = 1.0000000001;
    std::vector< std::pair<Ident, Ident> > boxIdPairResults;
    runTwoSpheresTest(stk::search::BOOST_RTREE, distanceBetweenSphereCenters, radiusOfOneHalf, boxIdPairResults);

    ASSERT_EQ(0u, boxIdPairResults.size());
}

TEST(Verification, NotQuiteEdgeOverlappingSpheres_OCTREE)
{
    double distanceBetweenSphereCenters = 1.0000000001;
    std::vector< std::pair<Ident, Ident> > boxIdPairResults;
    runTwoSpheresTest(stk::search::OCTREE, distanceBetweenSphereCenters, radiusOfOneHalf, boxIdPairResults);

    ASSERT_EQ(0u, boxIdPairResults.size());
}

template<class InnerBoundingBoxType, class OuterBoundingBoxType>
void runBoxOverlappingEightSurroundingBoxes(stk::search::SearchMethod searchMethod, const double radius, const unsigned numExpectedResults)
{
    int numProc = numProcessors();
    MPI_Comm comm = MPI_COMM_WORLD;
    int procId = stk::parallel_machine_rank(comm);

    std::vector< std::pair<OuterBoundingBoxType,Ident> > boxVector1;
    if(procId == 0)
    {
        boxVector1.push_back(generateBoundingVolume<OuterBoundingBoxType>(0, 0, 0, radius, 1, procId));
        boxVector1.push_back(generateBoundingVolume<OuterBoundingBoxType>(1, 0, 0, radius, 2, procId));
        boxVector1.push_back(generateBoundingVolume<OuterBoundingBoxType>(2, 0, 0, radius, 3, procId));
        boxVector1.push_back(generateBoundingVolume<OuterBoundingBoxType>(0, 1, 0, radius, 4, procId));
        //skip middle one
        boxVector1.push_back(generateBoundingVolume<OuterBoundingBoxType>(2, 1, 0, radius, 6, procId));
        boxVector1.push_back(generateBoundingVolume<OuterBoundingBoxType>(0, 2, 0, radius, 7, procId));
        boxVector1.push_back(generateBoundingVolume<OuterBoundingBoxType>(1, 2, 0, radius, 8, procId));
        boxVector1.push_back(generateBoundingVolume<OuterBoundingBoxType>(2, 2, 0, radius, 9, procId));
    }

    std::vector< std::pair<InnerBoundingBoxType,Ident> > boxVector2;
    if(procId == numProc-1)
    {
        boxVector2.push_back(generateBoundingVolume<InnerBoundingBoxType>(1, 1, 0, radius, 5, procId));
    }

    std::vector< std::pair<Ident, Ident> > boxIdPairResults;
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
    runBoxOverlappingEightSurroundingBoxes<Sphere,Sphere>(stk::search::BOOST_RTREE, radius, numExpectedResults);
}

TEST(Verification, SphereOverlappingEightSurroundingSpheres_OCTREE)
{
    const double radius = 0.708;
    const unsigned numExpectedResults = 8;
    runBoxOverlappingEightSurroundingBoxes<Sphere,Sphere>(stk::search::OCTREE, radius, numExpectedResults);
}

TEST(Verification, SphereOverlappingNoSurroundingPoints_BOOST_RTREE)
{
    const double radius = 0.99;
    const unsigned numExpectedResults = 0;
    runBoxOverlappingEightSurroundingBoxes<Sphere,Point>(stk::search::BOOST_RTREE, radius, numExpectedResults);
}

TEST(Verification, SphereOverlappingFourSurroundingPoints_BOOST_RTREE)
{
    const double radius = 1.41;
    const unsigned numExpectedResults = 4;
    runBoxOverlappingEightSurroundingBoxes<Sphere,Point>(stk::search::BOOST_RTREE, radius, numExpectedResults);
}

TEST(Verification, SphereOverlappingEightSurroundingPoints_BOOST_RTREE)
{
    const double radius = 1.42;
    const unsigned numExpectedResults = 8;
    runBoxOverlappingEightSurroundingBoxes<Sphere,Point>(stk::search::BOOST_RTREE, radius, numExpectedResults);
}

TEST(Verification, SphereOverlappingNoSurroundingPoints_OCTREE)
{
    const double radius = 0.99;
    const unsigned numExpectedResults = 0;
    runBoxOverlappingEightSurroundingBoxes<Sphere,Point>(stk::search::OCTREE, radius, numExpectedResults);
}

TEST(Verification, SphereOverlappingFourSurroundingPoints_OCTREE)
{
    const double radius = 1.41;
    const unsigned numExpectedResults = 4;
    runBoxOverlappingEightSurroundingBoxes<Sphere,Point>(stk::search::OCTREE, radius, numExpectedResults);
}

TEST(Verification, SphereOverlappingEightSurroundingPoints_OCTREE)
{
    const double radius = 1.42;
    const unsigned numExpectedResults = 8;
    runBoxOverlappingEightSurroundingBoxes<Sphere,Point>(stk::search::OCTREE, radius, numExpectedResults);
}

TEST(Verification, SphereOverlappingFourOfEightSurroundingSpheres_BOOST_RTREE)
{
    const double radius = 0.706;
    const unsigned numExpectedResults = 4;
    runBoxOverlappingEightSurroundingBoxes<Sphere,Sphere>(stk::search::BOOST_RTREE, radius, numExpectedResults);
}

TEST(Verification, SphereOverlappingFourOfEightSurroundingSpheres_OCTREE)
{
    const double radius = 0.706;
    const unsigned numExpectedResults = 4;
    runBoxOverlappingEightSurroundingBoxes<Sphere,Sphere>(stk::search::OCTREE, radius, numExpectedResults);
}

TEST(Verification, BoxOverlappingNoSurroundingPoints_BOOST_RTREE)
{
    const double radius = 0.99;
    const unsigned numExpectedResults = 0;
    runBoxOverlappingEightSurroundingBoxes<StkBox,Point>(stk::search::BOOST_RTREE, radius, numExpectedResults);
}

TEST(Verification, BoxOverlappingNoSurroundingPoints_OCTREE)
{
    const double radius = 0.99;
    const unsigned numExpectedResults = 0;
    runBoxOverlappingEightSurroundingBoxes<StkBox,Point>(stk::search::OCTREE, radius, numExpectedResults);
}

TEST(Verification, BoxOverlappingEightSurroundingPoints_BOOST_RTREE)
{
    const double radius = 1.01;
    const unsigned numExpectedResults = 8;
    runBoxOverlappingEightSurroundingBoxes<StkBox,Point>(stk::search::BOOST_RTREE, radius, numExpectedResults);
}

TEST(Verification, BoxOverlappingEightSurroundingPoints_OCTREE)
{
    const double radius = 1.01;
    const unsigned numExpectedResults = 8;
    runBoxOverlappingEightSurroundingBoxes<StkBox,Point>(stk::search::OCTREE, radius, numExpectedResults);
}

TEST(Verification, PointOverlappingNoSurroundingBoxes_BOOST_RTREE)
{
    const double radius = 0.99;
    const unsigned numExpectedResults = 0;
    runBoxOverlappingEightSurroundingBoxes<Point,StkBox>(stk::search::BOOST_RTREE, radius, numExpectedResults);
}

TEST(Verification, PointOverlappingEightSurroundingBoxes_BOOST_RTREE)
{
    const double radius = 1.01;
    const unsigned numExpectedResults = 8;
    runBoxOverlappingEightSurroundingBoxes<Point,StkBox>(stk::search::BOOST_RTREE, radius, numExpectedResults);
}

TEST(Verification, PointOverlappingNoSurroundingBoxes_OCTREE)
{
    const double radius = 0.99;
    const unsigned numExpectedResults = 0;
    runBoxOverlappingEightSurroundingBoxes<Point,StkBox>(stk::search::OCTREE, radius, numExpectedResults);
}

TEST(Verification, PointOverlappingEightSurroundingBoxes_OCTREE)
{
    const double radius = 1.01;
    const unsigned numExpectedResults = 8;
    runBoxOverlappingEightSurroundingBoxes<Point,StkBox>(stk::search::OCTREE, radius, numExpectedResults);
}

enum Axis
{
    xDim, yDim, zDim
};
template<class BoundingBoxType>
void runLineOfBoundingBoxes(stk::search::SearchMethod searchMethod, enum Axis axis)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int procId = stk::parallel_machine_rank(comm);

    const double radius = 0.708;
    const double distanceBetweenCenters = 1.0;

    const double paramCoord = procId * distanceBetweenCenters;

    std::vector< std::pair<BoundingBoxType, Ident> > boxVector1;
    std::vector< std::pair<BoundingBoxType, Ident> > boxVector2;
    if(procId % 2 == 0)
    {
        switch(axis)
        {
            case xDim:
                boxVector1.push_back(generateBoundingVolume<BoundingBoxType>(paramCoord, 0, 0, radius, 1, procId));
                break;
            case yDim:
                boxVector1.push_back(generateBoundingVolume<BoundingBoxType>(0, paramCoord, 0, radius, 1, procId));
                break;
            case zDim:
                boxVector1.push_back(generateBoundingVolume<BoundingBoxType>(0, 0, paramCoord, radius, 1, procId));
                break;
        }
    }
    else
    {
        switch(axis)
        {
            case xDim:
                boxVector2.push_back(generateBoundingVolume<BoundingBoxType>(paramCoord, 0, 0, radius, 1, procId));
                break;
            case yDim:
                boxVector2.push_back(generateBoundingVolume<BoundingBoxType>(0, paramCoord, 0, radius, 1, procId));
                break;
            case zDim:
                boxVector2.push_back(generateBoundingVolume<BoundingBoxType>(0, 0, paramCoord, radius, 1, procId));
                break;
        }
    }

    std::vector< std::pair<Ident, Ident> > boxIdPairResults;
    stk::search::coarse_search(boxVector1, boxVector2, searchMethod, comm, boxIdPairResults);

    int numProcs = stk::parallel_machine_size(comm);

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
    runLineOfBoundingBoxes<Sphere>(stk::search::BOOST_RTREE, xDim);
}

TEST(Verification, LineOfSpheres_OCTREE)
{
    runLineOfBoundingBoxes<Sphere>(stk::search::OCTREE, xDim);
}

TEST(Verification, LineOfBoxes_BOOST_RTREE)
{
    runLineOfBoundingBoxes<StkBox>(stk::search::BOOST_RTREE, yDim);
}

TEST(Verification, LineOfBoxes_OCTREE)
{
    runLineOfBoundingBoxes<StkBox>(stk::search::OCTREE, yDim);
}

TEST(Verification, LineOfSpheresZDimension_BOOST_RTREE)
{
    runLineOfBoundingBoxes<Sphere>(stk::search::BOOST_RTREE, zDim);
}

TEST(Verification, LineOfSpheresZDimension_OCTREE)
{
    runLineOfBoundingBoxes<Sphere>(stk::search::OCTREE, zDim);
}

}
