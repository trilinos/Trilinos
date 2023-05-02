// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
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
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
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

//-BEGIN
#include <stk_search/CoarseSearch.hpp>
#include <stk_search/BoundingBox.hpp>
#include <gtest/gtest.h>

namespace
{
typedef stk::search::Box<double> Box;
typedef stk::search::IdentProc<int, int> Id;
void assertPairInResults(Id a, Id b, const std::vector<std::pair<Id, Id> > &searchResults);
TEST(StkSearchHowTo, useKdtreeSearch)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int myProcId = stk::parallel_machine_rank(comm);
    std::vector<std::pair<Box, Id> > firstList, secondList;
    Box unitBox(Box::point_type(0, 0, 0), Box::point_type(1, 1, 1));
    Id firstId(0, myProcId);
    firstList.push_back(std::make_pair(unitBox, firstId));
    Id secondId(1, myProcId);
    secondList.push_back(std::make_pair(unitBox, secondId));

    std::vector<std::pair<Id, Id> > searchResults;
    stk::search::coarse_search(firstList, secondList, stk::search::KDTREE, comm, searchResults);

    int numProc = stk::parallel_machine_size(comm);
    for(int procId = 0; procId < numProc; procId++)
    {
        assertPairInResults(Id(0, myProcId), Id(1, procId), searchResults);
        assertPairInResults(Id(0, procId), Id(1, myProcId), searchResults);
    }
}
TEST(StkSearchHowTo, useSphereAndPointBoundingVolumes)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int myProcId = stk::parallel_machine_rank(comm);
    std::vector<std::pair<stk::search::Sphere<double>, Id> > firstList;
    stk::search::Point<double> center(0, 0, 0);
    const double radius = 0.5;
    stk::search::Sphere<double> unitSphere(center, radius);
    Id firstId(0, myProcId);
    firstList.push_back(std::make_pair(unitSphere, firstId));
    std::vector<std::pair<stk::search::Point<double>, Id> > secondList;
    stk::search::Point<double> point(0.1, 0.2, 0.3);
    Id secondId(1, myProcId);
    secondList.push_back(std::make_pair(point, secondId));

    std::vector<std::pair<Id, Id> > searchResults;
    stk::search::coarse_search(firstList, secondList, stk::search::KDTREE, comm, searchResults);

    int numProc = stk::parallel_machine_size(comm);
    for(int procId = 0; procId < numProc; procId++)
    {
        assertPairInResults(Id(0, myProcId), Id(1, procId), searchResults);
        assertPairInResults(Id(0, procId), Id(1, myProcId), searchResults);
    }
}
void assertPairInResults(Id a, Id b, const std::vector<std::pair<Id, Id> > &searchResults)
{
    std::pair<Id, Id> expectedIntersectionPair(a, b);
    std::vector<std::pair<Id, Id> >::const_iterator resultsIter =
            std::find(searchResults.begin(), searchResults.end(), expectedIntersectionPair);
    bool foundExpectedPairInResults = resultsIter != searchResults.end();
    ASSERT_TRUE(foundExpectedPairInResults);
}
}
//-END
