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

#include <gtest/gtest.h>
#include <stk_unit_test_utils/Search_UnitTestUtils.hpp>
#include <stk_unit_test_utils/MeshUtilsForBoundingVolumes.hpp>
#include <stk_unit_test_utils/getOption.h>

namespace
{

void printGoldResults(const FloatBoxIdentProcVector &domainBoxes, const std::vector< std::pair<Sphere, IdentProc> > &spheres)
{
    SearchResults boxIdPairResults;
    for (size_t i=0;i<domainBoxes.size();++i)
    {
        for (size_t j=0;j<spheres.size();++j)
        {
            if ( stk::search::intersects(domainBoxes[i].first, spheres[j].first) )
            {
                boxIdPairResults.push_back(std::make_pair(domainBoxes[i].second, spheres[j].second));
            }
        }
    }
    std::cerr << "Gold: Found " << boxIdPairResults.size() << " interactions.\n";
}

void printGoldResults(const FloatBoxIdentProcVector &domainBoxes, const FloatBoxIdentProcVector &spheres)
{
    SearchResults boxIdPairResults;
    for (size_t i=0;i<domainBoxes.size();++i)
    {
        for (size_t j=0;j<spheres.size();++j)
        {
          if ( stk::search::intersects(domainBoxes[i].first, spheres[j].first) )
            {
                boxIdPairResults.push_back(std::make_pair(domainBoxes[i].second, spheres[j].second));
            }
        }
    }
    std::cerr << "Gold: Found " << boxIdPairResults.size() << " interactions.\n";
}

struct Options
{
     std::string mSphereFile;
     std::string mVolumeFile;
     bool mCommunicateRangeBoxes;
     stk::search::SearchMethod mSearchMethod;
     bool mSpheresFirstThenBoxes;
     bool mTestToGetGoldResults;

     void checkForRequiredFile(const std::string &option, const std::string &file)
     {
         STK_ThrowRequireMsg(file != "NO_FILE_SPECIFIED", option << " required for this unit test.");
     }

     void setSphereFile()
     {
         std::string optionString = "-sphere";
         mSphereFile = stk::unit_test_util::get_option(optionString, "NO_FILE_SPECIFIED");
         checkForRequiredFile(optionString, mSphereFile);
     }

     void setVolumeFile()
     {
         std::string optionString = "-volume";
         mVolumeFile = stk::unit_test_util::get_option(optionString, "NO_FILE_SPECIFIED");
         checkForRequiredFile(optionString, mVolumeFile);
     }

     void setSearchMethod()
     {
         std::string optionString = "-method";
         std::string searchString = stk::unit_test_util::get_option(optionString, "gtk");
         if ( searchString == "kdtree" || searchString == "KDTREE" || searchString == "gtk") {
           mSearchMethod = stk::search::KDTREE;
           return;
         }
         if ( searchString == "MORTON_LBVH") {
           mSearchMethod = stk::search::MORTON_LBVH;
           return;
         }
         if ( searchString == "ARBORX") {
           mSearchMethod = stk::search::ARBORX;
           return;
         }

         STK_ThrowErrorMsg("unrecognized search method: "<<searchString);
     }

     void setRangeBoxCommunication()
     {
         std::string optionString = "-rangeBoxComm";
         mCommunicateRangeBoxes = true;
         if ( stk::unit_test_util::get_option(optionString, "yes") == "no" )
         {
             mCommunicateRangeBoxes = false;
         }
     }

     void setSphereBoxes()
     {
         std::string optionString = "-sb";
         mSpheresFirstThenBoxes = false;
         if ( stk::unit_test_util::get_option(optionString, "no" ) == "yes" )
         {
             mSpheresFirstThenBoxes = true;
         }
     }

     void setIfTestIsGoldTestRun()
     {
         mTestToGetGoldResults = false;
         std::string optionString = "-getGold";
         if ( stk::unit_test_util::get_option(optionString, "no") == "yes" )
         {
             mTestToGetGoldResults = true;
         }
     }

};

Options getOptionsForTest()
{
    Options local;
    local.setSphereFile();
    local.setVolumeFile();
    local.setSearchMethod();
    local.setRangeBoxCommunication();
    local.setSphereBoxes();
    local.setIfTestIsGoldTestRun();

    return local;
}

void printOptions(const Options& options)
{
    int procId = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &procId);
    if ( procId == 0 )
    {
        std::cerr << "Sphere file: " << options.mSphereFile << std::endl;
        std::cerr << "Volume file: " << options.mVolumeFile << std::endl;
        std::cerr << "Search Method: " << options.mSearchMethod << std::endl;
    }
}

TEST(NaluPerformance, BoxSphereIntersections)
{
    Options options = getOptionsForTest();

    MPI_Comm comm = MPI_COMM_WORLD;
    std::vector< std::pair<Sphere, IdentProc> > spheres;
    fillBoundingVolumesUsingNodesFromFile(comm, options.mSphereFile, spheres);

    FloatBoxIdentProcVector domainBoxes;
    fillBoxesUsingElementBlocksFromFile(comm, options.mVolumeFile, domainBoxes);

    SearchResults searchResults;

    double startTime = stk::wall_time();

    if ( options.mSpheresFirstThenBoxes )
    {
        stk::search::coarse_search(spheres, domainBoxes, options.mSearchMethod, comm, searchResults, options.mCommunicateRangeBoxes);
    }
    else
    {
        stk::search::coarse_search(domainBoxes, spheres, options.mSearchMethod, comm, searchResults, options.mCommunicateRangeBoxes);
    }

    double elapsedTime = stk::wall_time() - startTime;
    stk::unit_test_util::printPeformanceStats(elapsedTime, comm);

    if ( !options.mTestToGetGoldResults )
    {
        stk::unit_test_util::gatherResultstoProcZero(comm, searchResults);

        int procId = -1;
        MPI_Comm_rank(comm, &procId);
        if ( procId == 0 )
        {
            std::vector< std::pair<int,int> > globalIdMapping(searchResults.size());
            for (size_t i=0; i<searchResults.size(); i++)
            {
                globalIdMapping[i] = std::make_pair(searchResults[i].first.id(), searchResults[i].second.id());
            }
            std::sort(globalIdMapping.begin(), globalIdMapping.end());
            std::vector< std::pair<int,int> >::iterator iter_end = std::unique(globalIdMapping.begin(), globalIdMapping.end());
            globalIdMapping.erase(iter_end, globalIdMapping.end());

            size_t numInteractions = stk::unit_test_util::getGoldValueForTest();
            EXPECT_EQ(numInteractions, globalIdMapping.size());
        }
    }
    else
    {
        int numProcs=0;
        MPI_Comm_size(comm, &numProcs);
        if ( numProcs != 1 )
        {
            std::cerr << "Gold results are available only on serial runs.\n";
        }
        else
        {
            printGoldResults(domainBoxes, spheres);
        }
    }
}


TEST(NaluPerformance, BoxBoxIntersections)
{
    Options options = getOptionsForTest();
    printOptions(options);

    MPI_Comm comm = MPI_COMM_WORLD;

    FloatBoxIdentProcVector spheres;
    fillBoundingVolumesUsingNodesFromFile(comm, options.mSphereFile, spheres);

    FloatBoxIdentProcVector domainBoxes;
    fillBoxesUsingElementBlocksFromFile(comm, options.mVolumeFile, domainBoxes);

    SearchResults searchResults;

    double startTime = stk::wall_time();

    int procId = -1;
    MPI_Comm_rank(comm, &procId);

    if ( options.mSpheresFirstThenBoxes )
    {
        stk::search::coarse_search(spheres, domainBoxes, options.mSearchMethod, comm, searchResults);
    }
    else
    {
        stk::search::coarse_search(domainBoxes, spheres, options.mSearchMethod, comm, searchResults);
    }

    double elapsedTime = stk::wall_time() - startTime;
    stk::unit_test_util::printPeformanceStats(elapsedTime, comm);

    if ( options.mTestToGetGoldResults )
    {
        int numProcs=0;
        MPI_Comm_size(comm, &numProcs);
        if ( numProcs != 1 )
        {
            std::cerr << "Gold results are available only on serial runs.\n";
        }
        else
        {
            printGoldResults(domainBoxes, spheres);
        }
    }
    else
    {
        stk::unit_test_util::gatherResultstoProcZero(comm, searchResults);

        if ( procId == 0 )
        {
            std::vector< std::pair<int,int> > globalIdMapping(searchResults.size());
            for (size_t i=0; i<searchResults.size(); i++)
            {
                globalIdMapping[i] = std::make_pair(searchResults[i].first.id(), searchResults[i].second.id());
            }
            std::sort(globalIdMapping.begin(), globalIdMapping.end());
            std::vector< std::pair<int,int> >::iterator iter_end = std::unique(globalIdMapping.begin(), globalIdMapping.end());
            globalIdMapping.erase(iter_end, globalIdMapping.end());

            size_t numInteractions = stk::unit_test_util::getGoldValueForTest();
            EXPECT_EQ(numInteractions, globalIdMapping.size());
        }
    }
}
 
TEST(stkSearch, boxSphereIntersection)
{
    FloatBox box(0,0,0,1,1,1);
    Sphere sphere(Point(2,2,2), 0.5);
    EXPECT_FALSE(stk::search::intersects(box, sphere));
    EXPECT_FALSE(stk::search::intersects(sphere, box));
    Sphere sphere1(Point(1.1, 1.1, 1.1), 0.2);
    EXPECT_TRUE(stk::search::intersects(box, sphere1));
    EXPECT_TRUE(stk::search::intersects(sphere1, box));
    Sphere sphere2(Point(1.1, 1.1, 1.1), 0.17321);
    EXPECT_TRUE(stk::search::intersects(box, sphere2));
    EXPECT_TRUE(stk::search::intersects(sphere2, box));
    Sphere sphere3(Point(0.5, 0.5, 0.5), 1);
    EXPECT_TRUE(stk::search::intersects(box, sphere3));
    EXPECT_TRUE(stk::search::intersects(sphere3, box));
    Sphere sphere4(Point(0.5, 0.5, 0.5), 0.1);
    EXPECT_TRUE(stk::search::intersects(box, sphere4));
    EXPECT_TRUE(stk::search::intersects(sphere4, box));
    Sphere sphere5(Point(1.5, 1.5, 1.5), 0.5);
    EXPECT_FALSE(stk::search::intersects(box, sphere5));
    EXPECT_FALSE(stk::search::intersects(sphere5, box));
    FloatBox box1(1,1,1,2,2,2);
    EXPECT_TRUE(stk::search::intersects(box1,box));
}

}
