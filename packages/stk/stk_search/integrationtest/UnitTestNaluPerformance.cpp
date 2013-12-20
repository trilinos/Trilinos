#include <gtest/gtest.h>
#include <unit_tests/UnitTestUtils.hpp>
#include <unit_tests/MeshUtilsForBoundingVolumes.hpp>

namespace
{

void printGoldResults(const GtkBoxVector &domainBoxes, const std::vector< std::pair<Sphere, Ident> > &spheres)
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

TEST(NaluPerformance, BoxSphereIntersections)
{
    std::string sphereFile = getOption("-sphere", "junk");
    MPI_Comm comm = MPI_COMM_WORLD;
    std::vector< std::pair<Sphere, Ident> > spheres;
    fillBoundingVolumesUsingNodesFromFile(comm, sphereFile, spheres);

    std::string volumeFilename = getOption("-volume", "junk");
    GtkBoxVector domainBoxes;
    fillBoxesUsingElementBlocksFromFile(comm, volumeFilename, domainBoxes);

    SearchResults searchResults;
    stk::search::SearchMethod searchMethod = stk::search::OCTREE;
    if (getOption("-boost", "no") == "yes")
    {
        searchMethod = stk::search::BOOST_RTREE;
    }
    bool communicateRangeBoxInfo = true;
    if (getOption("-rangeBoxComm", "yes") == "no")
    {
        communicateRangeBoxInfo = false;
    }

    if ( getOption("-sb", "no" ) == "yes" )
    {
        stk::search::coarse_search(spheres, domainBoxes, searchMethod, comm, searchResults, communicateRangeBoxInfo);
    }
    else
    {
        stk::search::coarse_search(domainBoxes, spheres, searchMethod, comm, searchResults, communicateRangeBoxInfo);
    }

    if ( getOption("-getGold", "no") == "yes" )
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
        gatherResultstoProcZero(comm, searchResults);

        int procId=-1;
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

            size_t numInteractions = getGoldValueForTest();
            EXPECT_EQ(numInteractions, globalIdMapping.size());
        }
    }
}

TEST(stkSearch, boxSphereIntersection)
{
    GtkBox box(0,0,0,1,1,1);
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
}

}
