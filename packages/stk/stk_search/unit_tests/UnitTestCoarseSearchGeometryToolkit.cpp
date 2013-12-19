#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_search/Box.hpp>
#include <stk_search/BoundingBox.hpp>
#include <stk_search/IdentProc.hpp>
#include <stk_search/OctTreeOps.hpp>

namespace {

typedef stk::search::IdentProc<uint64_t, unsigned> Ident;
typedef stk::search::Point<double> Point;
typedef stk::search::Box<double> Box;

STKUNIT_UNIT_TEST(stk_search_oct_tree, checkCuts)
{
    typedef std::pair<Box, Ident> BoxIdent;
    typedef std::vector<BoxIdent> BoxVector;

    MPI_Comm comm = MPI_COMM_WORLD;
    int proc_id = -1;
    int num_procs = -1;
    MPI_Comm_rank(comm, &proc_id);
    MPI_Comm_size(comm, &num_procs);

    double offsetFromEdgeOfProcessorBoundary=0.1;
    double sizeOfDomainPerProcessor=1.0;
    double boxSize=0.8;
    ASSERT_TRUE(offsetFromEdgeOfProcessorBoundary<=sizeOfDomainPerProcessor-offsetFromEdgeOfProcessorBoundary);
    double min=offsetFromEdgeOfProcessorBoundary;
    double max=boxSize+offsetFromEdgeOfProcessorBoundary+1;

    if(num_procs >= 4)
    {
        Point min_corner, max_corner;

        min_corner[0] = offsetFromEdgeOfProcessorBoundary;
        min_corner[1] = offsetFromEdgeOfProcessorBoundary;
        min_corner[2] = offsetFromEdgeOfProcessorBoundary;
        max_corner[0] = boxSize+offsetFromEdgeOfProcessorBoundary;
        max_corner[1] = boxSize+offsetFromEdgeOfProcessorBoundary;
        max_corner[2] = boxSize+offsetFromEdgeOfProcessorBoundary;
        if (proc_id == 1)
        {
            min_corner[0] += sizeOfDomainPerProcessor;
            max_corner[0] += sizeOfDomainPerProcessor;
        }
        else if (proc_id == 3)
        {
            min_corner[1] += sizeOfDomainPerProcessor;
            max_corner[1] += sizeOfDomainPerProcessor;
        }
        else if (proc_id == 2)
        {
            min_corner[0] += sizeOfDomainPerProcessor;
            min_corner[1] += sizeOfDomainPerProcessor;
            max_corner[0] += sizeOfDomainPerProcessor;
            max_corner[1] += sizeOfDomainPerProcessor;
        }

        BoxVector local_domain;
        Ident domainBox(proc_id, proc_id);
        if ( proc_id < 4 )
        {
          local_domain.push_back(std::make_pair(Box(min_corner, max_corner), domainBox));
        }

        BoxVector local_range;
        local_range = local_domain;

        std::vector<float> global_box(6);
        stk::search::box_global_bounds(comm,
            local_domain.size(),
            &local_domain[0] ,
            local_range.size(),
            &local_range[0],
            &global_box[0]);

        float tolerance = 2*std::numeric_limits<float>::epsilon();
        EXPECT_NEAR(float(min), global_box[0], tolerance);
        EXPECT_NEAR(float(min), global_box[1], tolerance);
        EXPECT_NEAR(float(min), global_box[2], tolerance);
        EXPECT_NEAR(float(max), global_box[3], tolerance);
        EXPECT_NEAR(float(max), global_box[4], tolerance);
        EXPECT_NEAR(float(boxSize+offsetFromEdgeOfProcessorBoundary), global_box[5], tolerance);

        typedef std::map< stk::OctTreeKey, std::pair< std::list< BoxIdent >, std::list< BoxIdent > > > SearchTree ;

        SearchTree searchTree;

        bool local_violations = true;
        unsigned Dim = 3;
        stk::search::createSearchTree(&global_box[0], local_domain.size(), &local_domain[0],
                local_range.size(), &local_range[0], Dim, proc_id, local_violations,
                searchTree);

        const double tol = 0.001 ;
        std::vector< stk::OctTreeKey > cuts ;

        stk::search::oct_tree_partition( comm , searchTree , tol , cuts );

        const int tree_depth = 4;
        EXPECT_EQ(0u, stk::oct_tree_offset(tree_depth, cuts[0]));
        EXPECT_EQ(2u, stk::oct_tree_offset(tree_depth, cuts[1]));
        EXPECT_EQ(587u, stk::oct_tree_offset(tree_depth, cuts[2]));
        EXPECT_EQ(1172u, stk::oct_tree_offset(tree_depth, cuts[3]));

    }
    else
    {
        std::cerr << "WARNING: Test not setup for anything other than 4 processors; ran with " << num_procs << "." << std::endl;
    }
}

} //namespace
