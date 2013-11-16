#include <stk_search/CoarseSearch.hpp>
#include <stk_search/BoundingBox.hpp>
#include <stk_search/IdentProc.hpp>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <algorithm>
#include <vector>
#include <iterator>
#include <cstdlib>
#include <sstream>
#include <fstream>

#include <Geom_AxisAlignedBB.h>
#include <Geom_Search.h>
#include <search/ContactRangeSearch.h>
#include <search/ContactCommunication.h>

namespace std {
template <typename Key, typename Proc>
std::ostream & operator<<(std::ostream & out, std::pair<stk::search::ident::IdentProc<Key,Proc>,stk::search::ident::IdentProc<Key,Proc> > const& ident)
{
  return out << "[" << ident.first << ":" << ident.second << "]";
}
} // namespace std


namespace {

typedef stk::search::ident::IdentProc<uint64_t, unsigned> Ident;
typedef std::vector<std::pair<Ident,Ident> > SearchResults;

void checkSearchResults(const int proc_id, SearchResults &goldResults, SearchResults& searchResults)
{
    ASSERT_TRUE(searchResults.size() >= goldResults.size());

    Ident unUsedBox1(987654321, 1000000);
    Ident unUsedBox2(987654322, 1000000);
    size_t numResultsMatchingGoldResults = 0;

    for (size_t i = 0; i < searchResults.size(); i++ )
    {
//        if ( proc_id == 0 )
//        {
//            std::cerr << searchResults[i].first << "\t" << searchResults[i].second << std::endl;
//        }
        for (size_t j = 0; j < goldResults.size(); j++)
        {
            if ( searchResults[i] == goldResults[j] )
            {
                goldResults[j] = std::make_pair(unUsedBox1, unUsedBox2);
                numResultsMatchingGoldResults++;
            }
        }
    }
    EXPECT_EQ(goldResults.size(), numResultsMatchingGoldResults) << "proc id = " << proc_id << std::endl;
}

void testCoarseSearchAABBForAlgorithm(stk::search::FactoryOrder::Algorithm algorithm, MPI_Comm comm)
{
    stk::search::FactoryOrder searchParameters;
    searchParameters.m_communicator = comm;
    searchParameters.m_algorithm = algorithm;

    typedef stk::search::box::AxisAlignedBoundingBox<Ident, double, 3> Box;
    typedef std::vector<Box> BoxVector;

    int num_procs = stk::parallel_machine_size(comm);
    int proc_id   = stk::parallel_machine_rank(comm);

    double data[6];

    BoxVector local_domain, local_range;
    // what if identifier is NOT unique
    // x_min <= x_max
    // y_min <= y_max
    // z_min <= z_max

    data[0] = proc_id + 0.1; data[1] = 0.0; data[2] = 0.0;
    data[3] = proc_id + 0.9; data[4] = 1.0; data[5] = 1.0;

    Ident domainBox1(proc_id*4, proc_id);
    local_domain.push_back(Box(data, domainBox1));

    data[0] = proc_id + 0.1; data[1] = 2.0; data[2] = 0.0;
    data[3] = proc_id + 0.9; data[4] = 3.0; data[5] = 1.0;

    Ident domainBox2(proc_id*4+1, proc_id);
    local_domain.push_back(Box(data, domainBox2));

    data[0] = proc_id + 0.6; data[1] = 0.5; data[2] = 0.0;
    data[3] = proc_id + 1.4; data[4] = 1.5; data[5] = 1.0;

    Ident rangeBox1(proc_id*4+2, proc_id);
    local_range.push_back(Box(data, rangeBox1));

    data[0] = proc_id + 0.6; data[1] = 2.5; data[2] = 0.0;
    data[3] = proc_id + 1.4; data[4] = 3.5; data[5] = 1.0;

    Ident rangeBox2(proc_id*4+3, proc_id);
    local_range.push_back(Box(data, rangeBox2));

    SearchResults searchResults;

    stk::search::coarse_search(searchResults, local_range, local_domain, searchParameters);
    SearchResults goldResults;

    if (num_procs == 1) {
        goldResults.push_back(std::make_pair(domainBox1, rangeBox1));
        goldResults.push_back(std::make_pair(domainBox2, rangeBox2));

        checkSearchResults(proc_id, goldResults, searchResults);
    }
    else {
        if (proc_id == 0) {
            Ident domainBox1OnProcessor1(4,1);
            Ident domainBox2OnProcessor1(5,1);
            goldResults.push_back(std::make_pair(domainBox1, rangeBox1));
            goldResults.push_back(std::make_pair(domainBox1OnProcessor1, rangeBox1));
            goldResults.push_back(std::make_pair(domainBox2, rangeBox2));
            goldResults.push_back(std::make_pair(domainBox2OnProcessor1, rangeBox2));

            checkSearchResults(proc_id, goldResults, searchResults);
        }
        else if (proc_id == num_procs - 1) {
            Ident rangeBox1OnPreviousProcessor1((proc_id-1)*4 + 2, proc_id - 1);
            Ident rangeBox2OnPreviousProcessor1((proc_id-1)*4 + 3, proc_id - 1);

            goldResults.push_back(std::make_pair(domainBox1, rangeBox1OnPreviousProcessor1));
            goldResults.push_back(std::make_pair(domainBox1, rangeBox1));
            goldResults.push_back(std::make_pair(domainBox2, rangeBox2OnPreviousProcessor1));
            goldResults.push_back(std::make_pair(domainBox2, rangeBox2));

            checkSearchResults(proc_id, goldResults, searchResults);
        }
        else {
            Ident rangeBox1OnPreviousProcessor((proc_id-1)*4 + 2, proc_id - 1);
            Ident rangeBox2OnPreviousProcessor((proc_id-1)*4 + 3, proc_id - 1);
            Ident domainBox1OnNextProcessor((proc_id+1)*4,     proc_id + 1);
            Ident domainBox2OnNextProcessor((proc_id+1)*4 + 1, proc_id + 1);

            goldResults.push_back(std::make_pair(domainBox1, rangeBox1OnPreviousProcessor));
            goldResults.push_back(std::make_pair(domainBox1, rangeBox1));
            goldResults.push_back(std::make_pair(domainBox2, rangeBox2OnPreviousProcessor));
            goldResults.push_back(std::make_pair(domainBox2, rangeBox2));
            goldResults.push_back(std::make_pair(domainBox1OnNextProcessor, rangeBox1));
            goldResults.push_back(std::make_pair(domainBox2OnNextProcessor, rangeBox2));

            checkSearchResults(proc_id, goldResults, searchResults);
        }
    }
}

bool compare(const std::pair<geometry::AxisAlignedBB, geometry::AxisAlignedBB> &gold,
             const std::pair<geometry::AxisAlignedBB, geometry::AxisAlignedBB> &result)
{
    bool retVal = false;
    if ( gold.first.get_x_min() == result.first.get_x_min() &&
         gold.first.get_y_min() == result.first.get_y_min() &&
         gold.first.get_z_min() == result.first.get_z_min() )
    {
         if (gold.second.get_x_min() == result.second.get_x_min() &&
             gold.second.get_y_min() == result.second.get_y_min() &&
             gold.second.get_z_min() == result.second.get_z_min() )
        {
            retVal = true;
        }
    }
    return retVal;
}


void checkSearchResults(const int proc_id, const std::vector<std::pair<geometry::AxisAlignedBB, geometry::AxisAlignedBB> > &goldResults,
        const std::vector<std::pair<geometry::AxisAlignedBB, geometry::AxisAlignedBB> > &searchResults)
{
    ASSERT_EQ(goldResults.size(), searchResults.size());
    for (size_t i=0; i<goldResults.size();i++)
    {
        EXPECT_TRUE(compare(goldResults[i], searchResults[i])) << " failed for processor " << proc_id << " for interaction " << i << std::endl;
    }
}

void testCoarseSearchUsingGeometryToolkit(MPI_Comm comm)
{
    int num_procs = -1;
    int proc_id   = -1;
    MPI_Comm_rank(comm, &proc_id);
    MPI_Comm_size(comm, &num_procs);

    double data[6];

    std::vector<geometry::AxisAlignedBB> domainBoxes;
    data[0] = proc_id + 0.1; data[1] = 0.0; data[2] = 0.0;
    data[3] = proc_id + 0.9; data[4] = 1.0; data[5] = 1.0;
    domainBoxes.push_back(geometry::AxisAlignedBB(data[0], data[1], data[2], data[3], data[4], data[5]));

    data[0] = proc_id + 0.1; data[1] = 2.0; data[2] = 0.0;
    data[3] = proc_id + 0.9; data[4] = 3.0; data[5] = 1.0;
    domainBoxes.push_back(geometry::AxisAlignedBB(data[0], data[1], data[2], data[3], data[4], data[5]));

    std::vector<geometry::AxisAlignedBB> rangeBoxes;
    std::vector<int> procThatOwnsBox;

    data[0] = proc_id + 0.6; data[1] = 0.5; data[2] = 0.0;
    data[3] = proc_id + 1.4; data[4] = 1.5; data[5] = 1.0;
    rangeBoxes.push_back(geometry::AxisAlignedBB(data[0], data[1], data[2], data[3], data[4], data[5]));
    procThatOwnsBox.push_back(proc_id);

    data[0] = proc_id + 0.6; data[1] = 2.5; data[2] = 0.0;
    data[3] = proc_id + 1.4; data[4] = 3.5; data[5] = 1.0;
    rangeBoxes.push_back(geometry::AxisAlignedBB(data[0], data[1], data[2], data[3], data[4], data[5]));
    procThatOwnsBox.push_back(proc_id);

    std::vector<int> ghost_indices;
    std::vector<int> ghost_procs;
    ACME::BoxA_BoxB_Ghost(domainBoxes, rangeBoxes, comm, ghost_indices, ghost_procs);

    std::vector< std::vector<geometry::AxisAlignedBB> > send_list(num_procs);
    std::vector< std::vector<geometry::AxisAlignedBB> > recv_list(num_procs);

    for (size_t i=0;i<ghost_indices.size();i++)
    {
        send_list[ghost_procs[i]].push_back(rangeBoxes[ghost_indices[i]]);
    }

    ACME::Parallel_Data_Exchange(send_list, recv_list, comm);

    ASSERT_EQ((size_t)num_procs, recv_list.size());
    for (size_t i=0;i<recv_list.size();i++)
    {
        for (size_t j=0;j<recv_list[i].size();j++)
        {
            rangeBoxes.push_back(recv_list[i][j]);
            procThatOwnsBox.push_back(i);
        }
    }
    std::vector<int> interaction_list;
    std::vector<int> first_interaction;
    std::vector<int> last_interaction;

    geometry::BoxA_BoxB_Search(domainBoxes, rangeBoxes, interaction_list, first_interaction, last_interaction);

    EXPECT_EQ(domainBoxes.size(), first_interaction.size());
    EXPECT_EQ(domainBoxes.size(), last_interaction.size());

    std::vector<std::pair<geometry::AxisAlignedBB, geometry::AxisAlignedBB> > searchResults;

    std::stringstream os;
    os << "cubitFile_" << proc_id << ".jou";
    std::ofstream file(os.str().c_str());

    for (size_t i=0;i<domainBoxes.size();i++)
    {
        domainBoxes[i].PrintToCubit(file);
        for (int j=first_interaction[i];j<last_interaction[i];j++)
        {
            searchResults.push_back(std::make_pair(domainBoxes[i], rangeBoxes[interaction_list[j]]));
            rangeBoxes[interaction_list[j]].PrintToCubit(file);
        }
    }

    file.close();

    std::vector<std::pair<geometry::AxisAlignedBB, geometry::AxisAlignedBB> > goldResults;

    if (num_procs == 1)
    {
        goldResults.push_back(std::make_pair(domainBoxes[0], rangeBoxes[0]));
        goldResults.push_back(std::make_pair(domainBoxes[1], rangeBoxes[1]));
        checkSearchResults(proc_id, goldResults, searchResults);
    }
    else
    {
        if ( proc_id == 0 )
        {
            goldResults.push_back(std::make_pair(domainBoxes[0], rangeBoxes[0]));
            goldResults.push_back(std::make_pair(domainBoxes[1], rangeBoxes[1]));
            checkSearchResults(proc_id, goldResults, searchResults);
        }
        else if ( proc_id == num_procs - 1)
        {
            // Ordering is unfortunately important: for each domain box start with range box on this proc then the other procs
            goldResults.push_back(std::make_pair(domainBoxes[0], rangeBoxes[0]));
            goldResults.push_back(std::make_pair(domainBoxes[0], rangeBoxes[2]));
            goldResults.push_back(std::make_pair(domainBoxes[1], rangeBoxes[1]));
            goldResults.push_back(std::make_pair(domainBoxes[1], rangeBoxes[3]));
            checkSearchResults(proc_id, goldResults, searchResults);
        }
        else
        {
            goldResults.push_back(std::make_pair(domainBoxes[0], rangeBoxes[0]));
            goldResults.push_back(std::make_pair(domainBoxes[0], rangeBoxes[2]));
            goldResults.push_back(std::make_pair(domainBoxes[1], rangeBoxes[1]));
            goldResults.push_back(std::make_pair(domainBoxes[1], rangeBoxes[3]));
            checkSearchResults(proc_id, goldResults, searchResults);
        }
    }
    unlink(os.str().c_str()); // comment this out to view Cubit files for bounding boxes
}

//  axis aligned bounding box search

STKUNIT_UNIT_TEST(stk_search_not_boost, coarse_search_geometry_toolkit)
{
  testCoarseSearchUsingGeometryToolkit(MPI_COMM_WORLD);
}

STKUNIT_UNIT_TEST(stk_search_not_boost, coarse_search_3D_oct_tree)
{
  testCoarseSearchAABBForAlgorithm(stk::search::FactoryOrder::OCTREE, MPI_COMM_WORLD);
}

STKUNIT_UNIT_TEST(stk_search_not_boost, coarse_search_3D_bih_tree)
{
  testCoarseSearchAABBForAlgorithm(stk::search::FactoryOrder::BIHTREE, MPI_COMM_WORLD);
}


STKUNIT_UNIT_TEST(stk_search, coarse_search_3D_one_point)
{
  typedef stk::search::ident::IdentProc<uint64_t, unsigned> Ident;
  typedef stk::search::box::AxisAlignedBoundingBox<Ident, double, 3> Box;
  typedef std::vector<Box> BoxVector;
  typedef std::vector<std::pair<Ident,Ident> > SearchResults;

  stk::ParallelMachine comm = MPI_COMM_WORLD;
  //int num_procs = stk::parallel_machine_size(comm);
  int proc_id   = stk::parallel_machine_rank(comm);

  double data[6];

  BoxVector local_domain, local_range;
  // what if identifier is NOT unique
  // x_min <= x_max
  // y_min <= y_max
  // z_min <= z_max

  data[0] = 0.0; data[1] = 0.0; data[2] = 0.0;
  data[3] = 1.0; data[4] = 1.0; data[5] = 1.0;

  // One bounding box on processor 0 with the label:  0
  // All other processors have empty domain.
  Ident domainBox1(0, 0);
  if (proc_id == 0) {
    local_domain.push_back(Box(data, domainBox1));
  }

  data[0] = 0.5; data[1] = 0.5; data[2] = 0.5;
  data[3] = 0.5; data[4] = 0.5; data[5] = 0.5;

  // One range target on processor 0 with the label:  1
  // All other processors have empty range.
  Ident rangeBox1(1, 0);
  if (proc_id == 0) {
    local_range.push_back(Box(data, rangeBox1));
  }

  SearchResults searchResults;

  stk::search::FactoryOrder searchParameters;
  searchParameters.m_algorithm = stk::search::FactoryOrder::OCTREE;
  searchParameters.m_communicator = comm;
  stk::search::coarse_search(searchResults, local_range, local_domain, searchParameters);

  if (proc_id == 0) {
    STKUNIT_ASSERT_EQ(searchResults.size(), 1u);
    STKUNIT_EXPECT_EQ(searchResults[0], std::make_pair(domainBox1, rangeBox1));
  } else {
    STKUNIT_ASSERT_EQ(searchResults.size(), 0u);
  }
}

void printToCubit(std::ostream& out, const double* data)
{
    std::vector < geometry::Vec3d > corners;
    corners.push_back(geometry::Vec3d(data[0], data[1], data[2]));
    corners.push_back(geometry::Vec3d(data[3], data[1], data[2]));
    corners.push_back(geometry::Vec3d(data[3], data[4], data[2]));
    corners.push_back(geometry::Vec3d(data[0], data[4], data[2]));
    corners.push_back(geometry::Vec3d(data[0], data[1], data[5]));
    corners.push_back(geometry::Vec3d(data[3], data[1], data[5]));
    corners.push_back(geometry::Vec3d(data[3], data[4], data[5]));
    corners.push_back(geometry::Vec3d(data[0], data[4], data[5]));

    corners[0].PrintToCubit(out);
    out << "# {v0 = Id(\"vertex\")}" << std::endl;
    corners[1].PrintToCubit(out);
    out << "# {v1 = Id(\"vertex\")}" << std::endl;
    corners[2].PrintToCubit(out);
    out << "# {v2 = Id(\"vertex\")}" << std::endl;
    corners[3].PrintToCubit(out);
    out << "# {v3 = Id(\"vertex\")}" << std::endl;
    corners[4].PrintToCubit(out);
    out << "# {v4 = Id(\"vertex\")}" << std::endl;
    corners[5].PrintToCubit(out);
    out << "# {v5 = Id(\"vertex\")}" << std::endl;
    corners[6].PrintToCubit(out);
    out << "# {v6 = Id(\"vertex\")}" << std::endl;
    corners[7].PrintToCubit(out);
    out << "# {v7 = Id(\"vertex\")}" << std::endl;

    out << "create surface vertex {v0} {v1} {v2} {v3}" << std::endl;
    out << "create surface vertex {v0} {v4} {v5} {v1}" << std::endl;
    out << "create surface vertex {v1} {v5} {v6} {v2}" << std::endl;
    out << "create surface vertex {v3} {v2} {v6} {v7}" << std::endl;
    out << "create surface vertex {v0} {v3} {v7} {v4}" << std::endl;
    out << "create surface vertex {v4} {v7} {v6} {v5}" << std::endl;
}

STKUNIT_UNIT_TEST(stk_search_not_boost, checkCuts)
{
    typedef stk::search::box::AxisAlignedBoundingBox<Ident, double, 3> Box;
    typedef std::vector<Box> BoxVector;

    MPI_Comm comm = MPI_COMM_WORLD;
    int proc_id = -1;
    int num_procs = -1;
    MPI_Comm_rank(comm, &proc_id);
    MPI_Comm_size(comm, &num_procs);

    double offset=0.1;
    double size=1.0;
    double boxSize=0.1;
    ASSERT_TRUE(offset<=size-offset);
    double min=offset;
    double max=boxSize+offset+1;

    if(num_procs == 4)
    {
        double data[6];
        data[0] = offset;
        data[1] = offset;
        data[2] = offset;
        data[3] = boxSize+offset;
        data[4] = boxSize+offset;
        data[5] = boxSize+offset;
        if (proc_id == 1)
        {
            data[0] += 1.0;
            data[3] += 1.0;
        }
        else if (proc_id == 2)
        {
            data[1] += 1.0;
            data[4] += 1.0;
        }
        else if (proc_id == 3)
        {
            data[0] += 1.0;
            data[1] += 1.0;
            data[3] += 1.0;
            data[4] += 1.0;
        }

        std::stringstream os;
        os << "cubitFile_" << proc_id << ".jou";
        std::ofstream file(os.str().c_str());
        printToCubit(file, data);

        BoxVector local_domain;
        Ident domainBox(proc_id, proc_id);
        local_domain.push_back(Box(data, domainBox));

        BoxVector local_range;
        local_range = local_domain;

        std::vector<float> global_box(6);
        stk::search::box_global_bounds(comm,
            local_domain.size(),
            &local_domain[0] ,
            local_range.size(),
            &local_range[0],
            &global_box[0]);

        std::vector<double> globalBoxDouble(global_box.begin(), global_box.end());
        printToCubit(file, &globalBoxDouble[0]);

        float tolerance = 2*std::numeric_limits<float>::epsilon();
        EXPECT_NEAR(float(min), global_box[0], tolerance);
        EXPECT_NEAR(float(min), global_box[1], tolerance);
        EXPECT_NEAR(float(min), global_box[2], tolerance);
        EXPECT_NEAR(float(max), global_box[3], tolerance);
        EXPECT_NEAR(float(max), global_box[4], tolerance);
        EXPECT_NEAR(float(boxSize+offset), global_box[5], tolerance);

        typedef std::map< stk::OctTreeKey, std::pair< std::list< Box >, std::list< Box > > > SearchTree ;

        SearchTree searchTree;

        bool local_violations = true;
        unsigned Dim = 3;
        stk::search::createSearchTree(&global_box[0], local_domain.size(), &local_domain[0],
                local_range.size(), &local_range[0], Dim, proc_id, local_violations,
                searchTree);

        const int tree_depth = 4;
        for (int procCounter = 0; procCounter < 4; procCounter++)
        {
            if(proc_id == procCounter)
            {
                for(typename SearchTree::const_iterator i = searchTree.begin();
                        i != searchTree.end(); ++i)
                {
                    const stk::OctTreeKey & key = (*i).first;

                    const std::list<Box> & domain = (*i).second.first;
                    const std::list<Box> & range = (*i).second.second;

                    std::cerr << "\t depth = " << key.depth() << std::endl;
                    std::cerr << "\t ordinal = " << oct_tree_offset(tree_depth, key) << std::endl;
                    std::cerr << "\t num_d = " << domain.size() << std::endl;
                    std::cerr << "\t num_r = " << range.size() << std::endl;
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
        file.close();
        unlink(os.str().c_str()); // comment this out to view Cubit files for bounding boxes
    }
}


} //namespace
