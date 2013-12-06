#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_search/Box.hpp>
#include <stk_search/BoundingBox.hpp>
#include <stk_search/IdentProc.hpp>
#include <stk_search/OctTreeOps.hpp>

#include <Geom_AxisAlignedBB.h>
#include <Geom_Search.h>
#include <search/ContactRangeSearch.h>
#include <search/ContactCommunication.h>

namespace {

typedef stk::search::IdentProc<uint64_t, unsigned> Ident;
typedef stk::search::Point<double> Point;
typedef stk::search::Box<double> Box;

void printToCubit(std::ostream& out, Point const& min_corner, Point const& max_corner);

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

// void checkSearchResults(const int proc_id, SearchResults &goldResults, SearchResults& searchResults)
// {
//     ASSERT_TRUE(searchResults.size() >= goldResults.size());

//     Ident unUsedBox1(987654321, 1000000);
//     Ident unUsedBox2(987654322, 1000000);
//     size_t numResultsMatchingGoldResults = 0;

//     for (size_t i = 0; i < searchResults.size(); i++ )
//     {
// //        if ( proc_id == 0 )
// //        {
// //            std::cerr << searchResults[i].first << "\t" << searchResults[i].second << std::endl;
// //        }
//         for (size_t j = 0; j < goldResults.size(); j++)
//         {
//             if ( searchResults[i] == goldResults[j] )
//             {
//                 goldResults[j] = std::make_pair(unUsedBox1, unUsedBox2);
//                 numResultsMatchingGoldResults++;
//             }
//         }
//     }
//     EXPECT_EQ(goldResults.size(), numResultsMatchingGoldResults) << "proc id = " << proc_id << std::endl;
// }

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

  ASSERT_EQ(static_cast<size_t>(num_procs), recv_list.size());
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

void printToCubit(std::ostream& out, Point const& min_corner, Point const& max_corner)
{
    std::vector < geometry::Vec3d > corners;
    corners.push_back(geometry::Vec3d(min_corner[0], min_corner[1], min_corner[2]));
    corners.push_back(geometry::Vec3d(max_corner[0], min_corner[1], min_corner[2]));
    corners.push_back(geometry::Vec3d(max_corner[0], max_corner[1], min_corner[2]));
    corners.push_back(geometry::Vec3d(min_corner[0], max_corner[1], min_corner[2]));
    corners.push_back(geometry::Vec3d(min_corner[0], min_corner[1], max_corner[2]));
    corners.push_back(geometry::Vec3d(max_corner[0], min_corner[1], max_corner[2]));
    corners.push_back(geometry::Vec3d(max_corner[0], max_corner[1], max_corner[2]));
    corners.push_back(geometry::Vec3d(min_corner[0], max_corner[1], max_corner[2]));

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

        std::stringstream os;
        os << "cubitFile_" << proc_id << ".jou";
        std::ofstream file(os.str().c_str());
        printToCubit(file, min_corner, max_corner);

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

        {
          std::vector<double> globalBoxDouble(global_box.begin(), global_box.end());
          Point min_c(global_box[0], global_box[1], global_box[2]);
          Point max_c(global_box[2], global_box[3], global_box[4]);
          printToCubit(file, min_c, max_c);
        }

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

        const int tree_depth = 4;
        // unsigned maxOffsetForTreeDepthThree = stk::oct_tree_size(tree_depth-1);
        MPI_Barrier(MPI_COMM_WORLD);

        for (int procCounter = 0; procCounter < num_procs; procCounter++)
        {
            if(proc_id == procCounter)
            {
                std::cerr << "\t Nate:=====================================\n";
                std::cerr << "\t Nate: proc_id = " << procCounter << std::endl;
                for(SearchTree::const_iterator i = searchTree.begin(); i != searchTree.end(); ++i)
                {
                    const stk::OctTreeKey & key = (*i).first;

                    // EXPECT_EQ(0u, (stk::oct_tree_offset(tree_depth, key) - 1) % maxOffsetForTreeDepthThree);

                    const std::list<BoxIdent> & domain = (*i).second.first;
                    const std::list<BoxIdent> & range = (*i).second.second;
                    std::cerr << "\t Nate: depth = " << key.depth() << std::endl;
                    std::cerr << "\t Nate: ordinal = " << stk::oct_tree_offset(tree_depth, key) << std::endl;
                    std::cerr << "\t Nate: num_d = " << domain.size() << std::endl;
                    std::cerr << "\t Nate: num_r = " << range.size() << std::endl;
                    std::cerr << "\t Nate: key = " << key << std::endl;
                }
                std::cerr << "\t Nate:=====================================\n";
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
        file.close();
        unlink(os.str().c_str()); // comment this out to view Cubit files for bounding boxes

        const double tol = 0.001 ;

        std::vector< stk::OctTreeKey > cuts ;

        stk::search::oct_tree_partition( comm , searchTree , tol , cuts );

        if ( proc_id == 0 )
        {
            std::cerr << "Nate: For proc size of " << num_procs << std::endl;
            for (int i=0;i<num_procs;i++)
            {
              std::cerr << "Nate: cuts[" << i << "] = " << cuts[i] <<"\t with ordinal = " << stk::oct_tree_offset(tree_depth, cuts[i]) << std::endl;
            }
        }
    }
    else
    {
        std::cerr << "WARNING: Test not setup for anything other than 4 processors; ran with " << num_procs << "." << std::endl;
    }
}

STKUNIT_UNIT_TEST(stk_search, coarse_search_geometry_toolkit)
{
  testCoarseSearchUsingGeometryToolkit(MPI_COMM_WORLD);
}

} //namespace
