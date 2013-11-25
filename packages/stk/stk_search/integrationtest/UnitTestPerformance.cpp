#include <stk_search/CoarseSearch.hpp>
#include <stk_search/BoundingBox.hpp>
#include <stk_search/IdentProc.hpp>
#include <stk_util/util/memory_util.hpp>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/parallel/ParallelComm.hpp>

#include <mpi.h>
#include <gtest/gtest.h>
#include <vector>
#include <fstream>

#include <unit_tests/UnitTestUtils.hpp>

#include <Geom_AxisAlignedBB.h>
#include <Geom_Search.h>
#include <search/ContactRangeSearch.h>
#include <search/ContactCommunication.h>

#include <exodusMeshInterface.h>

extern int gl_argc;
extern char** gl_argv;

namespace
{

typedef stk::search::IdentProc<int,int> Ident;
typedef stk::search::Box<double> Box;
typedef std::vector< std::pair<Box,Ident> > BoxVector;
typedef std::vector<std::pair<Ident,Ident> > SearchResults;
typedef stk::search::Point<double> Point;

struct mybox
{
    double coordinates[6];
    void setCoordinates(double *inputCoord)
    {
        for (int i=0;i<6;i++)
        {
            coordinates[i] = inputCoord[i];
        }
    }
};

void printPeformanceStats(double elapsedTime, MPI_Comm comm);
void createBoundingBoxForElement(const sierra::Mesh::LocalNodeId *connectivity, const int numNodesPerElement,
        const std::vector<double> &coordinates, std::vector<double>& boxCoordinates);
void writeExodusFileUsingBoxes(const std::vector<mybox> & boxes, const std::string &filename);
void runStkSearchTest(stk::search::SearchMethod search);
void printSumOfResults(MPI_Comm comm, const size_t sizeResults);
void testGtkSearch(MPI_Comm comm, const std::vector<mybox>&domainBoxes, SearchResults& boxIdPairResults);
void fillBoxesUsingSidesetsFromFile(MPI_Comm comm, const std::string& filename, std::vector<mybox> &domainBoxes);
void writeExodusFileUsingBoxes(const std::vector<mybox>& boxes, const std::string &filename);
void testPerformanceOfAxisAlignedBoundingBoxes(stk::search::SearchMethod searchMethod, MPI_Comm comm);
std::string getOption(const std::string& option, const std::string defaultString = std::string("false"));
void testStkSearch(MPI_Comm comm, std::vector<mybox> &domainBoxes,
        stk::search::SearchMethod searchMethod, SearchResults boxIdPairResults);
void fillDomainBoxes(MPI_Comm comm, std::vector<mybox>& domainBoxes);

TEST(Performance, ofAxisAlignedBoundingBoxesUsingOctTree)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  stk::search::SearchMethod searchMethod = stk::search::OCTREE;
  testPerformanceOfAxisAlignedBoundingBoxes(searchMethod, comm);
}

TEST(Performance, ofAxisAlignedBoundingBoxesUsingBoostRtree)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  stk::search::SearchMethod searchMethod = stk::search::BOOST_RTREE;
  testPerformanceOfAxisAlignedBoundingBoxes(searchMethod, comm);
}

void testPerformanceOfAxisAlignedBoundingBoxes(stk::search::SearchMethod searchMethod, MPI_Comm comm)
{
    int proc = 0;
    MPI_Comm_rank(comm, &proc);
    int numProcs = 1;
    MPI_Comm_size(comm, &numProcs);

    size_t numColumnsPerProcessor = 1000;
    double boxSize = 1.0;

    std::vector<std::pair<Box, Ident> > smallBoxVector;
    std::vector<std::pair<Box, Ident> > bigBoxVector;

    std::vector<std::pair<Ident, Ident> > boxIdPairResults;

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

                smallBoxVector.push_back(generateBoundingVolume<Box>(centerX,
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
        bigBoxVector.push_back(generateBoundingVolume<Box>(centerX,
                centerY,
                centerZ,
                radius - boxSize / 2,
                id,
                proc));
    }

    double startTime = stk::wall_time();
    stk::search::coarse_search(smallBoxVector, bigBoxVector, searchMethod, comm, boxIdPairResults);
    double elapsedTime = stk::wall_time() - startTime;

    printPeformanceStats(elapsedTime, comm);

    size_t numExpectedResults = numColumnsPerProcessor * numColumnsPerProcessor;
    bool lastProcessorWithOddNumberOfProcs = numProcs % 2 != 0 && proc == numProcs - 1;
    if(lastProcessorWithOddNumberOfProcs)
    {
        numExpectedResults = 0;
    }

    EXPECT_EQ(numExpectedResults, boxIdPairResults.size());
}

////////////////////////////////////////////////////////////

TEST(Performance, stkSearchUsingBoost)
{
    runStkSearchTest(stk::search::BOOST_RTREE);
}

TEST(Performance, stkSearchUsingOcttree)
{
    runStkSearchTest(stk::search::OCTREE);
}

void runStkSearchTest(stk::search::SearchMethod searchMethod)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    std::vector<mybox> domainBoxes;
    fillDomainBoxes(comm, domainBoxes);

    SearchResults boxIdPairResults;
    testStkSearch(comm, domainBoxes, searchMethod, boxIdPairResults);
}

TEST(Performance, gtkSearch)
{
    MPI_Comm comm = MPI_COMM_WORLD;

    std::vector<mybox> domainBoxes;
    fillDomainBoxes(comm, domainBoxes);

    SearchResults boxIdPairResults;
    testGtkSearch(comm, domainBoxes, boxIdPairResults);
}

void testGtkSearch(MPI_Comm comm, const std::vector<mybox>&inputBoxes, SearchResults& searchResults)
{
    int num_procs = -1;
    int proc_id   = -1;
    MPI_Comm_rank(comm, &proc_id);
    MPI_Comm_size(comm, &num_procs);
    std::vector<int> procThatOwnsBox;

    std::vector<geometry::AxisAlignedBB> domainBoxes(inputBoxes.size());
    for(size_t i=0;i<inputBoxes.size();i++)
    {
        domainBoxes[i] = geometry::AxisAlignedBB(inputBoxes[i].coordinates[0],
                                                 inputBoxes[i].coordinates[1],
                                                 inputBoxes[i].coordinates[2],
                                                 inputBoxes[i].coordinates[3],
                                                 inputBoxes[i].coordinates[4],
                                                 inputBoxes[i].coordinates[5]
        );
        procThatOwnsBox.push_back(proc_id);
    }

    std::vector<geometry::AxisAlignedBB> rangeBoxes(domainBoxes);
    double startTime = stk::wall_time();
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

    double elapsedTime = stk::wall_time() - startTime;
    printPeformanceStats(elapsedTime, comm);

    EXPECT_EQ(domainBoxes.size(), first_interaction.size());
    EXPECT_EQ(domainBoxes.size(), last_interaction.size());

    // Ident box1, Ident box2
    for (size_t i=0;i<domainBoxes.size();i++)
    {
        Ident box1(i, procThatOwnsBox[i]);
        for (int j=first_interaction[i];j<last_interaction[i];j++)
        {
            Ident box2(interaction_list[j], procThatOwnsBox[interaction_list[j]]);
            searchResults.push_back(std::make_pair(box1, box2));
        }
    }

    printSumOfResults(comm, searchResults.size());
}

void testStkSearch(MPI_Comm comm, std::vector<mybox> &domainBoxes,
        stk::search::SearchMethod searchMethod, SearchResults boxIdPairResults)
{
    int procId=-1;
    MPI_Comm_rank(comm, &procId);

    int numProc=-1;
    MPI_Comm_size(comm, &numProc);

    BoxVector stkBoxes(domainBoxes.size());
    for (size_t i=0;i<domainBoxes.size();i++)
    {
        Point min(domainBoxes[i].coordinates[0], domainBoxes[i].coordinates[1], domainBoxes[i].coordinates[2]);
        Point max(domainBoxes[i].coordinates[3], domainBoxes[i].coordinates[4], domainBoxes[i].coordinates[5]);
        Ident domainBoxId(i, procId);
        stkBoxes[i] = std::make_pair(Box(min,max), domainBoxId);
    }

    double startTime = stk::wall_time();
    stk::search::coarse_search(stkBoxes, stkBoxes, searchMethod, comm, boxIdPairResults);
    double elapsedTime = stk::wall_time() - startTime;

    printPeformanceStats(elapsedTime, comm);
    printSumOfResults(comm, boxIdPairResults.size());

    int procIdDestination = 0;
    stk::CommAll gather(comm);
    for (int phase=0; phase<2; ++phase)
    {
        if ( procId != procIdDestination )
        {
            for (size_t j=0;j<boxIdPairResults.size();++j)
            {
                gather.send_buffer(procIdDestination).pack< std::pair<Ident, Ident> >(boxIdPairResults[j]);
            }
        }

        if (phase == 0) { //allocation phase
          gather.allocate_buffers( numProc / 4 );
        }
        else { // communication phase
          gather.communicate();
        }
    }

    if ( procId == procIdDestination )
    {
        for ( int p = 0 ; p < numProc ; ++p )
        {
            stk::CommBuffer &buf = gather.recv_buffer( p );
            while ( buf.remaining() )
            {
                std::pair<Ident, Ident> temp;
                buf.unpack< std::pair<Ident, Ident> >( temp );
                boxIdPairResults.push_back(temp);
            }
        }
        std::sort(boxIdPairResults.begin(), boxIdPairResults.end());
        SearchResults::iterator iter_end = std::unique(boxIdPairResults.begin(), boxIdPairResults.end());
        boxIdPairResults.erase(iter_end, boxIdPairResults.end());
        std::cerr << "Size of data on proc 0 after sorting: " << boxIdPairResults.size() << std::endl;
    }
}

TEST(Performance, getGoldResults)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int procId=-1;
    MPI_Comm_rank(comm, &procId);

    std::vector<mybox> domainBoxes;
    fillDomainBoxes(comm, domainBoxes);

    SearchResults boxIdPairResults;

    BoxVector stkBoxes(domainBoxes.size());
    for (size_t i=0;i<domainBoxes.size();i++)
    {
        Point min(domainBoxes[i].coordinates[0], domainBoxes[i].coordinates[1], domainBoxes[i].coordinates[2]);
        Point max(domainBoxes[i].coordinates[3], domainBoxes[i].coordinates[4], domainBoxes[i].coordinates[5]);
        Ident domainBoxId(i, procId);
        stkBoxes[i] = std::make_pair(Box(min,max), domainBoxId);
    }

    double startTime = stk::wall_time();
    for (size_t i=0;i<stkBoxes.size();++i)
    {
        for (size_t j=0;j<stkBoxes.size();++j)
        {
            if ( stk::search::intersects(stkBoxes[i].first, stkBoxes[j].first) )
            {
                boxIdPairResults.push_back(std::make_pair(stkBoxes[i].second, stkBoxes[j].second));
            }
        }
    }

    std::sort(boxIdPairResults.begin(), boxIdPairResults.end());
    SearchResults::iterator iter_end = std::unique(boxIdPairResults.begin(), boxIdPairResults.end());
    boxIdPairResults.erase(iter_end, boxIdPairResults.end());

    double elapsedTime = stk::wall_time() - startTime;
    printPeformanceStats(elapsedTime, comm);
    printSumOfResults(comm, boxIdPairResults.size());

    std::cerr << "Number of boxes: " << boxIdPairResults.size() << std::endl;
}

void createBoundingBoxesForElementsInElementBlocks(const int procId, const sierra::Mesh &mesh, const std::vector<double> &coordinates, BoxVector& domainBoxes)
{
    size_t numberBoundingBoxes = mesh.getNumberLocalElements();
    domainBoxes.resize(numberBoundingBoxes);

    sierra::Mesh::BlockIdVector blockIds;
    mesh.fillElementBlockIds(blockIds);

    size_t boxCounter = 0;

    std::vector<double> boxCoordinates(6);
    for (size_t elemBlockNum=0;elemBlockNum<mesh.getNumberElementBlocks();elemBlockNum++)
    {
        sierra::Mesh::LocalNodeIdVector connectivity;
        mesh.fillElementToLocalNodeConnectivityForBlock(blockIds[elemBlockNum], connectivity);
        int numNodesPerElement = mesh.getNumberNodesPerElement(blockIds[elemBlockNum]);
        size_t numElementsThisBlock = mesh.getNumberLocalElementsInBlock(blockIds[elemBlockNum]);
        for (size_t elemCounter=0;elemCounter<numElementsThisBlock;elemCounter++)
        {
            createBoundingBoxForElement(&connectivity[numNodesPerElement*elemCounter], numNodesPerElement, coordinates, boxCoordinates);
            Ident domainBoxId(boxCounter, procId);
            Point min(boxCoordinates[0], boxCoordinates[1], boxCoordinates[2]);
            Point max(boxCoordinates[3], boxCoordinates[4], boxCoordinates[5]);
            domainBoxes[boxCounter] = std::make_pair(Box(min, max), domainBoxId);
            boxCounter++;
        }
    }
}

void createBoundingBoxesForSidesInSidesets(const sierra::Mesh &mesh, const std::vector<double> &coordinates,
        std::vector<mybox>& domainBoxes)
{
    size_t numberBoundingBoxes = 0;
    sierra::Mesh::SideSetIdVector sidesetIds;
    mesh.fillSideSetIds(sidesetIds);

    for (size_t i=0;i<sidesetIds.size();i++)
    {
        numberBoundingBoxes += mesh.getNumberSidesInSideSet(sidesetIds[i]);
    }

    domainBoxes.resize(numberBoundingBoxes);

    size_t boxCounter = 0;

    std::vector<double> boxCoordinates(6);
    for (size_t ssetCounter=0;ssetCounter<sidesetIds.size();ssetCounter++)
    {
        sierra::Mesh::LocalNodeIdVector nodeIds;
        std::vector<int> numNodesPerFace;
        mesh.fillSideSetLocalNodeIds(sidesetIds[ssetCounter], nodeIds, numNodesPerFace);
        size_t offset=0;
        for (size_t i=0;i<numNodesPerFace.size();i++)
        {
            createBoundingBoxForElement(&nodeIds[offset], numNodesPerFace[i], coordinates, boxCoordinates);
            domainBoxes[boxCounter].setCoordinates(&boxCoordinates[0]);
            boxCounter++;
            offset += numNodesPerFace[i];
        }
    }

    ASSERT_EQ(boxCounter, numberBoundingBoxes);
}

void printSumOfResults(MPI_Comm comm, const size_t sizeResults)
{
    int procId=-1;
    MPI_Comm_rank(comm, &procId);
    int numResults=sizeResults;
    int sumOverAll=0;
    MPI_Allreduce(&numResults, &sumOverAll, 1, MPI_INT, MPI_SUM, comm);

    if (procId == 0 )
    {
        std::cerr << "Size of results: " << sumOverAll << std::endl;
    }
}

void fillBoxesUsingSidesetsFromFile(MPI_Comm comm, const std::string& filename, std::vector<mybox> &domainBoxes)
{
    sierra::ExodusMeshInterface mesh(filename, comm);

    std::vector<double> coordinates;
    mesh.fillCoordinates(coordinates);

    createBoundingBoxesForSidesInSidesets(mesh, coordinates, domainBoxes);
}

int openFileAndGetId(const int numBoxes, const int num_element_blocks, const std::string &filename)
{
    int CPU_word_size = sizeof(double);
    int IO_word_size = 8;
    int exoid = ex_create (filename.c_str(), EX_CLOBBER, &CPU_word_size, &IO_word_size);
    int num_dim = 3;
    int num_elements = numBoxes;
    int num_nodes_per_element = 8;
    int num_nodes = num_nodes_per_element*num_elements;

    int num_ns = 0, num_ss = 0;
    ex_put_init(exoid, "Boxes", num_dim, num_nodes, num_elements, num_element_blocks, num_ns, num_ss);
    return exoid;
}

void setHexCoordinates(const double &xmin, const double &ymin, const double &zmin, 
                       const double &xmax, const double &ymax, const double &zmax,
                       double* hexCoordinates)
{
//    int ordering[8] = { 4, 3, 2, 1, 8, 7, 6, 5 }; // one based!
    int ordering[8] = { 3, 2, 1, 0, 7, 6, 5, 4 }; 
   
    hexCoordinates[3*ordering[0]+0] = xmin; 
    hexCoordinates[3*ordering[0]+1] = ymin; 
    hexCoordinates[3*ordering[0]+2] = zmin; 

    hexCoordinates[3*ordering[1]+0] = xmax; 
    hexCoordinates[3*ordering[1]+1] = ymin; 
    hexCoordinates[3*ordering[1]+2] = zmin; 

    hexCoordinates[3*ordering[2]+0] = xmax; 
    hexCoordinates[3*ordering[2]+1] = ymin; 
    hexCoordinates[3*ordering[2]+2] = zmax; 

    hexCoordinates[3*ordering[3]+0] = xmin; 
    hexCoordinates[3*ordering[3]+1] = ymin; 
    hexCoordinates[3*ordering[3]+2] = zmax; 

    hexCoordinates[3*ordering[4]+0] = xmin; 
    hexCoordinates[3*ordering[4]+1] = ymax; 
    hexCoordinates[3*ordering[4]+2] = zmin; 

    hexCoordinates[3*ordering[5]+0] = xmax; 
    hexCoordinates[3*ordering[5]+1] = ymax; 
    hexCoordinates[3*ordering[5]+2] = zmin; 

    hexCoordinates[3*ordering[6]+0] = xmax; 
    hexCoordinates[3*ordering[6]+1] = ymax; 
    hexCoordinates[3*ordering[6]+2] = zmax; 

    hexCoordinates[3*ordering[7]+0] = xmin; 
    hexCoordinates[3*ordering[7]+1] = ymax; 
    hexCoordinates[3*ordering[7]+2] = zmax; 
}

void putCoordinatesInFile(const int exoid, const std::vector<mybox>& boxes)
{
    const int num_nodes_per_element = 8;
    const int spatialDim = 3;
    double *x = new double[num_nodes_per_element*boxes.size()];
    double *y = new double[num_nodes_per_element*boxes.size()];
    double *z = new double[num_nodes_per_element*boxes.size()];

    for (size_t i=0;i<boxes.size();i++)
    {
        double xmin=0, ymin=0, zmin=0;
        xmin = boxes[i].coordinates[0];
        ymin = boxes[i].coordinates[1];
        zmin = boxes[i].coordinates[2];

        double xmax=0, ymax=0, zmax=0;
        xmax = boxes[i].coordinates[3];
        ymax = boxes[i].coordinates[4];
        zmax = boxes[i].coordinates[5];
    
        double hexCoordinates[24];
        setHexCoordinates(xmin, ymin, zmin, xmax, ymax, zmax, &hexCoordinates[0]);
        unsigned offset = i*num_nodes_per_element;
        for (int j=0;j<num_nodes_per_element;j++)
        {
            x[offset+j] = hexCoordinates[spatialDim*j+0];            
            y[offset+j] = hexCoordinates[spatialDim*j+1];            
            z[offset+j] = hexCoordinates[spatialDim*j+2];            
        }
    }

    ex_put_coord(exoid, x, y, z);

    delete [] z;
    delete [] y;
    delete [] x;
}

void fillNumElementsPerBlock(const int num_elements, std::vector<int> &numElementsPerBlock)
{
    int numElementsPer=1;
    if ( num_elements < 100 )
    {
        numElementsPer = 1;
    }
    else if ( num_elements < 1000 )
    {
        numElementsPer=10;
    }
    else if ( num_elements < 10000 )
    {
        numElementsPer=100;
    }
    else
    {
        numElementsPer=1000;
    }

    for (int i=0;i<num_elements;i+=numElementsPer)
    {
        int numElementsThisBlock = (i+numElementsPer) < num_elements ? numElementsPer : num_elements-i;
        numElementsPerBlock.push_back(numElementsThisBlock);
    }
}

void writeExodusFileUsingBoxes(const std::vector<mybox>& boxes, const std::string &filename)
{
    if ( boxes.size() == 0 )
    {
        std::cerr << "Skipping writing of file. No boxes to write.\n";
    }

    const int num_nodes_per_elem = 8; 
    const int num_attr = 0;
    const unsigned num_elements = boxes.size();
    std::vector<int> numElementsPerBlock;
    fillNumElementsPerBlock(num_elements, numElementsPerBlock);
    const int num_blocks = numElementsPerBlock.size();
    const int exoid = openFileAndGetId(boxes.size(), num_blocks, filename);
    putCoordinatesInFile(exoid, boxes); 

    std::vector<int> connect(numElementsPerBlock[0]*num_nodes_per_elem);
    int ordering[8] = { 4, 3, 2, 1, 8, 7, 6, 5 }; // one based!
    unsigned offset = 0;
    for (int blockId=1;blockId<=num_blocks;blockId++)
    {
        const int num_elements_this_block = numElementsPerBlock[blockId-1];
        ex_put_elem_block(exoid, blockId, "HEX", num_elements_this_block, num_nodes_per_elem, num_attr);

        for (int j=0;j<num_nodes_per_elem*num_elements_this_block;j++)
        {
            connect[j] = ordering[j%num_nodes_per_elem]+offset+num_nodes_per_elem*(j/num_nodes_per_elem);
        }
        offset += num_elements_this_block*num_nodes_per_elem;

        ex_put_elem_conn(exoid, blockId, &connect[0]);
    }

    ex_close(exoid);
}

std::string getOption(const std::string& option, const std::string defaultString)
{
    std::string returnValue = defaultString;
    if ( gl_argv != 0 )
    {
        for (int i=0;i<gl_argc;i++)
        {
            std::string input_argv(gl_argv[i]);
            if ( option == input_argv )
            {
                if ( (i+1) < gl_argc )
                {
                    returnValue = std::string(gl_argv[i+1]);
                }
                break;
            }
        }
    }
    return returnValue;
}

void printPeformanceStats(double elapsedTime, MPI_Comm comm)
{
    long int maxHwm = 0, minHwm = 0;
    double avgHwm = 0;
    stk::get_memory_high_water_mark_across_processors(comm, maxHwm, minHwm, avgHwm);

    int proc=-1;
    MPI_Comm_rank(comm, &proc);

    int numProcs=0;
    MPI_Comm_size(comm, &numProcs);

    double minTime = 0, maxTime = 0, avgTime = 0;
    MPI_Allreduce(&elapsedTime, &maxTime, 1, MPI_DOUBLE, MPI_MAX, comm);
    MPI_Allreduce(&elapsedTime, &minTime, 1, MPI_DOUBLE, MPI_MIN, comm);
    double elapsedTimeDivided = elapsedTime/numProcs;
    MPI_Allreduce(&elapsedTimeDivided, &avgTime, 1, MPI_DOUBLE, MPI_SUM, comm);

    if (proc == 0)
    {
      double bytesInMegabyte = 1024*1024;
      std::cout << "Max time: "  << maxTime << ", Min time: " << minTime << ", Avg time: " << avgTime << std::endl;
      std::cout << std::setw(6) << std::fixed << std::setprecision(1) << "Max HWM: "<<double(maxHwm)/double(bytesInMegabyte)
        <<", Min HWM: "<<double(minHwm)/double(bytesInMegabyte)<<", Avg HWM: "<<avgHwm/bytesInMegabyte<<std::endl;
      std::cout<<"### Total Number of Steps Taken ###: 1"<<std::endl;
      std::cout<<"### Total Wall Clock Run Time Used ###: "<< maxTime <<std::endl;
    }
}

void createBoundingBoxForElement(const sierra::Mesh::LocalNodeId *connectivity, const int numNodesPerElement,
        const std::vector<double> &coordinates, std::vector<double>& boxCoordinates)
{
    int spatialDim = 3;
    double *minCoordinates = &boxCoordinates[0];
    double *maxCoordinates = &boxCoordinates[spatialDim];

    int firstNode=0;
    for (int j=0;j<spatialDim;j++)
    {
        minCoordinates[j] = coordinates[spatialDim*connectivity[firstNode]+j];
        maxCoordinates[j] = coordinates[spatialDim*connectivity[firstNode]+j];
    }

    for (int i=1;i<numNodesPerElement;i++)
    {
        sierra::Mesh::LocalNodeId nodeId = connectivity[i];
        for (int j=0;j<spatialDim;j++)
        {
           minCoordinates[j] = std::min(minCoordinates[j], coordinates[spatialDim*nodeId+j]);
           maxCoordinates[j] = std::max(maxCoordinates[j], coordinates[spatialDim*nodeId+j]);
        }
    }
    bool inflateBox = true;
    double percentInflation = 10;
    if ( inflateBox )
    {
        for (int i=0;i<spatialDim;i++)
        {
            double dist = maxCoordinates[i]-minCoordinates[i];
            if ( dist <= 1e-8 ) dist = 0.001;
            double inflation = dist*(0.5*percentInflation)/100.0;
            minCoordinates[i] -= inflation;
            maxCoordinates[i] += inflation;
        }
    }
}

void fillDomainBoxes(MPI_Comm comm, std::vector<mybox>& domainBoxes)
{
    std::string filename = getOption("-i", "input.exo");
    fillBoxesUsingSidesetsFromFile(comm, filename, domainBoxes);

    std::string exodusFilename = getOption("-o", "boxes.exo");
    writeExodusFileUsingBoxes(domainBoxes, exodusFilename);
}

}
