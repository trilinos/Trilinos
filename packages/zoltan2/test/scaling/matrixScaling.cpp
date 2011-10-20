// @HEADER
// ***********************************************************************
//
//         Zoltan2: Sandia Partitioning Ordering & Coloring Library
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
//
//  Create a large distributed matrix. Create a Zoltan2::GraphModel from
//  this matrix.  Observe memory usage and runtime.
//

#include <iostream>
#include <string>
#include <sstream>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_Time.hpp>
#include <Teuchos_CommandLineProcessor.hpp>

#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Map.hpp>

#include <Zoltan2_XpetraCrsMatrixInput.hpp>
#include <Zoltan2_GraphModel.hpp>

#include <MueLu_MatrixFactory.hpp>
#include <MueLu_GalleryParameters.hpp>

using namespace std;
using Teuchos::RCP;
using Teuchos::Comm;
using Teuchos::DefaultComm;
using Teuchos::CommandLineProcessor;
using Teuchos::REDUCE_SUM;
using Teuchos::REDUCE_MIN;
using Teuchos::REDUCE_MAX;
using Teuchos::reduceAll;

typedef float Scalar;
typedef int32_t LNO;
typedef int64_t GNO;

#define TEST_FAIL_AND_RETURN(comm, ok, s, code){ \
int gval, lval=( (ok) ? 0 : 1);       \
reduceAll<int,int>(comm, REDUCE_SUM, 1, &lval, &gval);\
if (gval){ \
  if ((comm).getRank() == 0){\
    std::cerr << "Error: " << s << std::endl;\
    std::cout << "FAIL" << std::endl;\
  } \
  return code;\
} \
}

#define START_TIMER { if (printTimingStats) timer.start(true); }
#define END_TIMER(s)   \
  if (printTimingStats){   \
    double elapsed = timer.stop(); \
    globalResults(s+std::string(" (s)"), elapsed, *comm); \
  }

#ifdef HAVE_MALLINFO
static LNO memBytes;
#define START_MEMCOUNT { \
  if (printMemoryStats) memBytes = Zoltan2::getAllocatedMemory();}
#define END_MEMCOUNT(s)  \
  if (printMemoryStats){   \
    memBytes = Zoltan2::getAllocatedMemory() - memBytes;   \
    globalResults(s+std::string(" (bytes)"), memBytes, *comm); \
  }
#else
#define START_MEMCOUNT
#define END_MEMCOUNT(s)
#endif

std::vector<std::string> resultLines;

template <typename T>
void globalResults(const std::string s, const T &lval, const Comm<int> &comm)
{
  T gval=0, min=0, max=0;
  reduceAll<int, T>(comm, REDUCE_SUM, 1, &lval, &gval);
  reduceAll<int, T>(comm, REDUCE_MIN, 1, &lval, &min);
  reduceAll<int, T>(comm, REDUCE_MAX, 1, &lval, &max);
  if (comm.getRank()==0){   \
    ostringstream oss;
    oss << s << ": min " << min;
    oss << ", max " << max;
    oss << ", avg " << gval/comm.getSize();
    oss << ", TOTAL " << gval << std::endl;
    resultLines.push_back(oss.str());
  }
}

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv);
  const RCP<const Comm<int> > &comm = DefaultComm<int>::getComm();
  Teuchos::Time timer("timer");
  int rank = comm->getRank();
  int fail=0;

  RCP<const Zoltan2::Environment> env =
    Teuchos::rcp(new Zoltan2::Environment);

  ///////////////////////////////////////////////////
  // options

  int xdim=10;  
  int ydim=10;
  int zdim=10;
  std::string matrixType("Laplace2D");
  int printMemoryStats=1;
  int printTimingStats=1;

  CommandLineProcessor clp(false, true);
  clp.setOption("xdim", &xdim, "X dimension of mesh that generates matrix.");
  clp.setOption("ydim", &ydim, "Y dimension of mesh that generates matrix.");
  clp.setOption("zdim", &zdim, "Z dimension of mesh that generates matrix.");
  clp.setOption("matrix", &matrixType, 
     "Matrix type: Laplace1D, Laplace2D, or Laplace3D");
  clp.setOption("memory", &printMemoryStats, "Print stats about memory usage.");
  clp.setOption("time", &printTimingStats, "Print timing stats.");
  int status = clp.parse(argc, argv) ;
  if (status != CommandLineProcessor::PARSE_SUCCESSFUL){
    if (status == CommandLineProcessor::PARSE_HELP_PRINTED) return 0;
    if (rank == 0) std::cout << "FAIL" << std::endl;
    return 1;
  }

  ///////////////////////////////////////////////////
  // Use Muelu to create a distributed matrix.

  typedef Tpetra::CrsMatrix<Scalar, LNO, GNO> tcrsMatrix_t;
  typedef Tpetra::Map<LNO, GNO> tMap_t;

  START_MEMCOUNT;
  Teuchos::CommandLineProcessor tclp;
  MueLu::Gallery::Parameters<GNO> params(tclp, xdim, ydim, zdim, matrixType);

  RCP<const tMap_t> map = Teuchos::rcp(new tMap_t(
      params.GetNumGlobalElements(), 0, comm));

  RCP<tcrsMatrix_t> matrix;

  try{
    START_TIMER;
    matrix = MueLu::Gallery::CreateCrsMatrix<Scalar, LNO, GNO,
      Tpetra::Map<LNO, GNO>, Tpetra::CrsMatrix<Scalar, LNO, GNO> >(
        params.GetMatrixType(), map, params.GetParameterList());
    END_TIMER("Build matrix");
  }
  catch (std::exception &e) { 
    std::cerr << rank << ": " << e.what();
    fail = 1;   // Probably not enough memory
  }

  END_MEMCOUNT("Matrix");

  TEST_FAIL_AND_RETURN(*comm, fail==0, "Creating a matrix", 1)

  ///////////////////////////////////////////////////
  // Create in input adapter for Zoltan2

  typedef Zoltan2::XpetraCrsMatrixInput<tcrsMatrix_t> inputAdapter_t;

  RCP<const inputAdapter_t> adapter;
  RCP<const tcrsMatrix_t> constMatrix = 
    Teuchos::rcp_const_cast<const tcrsMatrix_t>(matrix);

  START_MEMCOUNT
  try{
    START_TIMER;
    adapter = Teuchos::rcp(new inputAdapter_t(constMatrix));
    END_TIMER("Build input adapter");
  }
  catch(std::exception &e){
    std::cerr << rank << ": " << e.what();
    fail = 1;
  }
  END_MEMCOUNT("Input adapter");

  TEST_FAIL_AND_RETURN(*comm, fail==0, "Creating an input adapter", 1)

  ///////////////////////////////////////////////////
  // Create the graph model suitable for an algorithm

  typedef Zoltan2::GraphModel<inputAdapter_t> graphModel_t;

  RCP<graphModel_t> model;

  START_MEMCOUNT
  try{
    START_TIMER;
    model = Teuchos::rcp(new graphModel_t(adapter, comm, env));
    END_TIMER("Create graph model");
  }
  catch(std::exception &e){
    std::cerr << rank << ": " << e.what();
    fail = 1;
  }
  END_MEMCOUNT("Graph");

  TEST_FAIL_AND_RETURN(*comm, fail==0, "Creating a graph model", 1)

  if (rank == 0){
    std::vector<std::string>::iterator next = resultLines.begin();
    while (next != resultLines.end())
      std::cout << *next++;
    std::cout << "Graph: #vertices " << model->getGlobalNumVertices();
    std::cout << " #edges " << model->getGlobalNumEdges() << std::endl;
    std::cout << "PASS" << std::endl;
  }
}
