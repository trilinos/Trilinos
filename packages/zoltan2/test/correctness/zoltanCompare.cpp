// @HEADER
// ***********************************************************************
//                Copyright message goes here. 
// ***********************************************************************
// @HEADER

/*! \file rcbPerformance.cpp
    \brief A test that can do large scale problems and time them.
*/

#include <Zoltan2_TestHelpers.hpp>
#include <Zoltan2_BasicCoordinateInput.hpp>
#include <Zoltan2_PartitioningSolution.hpp>
#include <Zoltan2_PartitioningProblem.hpp>

using namespace std;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::Comm;

typedef Tpetra::MultiVector<scalar_t, lno_t, gno_t, node_t> tMVector_t;
typedef Zoltan2::BasicCoordinateInput<tMVector_t> inputAdapter_t;

int runRCP(const RCP<const Comm<int> > &comm,
  string fname, bool average_cuts, bool rectilinear_blocks)
{
  int rank = comm->getRank();
  int nprocs = comm->getSize();

  if (rank == 0){
    std::cout << "=================================" << std::endl;
    std::cout << "Number of processes: " << nprocs << std::endl;
    std::cout << fname << std::endl;
    if (average_cuts)
      std::cout << "average cuts is ON" << std::endl;
    if (rectilinear_blocks)
      std::cout << "rectilinear blocks is ON" << std::endl;
    std::cout << std::endl;
  }
  // Read in coordinates from file.

  std::string fullname(testDataFilePath+"/zoltan/"+fname);
  outputFlag_t flags;
  flags.set(OBJECT_COORDINATES);

  RCP<UserInputForTests> uinput;
  try{
    uinput = rcp(new UserInputForTests(fullname, comm, flags));
  }
  catch(...){
    if (rank == 0)
      std::cout << "FAIL: UserInputForTests" << std::endl;
    return 1;
  }

  RCP<tMVector_t> coords;
  try{
   coords = uinput->getCoordinates();
  }
  catch(...){
    if (rank == 0)
      std::cout << "FAIL: get coordinates" << std::endl;
    return 1;
  }

  // Create an input adapter.

  const gno_t *globalIds = 
    coords->getMap()->getNodeElementList().getRawPtr();

  size_t localCount = coords->getLocalLength();
  int dim = coords->getNumVectors();

  scalar_t *x=NULL, *y=NULL, *z=NULL;
  x = coords->getDataNonConst(0).getRawPtr();

  if (dim > 1){
    y = coords->getDataNonConst(1).getRawPtr();
    if (dim > 2)
      z = coords->getDataNonConst(2).getRawPtr();
  }

  RCP<inputAdapter_t> ia;

  try{
    ia = rcp(new inputAdapter_t(
      localCount, globalIds, x, y, z, 1, 1, 1));
  }
  catch(...){
    if (rank == 0)
      std::cout << "FAIL: input adapter" << std::endl;
    return 1;
  }

 // Parameters

  Teuchos::ParameterList params;
  params.set("timing_output_stream" , "std::cout");

  Teuchos::ParameterList &parParams = params.sublist("partitioning");
  parParams.set("algorithm", "rcb");
  double tolerance = 1.1;
  parParams.set("imbalance_tolerance", tolerance );

  Teuchos::ParameterList &geoParams = parParams.sublist("geometric");
  geoParams.set("bisection_num_test_cuts", 5);
  if (rectilinear_blocks)
    geoParams.set("rectilinear_blocks", "yes");
  if (average_cuts)
    geoParams.set("average_cuts", "yes");

  // Create the problem.

  RCP<Zoltan2::PartitioningProblem<inputAdapter_t> > problem;
  try{
    problem = rcp(new Zoltan2::PartitioningProblem<inputAdapter_t>(
      ia.getRawPtr(), &params));
  }
  catch(...){
    if (rank == 0)
      std::cout << "FAIL: problem" << std::endl;
    return 1;
  }

  try{
    problem->solve();
  }
  catch(...){
    if (rank == 0)
      std::cout << "FAIL: solve" << std::endl;
    return 1;
  }

  if (rank == 0)
    problem->getSolution().printMetrics(cout);

  return 0;
}
  
int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv);
  RCP<const Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
  int rank = comm->getRank();

  Teuchos::CommandLineProcessor cmdp (false, false);

  string inputFile("none"), average_cuts("no"), rectilinear_blocks("no");

  cmdp.setOption("inputFile", &inputFile, 
    "root of file name: \"grid20x19\" for \"grid20x19_coord.mtx\"", true);
  cmdp.setOption("average_cuts", &average_cuts, 
    "yes or no");
  cmdp.setOption("rectilinear_blocks", &rectilinear_blocks, 
    "yes or no");

  try{
    cmdp.parse(argc, argv);
  }
  catch(...){
    if (rank == 0)
      std::cout << "FAIL: arguments" << std::endl;
    return 1;
  }

  if (inputFile == string("none"))
    return 0;

  bool ac = false, rb = false;
  if (average_cuts == string("yes"))
    ac = true;
  if (rectilinear_blocks == string("yes"))
    rb = true;

  string fname(inputFile+".mtx");

  int fail = runRCP(comm, fname, ac, rb);

  if (rank == 0 && !fail)
    std::cout << "PASS" << std::endl;

  return 0;
}
