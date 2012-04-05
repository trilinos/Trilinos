// @HEADER
// ***********************************************************************
//                Copyright message goes here. 
// ***********************************************************************
// @HEADER

/*! \file rcb.cpp
    \brief An example of partitioning coordinates with RCB.
*/

#include <Zoltan2_TestHelpers.hpp>
#include <Zoltan2_BasicCoordinateInput.hpp>
#include <Zoltan2_Partitioning.hpp>

using namespace std;

/*! \example rcb.cpp
    An example of the use of the RCB algorithm to partition coordinate data.
    \todo error handling
    \todo write some examples that don't use teuchos
    \todo check the solution, visualize it somehow
*/

int main(int argc, char *argv[])
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  ///////////////////////////////////////////////////////////////////////
  // Initialize MPI

  Teuchos::GlobalMPISession session(&argc, &argv);
  RCP<const Teuchos::Comm<int> > tcomm = Teuchos::DefaultComm<int>::getComm();
  int rank = tcomm->getRank();
  int nprocs = tcomm->getSize();

  ///////////////////////////////////////////////////////////////////////
  // Read coordinates from a file.

  std::string fname(testDataFilePath+"/USAir97.mtx");
  outputFlag_t flags(OBJECT_COORDINATES);

  typedef Tpetra::MultiVector<scalar_t, lno_t, gno_t, node_t> tMVector_t;

  UserInputForTests uinput(fname, tcomm, flags);

  RCP<tMVector_t> coords = uinput->getCoordinates();

  size_t localCount = coords->getLocalLength();
  size_t globalCount = coords->getGlobalLength();
  int dim = coords->getNumVectors();

  scalar_t *x=NULL, *y=NULL, *z=NULL;
  x = coords->getDataNonConst(0).getRawPtr();

  if (dim > 1){
    y = coords->getDataNonConst(1).getRawPtr();
    if (dim > 2)
      z = coords->getDataNonConst(2).getRawPtr();
  }

  gno_t *globalIds = coords->getMap()->getNodeElementList()->getRawPtr();
   
  ///////////////////////////////////////////////////////////////////////
  // Create a Zoltan2 input adapter for this geometry.

  typedef Zoltan2::BasicCoordinateInput<tMVector_t> inputAdapter_t;

  inputAdapter_t ia(localCount, globalIds, x, y, z, 1, 1, 1);
   
  ///////////////////////////////////////////////////////////////////////
  // Create parameters for an RCB problem

  Teuchos::ParameterList params("test params");
  myParams.set("debug_level", "detailed_status");
  myParams.set("debug_procs", "0");
  myParams.set("error_check_level", "debug_mode_assertions");

  Teuchos::ParameterList &parParams = params.sublist("partitioning");
  parParams.set("algorithm", "rcb");
  parParams.set("imbalance_tolerance", 1.1);
  parParams.set("num_global_parts", nprocs);
  parParams.set("bisection_num_test_cuts", 7);
   
  ///////////////////////////////////////////////////////////////////////
  // Create a Zoltan2 partitioning problem

  Zoltan2::PartitioningProblem<inputAdapter_t> problem(&ia, &params, 
    MPI_COMM_WORLD);
   
  ///////////////////////////////////////////////////////////////////////
  // Solve the problem

  problem.solve();
   
  ///////////////////////////////////////////////////////////////////////
  // Obtain the solution

  Zoltan2::PartitioningSolution<inputAdapter_t> &solution = 
    problem.getSolution();
   
  ///////////////////////////////////////////////////////////////////////
  // Check the solution.

  if (rank == 0)
    std::cout << "PASS" << std::endl;
}

