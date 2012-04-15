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
using Zoltan2::MetricValues;

/*! \example rcb.cpp
    An example of the use of the RCB algorithm to partition coordinate data.
    \todo error handling
    \todo write some examples that don't use teuchos
    \todo check the solution, visualize it somehow
    \todo this has become a test - move to test dir and create an example
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
  outputFlag_t flags;
  flags.set(OBJECT_COORDINATES);

  typedef Tpetra::MultiVector<scalar_t, lno_t, gno_t, node_t> tMVector_t;

  UserInputForTests uinput(fname, tcomm, flags);

  RCP<tMVector_t> coords = uinput.getCoordinates();

  size_t localCount = coords->getLocalLength();
  //size_t globalCount = coords->getGlobalLength();
  int dim = coords->getNumVectors();

  scalar_t *x=NULL, *y=NULL, *z=NULL;
  x = coords->getDataNonConst(0).getRawPtr();

  if (dim > 1){
    y = coords->getDataNonConst(1).getRawPtr();
    if (dim > 2)
      z = coords->getDataNonConst(2).getRawPtr();
  }

  const gno_t *globalIds = coords->getMap()->getNodeElementList().getRawPtr();
   
  ///////////////////////////////////////////////////////////////////////
  // Create a Zoltan2 input adapter for this geometry.

  typedef Zoltan2::BasicCoordinateInput<tMVector_t> inputAdapter_t;
  inputAdapter_t ia(localCount, globalIds, x, y, z, 1, 1, 1);
   
  ///////////////////////////////////////////////////////////////////////
  // Create parameters for an RCB problem

  Teuchos::ParameterList params("test params");
  params.set("debug_level", "detailed_status");
  params.set("debug_procs", "0");
  params.set("error_check_level", "debug_mode_assertions");

  Teuchos::ParameterList &parParams = params.sublist("partitioning");
  parParams.set("algorithm", "rcb");
  parParams.set("imbalance_tolerance", 1.1);
  parParams.set("num_global_parts", nprocs);

  Teuchos::ParameterList &geoParams = parParams.sublist("geometric");
  geoParams.set("bisection_num_test_cuts", 7);

  ///////////////////////////////////////////////////////////////////////
  // Create a Zoltan2 partitioning problem

#ifdef ZOLTAN2_HAVE_MPI
  Zoltan2::PartitioningProblem<inputAdapter_t> problem(&ia, &params,
    MPI_COMM_WORLD);
#else
  Zoltan2::PartitioningProblem<inputAdapter_t> problem(&ia, &params);
#endif
   
  ///////////////////////////////////////////////////////////////////////
  // Solve the problem

  problem.solve();
   
  ///////////////////////////////////////////////////////////////////////
  // Obtain the solution

  const Zoltan2::PartitioningSolution<inputAdapter_t> &solution = 
    problem.getSolution();
   
  ///////////////////////////////////////////////////////////////////////
  // Check the solution.

  const ArrayRCP<MetricValues<scalar_t> > & metrics1 =
    solution.getMetrics();

  if (rank == 0)
    Zoltan2::printMetrics<scalar_t>(cout, nprocs, nprocs, nprocs, 
      metrics1.view(0,metrics1.size()));

  ///////////////////////////////////////////////////////////////////////
  // Test serial partitioning.

  if (rank==0){

    cout << "Serial partitioning: " << 2*nprocs << " parts." << endl;
    parParams.set("num_global_parts", 2*nprocs);

    Zoltan2::PartitioningProblem<inputAdapter_t> serialProblem(
      &ia, &params, MPI_COMM_SELF);
 
    serialProblem.solve();

    const Zoltan2::PartitioningSolution<inputAdapter_t> &serialSolution = 
      serialProblem.getSolution();

    const ArrayRCP<MetricValues<scalar_t> > & metrics =
      serialSolution.getMetrics();

    Zoltan2::printMetrics<scalar_t>(cout, 2*nprocs, 2*nprocs, 2*nprocs, 
      metrics.view(0, metrics.size()));
  }

  if (rank == 0)
    std::cout << "PASS" << std::endl;
}

