// @HEADER
// ***********************************************************************
//                Copyright message goes here. 
// ***********************************************************************
// @HEADER

/*! \file rcb.cpp
    \brief An example of partitioning coordinates with RCB.
*/

#include <Zoltan2_TestHelpers.hpp>     // for UserInputForTests
#include <Zoltan2_Partitioning.hpp>
#include <Zoltan2_BasicCoordinateInput.hpp>
#include <vector>

using namespace std;
using std::vector;
using Zoltan2::MetricValues;

/*! \example rcb.cpp
    An example of the use of the RCB algorithm to partition coordinate data.
*/

int main(int argc, char *argv[])
{
  using Teuchos::RCP;
  using Teuchos::rcp;

#ifdef HAVE_ZOLTAN2_MPI                   
  MPI_Init(&argc, &argv);
  int rank, nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
  int rank=0, nprocs=1;
#endif

  ///////////////////////////////////////////////////////////////////////
  // Read coordinates from a file.

  std::string fname(testDataFilePath+"/USAir97.mtx");
  outputFlag_t flags;
  flags.set(OBJECT_COORDINATES);

  typedef Tpetra::MultiVector<scalar_t, lno_t, gno_t, node_t> tMVector_t;

  Teuchos::RCP<const Teuchos::Comm<int> > tcomm =
    Teuchos::DefaultComm<int>::getComm();
  UserInputForTests uinput(fname, tcomm, flags);

  RCP<tMVector_t> coords = uinput.getCoordinates();

  size_t localCount = coords->getLocalLength();
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
  ///////////////////////////////////////////////////////////////////////
  // A simple problem with no weights.
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////

  // Create a Zoltan2 input adapter for this geometry.

  typedef Zoltan2::BasicCoordinateInput<tMVector_t> inputAdapter_t;
  inputAdapter_t ia1(localCount, globalIds, x, y, z, 1, 1, 1);

  // Create a Zoltan2 partitioning problem

#ifdef HAVE_ZOLTAN2_MPI                   
  Zoltan2::PartitioningProblem<inputAdapter_t> problem1(&ia1, &params, 
    MPI_COMM_WORLD);
#else
  Zoltan2::PartitioningProblem<inputAdapter_t> problem1(&ia1, &params);
#endif
   
  // Solve the problem

  problem1.solve();
   
  // Obtain the solution

  const Zoltan2::PartitioningSolution<inputAdapter_t> &solution1 = 
    problem1.getSolution();
   
  // Check the solution.

  const ArrayRCP<MetricValues<scalar_t> > & metrics1 =
    solution1.getMetrics();

  if (rank == 0)
    Zoltan2::printMetrics<scalar_t>(cout, nprocs, nprocs, nprocs, 
      metrics1.view(0,metrics1.size()));

  if (rank == 0)
    if (solution1.getImbalance() < 1.1)   // object imbalance
      std::cout << "PASS" << std::endl;
    else
      std::cout << "FAIL" << std::endl;
   
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  // Try a problem with weights (1 dimension)
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////

  scalar_t *weights = new scalar_t [localCount];
  for (int i=0; i < localCount; i++){
    weights[i] = 1.0 + rank / nprocs;
  }

  // Create a Zoltan2 input adapter that includes weights.

  vector<const scalar_t *>coordVec(2);
  vector<int> coordStrides(2);

  coordVec[0] = x; coordStrides[0] = 1;
  coordVec[1] = y; coordStrides[1] = 1;

  vector<const scalar_t *>weightVec(1);
  vector<int> weightStrides(1);

  weightVec[0] = weights; weightStrides[0] = 1;

  inputAdapter_t ia2(
    localCount, globalIds,  
    coordVec, coordStrides, 
    weightVec, weightStrides);


  // Create a Zoltan2 partitioning problem

#ifdef HAVE_ZOLTAN2_MPI                   
  Zoltan2::PartitioningProblem<inputAdapter_t> problem2(
    &ia2, &params, MPI_COMM_WORLD);
#else
  Zoltan2::PartitioningProblem<inputAdapter_t> problem2(&ia2, &params);
#endif

  // Solve the problem

  problem2.solve();

  ///////////////////////////////////////////////////////////////////////
  // Obtain the solution

  const Zoltan2::PartitioningSolution<inputAdapter_t> &solution2 =
    problem2.getSolution();

  ///////////////////////////////////////////////////////////////////////
  // Check the solution.

  const ArrayRCP<MetricValues<scalar_t> > & metrics2 =
    solution2.getMetrics();

  if (rank == 0)
    Zoltan2::printMetrics<scalar_t>(cout, nprocs, nprocs, nprocs,
      metrics2.view(0,metrics2.size()));

  if (rank == 0)
    if (solution2.getImbalance() < 1.1)   // weight imbalance
      std::cout << "PASS" << std::endl;
    else
      std::cout << "FAIL" << std::endl;

  if (localCount > 0)
    delete [] weights;

  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  // Try a problem with multiple weights.
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////

#if 0
  // Add to the parameters the multicriteria objective.

  parParams.set("objective", "multicriteria_minimize_total_weight");

  // Create the new weights.

  weights = new scalar_t [localCount*3];
  srand(555);

  for (int i=0; i < localCount*3; i+=3){
    weights[i] = 1.0 + rank / nprocs;      // weight dimension 1
    weights[i+1] = rank<nprocs/2 ? 1 : 2;  // weight dimension 2
    weights[i+2] = rand() +.5;             // weight dimension 3
  }

  // Create a Zoltan2 input adapter with these weights.

  weightVec.resize(3);
  weightStrides.resize(3);

  weightVec[0] = weights;   weightStrides[0] = 3;
  weightVec[1] = weights+1; weightStrides[1] = 3;
  weightVec[2] = weights+2; weightStrides[2] = 3;

  inputAdapter_t ia3(
    localCount, globalIds,  
    coordVec, coordStrides, 
    weightVec, weightStrides);

  // Create a Zoltan2 partitioning problem.

#ifdef HAVE_ZOLTAN2_MPI                   
  Zoltan2::PartitioningProblem<inputAdapter_t> problem3(
    &ia3, &params, MPI_COMM_WORLD);
#else
  Zoltan2::PartitioningProblem<inputAdapter_t> problem3(&ia3, &params);
#endif

  // Solve the problem

  problem3.solve();

  ///////////////////////////////////////////////////////////////////////
  // Obtain the solution

  const Zoltan2::PartitioningSolution<inputAdapter_t> &solution3 =
    problem3.getSolution();

  ///////////////////////////////////////////////////////////////////////
  // Check the solution.

  const ArrayRCP<MetricValues<scalar_t> > & metrics3 =
    solution2.getMetrics();

  if (rank == 0)
    Zoltan2::printMetrics<scalar_t>(cout, nprocs, nprocs, nprocs,
      metrics3.view(0,metrics3.size()));

  if (rank == 0)
    if (solution3.getImbalance() < 1.1) 
      std::cout << "PASS" << std::endl;
    else
      std::cout << "FAIL" << std::endl;

  if (localCount)
    delete [] weights;

  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  // Using part sizes, ask for some parts to be empty.
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////

  // Change the number of parts to twice the number of processes to
  // ensure that we have more than one global part.

  parParams.set("num_global_parts", nprocs*2);
  parParams.set("num_local_parts", 2);

  // Using the initial problem that did not have any weights, reset
  // parameter list, and give it some part sizes.

  problem1.resetParameterList(&params);

  zoltan2_partId_t *partIds = new zoltan2_partId_t [2];
  scalar_t *partSizes = new scalar_t [2];

  partIds[0] = rank*2;    partSizes[0] = 0;
  partIds[1] = rank*2+1;  partSizes[1] = 1;

  problem1.setPartSizes(2, partIds, partSizes);

  // Solve the problem.  The flag "false" indicates that we
  // have not changed the input data, which allows the problem
  // so skip some work when re-solving.

  problem1.solve(false);

  // Obtain the solution

  const Zoltan2::PartitioningSolution<inputAdapter_t> &solution4 =
    problem1.getSolution();

  // Check the solution.

  const ArrayRCP<MetricValues<scalar_t> > & metrics4 =
    solution4.getMetrics();

  if (rank == 0)
    Zoltan2::printMetrics<scalar_t>(cout, nprocs, nprocs, nprocs,
      metrics4.view(0,metrics4.size()));

  if (rank == 0)
    if (solution4.getImbalance() < 1.1)   // object imbalance
      std::cout << "PASS" << std::endl;
    else
      std::cout << "FAIL" << std::endl;

  delete [] partIds;
  delete [] partSizes;
#endif
}

