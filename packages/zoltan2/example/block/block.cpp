// @HEADER
// ***********************************************************************
//                Copyright message goes here. 
// ***********************************************************************
// @HEADER

/*! \file block.cpp
    \brief An example of partitioning global ids with Block.
*/

#include <Zoltan2_TestHelpers.hpp>
#include <Zoltan2_BasicIdentifierInput.hpp>
#include <Zoltan2_Partitioning.hpp>

using namespace std;

/*! \example block.cpp
    An example of the use of the Block algorithm to partition data.
    \todo error handling
    \todo write some examples that don't use teuchos
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
  //size_t globalCount = coords->getGlobalLength();

  Teuchos::ArrayView<const gno_t> gnoList = 
    coords->getMap()->getNodeElementList();

  const gno_t *globalIds = gnoList.getRawPtr();
   
  ///////////////////////////////////////////////////////////////////////
  // Create a Zoltan2 input adapter with no weights

  typedef Zoltan2::BasicIdentifierInput<tMVector_t> inputAdapter_t;

  std::vector<const scalar_t *> noWeights;
  std::vector<int> noStrides;

  inputAdapter_t ia(localCount, globalIds, noWeights, noStrides);
   
  ///////////////////////////////////////////////////////////////////////
  // Create parameters for an Block problem

  Teuchos::ParameterList params("test params");
  params.set("debug_level", "detailed_status");
  params.set("debug_procs", "0");
  params.set("error_check_level", "debug_mode_assertions");

  Teuchos::ParameterList &parParams = params.sublist("partitioning");
  parParams.set("algorithm", "block");
  parParams.set("imbalance_tolerance", 1.1);
  parParams.set("num_global_parts", nprocs);
   
  ///////////////////////////////////////////////////////////////////////
  // Create a Zoltan2 partitioning problem

#ifdef HAVE_ZOLTAN2_MPI
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

  const ArrayRCP<Zoltan2::MetricValues<scalar_t> > & metrics1 =
    solution.getMetrics();

  if (rank == 0)
    Zoltan2::printMetrics<scalar_t>(cout, nprocs, nprocs, nprocs,
      metrics1.view(0,metrics1.size()));

   
  if (rank == 0)
    std::cout << "PASS" << std::endl;
}

