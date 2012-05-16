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
typedef Tpetra::Map<lno_t, gno_t, node_t> tMap_t;
typedef Zoltan2::BasicCoordinateInput<tMVector_t> inputAdapter_t;

UserInputForTests *uinput=NULL;

const RCP<const tMap_t> & getMeshCoordinates(
    const RCP<const Teuchos::Comm<int> > & comm,
    int xdim, int ydim, int zdim,
    ArrayView<ArrayRCP<const scalar_t> > coords)
{
  int dim = 3;
  if (zdim == 0){
    dim--;
    if (ydim  == 0)
      dim--;
  }

  if (dim == 3)
    uinput = new UserInputForTests(xdim, ydim, zdim, 
    string("Laplace3D"), comm, true);
  else if (dim == 2)
    uinput = new UserInputForTests(xdim, ydim, zdim, 
    string("Laplace2D"), comm, true);
  else
    uinput = new UserInputForTests(xdim, ydim, zdim,
    string("Identity"), comm, true);

  RCP<tMVector_t> meshCoords = uinput->getCoordinates();

  for (int i=0; i < dim; i++)
    coords[i] = meshCoords->getData(i);

  return meshCoords->getMap();
}


int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv);
  RCP<const Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
  int rank = comm->getRank();

  // For now the only command line option is the global number
  // of coordinates.
  // Eventually when our parameters are xml-ized we will read
  // in the Zoltan2 parameters and the test parameters from an
  // xml file.

  gno_t numGlobalCoords = 1000;
  if (argc > 1)
    numGlobalCoords = atol(argv[1]);

  double k = log(numGlobalCoords) / 3;
  double xdimf = exp(k) + 0.5;
  int xdim = static_cast<int>(floor(xdimf));
  int ydim = xdim;
  int zdim = numGlobalCoords / (xdim*ydim);

  // Obtain coordinates.

  Array<ArrayRCP<const scalar_t> > coordinates(3);

  const RCP<const tMap_t> &map = getMeshCoordinates(comm,
    xdim, ydim, zdim, coordinates.view(0,3));

  // Create an input adapter.

  const gno_t *globalIds = map->getNodeElementList().getRawPtr();
  size_t localCount = coordinates[0].size();

  inputAdapter_t ia(localCount, globalIds, 
    coordinates[0].getRawPtr(),
    coordinates[1].getRawPtr(),
    coordinates[2].getRawPtr(),
    1,1,1);

  // Parameters

  Teuchos::ParameterList params;
  params.set("timing_output_stream" , "std::cout");

  Teuchos::ParameterList &parParams = params.sublist("partitioning");
  parParams.set("algorithm", "rcb");
  double tolerance = 1.1;
  parParams.set("imbalance_tolerance", tolerance );

  Teuchos::ParameterList &geoParams = parParams.sublist("geometric");
  geoParams.set("bisection_num_test_cuts", 20);
  geoParams.set("rectilinear_blocks", "no");
  geoParams.set("average_cuts", "no");

  // Create a problem, solve it, and display the quality.

#ifdef HAVE_ZOLTAN2_MPI
  Zoltan2::PartitioningProblem<inputAdapter_t> problem(&ia, &params,
    MPI_COMM_WORLD);
#else
  Zoltan2::PartitioningProblem<inputAdapter_t> problem(&ia, &params);
#endif

  problem.solve();

  if (uinput){
    delete uinput;
    uinput = NULL;
  }

  comm->barrier();

  problem.printTimers();

  comm->barrier();

  problem.getSolution().printMetrics(std::cout);

  if (rank == 0)
    std::cout << "PASS" << std::endl;
}
