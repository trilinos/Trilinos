// @HEADER
// ***********************************************************************
//                Copyright message goes here. 
// ***********************************************************************
// @HEADER

/*! \file rcb.cpp
    \brief An example of partitioning coordinates with RCB.
    \todo add more cases to this test.
*/

#include <Zoltan2_TestHelpers.hpp>
#include <Zoltan2_BasicCoordinateInput.hpp>
#include <Zoltan2_XpetraMultiVectorInput.hpp>
#include <Zoltan2_PartitioningSolution.hpp>
#include <Zoltan2_PartitioningProblem.hpp>

using namespace std;
using Teuchos::RCP;
using Teuchos::rcp;

typedef Tpetra::MultiVector<scalar_t, lno_t, gno_t, node_t> tMVector_t;
typedef Zoltan2::BasicUserTypes<scalar_t, gno_t, lno_t, gno_t> myTypes_t;


/*! \test rcbTest.cpp
    An example of the use of the RCB algorithm to partition coordinate data.
    \todo error handling
    \todo write some examples that don't use teuchos
    \todo check the solution, visualize it somehow
*/

void testFromDataFile(const RCP<const Teuchos::Comm<int> > & comm)
{
  std::string fname("simple");
  UserInputForTests uinput(testDataFilePath, fname, comm, true);

  RCP<tMVector_t> coords = uinput.getCoordinates();

  RCP<const tMVector_t> coordsConst = rcp_const_cast<const tMVector_t>(coords);

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
   
#if 0
  typedef Zoltan2::BasicCoordinateInput<tMVector_t> inputAdapter_t;
  inputAdapter_t ia(localCount, globalIds, x, y, z, 1, 1, 1);
#else
  typedef Zoltan2::XpetraMultiVectorInput<tMVector_t> inputAdapter_t;
  inputAdapter_t ia(coordsConst);
#endif
   
  Teuchos::ParameterList params("test params");
  Teuchos::ParameterList &parParams = params.sublist("partitioning");
  Teuchos::ParameterList &geoParams = parParams.sublist("geometric");
  geoParams.set("bisection_num_test_cuts", 7);

#ifdef HAVE_ZOLTAN2_MPI                   
  Zoltan2::PartitioningProblem<inputAdapter_t> problem(&ia, &params,
    MPI_COMM_WORLD);
#else
  Zoltan2::PartitioningProblem<inputAdapter_t> problem(&ia, &params);
#endif

  problem.solve();

  const Zoltan2::PartitioningSolution<inputAdapter_t> &solution = 
    problem.getSolution();

  if (comm->getRank() == 0)
    problem.printMetrics(cout);
}

void serialTest()
{
  int numParts = 8;
  int numCoords = 1000;

  cout << "Serial partitioning: " << numParts << " parts." << endl;

  gno_t *ids = new gno_t [numCoords];
  if (!ids)
    throw std::bad_alloc();
  for (lno_t i=0; i < numCoords; i++)
    ids[i] = i;
  ArrayRCP<gno_t> globalIds(ids, 0, numCoords, true);

  Array<ArrayRCP<scalar_t> > randomCoords(3);
  UserInputForTests::getRandomData(555, numCoords, 0, 10, 
    randomCoords.view(0,3));

  typedef Zoltan2::BasicCoordinateInput<myTypes_t> inputAdapter_t;

  inputAdapter_t ia(numCoords, ids, 
    randomCoords[0].getRawPtr(), randomCoords[1].getRawPtr(),
     randomCoords[2].getRawPtr(), 1,1,1);

  Teuchos::ParameterList params("test params");
  params.set("debug_level", "basic_status");

  Teuchos::ParameterList &parParams = params.sublist("partitioning");
  parParams.set("num_global_parts", numParts);
  parParams.set("algorithm", "rcb");
  parParams.set("imbalance_tolerance", 1.1);

  Teuchos::ParameterList &geoParams = parParams.sublist("geometric");
  geoParams.set("bisection_num_test_cuts", 7);

#ifdef HAVE_ZOLTAN2_MPI                   
  Zoltan2::PartitioningProblem<inputAdapter_t> serialProblem(
    &ia, &params, MPI_COMM_SELF);
#else
  Zoltan2::PartitioningProblem<inputAdapter_t> serialProblem(&ia, &params);
#endif
 
  serialProblem.solve();

  const Zoltan2::PartitioningSolution<inputAdapter_t> &serialSolution = 
    serialProblem.getSolution();

  serialProblem.printMetrics(cout);
}

void meshCoordinatesTest(const RCP<const Teuchos::Comm<int> > & comm)
{
  int xdim = 40;
  int ydim = 60;
  int zdim = 20;
  UserInputForTests uinput(xdim, ydim, zdim, string("Laplace3D"), comm, true);

  RCP<tMVector_t> coords = uinput.getCoordinates();

  size_t localCount = coords->getLocalLength();

  scalar_t *x=NULL, *y=NULL, *z=NULL;
  x = coords->getDataNonConst(0).getRawPtr();
  y = coords->getDataNonConst(1).getRawPtr();
  z = coords->getDataNonConst(2).getRawPtr();

  const gno_t *globalIds = coords->getMap()->getNodeElementList().getRawPtr();
  typedef Zoltan2::BasicCoordinateInput<tMVector_t> inputAdapter_t;

  inputAdapter_t ia(localCount, globalIds, x, y, z, 1, 1, 1);

  Teuchos::ParameterList params("test params");
  Teuchos::ParameterList &parParams = params.sublist("partitioning");
  Teuchos::ParameterList &geoParams = parParams.sublist("geometric");
  geoParams.set("bisection_num_test_cuts", 7);
  geoParams.set("rectilinear_blocks", "yes");

#ifdef HAVE_ZOLTAN2_MPI
  Zoltan2::PartitioningProblem<inputAdapter_t> problem(&ia, &params,
    MPI_COMM_WORLD);
#else
  Zoltan2::PartitioningProblem<inputAdapter_t> problem(&ia, &params);
#endif

  problem.solve();

  const Zoltan2::PartitioningSolution<inputAdapter_t> &solution =
    problem.getSolution();

  if (comm->getRank()  == 0)
    problem.printMetrics(cout);
}

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv);
  RCP<const Teuchos::Comm<int> > tcomm = Teuchos::DefaultComm<int>::getComm();
  int rank = tcomm->getRank();

  meshCoordinatesTest(tcomm);

  testFromDataFile(tcomm);

  if (rank == 0)
    serialTest();

  if (rank == 0)
    std::cout << "PASS" << std::endl;
}
