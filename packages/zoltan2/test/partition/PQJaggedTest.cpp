// @HEADER
// ***********************************************************************
//                Copyright message goes here. 
// ***********************************************************************
// @HEADER

/*! \file PQJagged.cpp
    \brief An example of partitioning coordinates with PQJagged.
    \todo add more cases to this test.
*/

#include <Zoltan2_TestHelpers.hpp>
#include <Zoltan2_BasicCoordinateInput.hpp>
#include <Zoltan2_XpetraMultiVectorInput.hpp>
#include <Zoltan2_PartitioningSolution.hpp>
#include <Zoltan2_PartitioningProblem.hpp>
#include <string>
using namespace std;
using Teuchos::RCP;
using Teuchos::rcp;

typedef Tpetra::MultiVector<scalar_t, lno_t, gno_t, node_t> tMVector_t;
typedef Zoltan2::BasicUserTypes<scalar_t, gno_t, lno_t, gno_t> myTypes_t;


/*! \test PQJaggedTest.cpp
    An example of the use of the PQJagged algorithm to partition coordinate data.
    \todo error handling
    \todo write some examples that don't use teuchos
    \todo check the solution, visualize it somehow
*/

void testFromDataFile(const RCP<const Teuchos::Comm<int> > & comm, int numParts, float imbalance, std::string fname)
{
  //std::string fname("simple");
	cout << "running " << fname << endl;
  UserInputForTests uinput(testDataFilePath, fname, comm, true);

  cout << "n" << endl;
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

  parParams.set("num_global_parts", numParts);
  parParams.set("algorithm", "PQJagged");
  parParams.set("imbalance_tolerance", double(imbalance));

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

  //if (comm->getRank() == 0)
  //  solution.printMetrics(cout);

  cout << "testFromDataFile is done " << endl;

}

void serialTest(int numParts, int numCoords, float imbalance)
{
  //int numCoords = 1000;


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
  parParams.set("algorithm", "PQJagged");
  parParams.set("imbalance_tolerance", double(imbalance));

  //string algorithmss("dd");
  //bool isSets;
  //getParameterValue(parParams, "partitioning", "algorithm", isSets, algorithmss);
  //cout << "algo:" << algorithmss << endl;

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

  //serialSolution.printMetrics(cout);


}

void meshCoordinatesTest(const RCP<const Teuchos::Comm<int> > & comm)
{
  int xdim = 80;
  int ydim = 60;
  int zdim = 40;

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
  parParams.set("algorithm", "PQJagged");

  //parParams.set("algorithm", "rcb");
  Teuchos::ParameterList &geoParams = parParams.sublist("geometric");
  geoParams.set("bisection_num_test_cuts", 7);
  geoParams.set("rectilinear_blocks", "yes");

  parParams.set("num_global_parts", 10);

#ifdef HAVE_ZOLTAN2_MPI

  Zoltan2::PartitioningProblem<inputAdapter_t> problem(&ia, &params,
    MPI_COMM_WORLD);
#else
  Zoltan2::PartitioningProblem<inputAdapter_t> problem(&ia, &params);
#endif

  problem.solve();


  const Zoltan2::PartitioningSolution<inputAdapter_t> &solution =
    problem.getSolution();

  //if (comm->getRank()  == 0)
  //  solution.printMetrics(cout);
}



void meshCoordinatesTest2(const RCP<const Teuchos::Comm<int> > & comm, string pqParts, int numCoords, float imbalance)
{


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

	  int np = comm->getSize();
	  inputAdapter_t ia(numCoords/np, ids,
	    randomCoords[0].getRawPtr(), randomCoords[1].getRawPtr(),
	     randomCoords[2].getRawPtr(), 1,1,1);

  Teuchos::ParameterList params("test params");
  Teuchos::ParameterList &parParams = params.sublist("partitioning");
  parParams.set("algorithm", "PQJagged");

  params.set("pqParts", pqParts);
  //parParams.set("algorithm", "rcb");
  Teuchos::ParameterList &geoParams = parParams.sublist("geometric");
  geoParams.set("bisection_num_test_cuts", 7);
  geoParams.set("rectilinear_blocks", "yes");

  //parParams.set("num_global_parts", numParts);
  parParams.set("imbalance_tolerance", double(imbalance));

#ifdef HAVE_ZOLTAN2_MPI

  Zoltan2::PartitioningProblem<inputAdapter_t> problem(&ia, &params,
    MPI_COMM_WORLD);
#else
  Zoltan2::PartitioningProblem<inputAdapter_t> problem(&ia, &params);
#endif

  problem.solve();


  const Zoltan2::PartitioningSolution<inputAdapter_t> &solution =
    problem.getSolution();

  //if (comm->getRank()  == 0)
  //  solution.printMetrics(cout);
}


string convert_to_string(char *args){
    string tmp = "";
    for(int i = 0; args[i] != 0; i++)
        tmp += args[i];
    return tmp;
}
bool getArgumentValue(string &argumentid, float &argumentValue, string argumentline){
    stringstream stream(stringstream::in | stringstream::out);
    stream << argumentline;
    getline(stream, argumentid, '=');
    if (stream.eof()){
        return false;
    }
    stream >> argumentValue;
    return true;
}



int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv);
  RCP<const Teuchos::Comm<int> > tcomm = Teuchos::DefaultComm<int>::getComm();
  int rank = tcomm->getRank();


  int numParts = 10; float imbalance = 0.03;
  int numCoords = 1000;
  string pqParts = "";
  int opt = 0;
  std::string fname = "simple";
  for(int i = 0; i < argc; ++i){
	  string tmp = convert_to_string(argv[i]);
	  string identifier = "";
	  int value = -1; float fval = -1;
	  if(!getArgumentValue(identifier, fval, tmp)) continue;
	  value = int (fval);

	  if(identifier == "P"){
		  /*
		  if(value > 0){
			  numParts=value;
		  } else {
			  cout << "Invalid argument at " << tmp;
			  exit(1);
		  }
		  */
		  stringstream stream(stringstream::in | stringstream::out);

		  stream << tmp;

		  getline(stream, fname, '=');
		  stream >> pqParts;
	  } else if(identifier == "D"){
		  if(value > 0){
			  numCoords=value;
		  } else {
			  cout << "Invalid argument at " << tmp;
			  exit(1);
		  }
	  }else if(identifier == "I"){
		  if(fval > 0){
			  imbalance=fval;
			  cout << imbalance << endl;
		  } else {
			  cout << "Invalid argument at " << tmp;
			  exit(1);
		  }
	  } else if(identifier == "F"){
		  stringstream stream(stringstream::in | stringstream::out);
		  stream << tmp;
		  getline(stream, fname, '=');

		  stream >> fname;
	  }
	  else if(identifier == "O"){
	  		  if(value > 0 && value < 3){
	  			  opt = value;
	  		  } else {
	  			  cout << "Invalid argument at " << tmp;
	  			  exit(1);
	  		  }
	  	  }
	  else {
		  cout << "Invalid argument at " << tmp;
		  exit(1);

	  }

  }

  //MPI_Init(&argc, &argv);
  //int myrank; int numprocs;
  //MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  //MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  //meshCoordinatesTest2(tcomm,numParts, numCoords, imbalance);
  //meshCoordinatesTest(tcomm);

  switch (opt){
  case 0:
	  testFromDataFile(tcomm,numParts, imbalance,fname);
	  break;
  case 1:
	  meshCoordinatesTest2(tcomm,pqParts, numCoords, imbalance);
	  break;
  case 2:
	  serialTest(numParts, numCoords, imbalance);
	  break;
  }
  //if (rank == 0)
  //  serialTest(numParts, numCoords, imbalance);

  if (rank == 0)
    std::cout << "PASS" << std::endl;
}
