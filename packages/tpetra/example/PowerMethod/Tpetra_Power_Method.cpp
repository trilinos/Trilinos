#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_Array.hpp>

#include "Tpetra_Power_Method.hpp"

#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Version.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_CrsMatrix.hpp"

int main(int argc, char *argv[]) {
  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);

  //
  // Specify types used in this example
  // 
  typedef double Scalar;
  typedef Teuchos::ScalarTraits<Scalar>::magnitudeType Magnitude;
  typedef int Ordinal;
  typedef Tpetra::DefaultPlatform::DefaultPlatformType           Platform;
  typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType Node;
  using Tpetra::global_size_t;

  // 
  // Get the default communicator and node
  //
  Platform &platform = Tpetra::DefaultPlatform::getDefaultPlatform();
  Teuchos::RCP<const Teuchos::Comm<int> > comm = platform.getComm();
  Teuchos::RCP<Node>             node = platform.getNode();
  const int myRank = comm->getRank();

  //
  // Get example parameters from command-line processor
  //  
  bool printMatrix = false;
  bool verbose = (myRank==0);
  int niters = 100;
  int numGlobalElements = 100;
  Magnitude tolerance = 1.0e-2;
  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("numGlobalElements",&numGlobalElements,"Global problem size.");
  cmdp.setOption("tolerance",&tolerance,"Relative residual tolerance used for solver.");
  cmdp.setOption("iterations",&niters,"Maximum number of iterations.");
  cmdp.setOption("printMatrix","noPrintMatrix",&printMatrix,"Print the full matrix after reading it.");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }

  // 
  // Say hello, print some communicator info
  //
  if (verbose) {
    std::cout << "\n" << Tpetra::version() << std::endl << std::endl;
  }
  std::cout << "Comm info: " << *comm;

  //
  // Construct the problem
  //
  // Construct a Map that puts approximately the same number of equations on each processor.
  Teuchos::RCP<const Tpetra::Map<Ordinal> > map = Tpetra::createUniformContigMap<Ordinal,Ordinal>(numGlobalElements, comm);
  // Get update list and number of local equations from newly created map.
  const size_t numMyElements = map->getNodeNumElements();
  Teuchos::ArrayView<const Ordinal> myGlobalElements = map->getNodeElementList();
  // Create an OTeger vector NumNz that is used to build the Petra Matrix.
  // NumNz[i] is the Number of OFF-DIAGONAL term for the ith global equation 
  // on this processor
  Teuchos::ArrayRCP<size_t> NumNz = Teuchos::arcp<size_t>(numMyElements);
  // We are building a tridiagonal matrix where each row has (-1 2 -1)
  // So we need 2 off-diagonal terms (except for the first and last equation)
  for (size_t i=0; i < numMyElements; ++i) {
    if (myGlobalElements[i] == 0 || myGlobalElements[i] == numGlobalElements-1) {
      // boundary
      NumNz[i] = 2;
    }
    else {
      NumNz[i] = 3;
    }
  }
  // Create a Tpetra::Matrix using the Map, with a static allocation dictated by NumNz
  Teuchos::RCP< Tpetra::CrsMatrix<Scalar,Ordinal> > A;
  A = Teuchos::rcp( new Tpetra::CrsMatrix<Scalar,Ordinal>(map, NumNz, Tpetra::StaticProfile) );
  // We are done with NumNZ
  NumNz = Teuchos::null;
  // Add rows one-at-a-time
  // Off diagonal values will always be -1
  const Scalar two    = static_cast<Scalar>( 2.0);
  const Scalar negOne = static_cast<Scalar>(-1.0);
  for (size_t i=0; i<numMyElements; i++) {
    if (myGlobalElements[i] == 0) {
      A->insertGlobalValues( myGlobalElements[i],
                             Teuchos::tuple<Ordinal>( myGlobalElements[i], myGlobalElements[i]+1 ),
                             Teuchos::tuple<Scalar> ( two, negOne ) );
    }
    else if (myGlobalElements[i] == numGlobalElements-1) {
      A->insertGlobalValues( myGlobalElements[i],
                             Teuchos::tuple<Ordinal>( myGlobalElements[i]-1, myGlobalElements[i] ),
                             Teuchos::tuple<Scalar> ( negOne, two ) );
    }
    else {
      A->insertGlobalValues( myGlobalElements[i],
                             Teuchos::tuple<Ordinal>( myGlobalElements[i]-1, myGlobalElements[i], myGlobalElements[i]+1 ),
                             Teuchos::tuple<Scalar> ( negOne, two, negOne ) );
    }
  }
  // Finish up
  A->fillComplete(Tpetra::DoOptimizeStorage);
  if (printMatrix) {
    Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    A->describe(*fos, Teuchos::VERB_EXTREME);
  }
  else if (verbose) {
    std::cout << std::endl << A->description() << std::endl << std::endl;
  }

  //
  // Iterate
  //
  Scalar lambda; (void)lambda;
  lambda = TpetraExamples::powerMethod<Scalar,Ordinal>(A, niters, tolerance, verbose);

  // Increase diagonal dominance
  if (verbose) {
    std::cout << "\nIncreasing magnitude of first diagonal term, solving again\n"
              << std::endl;
  }

  //A->resumeFill();
  if (A->getRowMap()->isNodeGlobalElement(0)) {
    // get a copy of the row with with global index 0
    // modify the diagonal entry of that row
    // submit the modified values to the matrix
    const Ordinal ID = 0;
    size_t numVals = A->getNumEntriesInGlobalRow(ID);
    Teuchos::Array<Scalar>  rowvals(numVals);
    Teuchos::Array<Ordinal> rowinds(numVals);
    A->getGlobalRowCopy(ID, rowinds, rowvals, numVals);       // Get A(0,:)
    for (size_t i=0; i<numVals; i++) {
      if (rowinds[i] == ID) {
        // we have found the diagonal; modify it and break the loop
        rowvals[i] *= 10.0;
        break;
      }
    }
    A->replaceGlobalValues(ID, rowinds(), rowvals());
  }
  //A->fillComplete();

  // Iterate (again)
  lambda = TpetraExamples::powerMethod<Scalar,Ordinal>(A, niters, tolerance, verbose);  

  if (verbose) {
    std::cout << "\nEnd Result: TEST PASSED" << std::endl;
  }
  return 0;
}
