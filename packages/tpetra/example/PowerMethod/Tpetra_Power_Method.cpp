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

#include <algorithm>
#include <functional>

int main(int argc, char *argv[]) {
  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);

  //
  // Specify types used in this example
  // 
  typedef double                                                  Scalar;
  typedef Teuchos::ScalarTraits<Scalar>::magnitudeType            Magnitude;
  typedef int                                                     Ordinal;
  typedef Tpetra::DefaultPlatform::DefaultPlatformType            Platform;
  typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType  Node;
  typedef Tpetra::Map<Ordinal,Ordinal,Node>                       Map;
  typedef Tpetra::CrsMatrix<Scalar,Ordinal,Ordinal,Node>          CrsMatrix;
  using Teuchos::RCP;
  using Teuchos::tuple;

  // 
  // Get the default communicator and node
  //
  Platform &platform = Tpetra::DefaultPlatform::getDefaultPlatform();
  RCP<const Teuchos::Comm<int> > comm = platform.getComm();
  RCP<Node>                      node = platform.getNode();
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
  RCP<const Map> map = Tpetra::createUniformContigMap<Ordinal,Ordinal>(numGlobalElements, comm);
  // Get update list and number of local equations from newly created map.
  const size_t numMyElements = map->getNodeNumElements();
  Teuchos::ArrayView<const Ordinal> myGlobalElements = map->getNodeElementList();
  // Create a CrsMatrix using the map, with a dynamic allocation of 3 entries per row
  RCP<CrsMatrix> A = Tpetra::createCrsMatrix<Scalar>(map,3);
  // Add rows one-at-a-time
  for (size_t i=0; i<numMyElements; i++) {
    if (myGlobalElements[i] == 0) {
      A->insertGlobalValues( myGlobalElements[i],
                             tuple<Ordinal>( myGlobalElements[i], myGlobalElements[i]+1 ),
                             tuple<Scalar> ( 2.0, -1.0 ) );
    }
    else if (myGlobalElements[i] == numGlobalElements-1) {
      A->insertGlobalValues( myGlobalElements[i],
                             tuple<Ordinal>( myGlobalElements[i]-1, myGlobalElements[i] ),
                             tuple<Scalar> ( -1.0, 2.0 ) );
    }
    else {
      A->insertGlobalValues( myGlobalElements[i],
                             tuple<Ordinal>( myGlobalElements[i]-1, myGlobalElements[i], myGlobalElements[i]+1 ),
                             tuple<Scalar> ( -1.0, 2.0, -1.0 ) );
    }
  }
  // Complete the fill, ask that storage be reallocated and optimized
  A->fillComplete(Tpetra::DoOptimizeStorage);
  if (printMatrix) {
    RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    A->describe(*fos, Teuchos::VERB_EXTREME);
  }
  else if (verbose) {
    std::cout << std::endl << A->description() << std::endl << std::endl;
  }

  //
  // Iterate
  //
  TpetraExamples::powerMethod<Scalar,Ordinal>(A, niters, tolerance, verbose);

  // Increase diagonal dominance
  if (verbose) {
    std::cout << "\nIncreasing magnitude of first diagonal term, solving again\n"
              << std::endl;
  }

  A->resumeFill();
  if (A->getRowMap()->isNodeGlobalElement(0)) {
    // get a copy of the row with with global index 0
    // modify the diagonal entry of that row
    // submit the modified values to the matrix
    size_t numVals = A->getNumEntriesInGlobalRow(0);
    Teuchos::Array<Scalar>  rowvals(numVals);
    Teuchos::Array<Ordinal> rowinds(numVals);
    A->getGlobalRowCopy(0, rowinds, rowvals, numVals);
    // replace the diagonal with 10.0
    std::replace_if( rowvals.begin(), rowvals.end(), 
                     std::bind2nd(std::equal_to<Ordinal>(), 0), 
                     10.0 
                   );
    A->replaceGlobalValues(0, rowinds(), rowvals());
  }
  A->fillComplete();

  // Iterate again
  TpetraExamples::powerMethod<Scalar,Ordinal>(A, niters, tolerance, verbose);  

  if (verbose) {
    std::cout << "\nEnd Result: TEST PASSED" << std::endl;
  }
  return 0;
}
