#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_CommandLineProcessor.hpp>

#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_CrsMatrix.hpp>

// I/O for Harwell-Boeing files
#define HIDE_TPETRA_INOUT_IMPLEMENTATIONS
#include <Tpetra_MatrixIO.hpp>

#include <Amesos2.hpp>
#include <Amesos2_Version.hpp>


int main(int argc, char *argv[]) {
  Teuchos::GlobalMPISession mpiSession(&argc,&argv);
  typedef double Scalar;
  typedef Teuchos::ScalarTraits<Scalar>::magnitudeType Magnitude;
  typedef int Ordinal;

  typedef double Scalar;
  typedef int LO;
  typedef int GO;
  typedef Tpetra::DefaultPlatform::DefaultPlatformType           Platform;
  typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType Node;

  typedef Tpetra::CrsMatrix<Scalar,LO,GO,Node> MAT;
  typedef Tpetra::MultiVector<Scalar,LO,GO,Node> MV;

  using Tpetra::global_size_t;
  using Teuchos::tuple;
  using Teuchos::RCP;
  using Teuchos::rcp;


  // 
  // Get the default communicator
  //
  Platform &platform = Tpetra::DefaultPlatform::getDefaultPlatform();
  Teuchos::RCP<const Teuchos::Comm<int> > comm = platform.getComm();
  Teuchos::RCP<Node>             node = platform.getNode();
  int myRank  = comm->getRank();

  Teuchos::oblackholestream blackhole;
  std::ostream &out = ( myRank == 0 ? std::cout : blackhole );
  RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));

  bool printMatrix   = false;
  bool printSolution = false;
  bool printTiming   = false;
  bool verbose = (myRank==0);
  std::string filename("bcsstk14.hb");
  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("filename",&filename,"Filename for Harwell-Boeing test matrix.");
  cmdp.setOption("print_matrix","no_print_matrix",&printMatrix,"Print the full matrix after reading it.");
  cmdp.setOption("print_solution","no_print_solution",&printSolution,"Print solution vector after solve.");
  cmdp.setOption("print_timing","no_print_timing",&printTiming,"Print solver timing statistics");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }

  // Say hello
  out << Amesos::version() << std::endl << std::endl;

  const size_t numVectors = 1;

//  RCP<MAT> A = rcp( new MAT(map,3) ); // max of three entries in a row

  RCP<MAT> A;
  Tpetra::Utils::readHBMatrix(filename,comm,node,A);
  if (printMatrix) {
    A->describe(*fos, Teuchos::VERB_EXTREME);
  }
  else if (verbose) {
    std::cout << std::endl << A->description() << std::endl << std::endl;
  }

  // create a Map
  global_size_t nrows = A->getGlobalNumRows();
  RCP<Tpetra::Map<LO,GO,Node> > map = rcp( new Tpetra::Map<LO,GO,Node>(nrows,0,comm) );

  // Create random X
  RCP<MV> X = rcp(new MV(map,numVectors));
  X->randomize();

  /* Create B
   *
   * Use RHS:
   *
   *  [[10]
   *   [10]
   *   [10]
   *   [10]
   *   [10]
   *   [10]]
   */
  RCP<MV> B = rcp(new MV(map,numVectors));
  B->putScalar(10);

  // Constructor from Factory
  RCP<Amesos::SolverBase> solver = Amesos::Factory<MAT,MV>::create("Superlu",A,X,B);

  solver->symbolicFactorization().numericFactorization().solve();

  // change one of the matrix values and re-solve.
  //
  // Replace the lowest column index and lowest row index entry with "20"
  A->replaceGlobalValues(
    Teuchos::as<GO>(A->getRowMap()->getMinAllGlobalIndex()),
    tuple<GO>(A->getColMap()->getMinAllGlobalIndex()),
    tuple<Scalar>(20));

  solver->numericFactorization().solve();

  // change the RHS vector and re-solve.
  B->replaceGlobalValue(numVectors/3,0,7);
  B->replaceGlobalValue(numVectors/2,0,15);

  solver->solve();

  if( printSolution ){
    // Print the solution
    X->describe(*fos,Teuchos::VERB_EXTREME);
  }

  if( printTiming ){
    // Print some timing statistics
    solver->printTiming(*fos);
  }

  // We are done.
  return 0;
}
