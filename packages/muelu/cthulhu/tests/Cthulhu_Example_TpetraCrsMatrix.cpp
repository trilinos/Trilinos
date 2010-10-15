#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_Array.hpp>

#include "Tpetra_DefaultPlatform.hpp"
// #include "Tpetra_Version.hpp"
// #include "Tpetra_Map.hpp"
// #include "Tpetra_MultiVector.hpp"
// #include "Tpetra_Vector.hpp"
// #include "Tpetra_CrsMatrix.hpp"

#include "Cthulhu_DefaultPlatform.hpp"
#include "Cthulhu_TpetraCrsMatrix.hpp"

int main(int argc, char *argv[]) {
  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);

  //
  // Specify types used in this example
  // 
  typedef double Scalar;
  typedef Teuchos::ScalarTraits<Scalar>::magnitudeType Magnitude;
  typedef int Ordinal;
  typedef Cthulhu::DefaultPlatform::DefaultPlatformType           Platform;
  typedef Cthulhu::DefaultPlatform::DefaultPlatformType::NodeType Node;
  using Cthulhu::global_size_t;

  // 
  // Get the default communicator and node
  //
  Platform &platform = Cthulhu::DefaultPlatform::getDefaultPlatform();
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

  std::cout << "Comm info: " << *comm;

  return 0;
}
