#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_TestForException.hpp>

#include <Kokkos_TBBNode.hpp>

// I/O for Harwell-Boeing files
#include "Tpetra_MatrixIO.hpp"

#include "Tpetra_Version.hpp"
#include "CRSTiming.hpp"

int main(int argc, char *argv[]) {
  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv);

  //
  // Get example parameters from command-line processor
  //  
  int numThreads = -1;
  std::string filename("bcsstk14.hb");
  int verbose = 1;
  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("num-threads",&numThreads,"Number of threads.");
  cmdp.setOption("verbose",&verbose,"Verbose (zero for silent).");
  cmdp.setOption("filename",&filename,"Filename for Harwell-Boeing test matrix.");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }

  // 
  // Say hello, print some communicator info
  //
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Teuchos::createMpiComm<int>(Teuchos::opaqueWrapper<MPI_Comm>(MPI_COMM_WORLD));
  if (comm->getRank() == 0) {
    std::cout << "\n" << Tpetra::version() << std::endl << std::endl;
    std::cout << argv[0] << filename << std::endl;
    std::cout << "Comm info: " << *comm;
  }

  typedef Kokkos::TBBNode Node;
  Teuchos::ParameterList params;
  params.set<int>("Num Threads",numThreads);
  params.set<int>("Verbose",verbose);
  Teuchos::RCP<Node> node = Teuchos::rcp(new Node(params));

  if (comm->getRank() == 0) {
    typedef Kokkos::DefaultKernels<double,int,Node>::SparseOps DSM;
#ifndef HAVE_KOKKOS_NO_FIRST_TOUCH_MATVEC_ALLOCATION
    std::cout << "Using Kokkos first-touch matrix objects." << std::endl;
    Kokkos::CrsMatrix<double,int,Node,DSM> *mat;
    // this will fail to compile if the above macro doesn't appropriately correspond with the reality of the CrsMatrix inheritance
    TEST_FOR_EXCEPT( (static_cast<Kokkos::FirstTouchHostCrsMatrix<double,int,Node,DSM> *>(mat) != 0) );
#else
    std::cout << "Not using Kokkos first-touch matrix objects." << std::endl;
    // this will fail to compile if the above macro doesn't appropriately correspond with the reality of the CrsMatrix inheritance
    TEST_FOR_EXCEPT( (static_cast<   Kokkos::CrsMatrixHostCompute<double,int,Node,DSM> *>(mat) != 0) );
#endif
  }

  //
  // Read Tpetra::CrsMatrix from file
  //
  Teuchos::RCP< Tpetra::CrsMatrix<double,int,int,Node> > A;
  Tpetra::Utils::readHBMatrix(filename,comm,node,A);
  if (verbose) {
    std::cout << std::endl << A->description() << std::endl << std::endl;
  }

  CRSTiming<Node>(A);

  return 0;
}
