#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>

#include <examples/Kokkos_DummySparseKernelClass.hpp>
#include <Kokkos_DefaultNode.hpp>

#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_CrsMatrix.hpp"

/** \file LocalMatOpsExample.cpp
    \brief A file testing the build of the KokkosExamples::DummySparseKernel and illustrating a custom sparse mat-vec with Tpetra::CrsMatrix.
 */

int main(int argc, char *argv[]) {
  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);

  typedef Tpetra::DefaultPlatform::DefaultPlatformType           Platform;
  typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType Node;
  typedef KokkosExamples::DummySparseKernel<Node>                SparseOps;
  typedef Tpetra::Map<int,int,Node>                              Map;
  typedef Tpetra::CrsMatrix<float,int,int,Node,SparseOps>        Matrix;
  typedef Tpetra::MultiVector<float,int,int,Node>                MultiVector;

  // 
  // Get the default communicator and node
  //
  Platform &platform = Tpetra::DefaultPlatform::getDefaultPlatform();
  Teuchos::RCP<const Teuchos::Comm<int> > comm = platform.getComm();
  Teuchos::RCP<Node>             node = platform.getNode();

  std::cout << "Note, this class doesn't actually do anything. We are only testing that it compiles." << std::endl;

  // create the matrix with a custom sparse mat-vec kernel
  const size_t numRows = 5;
  const Tpetra::global_size_t numGlobalRows = numRows*comm->getSize();
  Teuchos::RCP<const Map> map = Tpetra::createUniformContigMapWithNode<int,int,Node>(numGlobalRows, comm, node);
  Teuchos::RCP<Matrix> matrix = Teuchos::rcp( new Matrix(map,1,Tpetra::DynamicProfile) );
  matrix->fillComplete();

  Teuchos::RCP<MultiVector> X = Tpetra::createMultiVector<float>(map, 2),
                            Y = Tpetra::createMultiVector<float>(map, 2);
  matrix->apply(*X, *Y);

  std::cout << "End Result: TEST PASSED" << std::endl;
  return 0;
}
