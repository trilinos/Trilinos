#include "Kokkos_DummySparseKernelClass.hpp"
#include <Kokkos_MultiVector.hpp>
#include <Kokkos_DefaultNode.hpp>

/** \file DummySparseKernelDriver.cpp
    \brief A file testing the build of the DummySparseKernel class and illustrating its usage.
 */

int main() {

  typedef Kokkos::DefaultNode::DefaultNodeType                 Node;
  typedef KokkosExamples::DummySparseKernel<Node>         SparseOps;
  typedef Kokkos::CrsGraph <       int,Node,SparseOps>        Graph;
  typedef Kokkos::CrsMatrix<double,int,Node,SparseOps>    DoubleMat;
  typedef Kokkos::CrsMatrix< float,int,Node,SparseOps>     FloatMat;
  typedef Kokkos::MultiVector<double,Node>                DoubleVec;
  typedef Kokkos::MultiVector<float,Node>                  FloatVec;

  std::cout << "Note, this class doesn't actually do anything. We are only testing that it compiles." << std::endl;

  // get a pointer to the default node
  Teuchos::RCP<Node> node = Kokkos::DefaultNode::getDefaultNode();

  // create the graph G
  const size_t numRows = 5;
  Graph G(numRows,node);

  // create a double-valued matrix dM using the graph G
  DoubleMat dM(G);
  // create a double-valued sparse kernel using the rebind functionality
  SparseOps::rebind<double>::other doubleKernel(node);
  // initialize it with G and dM
  doubleKernel.initializeStructure(G);
  doubleKernel.initializeValues(dM);
  // create double-valued vectors and initialize them
  DoubleVec dx(node), dy(node);
  // test the sparse kernel operator interfaces
  doubleKernel.multiply( Teuchos::NO_TRANS, 1.0, dx, dy);
  doubleKernel.multiply( Teuchos::NO_TRANS, 1.0, dx, 1.0, dy);
  doubleKernel.solve( Teuchos::NO_TRANS, Teuchos::UPPER_TRI, Teuchos::UNIT_DIAG, dy, dx);

  // create a float-valued matrix fM using the graph G
  FloatMat fM(G);
  // create a double-valued sparse kernel using the rebind functionality
  SparseOps::rebind<float>::other floatKernel(node);
  // initialize it with G and fM
  floatKernel.initializeStructure(G);
  floatKernel.initializeValues(fM);
  // create float-valued vectors and initialize them
  FloatVec fx(node), fy(node);
  // test the sparse kernel operator interfaces
  floatKernel.multiply( Teuchos::NO_TRANS, 1.0f, fx, fy);
  floatKernel.multiply( Teuchos::NO_TRANS, 1.0f, fx, 1.0f, fy);
  floatKernel.solve( Teuchos::NO_TRANS, Teuchos::UPPER_TRI, Teuchos::UNIT_DIAG, fy, fx);

  std::cout << "End Result: TEST PASSED" << std::endl;
  return 0;
}
