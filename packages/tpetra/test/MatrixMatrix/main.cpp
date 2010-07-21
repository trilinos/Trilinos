#include "Tpetra_MatrixMatrix.hpp"
#include "Tpetra_MatrixIO.hpp"
#include "Kokkos_DefaultNode.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_VerboseObject.hpp"

int main(int argc, char* argv[]) {

  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
  Teuchos::RCP<Kokkos::DefaultNode::DefaultNodeType> node = Kokkos::DefaultNode::getDefaultNode();
  Teuchos::RCP<Tpetra::CrsMatrix<double, int> > A = Teuchos::null;
  Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();
  
  Tpetra::Utils::readHBMatrix("A2.mat", comm, node, A);

  A->describe(*out, Teuchos::VERB_EXTREME);

  return(0);
}
