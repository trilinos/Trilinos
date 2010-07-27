#include "Tpetra_MatrixMatrix.hpp"
#include "Tpetra_MatrixIO.hpp"
#include "Kokkos_DefaultNode.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_VerboseObject.hpp"

int main(int argc, char* argv[]) {
  typedef Tpetra::CrsMatrix<double, int> MyMatrix;
  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
  Teuchos::RCP<Kokkos::DefaultNode::DefaultNodeType> node = Kokkos::DefaultNode::getDefaultNode();
  Teuchos::RCP<MyMatrix> A = Teuchos::null;
  Teuchos::RCP<MyMatrix> B = Teuchos::null;
  Teuchos::RCP<MyMatrix> C = Teuchos::null;
  Teuchos::RCP<MyMatrix> Ccheck = Teuchos::null;
  Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();
  
  Tpetra::Utils::readHBMatrix("A.hb", comm, node, A);
  Tpetra::Utils::readHBMatrix("B.hb", comm, node, B);
  Tpetra::Utils::readHBMatrix("C.hb", comm, node, Ccheck);

  C = rcp(new MyMatrix(A->getRowMap(), A->getGlobalNumRows()));


  Teuchos::RCP<const MyMatrix> Aconst = A;
  Teuchos::RCP<const MyMatrix> Bconst = B;


  Tpetra::MatrixMatrix::Multiply(Aconst, false, Bconst, false, C, true);

  //Ccheck->describe(*out, Teuchos::VERB_EXTREME);
  Teuchos::RCP<Tpetra::Map<int, int> > cTargetMap;
  if(comm->getRank() == 0){
    size_t numLocals = A->getGlobalNumRows();
    cTargetMap= Teuchos::rcp(new Tpetra::Map<int, int>(A->getGlobalNumRows(), numLocals, A->getRowMap()->getIndexBase(), A->getComm()));
  }
  else{
    cTargetMap = Teuchos::rcp(new Tpetra::Map<int, int>(A->getGlobalNumRows(), 0, A->getRowMap()->getIndexBase(), A->getComm()));
  }
  Teuchos::RCP<MyMatrix> importedC = Teuchos::rcp(new MyMatrix(cTargetMap, C->getGlobalMaxNumRowEntries()));
  Tpetra::Import<int,int> importer(C->getRowMap(), cTargetMap);
  importedC->doImport(*C, importer, Tpetra::ADD);
  importedC->fillComplete();
  importedC->describe(*out, Teuchos::VERB_EXTREME);

  	


  //C->describe(*out, Teuchos::VERB_EXTREME);

  return(0);
}
