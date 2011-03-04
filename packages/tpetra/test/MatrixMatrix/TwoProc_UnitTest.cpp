#include "Tpetra_MatrixMatrix.hpp"
#include "Tpetra_MatrixIO.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_UnitTestHarness.hpp"

namespace Tpetra{



TEUCHOS_UNIT_TEST(Tpetra_MatMat, TwoProcTest){

  typedef Map<int>                                       Map;
  typedef CrsMatrix<double,int>                          CrsMatrix;
  typedef CrsMatrix::mat_vec_type                       MatVec;
  typedef CrsMatrix::node_type                          DNode;
   
  RCP<const Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  int thisproc = comm->getRank();
  int numprocs = comm->getSize();

  out << "num procs " << numprocs << std::endl;

  if(numprocs <2){ return; }



  // set up a row-std::map with 2 global elements,
  // 1 on each proc.
  const int numGlobalRows = 2;
  ArrayRCP<int> myrow(1);
  if (thisproc == 1) myrow[0] = 7;
  else               myrow[0] = 3;
  RCP<const Map> rowmap = rcp(new Map(numGlobalRows, myrow(), 0, comm));

  //set up a domain-std::map with columns 0 - 4 on proc 0,
  //and columns 5 - 9 on proc 1.
  const int numGlobalCols = 10;
  const int numMyCols = 5;
  ArrayRCP<int> mycols(numGlobalCols);
  for(int i=0; i<numGlobalCols; ++i) {
    mycols[i] = i;
  }
  RCP<const Map> domainmap = 
    rcp(new Map(numGlobalCols, mycols(thisproc*numMyCols,numMyCols), 0, comm));

  // now create matrices A, B and C with rowmap; the second argument is just the suggested allocation size
  RCP<CrsMatrix> A = rcp(new CrsMatrix(rowmap, 1));
  A->setObjectLabel("Factor Matrix A");
  RCP<CrsMatrix> C = rcp(new CrsMatrix(rowmap, 1));
  C->setObjectLabel("Product matrix C");

  ArrayRCP<double> coefs(numGlobalCols);
  for(int i=0; i<numGlobalCols; ++i) {
    coefs[i] = 1.0*i;
  }

  A->insertGlobalValues(
    myrow[0], 
    mycols(thisproc*numMyCols, numMyCols), 
    coefs(thisproc*numMyCols, numMyCols));

  A->fillComplete(domainmap, rowmap);
  // A->describe(*out, Teuchos::VERB_EXTREME);

  MatrixMatrix::Multiply(*A, false, *A, true, *C);
  C->describe(out, Teuchos::VERB_EXTREME);


  TEST_EQUALITY(C->getGlobalNumEntries(), 2);


}


} //namespace Tpetra

