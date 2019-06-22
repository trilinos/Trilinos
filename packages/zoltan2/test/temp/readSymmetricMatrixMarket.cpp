
#include "Tpetra_Core.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "MatrixMarket_Tpetra.hpp"

int main(int narg, char **arg)
{

  Tpetra::ScopeGuard scope(&narg, &arg);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  std::cout << "Usage:  a.out generalMatrixMarketFile symmetricMatrixMarketFile"
            << std::endl;

  typedef Tpetra::CrsMatrix<> tcrsMatrix_t;
  typedef Tpetra::MatrixMarket::Reader<tcrsMatrix_t> reader_t;

  // read general matrix
  std::cout << comm->getRank() << "Reading general matrix " << std::endl;
  char *filename = arg[1];
  Teuchos::RCP<tcrsMatrix_t> genMatrix = 
           reader_t::readSparseFile(filename, comm, true, false, true);
  std::cout << comm->getRank() << "SUCCESS Read general matrix " << std::endl;


  // read symmetric matrix
  std::cout << comm->getRank() << "Reading general matrix " << std::endl;
  filename = arg[2];
  Teuchos::RCP<tcrsMatrix_t> symMatrix = 
           reader_t::readSparseFile(filename, comm, true, false, true);
  std::cout << comm->getRank() << "SUCCESS Read symmetric matrix " << std::endl;


  // write symmetric matrix to general file for inspection
  std::cout << comm->getRank() << "Writing symmetric matrix to general file " 
            << std::endl;
  std::string outputFile("out.mtx");
  Tpetra::MatrixMarket::Writer<tcrsMatrix_t>::writeSparseFile(outputFile,
                                                              symMatrix, true);
  std::cout << comm->getRank() << "SUCCESS Wrote symmetric matrix" 
            << std::endl;

  // check some metrics
  if (comm->getRank() == 0) {
    int nfail = 0;
    if (symMatrix->getGlobalNumRows() != genMatrix->getGlobalNumRows()) {
      std::cout << "FAIL:  globalNumRows " << symMatrix->getGlobalNumRows()
                << " != " << genMatrix->getGlobalNumRows() << std::endl;
      nfail++;
    }
    if (symMatrix->getGlobalNumCols() != genMatrix->getGlobalNumCols()) {
      std::cout << "FAIL:  globalNumCols " << symMatrix->getGlobalNumCols()
                << " != " << genMatrix->getGlobalNumCols() << std::endl;
      nfail++;
    }
    if (symMatrix->getGlobalNumEntries() != genMatrix->getGlobalNumEntries()) {
      std::cout << "FAIL:  globalNumEntries " 
                << symMatrix->getGlobalNumEntries() << " != " 
                << genMatrix->getGlobalNumEntries() << std::endl;
      nfail++;
    }

    if (nfail == 0) std::cout << "PASS" << std::endl;
  }

  return 0;
}
