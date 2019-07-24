
#include "Tpetra_Core.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "MatrixMarket_Tpetra.hpp"

typedef Tpetra::Map<> map_t;
typedef Tpetra::CrsMatrix<> tcrsMatrix_t;
typedef typename tcrsMatrix_t::global_ordinal_type gno_t;
typedef typename tcrsMatrix_t::local_ordinal_type lno_t;
typedef typename tcrsMatrix_t::scalar_type scalar_t;

typedef Tpetra::MatrixMarket::Reader<tcrsMatrix_t> reader_t;

void printMetrics(const Teuchos::RCP<tcrsMatrix_t> &mat, const char *str)
{
  size_t ninternal = 0;
  size_t nboundary = 0;

  size_t numRows = mat->getNodeNumRows();
  const Teuchos::RCP<const map_t> rowMap = mat->getRowMap();
  
  for (size_t i = 0; i < numRows; i++) {
    Teuchos::ArrayView<const lno_t> indices;
    Teuchos::ArrayView<const scalar_t> values;
    gno_t gno = rowMap->getGlobalElement(i);
    mat->getLocalRowView(i, indices, values);
    size_t nnz = indices.size();
    bool internal = true;
    for (size_t j = 0; j < nnz; j++) {
      gno_t jidx = rowMap->getGlobalElement(indices[j]);
      if (jidx == Teuchos::OrdinalTraits<lno_t>::invalid()) {
        internal = false;
        break;
      }
    }
    if (internal) ninternal++;
    else nboundary++;
  }

  std::cout << rowMap->getComm()->getRank() << " " << str
            << " Number of rows " << numRows << " of " << rowMap->getGlobalNumElements() << "\n"
            << rowMap->getComm()->getRank() << " " << str
            << " Number of nonzeros " << mat->getNodeNumEntries() << "\n"
            << rowMap->getComm()->getRank() << " " << str
            << " Number of internal " << ninternal << "\n" 
            << rowMap->getComm()->getRank() << " " << str
            << " Number of boundary " << nboundary  << "\n" 
            << std::endl;
}


int main(int narg, char **arg)
{

  Tpetra::ScopeGuard scope(&narg, &arg);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  std::cout << "Usage:  a.out generalMatrixMarketFile symmetricMatrixMarketFile"
            << std::endl;


  // read general matrix
  std::cout << comm->getRank() << "Reading general matrix " << std::endl;
  char *filename = arg[1];
  Teuchos::RCP<tcrsMatrix_t> genMatrix = 
           reader_t::readSparseFile(filename, comm, true, false, true);
  std::cout << comm->getRank() << "SUCCESS Read general matrix " << std::endl;

  printMetrics(genMatrix, "general");

  // read symmetric matrix
  std::cout << comm->getRank() << "Reading general matrix " << std::endl;
  filename = arg[2];
  Teuchos::RCP<tcrsMatrix_t> symMatrix = 
           reader_t::readSparseFile(filename, comm, true, false, true);
  std::cout << comm->getRank() << "SUCCESS Read symmetric matrix " << std::endl;

  printMetrics(symMatrix, "symmetric");

#ifdef WRITE_IT_OUT
  // write symmetric matrix to general file for inspection
  std::cout << comm->getRank() << "Writing symmetric matrix to general file " 
            << std::endl;
  std::string outputFile("out.mtx");
  Tpetra::MatrixMarket::Writer<tcrsMatrix_t>::writeSparseFile(outputFile,
                                                              symMatrix, true);
  std::cout << comm->getRank() << "SUCCESS Wrote symmetric matrix" 
            << std::endl;
#endif

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
