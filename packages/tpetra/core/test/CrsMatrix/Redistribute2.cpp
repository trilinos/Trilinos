
#include "Tpetra_Core.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include <vector>

class TestRedistribute 
{
public:
  using map_t = Tpetra::Map<>;
  using GO = typename map_t::global_ordinal_type;
  using LO = typename map_t::local_ordinal_type;
  using matrix_t = Tpetra::CrsMatrix<GO>;
  using vector_t = Tpetra::Vector<GO>;

  TestRedistribute() { }

  int run() {

    auto comm = Tpetra::getDefaultComm();
    int me = comm->getRank(); 
    int np = comm->getSize();
    int npM1 = np - 1;

    int nrows = 3 * np * (np - 1);
    int myRowsNp = nrows / np;
    int myRowsNpM1 = (me < npM1 ? (nrows/npM1) : 0);
  
    Tpetra::global_size_t dummy = 
            Teuchos::OrdinalTraits<Tpetra::global_size_t>:: invalid();

    // Map with rows across all processors
    Teuchos::RCP<const map_t> mapNp = 
             rcp(new map_t(dummy, myRowsNp, 0, comm));

    // Map with rows across np-1 processors
    Teuchos::RCP<const map_t> mapNpM1 = 
             rcp(new map_t(dummy, myRowsNpM1, 0, comm));

    // Matrix across all processors
    Teuchos::RCP<matrix_t> Amat = rcp(new matrix_t(mapNp, 5));
    
    Teuchos::Array<GO> cols(5);
    Teuchos::Array<GO> vals(5);
    for (int i = 0; i < myRowsNp; i++) {
      GO gid = mapNp->getGlobalElement(i);
      for (int j = 0; j < 5; j++) vals[j] = gid;
      int nz = 0;
      cols[nz++] = gid;
      if (gid+1 < nrows) cols[nz++] = gid+1;
      if (gid+2 < nrows) cols[nz++] = gid+2;
      if (gid-1 >= 0) cols[nz++] = gid-1;
      if (gid-2 >= 0) cols[nz++] = gid-2;
      Amat->insertGlobalValues(gid, cols(0,nz), vals(0,nz));
    }
    Amat->fillComplete();
  
    Teuchos::FancyOStream foo(Teuchos::rcp(&std::cout,false));
    std::cout << me << " AMAT" << std::endl;
    Amat->describe(foo, Teuchos::VERB_EXTREME);

    // Build redistributed matrix 

    Teuchos::RCP<matrix_t> Rmat; 
    redistributeMatrix(myRowsNpM1, Amat, Rmat);

    std::cout << me << " RMAT" << std::endl;
    Rmat->describe(foo, Teuchos::VERB_EXTREME);

    // Test Rmat for correctness
    int ierr = 0;
    for (size_t i = 0; i < Rmat->getRowMap()->getNodeNumElements(); i++) {
      GO gid = Rmat->getRowMap()->getGlobalElement(i);
      Teuchos::ArrayView<const LO> rcols;
      Teuchos::ArrayView<const GO> rvals;
      Rmat->getLocalRowView(i, rcols, rvals);
  
      // Check each nonzero for existence and correct value
      for (GO j = gid-2; j <= gid+2; j++) {
        if (j >= 0 && j < nrows) {
          bool foundCol = false;
          for (int k = 0; k < rcols.size(); k++) {
            if (j == Rmat->getColMap()->getGlobalElement(rcols[k])) {
              foundCol = true;
              if (rvals[k] != gid) {
                ierr++;
                std::cout << me << " error Rmat[" << gid << "," << j << "] != "
                          << gid << std::endl;
              }
            }
          }
          if (!foundCol) {
            ierr++;
            std::cout << me << " nonzero Rmat[" << gid << "," << j << "] "
                      << "not found" << std::endl;
          }
        }
      }
    }
    return ierr;
  }

private:

  void determineCompressedGIDs(Teuchos::RCP<const matrix_t> A,
                               GO &baseGID,
                               std::vector<GO> & compressedColGIDs)
  {
    Teuchos::RCP<const map_t> rowMap = A->getRowMap();
    Teuchos::RCP<const map_t> colMap = A->getColMap();
    Teuchos::RCP<const Teuchos::Comm<int> > TComm = A->getComm();

    GO numMyRows = rowMap->getNodeNumElements();
    GO numRowsScanSum;
    Teuchos::scan<int, GO> (*TComm, Teuchos::REDUCE_SUM, 1, &numMyRows,
                            &numRowsScanSum);
    baseGID = numRowsScanSum - numMyRows;

    std::vector<GO> globalIDs(numMyRows);
    for (GO i=0; i<numMyRows; i++) globalIDs[i] = baseGID + i;

    vector_t rowVec(rowMap, Teuchos::ArrayView<GO>(globalIDs));
    vector_t colVec(colMap);

    Tpetra::Import<LO,GO> importer(rowMap, colMap);
    colVec.doImport(rowVec, importer, Tpetra::INSERT);

    GO numCols = colMap->getNodeNumElements();
    compressedColGIDs.resize(numCols);

    Teuchos::ArrayRCP<const GO> colGIDs = colVec.getData();
    for (GO i=0; i<numCols; i++) {
        compressedColGIDs[i] = colGIDs[i];
    }
  }

  void redistributeMatrix(GO m_numRowsMe,
                          Teuchos::RCP<const matrix_t> A,
                          Teuchos::RCP<matrix_t> &redistributedA)
  {
    Tpetra::global_size_t IGO = 
            Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();

    Teuchos::RCP<const map_t> rowMap = A->getRowMap();

    Teuchos::RCP<const map_t> rowMapSource = 
        rcp( new map_t(IGO, rowMap->getNodeNumElements(), 0, A->getComm()));

    Teuchos::RCP<const map_t> rowMapTarget = 
        rcp( new map_t(IGO, m_numRowsMe, 0, A->getComm()) );

    Teuchos::RCP<vector_t> m_rhs = rcp( new vector_t(rowMapSource, 1) );
    Teuchos::RCP<vector_t> m_rhsGrid = rcp( new vector_t(rowMapTarget, 1) );
    Teuchos::RCP<Tpetra::Import<LO,GO> > m_importer = 
            rcp( new Tpetra::Import<LO,GO>(rowMapSource, rowMapTarget) );

    Tpetra::Vector<int> countSource(rowMapSource), countTarget(rowMapTarget);
    auto countValues = countSource.getDataNonConst();
    for (size_t i=0; i < rowMapSource->getNodeNumElements(); i++) {
      countValues[i] = A->getNumEntriesInLocalRow(i);
    }
    countTarget.doImport(countSource, *m_importer, Tpetra::INSERT);

    Teuchos::Array<size_t> nnzTarget(rowMapTarget->getNodeNumElements());
    auto dataTarget = countTarget.getData();
    for (size_t i = 0; i < rowMapTarget->getNodeNumElements(); i++)
      nnzTarget[i] = dataTarget[i];

    GO baseGID;
    std::vector<GO> compressedColGIDs;
    this->determineCompressedGIDs(A, baseGID, compressedColGIDs);

    redistributedA = rcp(new matrix_t(rowMapTarget, nnzTarget()));
    Teuchos::RCP<const map_t> colMap = A->getColMap();

    std::vector<GO> globalCols;
    for(size_t i=0; i < rowMap->getNodeNumElements(); i++) {

      Teuchos::ArrayView<const LO> Indices;
      Teuchos::ArrayView<const GO> Values;
      A->getLocalRowView(i, Indices, Values);

      const Teuchos::Ordinal numIndices = Indices.size();
      globalCols.resize(numIndices);

      for(auto j = 0; j < numIndices; j++) {
        globalCols[j] = compressedColGIDs[Indices[j]];
      }

      GO globalRow = rowMapSource->getGlobalElement(i);
      redistributedA->insertGlobalValues(
                                  globalRow,
                                  Teuchos::ArrayView<const GO>(globalCols),
                                  Values);
    }
    redistributedA->fillComplete();
  }
};

int main(int narg, char** arg)
{
  Tpetra::ScopeGuard ts(&narg, &arg);
  auto comm = Tpetra::getDefaultComm();
  int np = comm->getSize();
  
  // This test requires at least two processors
  if (np == 1) {
    std::cout << "This test requires at least two processors" << std::endl;
    std::cout << "TEST PASSED" << std::endl;
    return 0;
  }

  TestRedistribute test;
  int ierr = test.run();

  int gerr;
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &ierr, &gerr);
  if (gerr) std::cout << "TEST FAILED with " << gerr << " errors" << std::endl;
  else std::cout << "TEST PASSED" << std::endl;

  return gerr;
}
