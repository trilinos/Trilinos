
#include "Tpetra_Core.hpp"
#include "Tpetra_CrsMatrix.hpp"

int main(int narg, char** arg)
{
  Tpetra::ScopeGuard ts(&narg, &arg);
  auto comm = Tpetra::getDefaultComm();
  int me = comm->getRank(); 
  int np = comm->getSize();
  int npM1 = np - 1;
  
  // This test requires at least two processors
  if (np == 1) {
    std::cout << "PASSED" << std::endl;
    return 0;
  }

  int nrows = 3 * np * (np - 1);
  int myRowsNp = nrows / np;
  int myRowsNpM1 = (me < npM1 ? (nrows/npM1) : 0);
  
  using map_t = Tpetra::Map<>;
  using gno_t = typename map_t::global_ordinal_type;
  using lno_t = typename map_t::local_ordinal_type;
  using matrix_t = Tpetra::CrsMatrix<gno_t>;
  using vector_t = Tpetra::Vector<gno_t>;

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
  
  Teuchos::Array<gno_t> cols(5);
  Teuchos::Array<gno_t> vals(5);
  for (int i = 0; i < myRowsNp; i++) {
    gno_t gid = mapNp->getGlobalElement(i);
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

  // Same matrix across np-1 processors:  
  // insert nonzeros on np processors; 
  // let fillComplete migrate according to mapNpM1

  // Get exact nnz per row
  vector_t countNp(mapNp);
  auto dataNp = countNp.getDataNonConst();
  for (size_t i = 0; i < mapNp->getNodeNumElements(); i++) {
    dataNp[i] = Amat->getNumEntriesInLocalRow(i);
  }
  vector_t countNpM1(mapNpM1);
  Tpetra::Import<lno_t, gno_t> importer(mapNp, mapNpM1);
  countNpM1.doImport(countNp, importer, Tpetra::INSERT);

  Teuchos::Array<size_t> nnzNpM1(mapNpM1->getNodeNumElements());
  auto dataNpM1 = countNpM1.getData();
  for (size_t i = 0; i < mapNpM1->getNodeNumElements(); i++) 
    nnzNpM1[i] = dataNpM1[i];
  
  // Create Redistributed matrix
  Teuchos::RCP<matrix_t> Rmat = rcp(new matrix_t(mapNpM1, nnzNpM1()));

  for (int i = 0; i < myRowsNp; i++) {
    gno_t gid = mapNp->getGlobalElement(i);
    for (int j = 0; j < 5; j++) vals[j] = gid;
    int nz = 0;
    cols[nz++] = gid;
    if (gid+1 < nrows) cols[nz++] = gid+1;
    if (gid+2 < nrows) cols[nz++] = gid+2;
    if (gid-1 >= 0) cols[nz++] = gid-1;
    if (gid-2 >= 0) cols[nz++] = gid-2;
    Rmat->insertGlobalValues(gid, cols(0,nz), vals(0,nz));
  }
  Rmat->fillComplete();
  std::cout << me << " RMAT" << std::endl;
  Rmat->describe(foo, Teuchos::VERB_EXTREME);

  // Test Rmat for correctness
  int ierr = 0;
  for (size_t i = 0; i < Rmat->getRowMap()->getNodeNumElements(); i++) {
    gno_t gid = Rmat->getRowMap()->getGlobalElement(i);
    Teuchos::ArrayView<const lno_t> rcols;
    Teuchos::ArrayView<const gno_t> rvals;
    Rmat->getLocalRowView(i, rcols, rvals);

    // Check each nonzero for existence and correct value
    for (gno_t j = gid-2; j <= gid+2; j++) {
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

  int gerr;
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &ierr, &gerr);
  if (gerr) std::cout << "FAILED with " << gerr << " errors" << std::endl;
  else std::cout << "PASSED" << std::endl;

  return gerr;
}
