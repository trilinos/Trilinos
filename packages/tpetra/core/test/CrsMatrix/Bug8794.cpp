
#include "Tpetra_Core.hpp"
#include "Tpetra_CrsMatrix.hpp"

int main(int narg, char** arg)
{
  Tpetra::ScopeGuard ts(&narg, &arg);

  using map_t = Tpetra::Map<>;
  using gno_t = typename map_t::global_ordinal_type;
  using lno_t = typename map_t::local_ordinal_type;
  using node_t = typename map_t::node_type;
  using scalar_t = double;
  using matrix_t = Tpetra::CrsMatrix<scalar_t>;
  using vector_t = Tpetra::Vector<scalar_t>;

  auto comm = Tpetra::getDefaultComm();
  int me = comm->getRank();
  int np = comm->getSize();

  int nrows = 10001;
  int divisor = 100;
  int maxNzPerRow = nrows / divisor + 1;

  // Map with rows across all processors
  Teuchos::RCP<const map_t> map = rcp(new map_t(nrows, 0, comm));

  // Vectors for SpMV and expected values
  vector_t expected(map);
  vector_t x(map);
  vector_t y(map);

  // Build matrix distributed across np-1 processors:  
  // insert nonzeros on np processors; 
  // let fillComplete migrate according to mapNpM1

  Teuchos::Array<gno_t> cols(maxNzPerRow);
  Teuchos::Array<scalar_t> vals(maxNzPerRow, 1.);

  matrix_t Amat(map, maxNzPerRow);

  // Initialize matrix and expected value of SpMV product
  {
    expected.putScalar(0.);
    expected.clear_sync_state();
    expected.modify_host();
    auto expectedData = expected.getDataNonConst();
    for (size_t i = 0; i < map->getNodeNumElements(); i++) {

      gno_t gid = map->getGlobalElement(i);
      int nz = 0;

      if (gid % (divisor+1) == 1) {  // dense row
        for (int j = 0; j < nrows; j+=divisor) {
          gno_t tmp = (gid + j) % nrows;
          cols[nz++] = tmp;
std::cout << me << " of " << np << " (i,j) " << gid << " " << tmp << std::endl;
          expectedData[i] += tmp;
        }
      }
      else {  // sparse row
        cols[nz++] = gid;
        expectedData[i] += scalar_t(gid);
        if (gid+1 < nrows) { cols[nz++] = gid+1; expectedData[i] += gid+1; }
        if (gid+2 < nrows) { cols[nz++] = gid+2; expectedData[i] += gid+2; }
        if (gid-1 >= 0)    { cols[nz++] = gid-1; expectedData[i] += gid-1; }
        if (gid-2 >= 0)    { cols[nz++] = gid-2; expectedData[i] += gid-2; }
      }
  
      Amat.insertGlobalValues(gid, cols(0,nz), vals(0,nz));
    }
  
    Amat.fillComplete();
    std::cout << me << " of " << np << ": \n"
              << "  nrows     " << Amat.getNodeNumRows() << "\n"
              << "  nnz       " << Amat.getNodeNumEntries() << "\n"
              << "  maxPerRow " << Amat.getNodeMaxNumRowEntries() << "\n"
              << std::endl;

    Teuchos::FancyOStream foo(Teuchos::rcp(&std::cout,false));
    Amat.describe(foo, Teuchos::VERB_EXTREME);
  }

  // Initialize domain vector for SpMV
  {
    x.clear_sync_state();
    x.modify_host();
    auto xData = x.getDataNonConst();
    for (int i = 0; i < map->getNodeNumElements(); i++) 
      xData[i] = map->getGlobalElement(i);
  }

  Amat.apply(x, y);

  // Test product for correctness
  int ierr = 0;
  {
    expected.sync_host();
    y.sync_host();
    auto expectedData = expected.getData();
    auto yData = y.getData();

    for (size_t i = 0; i < map->getNodeNumElements(); i++) {
      if (yData[i] != expectedData[i]) {
        std::cout << me << " of " << np << ": y[" << map->getGlobalElement(i)
                  << "] " << yData[i] << " != " << expectedData[i]
                  << " expected" << std::endl;
        ierr++;
      }
    }
  }

  int gerr;
  Teuchos::reduceAll<int,int>(*comm, Teuchos::REDUCE_SUM, 1, &ierr, &gerr);
  if (gerr) std::cout << "TEST FAILED with " << gerr << " errors" << std::endl;
  else std::cout << "TEST PASSED" << std::endl;

  return gerr;
}
