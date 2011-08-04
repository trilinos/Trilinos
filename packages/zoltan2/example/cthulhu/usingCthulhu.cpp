#include <string>
#include "Zoltan2_config.h"

// Cthulhu's build usually defines these
#define HAVE_CTHULHU_EPETRA
#define HAVE_CTHULHU_TPETRA

#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_RCP.hpp>
#include <Epetra_CrsMatrix.h>
#include <Tpetra_CrsMatrix.hpp>
#include <Cthulhu_EpetraCrsMatrix.hpp>
#include <Cthulhu_TpetraCrsMatrix.hpp>


using namespace std;

typedef int LO;
typedef long GO;
typedef double Scalar;

int globalNumRows = 80;
size_t maxRowSize = 8;

template<typename Scalar, typename LO, typename GO>
  void populate_crsmatrix(Cthulhu::CrsMatrix<Scalar, LO, GO> &matrix);
template<typename Scalar, typename LO, typename GO>
  void print_crsmatrix(Cthulhu::CrsMatrix<Scalar, LO, GO> &matrix, std::string caption);

int main(int argc, char *argv[])
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
  int rank = comm->getRank();
  int nprocs = comm->getSize();

  int base=0;

  // Feature #1:
  // Create a Cthulhu::EpetraMatrix and a Cthulhu::TpetraMatrix and
  //   perform operations on them using the same interface.

  // Epetra_CrsMatrix

  Cthulhu::EpetraMap emap(globalNumRows, base, comm);
  Cthulhu::EpetraCrsMatrix emtx(emap, maxRowSize);

  // Tpetra_CrsMatrix

  Cthulhu::TpetraMap<LO, GO> tmap(globalNumRows, base, comm);
  Cthulhu::TpetraCrsMatrix<Scalar, LO, GO> tmtx(tmap, maxRowSize);

  // Populate both matrices using the same interface.

  populate_crsmatrix<double, int, int>(emtx);
  populate_crsmatrix<Scalar, LO, GO>(tmtx);

  // Print both as a check

  print_crsmatrix<double, int, int>(emtx);
  print_crsmatrix<Scalar, LO, GO>(tmtx);


  // Feature #2:
  // Create an Epetra::CrsMatrix and a Tpetra::CrsMatrix and
  //  typecast them as Cthulhu objects and then use the
  //  same interface to work with them.

}

template<typename Scalar, typename LO, typename GO>
  void print_crsmatrix(Cthulhu::CrsMatrix<Scalar, LO, GO> &matrix, std::string caption)
{
  const Teuchos::RCP<const Cthulhu::Map<LO, GO> > map = matrix.getRowMap();
  const Teuchos::RCP<const Teuchos::Comm<int> > comm = map.getComm();

  int rank = comm->getRank();
  int nprocs = comm->getSize();

  size_t nGlobalCols = matrix.getGlobalNumCols(); 
  size_t nGlobalRows = matrix.getGlobalNumRows();
  size_t nrows = matrix.getNodeNumRows();

  if (nGlobalRows > 200){
    if (rank==0){
      std::cout << "Skip printing matrix because it has " << nGlobalRows << " rows." << std::endl;
    }
    return;
    
  }

  comm.barrier();

  GO base = map.getIndexBase();
  GO localBase = map.getMinLocalIndex();
  size_t maxEntries = matrix.getNodeMaxNumRowEntries();

  Teuchos::ArrayRCP<LO> colLID(maxEntries);
  Teuchos::ArrayRCP<Scalar> nonZero(maxEntries);

  std::string blankRow("  ",nGlobalCols);

  for (int p=0; p < nproc; p++){
    if (p==rank){
      if ((p == 0) && (caption.size() > 0)) std::cout << caption << std::endl;

      for (LO i=0; i < nrows; i++){
        size_t numnz = matrix.getNumEntriesInLocalRow(i);
        matrix.getLocalRowCopy(i, colLID.view(0, numnz), nonZero.view(0, numnz));

        std::string rowtext(blankRow);
        for (LO j=0; j < numnz; j++){
          rowtext.replace(colLID[j]*2, 2, " 1");
        }

        std::cout << std::setw(3) << localBase + i << " " << rowtext << std::endl;
      }
      std::cout.flush();
    }
    comm.barrier();
  }
}

template<typename Scalar, typename LO, typename GO>
  void populate_crsmatrix(Cthulhu::CrsMatrix<Scalar, LO, GO> &matrix)
{
  const Teuchos::RCP<const Cthulhu::Map<LO, GO> > map = matrix.getRowMap();
  const Teuchos::RCP<const Cthulhu::Map<LO, GO> > columnMap = matrix.getColumnMap();

  GO base = map.getIndexBase();
  GO from = map.getMinGlobalIndex();
  GO to = map.getMaxGlobalIndex();
  size_t globalNumCols = matrix.getGlobalNumCols();
  GO lastCol = base + globalNumCols - 1;
  LO numRows = to - from + 1;   // we know they're contiguous

  GO *column = new GO [numRows];
  Scalar *val = new Scalar[maxRowSize];

  for (int i=0; i < maxRowSize; i++){
    val[i] = 1.0;    // all nonzeros are 1.0
  }

  Teuchos::ArrayRCP<GO> colIds(column, 0, numRows, true);
  Teuchos::ArrayRCP<Scalar> nonZeros(val, 0, maxRowSize, true);

  for (GO rowId = from; rowId <= to; rowId++){
    GO fromColId = rowId - maxRowSize/2 + 1;
    GO toColId = rowId + maxRowSize/2 - 1;

    fromColId = (fromColId < base) ? base : fromColId;
    toColId = (toColId > lastCol) ? lastCol : toColId;

    GO *ids = column;
    size_t numCols = 0;
    for (GO col=fromColId; col <= toColId; col++){
      *ids++ = col;
      numCols++;
    }

    matrix.insertGlobalValues(rowId, colIds.view(0, numCols), nonZeros.view(0, numCols);
  }

  matrix.fillComplete(map, columnMap);
}
