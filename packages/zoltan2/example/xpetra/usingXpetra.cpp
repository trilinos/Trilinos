#include <string>
#include "Zoltan2_config.h"
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_RCP.hpp>
#include <Epetra_MpiComm.h>
#include <Epetra_SerialComm.h>
#include <Epetra_CrsMatrix.h>
#include <Tpetra_CrsMatrix.hpp>
#include <MueLu_config.hpp>   // For HAVE_MUELU_TPETRA HAVE_MUELU_EPETRA
#include <Xpetra_EpetraCrsMatrix.hpp>
#include <Xpetra_TpetraCrsMatrix.hpp>

using namespace std;

typedef int LO;
typedef long GO;
typedef double Scalar;

GO globalNumRows = 80;
LO maxRowSize = 8;

template<typename Scalar, typename LO, typename GO>
  void print_crsmatrix(const Xpetra::CrsMatrix<Scalar, LO, GO> &matrix, std::string caption)
{
  const Teuchos::RCP<const Xpetra::Map<LO, GO> > map = matrix.getRowMap();
  const Teuchos::RCP<const Xpetra::Map<LO, GO> > colmap = matrix.getColMap();
  const Teuchos::RCP<const Teuchos::Comm<int> > comm = map->getComm();

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

  comm->barrier();

  size_t maxEntries = matrix.getNodeMaxNumRowEntries();

  Teuchos::ArrayRCP<LO> colLID(maxEntries);
  Teuchos::ArrayRCP<Scalar> nonZero(maxEntries);

  std::string blankRow(nGlobalCols,' ');

  for (int p=0; p < nprocs; p++){
    if (p==rank){
      if ((p == 0) && (caption.size() > 0)) std::cout << caption << std::endl;

      for (size_t i=0; i < nrows; i++){
        size_t numnz = matrix.getNumEntriesInLocalRow(i);
        GO rowGID = map->getGlobalElement(i);

        matrix.getLocalRowCopy(i, colLID.view(0, numnz), nonZero.view(0, numnz), numnz);

        std::string rowtext(blankRow);
        for (size_t j=0; j < numnz; j++){
          rowtext.replace(colmap->getGlobalElement(colLID[j]), 1, "1");
        }

        std::cout << std::setw(3) << rowGID << " " << rowtext.c_str() << std::endl;
      }
      std::cout.flush();
    }
    comm->barrier();
  }
}

template<typename Scalar, typename LO, typename GO>
  void populate_crsmatrix(Xpetra::CrsMatrix<Scalar, LO, GO> &matrix)
{
  const Teuchos::RCP<const Xpetra::Map<LO, GO> > map = matrix.getRowMap();

  size_t nrows = map->getNodeNumElements();
  GO base = map->getIndexBase();
  GO final = base + map->getGlobalNumElements() - 1;
  LO halfRow = maxRowSize/2;
  Teuchos::ArrayRCP<GO> colIds(maxRowSize);
  Teuchos::ArrayRCP<Scalar> nonZeros(maxRowSize);

  for (LO i=0; i < maxRowSize; i++){
    nonZeros[i] = 1;
  }
  size_t len;
  for (size_t i=0; i < nrows; i++){
    GO rowGID = map->getGlobalElement(i);
    len = 0;
    for (GO j=rowGID-halfRow+1; j < rowGID+halfRow; j++){
      if ((j < base) || (j > final)) continue;
      colIds[len++] = j;
    }
    
    matrix.insertGlobalValues(rowGID, colIds.view(0, len), nonZeros.view(0, len));
  }

  matrix.fillComplete();
}

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  Teuchos::GlobalMPISession session(&argc, &argv);
  Epetra_MpiComm commObj(MPI_COMM_WORLD);
#else
  Epetra_SerialComm commObj();
#endif


  Teuchos::RCP<const Teuchos::Comm<int> > comm =
    Teuchos::DefaultComm<int>::getComm();

  int base=0;

  // Feature #1:
  // Create a Xpetra::EpetraCrsMatrix and a Xpetra::TpetraCrsMatrix and
  //   perform operations on them using the same interface.

  // Epetra_CrsMatrix

  Teuchos::RCP<const Xpetra::EpetraMap>  emap = 
    Teuchos::rcp(new Xpetra::EpetraMap(globalNumRows, base, comm));

  Xpetra::EpetraCrsMatrix emtx(emap, maxRowSize);

  // Tpetra_CrsMatrix

  Teuchos::RCP<const Xpetra::TpetraMap<LO, GO> >  tmap = 
    Teuchos::rcp(new Xpetra::TpetraMap<LO, GO>(globalNumRows, base, comm));

  Xpetra::TpetraCrsMatrix<Scalar, LO, GO> tmtx(tmap, maxRowSize);

  // Populate both matrices using the same interface.

  populate_crsmatrix<double, int, int>(emtx);
  populate_crsmatrix<Scalar, LO, GO>(tmtx);

  // Print both as a check

  print_crsmatrix<double, int, int>(emtx, "Xpetra::EpetraCrsMatrix");
  print_crsmatrix<Scalar, LO, GO>(tmtx, "Xpetra::TpetraCrsMatrix");

  // Feature #2:
  // Create an Epetra::CrsMatrix and a Tpetra::CrsMatrix and
  //  interact with them as const Xpetra objects.

  int rank = comm->getRank();

  Teuchos::ArrayRCP<double> double_values(10);  
  Teuchos::ArrayRCP<Scalar> scalar_values(10);  
  Teuchos::ArrayRCP<int> int_indices(10);  
  Teuchos::ArrayRCP<GO> go_indices(10);  

  for (int j=0; j < 10; j++){
     double_values[j] = double(1.0);
     scalar_values[j] = Scalar(1.0);
  }

  int nrows = (rank+1) * 2 + 10;
  int nGlobalRows=0;

  Teuchos::reduceAll<int, int>(*comm, Teuchos::REDUCE_SUM, 1, &nrows, &nGlobalRows);

  Epetra_Map map(nGlobalRows, nrows, base, commObj);
  Teuchos::RCP<Epetra_CrsMatrix> emtxPtr = 
    Teuchos::rcp<Epetra_CrsMatrix>(new Epetra_CrsMatrix(Copy, map, 10));

  Teuchos::RCP<Tpetra::Map<LO, GO> > tmapPtr = 
    Teuchos::rcp<Tpetra::Map<LO, GO> >(
      new Tpetra::Map<LO, GO>(nGlobalRows, nrows, base, comm));
  Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LO, GO> > tmtxPtr = 
    Teuchos::rcp<Tpetra::CrsMatrix<Scalar, LO, GO> >(
      new Tpetra::CrsMatrix<Scalar, LO, GO>(tmapPtr, 10));

  for (int i=0; i < nrows; i++){
    int first = map.GID(i)-4;
    int last = map.GID(i)+6;
    first = (first < base) ? base : first;
    last = (last > base+nGlobalRows) ? base+nGlobalRows : last;
    int numnz = last - first;      
    for (int j=0; j < numnz; j++,first++){
      int_indices[j] = first;
      go_indices[j] = GO(first);
    }
    emtxPtr->InsertGlobalValues(map.GID(i), numnz, double_values.get(), int_indices.get());

    tmtxPtr->insertGlobalValues(tmapPtr->getGlobalElement(i), 
      go_indices.view(0,numnz), scalar_values.view(0,numnz));
  }

  emtxPtr->FillComplete();
  tmtxPtr->fillComplete();

  const Xpetra::EpetraCrsMatrix constEMtx(emtxPtr);
  const Xpetra::TpetraCrsMatrix<Scalar, LO, GO>  constTMtx(tmtxPtr);

  print_crsmatrix(constEMtx, "The Epetra_CrsMatrix derived object");
  print_crsmatrix(constTMtx, "The Tpetra::CrsMatrix derived object");

  comm->barrier();

  std::cout << "PASS" << std::endl;

  return 0;
}

