#include "Epetra_SerialComm.h"
#ifdef HAVE_MPI
#  include "Epetra_MpiComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_as.hpp"


namespace {


//
// Helper code and declarations
//

using Teuchos::as;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::Array;


const int g_localDim = 4; // ToDo: Make variable!


RCP<const Epetra_Comm> getEpetraComm()
{
#ifdef HAVE_MPI
  return rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
  return rcp(new Epetra_SerialComm());
#endif
}


RCP<Epetra_CrsMatrix> getEpetraMatrix(int numRows, int numCols, double shift=0.0) 
{

  const RCP<const Epetra_Comm> comm = getEpetraComm();

  const Epetra_Map rowMap(numRows, 0, *comm);
  const Epetra_Map domainMap(numCols, numCols, 0, *comm);
 
  const RCP<Epetra_CrsMatrix> epetraCrsM =
    rcp(new Epetra_CrsMatrix(Copy, rowMap, numCols));

  Array<double> rowEntries(numCols);
  Array<int> columnIndices(numCols);
  for (int j = 0; j < numCols; ++j) {
    columnIndices[j] = j;
  }

  const int numLocalRows = rowMap.NumMyElements();

  for (int i = 0; i < numLocalRows; ++i) {
    
    for (int j = 0; j < numCols; ++j) {
      rowEntries[j] = as<double>(i+1) + as<double>(j+1) / 10 + shift;
    }

    epetraCrsM->InsertMyValues( i, numCols, &rowEntries[0], &columnIndices[0] );

  }

  epetraCrsM->FillComplete(domainMap, rowMap);

  return epetraCrsM;
}


} // namespace
