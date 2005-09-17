#ifndef GALERI_RIS_H
#define GALERI_RIS_H

#include "Galeri_Exception.h"
#include "Epetra_Comm.h"
#include "Epetra_BlockMap.h"
#include "Epetra_CrsMatrix.h"

namespace Galeri {
namespace Matrices {

inline
Epetra_CrsMatrix* Ris(const Epetra_Map* Map)
{
  // this is actually a dense matrix, stored into Crs format
  int NumGlobalElements = Map->NumGlobalElements();
  int NumMyElements     = Map->NumMyElements();
  int* MyGlobalElements = Map->MyGlobalElements();

  Epetra_CrsMatrix* Matrix = new Epetra_CrsMatrix(Copy, *Map, 
                                                  NumGlobalElements);

  vector<double> Values(NumGlobalElements);
  vector<int>    Indices(NumGlobalElements);

  for (int i = 0 ; i < NumGlobalElements ; ++i) Indices[i] = i;
  
  for (int i = 0 ; i < NumMyElements ; ++i) 
  {
    int iGlobal = MyGlobalElements[i];
    for (int jGlobal = 0 ; jGlobal < NumGlobalElements ; ++jGlobal) 
    {
      Values[jGlobal] = 0.5 / (NumGlobalElements - (iGlobal+1) - (jGlobal + 1) + 1.5);
    }
    Matrix->InsertGlobalValues(MyGlobalElements[i], NumGlobalElements, 
                               &Values[0], &Indices[0]);
  }
    
  // Finish up, trasforming the matrix entries into local numbering,
  // to optimize data transfert during matrix-vector products
  Matrix->FillComplete();
  Matrix->OptimizeStorage();

  return(Matrix);
}

} // namespace Matrices
} // namespace Galeri
#endif
