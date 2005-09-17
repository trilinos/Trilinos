#ifndef GALERI_DIAGMATRIX_H
#define GALERI_DIAGMATRIX_H

#include "Galeri_Exception.h"
#include "Epetra_Comm.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"

namespace Galeri {
namespace Matrices {

inline
Epetra_CrsMatrix* Diag(const Epetra_Map* Map, double Value)
{
  Epetra_CrsMatrix* Matrix = new Epetra_CrsMatrix(Copy, *Map,  1);

  int NumMyElements = Map->NumMyElements();
  int* MyGlobalElements = Map->MyGlobalElements();

  for (int i = 0 ; i < NumMyElements ; ++i) 
  {
    int Indices = MyGlobalElements[i];

    Matrix->InsertGlobalValues(MyGlobalElements[i], 1, &Value, &Indices);
  }

  Matrix->FillComplete();
  Matrix->OptimizeStorage();

  return(Matrix);
}

Epetra_CrsMatrix* Diag(const Epetra_Map* Map, Epetra_Vector& Vector)
{
  Epetra_CrsMatrix* Matrix = new Epetra_CrsMatrix(Copy, *Map,  1);

  if (!Map->SameAs(Vector.Map()))
    throw(Exception(__FILE__, __LINE__,
                    "Vector.Map() is not equivalent to input Map"));

  int NumMyElements = Map->NumMyElements();
  int* MyGlobalElements = Map->MyGlobalElements();

  for (int i = 0 ; i < NumMyElements ; ++i) 
  {
    int Indices = MyGlobalElements[i];
    double Value = Vector[i];

    Matrix->InsertGlobalValues(MyGlobalElements[i], 1, &Value, &Indices);
  }

  Matrix->FillComplete();
  Matrix->OptimizeStorage();

  return(Matrix);
}

} // namespace Matrices
} // namespace Galeri
#endif
