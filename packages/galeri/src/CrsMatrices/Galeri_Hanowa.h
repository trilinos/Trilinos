#ifndef GALERI_HANOWA_H
#define GALERI_HANOWA_H

#include "Galeri_Exception.h"
#include "Epetra_Comm.h"
#include "Epetra_BlockMap.h"
#include "Epetra_CrsMatrix.h"

namespace Galeri {
namespace Matrices {

inline
Epetra_CrsMatrix* Hanowa(const Epetra_Map* Map, const double value)
{
  int NumMyElements     = Map->NumMyElements();
  int NumGlobalElements = Map->NumGlobalElements();
  int* MyGlobalElements = Map->MyGlobalElements();

  if (NumGlobalElements % 2) 
    throw(Exception(__FILE__, __LINE__,
                    "`hanowa' matrix requires a even number of points"));

  Epetra_CrsMatrix* Matrix = new Epetra_CrsMatrix(Copy, *Map, 2);

  double Values[2];
  int    Indices[2];

  int half = NumGlobalElements / 2;
  
  for (int i = 0 ; i < NumMyElements ; ++i) 
  {
    int NumEntries = 2;
    int Global = MyGlobalElements[i];
    Indices[0] = Global;
    if (Global < half) Indices[1] = Global + half;
    else               Indices[1] = Global - half;
    Values[0] = value;
    // convert from C style to FORTRAN style
    if (Global < half) Values[1] = (double) -Global - 1;
    else               Values[1] = (double)  (Global - half) + 1;

    Matrix->InsertGlobalValues(MyGlobalElements[i], NumEntries, 
                               Values, Indices);
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
