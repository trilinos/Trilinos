#ifndef GALERI_STRETCHED2DMATRIX_H
#define GALERI_STRETCHED2DMATRIX_H

#include "Galeri_Exception.h"
#include "Galeri_Utils.h"
#include "Epetra_Comm.h"
#include "Epetra_BlockMap.h"
#include "Epetra_CrsMatrix.h"

namespace Galeri {
namespace Matrices {

  //def for espilon muyst be 1e-5
inline
Epetra_CrsMatrix* 
Stretched2D(const Epetra_Map* Map, const int nx, const int ny,
            const double epsilon)
{
  Epetra_CrsMatrix* Matrix = new Epetra_CrsMatrix(Copy, *Map,  9);

  int NumMyElements = Map->NumMyElements();
  int* MyGlobalElements = Map->MyGlobalElements();

  int left, right, lower, upper;
    
  double Values[9];
  int    Indices[9];

  for (int i = 0 ; i < NumMyElements ; ++i) 
  {
    int NumEntries = 0;
    GetNeighboursCartesian2d(MyGlobalElements[i], nx, ny, 
			     left, right, lower, upper);

    if (left != -1) 
    {
      Indices[NumEntries] = left;
      Values[NumEntries] = 2.0 - epsilon;
      ++NumEntries;
    }
    if (right != -1) 
    {
      Indices[NumEntries] = right;
      Values[NumEntries] = 2.0 - epsilon;
      ++NumEntries;
    }
    if (lower != -1) 
    {
      Indices[NumEntries] = lower;
      Values[NumEntries] = -4.0 + epsilon;
      ++NumEntries;
    }
    if (upper != -1) 
    {
      Indices[NumEntries] = upper;
      Values[NumEntries] = -4.0 + epsilon;
      ++NumEntries;
    }
    if (left != -1 && lower != -1) 
    {
      Indices[NumEntries] = lower - 1;
      Values[NumEntries] = -1.0;
      ++NumEntries;
    }
    if (right != -1 && lower != -1) 
    {
      Indices[NumEntries] = lower + 1;
      Values[NumEntries] = -1.0;
      ++NumEntries;
    }
    if (left != -1 && upper != -1) 
    {
      Indices[NumEntries] = upper - 1;
      Values[NumEntries] = -1.0;
      ++NumEntries;
    }
    if (right != -1 && upper != -1) 
    {
      Indices[NumEntries] = upper + 1;
      Values[NumEntries] = -1.0;
      ++NumEntries;
    }
    
    Indices[NumEntries] = MyGlobalElements[i];
    Values[NumEntries] = 8.0;
    ++NumEntries;

    Matrix->InsertGlobalValues(MyGlobalElements[i], NumEntries, 
                               Values, Indices);
  }

  Matrix->FillComplete();
  Matrix->OptimizeStorage();

  return(Matrix);
}

} // namespace Matrices
} // namespace Galeri
#endif
