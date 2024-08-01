// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  if(Map->GlobalIndicesInt()) {
    int* MyGlobalElements = Map->MyGlobalElements();

    for (int i = 0 ; i < NumMyElements ; ++i) 
    {
      int Indices = MyGlobalElements[i];

      Matrix->InsertGlobalValues(MyGlobalElements[i], 1, &Value, &Indices);
    }
  }
  else
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  if(Map->GlobalIndicesLongLong()) {
    long long* MyGlobalElements = Map->MyGlobalElements64();

    for (int i = 0 ; i < NumMyElements ; ++i) 
    {
      long long Indices = MyGlobalElements[i];

      Matrix->InsertGlobalValues(MyGlobalElements[i], 1, &Value, &Indices);
    }
  }
  else
#endif
    throw "Galeri::Matrices::Diag: GlobalIndices type unknown";

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
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  if(Map->GlobalIndicesInt()) {
    int* MyGlobalElements = Map->MyGlobalElements();

    for (int i = 0 ; i < NumMyElements ; ++i) 
    {
      int Indices = MyGlobalElements[i];
      double Value = Vector[i];

      Matrix->InsertGlobalValues(MyGlobalElements[i], 1, &Value, &Indices);
    }
  }
  else
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  if(Map->GlobalIndicesLongLong()) {
    long long* MyGlobalElements = Map->MyGlobalElements64();

    for (int i = 0 ; i < NumMyElements ; ++i) 
    {
      long long Indices = MyGlobalElements[i];
      double Value = Vector[i];

      Matrix->InsertGlobalValues(MyGlobalElements[i], 1, &Value, &Indices);
    }
  }
  else
#endif
    throw "Galeri::Matrices::Diag: GlobalIndices type unknown";

  Matrix->FillComplete();
  Matrix->OptimizeStorage();

  return(Matrix);
}

} // namespace Matrices
} // namespace Galeri
#endif
