// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef GALERI_TRIDIAG_H
#define GALERI_TRIDIAG_H

#include "Galeri_Exception.h"
#include "Epetra_Comm.h"
#include "Epetra_BlockMap.h"
#include "Epetra_CrsMatrix.h"

namespace Galeri {
namespace Matrices {

template<typename int_type>
inline
Epetra_CrsMatrix* 
Tridiag(const Epetra_Map* Map, const double a, const double b, const double c)
{
  Epetra_CrsMatrix* Matrix = new Epetra_CrsMatrix(Copy, *Map,  3);

  int_type NumGlobalElements = (int_type) Map->NumGlobalElements64();
  int NumMyElements = Map->NumMyElements();
  int_type* MyGlobalElements = 0;
  Map->MyGlobalElementsPtr(MyGlobalElements);

  std::vector<double> Values(2);
  std::vector<int_type> Indices(2);
  int NumEntries;

  for (int i = 0 ; i < NumMyElements ; ++i) 
  {
    if (MyGlobalElements[i] == 0) 
    {
      // off-diagonal for first row
      Indices[0] = 1;
      NumEntries = 1;
      Values[0] = c;
    } 
    else if (MyGlobalElements[i] == NumGlobalElements - 1) 
    {
      // off-diagonal for last row
      Indices[0] = NumGlobalElements - 2;
      NumEntries = 1;
      Values[0] = b;
    } 
    else 
    {
      // off-diagonal for internal row
      Indices[0] = MyGlobalElements[i] - 1;
      Values[1] = b;
      Indices[1] = MyGlobalElements[i] + 1;
      Values[0] = c;
      NumEntries = 2;
    }
    Matrix->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[0], 
                               &Indices[0]);
    // Put in the diagonal entry
    Values[0] = a;
    Matrix->InsertGlobalValues(MyGlobalElements[i], 1, &Values[0], 
                               MyGlobalElements + i);
  }
  
  // Finish up, trasforming the matrix entries into local numbering,
  // to optimize data transfert during matrix-vector products
  Matrix->FillComplete();
  Matrix->OptimizeStorage();

  return(Matrix);
}

Epetra_CrsMatrix* 
Tridiag(const Epetra_Map* Map, const double a, const double b, const double c)
{
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  if(Map->GlobalIndicesInt()) {
	  return Tridiag<int>(Map, a, b, c);
  }
  else
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  if(Map->GlobalIndicesLongLong()) {
	  return Tridiag<long long>(Map, a, b, c);
  }
  else
#endif
    throw "Galeri::Matrices::Tridiag: GlobalIndices type unknown";
}

} // namespace Matrices
} // namespace Galeri
#endif
