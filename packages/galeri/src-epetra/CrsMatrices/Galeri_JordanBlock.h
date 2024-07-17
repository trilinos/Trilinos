// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef GALERI_JORDANBLOCK_H
#define GALERI_JORDANBLOCK_H

#include "Galeri_Exception.h"
#include "Epetra_Comm.h"
#include "Epetra_BlockMap.h"
#include "Epetra_CrsMatrix.h"

namespace Galeri {
namespace Matrices {

template<typename int_type>
inline
Epetra_CrsMatrix* JordanBlock(const Epetra_Map* Map, const double value)
{
  // this is actually a dense matrix, stored into Crs format
  int NumMyElements     = Map->NumMyElements();
  int_type NumGlobalElements = (int_type) Map->NumGlobalElements64();
  int_type* MyGlobalElements = 0;
  Map->MyGlobalElementsPtr(MyGlobalElements);

  Epetra_CrsMatrix* Matrix = new Epetra_CrsMatrix(Copy, *Map, 2);

  int_type Indices[2];
  double Values[2];
  
  for (int i = 0 ; i < NumMyElements ; ++i) 
  {
    int NumEntries = 0;
    if (MyGlobalElements[i] != NumGlobalElements - 1) 
    {
      Indices[NumEntries] = MyGlobalElements[i] + 1;
      Values[NumEntries] = 1.0;
      NumEntries++;
    }
    // diagonal contribution
    Indices[NumEntries] = MyGlobalElements[i];
    Values[NumEntries] = value;
    NumEntries++;

    Matrix->InsertGlobalValues(MyGlobalElements[i], NumEntries,
                               Values, Indices);
  }

  // Finish up, trasforming the matrix entries into local numbering,
  // to optimize data transfert during matrix-vector products
  Matrix->FillComplete();
  Matrix->OptimizeStorage();

  return(Matrix);
}

Epetra_CrsMatrix* JordanBlock(const Epetra_Map* Map, const double value)
{
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  if(Map->GlobalIndicesInt()) {
	  return JordanBlock<int>(Map, value);
  }
  else
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  if(Map->GlobalIndicesLongLong()) {
	  return JordanBlock<long long>(Map, value);
  }
  else
#endif
    throw "Galeri::Matrices::JordanBlock: GlobalIndices type unknown";
}

} // namespace Matrices
} // namespace Galeri
#endif
