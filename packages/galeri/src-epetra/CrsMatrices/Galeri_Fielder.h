// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef GALERI_FIEDLER_H
#define GALERI_FIEDLER_H

#include "Galeri_Exception.h"
#include "Epetra_Comm.h"
#include "Epetra_BlockMap.h"
#include "Epetra_CrsMatrix.h"

namespace Galeri {
namespace Matrices {

template<typename int_type>
inline
Epetra_CrsMatrix* Fielder(const Epetra_Map* Map)
{
  // this is actually a dense matrix, stored into Crs format
  int_type NumGlobalElements = (int_type) Map->NumGlobalElements64();
  int NumMyElements     = Map->NumMyElements();
  int_type* MyGlobalElements = 0;
  Map->MyGlobalElementsPtr(MyGlobalElements);

  Epetra_CrsMatrix* Matrix = new Epetra_CrsMatrix(Copy, *Map, NumGlobalElements);

  std::vector<double> Values(NumGlobalElements);
  std::vector<int_type> Indices(NumGlobalElements);

  for (int i = 0 ; i < NumMyElements ; ++i) 
  {
    int_type NumEntries = NumGlobalElements;
    int_type iGlobal = MyGlobalElements[i];
    
    for (int_type j = 0 ; j < NumGlobalElements ; ++j) 
    {
      Indices[j] = j;
      Values[j]  = (double)abs((double) (iGlobal - j));
    }

    Matrix->InsertGlobalValues(MyGlobalElements[i], NumEntries, 
                               &Values[0], &Indices[0]);
  }

  // Finish up, trasforming the matrix entries into local numbering,
  // to optimize data transfert during matrix-vector products
  Matrix->FillComplete();
  Matrix->OptimizeStorage();

  return(Matrix);
}

Epetra_CrsMatrix* Fielder(const Epetra_Map* Map)
{
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  if(Map->GlobalIndicesInt()) {
	  return Fielder<int>(Map);
  }
  else
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  if(Map->GlobalIndicesLongLong()) {
	  return Fielder<long long>(Map);
  }
  else
#endif
    throw "Galeri::Matrices::Fielder: GlobalIndices type unknown";
}

} // namespace Matrices
} // namespace Galeri
#endif
