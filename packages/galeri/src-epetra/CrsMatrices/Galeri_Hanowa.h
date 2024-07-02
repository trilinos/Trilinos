// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef GALERI_HANOWA_H
#define GALERI_HANOWA_H

#include "Galeri_Exception.h"
#include "Epetra_Comm.h"
#include "Epetra_BlockMap.h"
#include "Epetra_CrsMatrix.h"

namespace Galeri {
namespace Matrices {

template<typename int_type>
inline
Epetra_CrsMatrix* Hanowa(const Epetra_Map* Map, const double value)
{
  int NumMyElements     = Map->NumMyElements();
  int_type NumGlobalElements = (int_type) Map->NumGlobalElements64();
  int_type* MyGlobalElements = 0;
  Map->MyGlobalElementsPtr(MyGlobalElements);

  if (NumGlobalElements % 2) 
    throw(Exception(__FILE__, __LINE__,
                    "`hanowa' matrix requires a even number of points"));

  Epetra_CrsMatrix* Matrix = new Epetra_CrsMatrix(Copy, *Map, 2);

  double Values[2];
  int_type Indices[2];

  int_type half = NumGlobalElements / 2;
  
  for (int i = 0 ; i < NumMyElements ; ++i) 
  {
    int NumEntries = 2;
    int_type Global = MyGlobalElements[i];
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

Epetra_CrsMatrix* Hanowa(const Epetra_Map* Map, const double value)
{
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  if(Map->GlobalIndicesInt()) {
	  return Hanowa<int>(Map, value);
  }
  else
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  if(Map->GlobalIndicesLongLong()) {
	  return Hanowa<long long>(Map, value);
  }
  else
#endif
    throw "Galeri::Matrices::Hanowa: GlobalIndices type unknown";
}

} // namespace Matrices
} // namespace Galeri
#endif
