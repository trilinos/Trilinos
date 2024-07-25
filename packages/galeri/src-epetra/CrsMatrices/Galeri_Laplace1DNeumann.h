// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#ifndef GALERI_LAPLACE1DNEUMANNMATRIX_H
#define GALERI_LAPLACE1DNEUMANNMATRIX_H

#include "Galeri_Exception.h"
#include "Epetra_Comm.h"
#include "Epetra_BlockMap.h"
#include "Epetra_CrsMatrix.h"

namespace Galeri {
namespace Matrices {

template<typename int_type>
inline
Epetra_CrsMatrix* Laplace1DNeumann(const Epetra_Map* Map)
{
  int NumMyElements     = Map->NumMyElements();
  int_type NumGlobalElements = (int_type) Map->NumGlobalElements64();
  int_type* MyGlobalElements = 0;
  Map->MyGlobalElementsPtr(MyGlobalElements);

  Epetra_CrsMatrix* Matrix = new Epetra_CrsMatrix(Copy, *Map, 3);

  std::vector<double> Values(2);
  std::vector<int_type> Indices(2);
  int NumEntries;

  for (int i = 0 ; i < NumMyElements ; ++i) 
  {
    if (MyGlobalElements[i] == 0) 
    {
      Indices[0] = 1;
      NumEntries = 1;
      Values[0] = -1.0;
    } 
    else if (MyGlobalElements[i] == NumGlobalElements - 1) 
    {
      Indices[0] = NumGlobalElements - 2;
      NumEntries = 1;
      Values[0] = -1.0;
    } else {
      Indices[0] = MyGlobalElements[i] - 1;
      Values[1] = -1.0;
      Indices[1] = MyGlobalElements[i] + 1;
      Values[0] = -1.0;
      NumEntries = 2;
    }

    Matrix->InsertGlobalValues(MyGlobalElements[i], NumEntries, 
                               &Values[0], &Indices[0]);

    // Put in the diagonal entry
    if (MyGlobalElements[i] == 0 || 
        (MyGlobalElements[i] == NumGlobalElements - 1))
      Values[0] = 1.0;
    else
      Values[0] = 2.0;
    
    Matrix->InsertGlobalValues(MyGlobalElements[i], 1, 
                               &Values[0], MyGlobalElements + i);
  }
  
  // Finish up, trasforming the matrix entries into local numbering,
  // to optimize data transfert during matrix-vector products
  Matrix->FillComplete();
  Matrix->OptimizeStorage();

  return(Matrix);
}

inline
Epetra_CrsMatrix* Laplace1DNeumann(const Epetra_Map* Map)
{
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  if(Map->GlobalIndicesInt()) {
	  return Laplace1DNeumann<int>(Map);
  }
  else
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  if(Map->GlobalIndicesLongLong()) {
	  return Laplace1DNeumann<long long>(Map);
  }
  else
#endif
    throw "Galeri::Matrices::Laplace1DNeumann: GlobalIndices type unknown";
}

} // namespace Matrices
} // namespace Galeri
#endif
