// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef GALERI_CROSS3D_H
#define GALERI_CROSS3D_H

#include "Galeri_Exception.h"
#include "Galeri_Utils.h"
#include "Epetra_Comm.h"
#include "Epetra_BlockMap.h"
#include "Epetra_CrsMatrix.h"

namespace Galeri {
namespace Matrices {

template<typename int_type>
inline
Epetra_CrsMatrix* 
Cross3D(const Epetra_Map* Map, 
        const int nx, const int ny, const int nz,
        const double a, const double b, const double c,
        const double d, const double e, const double f,
        const double g)
{
  Epetra_CrsMatrix* Matrix = new Epetra_CrsMatrix(Copy, *Map,  7);

  int NumMyElements = Map->NumMyElements();
  int_type* MyGlobalElements = 0;
  Map->MyGlobalElementsPtr(MyGlobalElements);

  int left, right, lower, upper, below, above;
  std::vector<double> Values(6);
  std::vector<int_type> Indices(6);

  //    e
  //  b a c
  //    d
  // + f below and g above
  
  for (int i = 0 ; i < NumMyElements ; ++i) 
  {
    int NumEntries = 0;
    GetNeighboursCartesian3d(MyGlobalElements[i], nx, ny, nz,
			     left, right, lower, upper, below, above);

    if (left != -1) 
    {
      Indices[NumEntries] = left;
      Values[NumEntries] = b;
      ++NumEntries;
    }
    if (right != -1) 
    {
      Indices[NumEntries] = right;
      Values[NumEntries] = c;
      ++NumEntries;
    }
    if (lower != -1) 
    {
      Indices[NumEntries] = lower;
      Values[NumEntries] = d;
      ++NumEntries;
    }
    if (upper != -1) 
    {
      Indices[NumEntries] = upper;
      Values[NumEntries] = e;
      ++NumEntries;
    }
    if (below != -1) 
    {
      Indices[NumEntries] = below;
      Values[NumEntries] = f;
      ++NumEntries;
    }
    if (above != -1) 
    {
      Indices[NumEntries] = above;
      Values[NumEntries] = g;
      ++NumEntries;
    }
    // put the off-diagonal entries
    Matrix->InsertGlobalValues(MyGlobalElements[i], NumEntries, 
                               &Values[0], &Indices[0]);
    // Put in the diagonal entry
    double diag = a;
	
    Matrix->InsertGlobalValues(MyGlobalElements[i], 1, 
                               &diag, MyGlobalElements + i);
  }
  Matrix->FillComplete();
  Matrix->OptimizeStorage();

  return(Matrix);
}

inline
Epetra_CrsMatrix* 
Cross3D(const Epetra_Map* Map, 
        const int nx, const int ny, const int nz,
        const double a, const double b, const double c,
        const double d, const double e, const double f,
        const double g)
{
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  if(Map->GlobalIndicesInt()) {
	  return Cross3D<int>(Map, nx, ny, nz, a, b, c, d, e, f, g);
  }
  else
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  if(Map->GlobalIndicesLongLong()) {
	  return Cross3D<long long>(Map, nx, ny, nz, a, b, c, d, e, f, g);
  }
  else
#endif
    throw "Galeri::Matrices::Cross3D: GlobalIndices type unknown";
}

} // namespace Matrices
} // namespace Galeri
#endif
