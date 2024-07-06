// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef GALERI_BIGCROSS2D_H
#define GALERI_BIGCROSS2D_H

#include "Galeri_Exception.h"
#include "Galeri_Utils.h"
#include "Epetra_Comm.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"

namespace Galeri {
namespace Matrices {

template<typename int_type>
inline
Epetra_CrsMatrix* 
BigCross2D(const Epetra_Map* Map, const int nx, const int ny,
           const double a, 
           const double b, const double c, const double d, const double e,
           const double bb, const double cc, const double dd, const double ee)
{
  Epetra_CrsMatrix* Matrix = new Epetra_CrsMatrix(Copy, *Map,  9);

  int NumMyElements = Map->NumMyElements();
  int_type* MyGlobalElements = 0;
  Map->MyGlobalElementsPtr(MyGlobalElements);

  int left, right, lower, upper;
  int left2, right2, lower2, upper2;

  double Values[9];
  int_type Indices[9];

  //       ee
  //       e
  //  bb b a c cc
  //       d
  //       dd
  
  for (int i = 0 ; i < NumMyElements ; ++i) 
  {
    int NumEntries = 0;
    GetNeighboursCartesian2d(MyGlobalElements[i], nx, ny, 
			     left, right, lower, upper,
                             left2, right2, lower2, upper2);
    // internal small cross
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
    // external elements of the cross
    if (left2 != -1) 
    {
      Indices[NumEntries] = left2;
      Values[NumEntries] = bb;
      ++NumEntries;
    }
    if (right2 != -1) 
    {
      Indices[NumEntries] = right2;
      Values[NumEntries] = cc;
      ++NumEntries;
    }
    if (lower2 != -1) 
    {
      Indices[NumEntries] = lower2;
      Values[NumEntries] = dd;
      ++NumEntries;
    }
    if (upper2 != -1) 
    {
      Indices[NumEntries] = upper2;
      Values[NumEntries] = ee;
      ++NumEntries;
    }
    
    // diagonal element
    Indices[NumEntries] = MyGlobalElements[i];
    Values[NumEntries] = a;
    ++NumEntries;

    Matrix->InsertGlobalValues(MyGlobalElements[i], NumEntries, 
                               &Values[0], &Indices[0]);
  }

  Matrix->FillComplete();
  Matrix->OptimizeStorage();

  return(Matrix);
}

inline
Epetra_CrsMatrix* 
BigCross2D(const Epetra_Map* Map, const int nx, const int ny,
           const double a, 
           const double b, const double c, const double d, const double e,
           const double bb, const double cc, const double dd, const double ee)
{
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  if(Map->GlobalIndicesInt()) {
	  return BigCross2D<int>(Map, nx, ny, a, b, c, d, e, bb, cc, dd, ee);
  }
  else
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  if(Map->GlobalIndicesLongLong()) {
	  return BigCross2D<long long>(Map, nx, ny, a, b, c, d, e, bb, cc, dd, ee);
  }
  else
#endif
    throw "Galeri::Matrices::BigCross2D: GlobalIndices type unknown";
}

} // namespace Matrices
} // namespace Galeri
#endif
