// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef GALERI_BENTPIPE2D_H
#define GALERI_BENTPIPE2D_H

#include "Galeri_Exception.h"
#include "Galeri_Cross2D.h"
#include "Galeri_Utils.h"
#include "Epetra_Comm.h"
#include "Epetra_BlockMap.h"
#include "Epetra_CrsMatrix.h"

namespace Galeri {
namespace Matrices {

template<typename int_type>
inline
Epetra_CrsMatrix* 
BentPipe2D(const Epetra_Map* Map, const int nx, const int ny,
           const double lx, const double ly,
           const double conv, const double diff)
{
  int NumMyElements = Map->NumMyElements();
  int_type* MyGlobalElements = 0;
  Map->MyGlobalElementsPtr(MyGlobalElements);

  Epetra_Vector A(*Map);
  Epetra_Vector B(*Map);
  Epetra_Vector C(*Map);
  Epetra_Vector D(*Map);
  Epetra_Vector E(*Map);
  
  A.PutScalar(0.0);
  B.PutScalar(0.0);
  C.PutScalar(0.0);
  D.PutScalar(0.0);
  E.PutScalar(0.0);
  
  double hx = lx / (nx + 1);
  double hy = ly / (ny + 1);

  for (int i = 0 ; i < NumMyElements ; ++i) 
  {
    int ix, iy;
    ix = (MyGlobalElements[i]) % nx;
    iy = (MyGlobalElements[i] - ix) / nx;
    double x = hx * (ix + 1);
    double y = hy * (iy + 1);
    double ConvX =  conv * 2 * x * (0.5 * x - 1.) * (1. - 2 * y) / hx;
    double ConvY = -conv * 4 * y * (y - 1.) * (1. - x) / hy;

    // convection part
    
    if (ConvX < 0) 
    {
      C[i] += ConvX;
      A[i] -= ConvX;
    } 
    else 
    {
      B[i] -= ConvX;
      A[i] += ConvX;
    }

    if (ConvY < 0) 
    {
      E[i] += ConvY;
      A[i] -= ConvY;
    } 
    else 
    {
      D[i] -= ConvY;
      A[i] += ConvY;
    }

    // add diffusion part
    A[i] += diff * 2. / (hx * hx) + diff * 2. / (hy * hy);
    B[i] -= diff / (hx * hx);
    C[i] -= diff / (hx * hx);
    D[i] -= diff / (hy * hy);
    E[i] -= diff / (hy * hy);
  }

  return(Matrices::Cross2D(Map, nx, ny, A, B, C, D, E));
}

inline
Epetra_CrsMatrix* 
BentPipe2D(const Epetra_Map* Map, const int nx, const int ny,
           const double lx, const double ly,
           const double conv, const double diff)
{
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  if(Map->GlobalIndicesInt()) {
	  return BentPipe2D<int>(Map, nx, ny, lx, ly, conv, diff);
  }
  else
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  if(Map->GlobalIndicesLongLong()) {
	  return BentPipe2D<long long>(Map, nx, ny, lx, ly, conv, diff);
  }
  else
#endif
    throw "Galeri::Matrices::BentPipe2D: GlobalIndices type unknown";
}

} // namespace Matrices
} // namespace Galeri
#endif
