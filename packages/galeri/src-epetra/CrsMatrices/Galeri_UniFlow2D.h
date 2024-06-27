// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef GALERI_UNIFLOW_H
#define GALERI_UNIFLOW_H

#include "Galeri_Exception.h"
#include "Galeri_Utils.h"
#include "Galeri_Cross2D.h"
#include "Epetra_Comm.h"
#include "Epetra_BlockMap.h"
#include "Epetra_CrsMatrix.h"

namespace Galeri {
namespace Matrices {

  // def for conv = 1;
  // def for diff = 1e-5
  // for alpha = 0
template<typename int_type>
inline
Epetra_CrsMatrix* UniFlow2D(const Epetra_Map* Map, 
                            const int nx, const int ny,
                            const double lx, const double ly,
                            const double conv, const double diff,
                            const double alpha)
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
    double ConvX = conv * cos(alpha) / hx;
    double ConvY = conv * sin(alpha) / hy;

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
    A[i] += diff *2. / (hx * hx) + diff * 2. /(hy * hy);
    B[i] -= diff / (hx * hx);
    C[i] -= diff / (hx * hx);
    D[i] -= diff / (hy * hy);
    E[i] -= diff / (hy * hy);
  }

  return(Cross2D(Map, nx, ny, A, B, C, D, E));
}

inline
Epetra_CrsMatrix* UniFlow2D(const Epetra_Map* Map, 
                            const int nx, const int ny,
                            const double lx, const double ly,
                            const double conv, const double diff,
                            const double alpha)
{
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  if(Map->GlobalIndicesInt()) {
	  return UniFlow2D<int>(Map, nx, ny, lx, ly, conv, diff, alpha);
  }
  else
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  if(Map->GlobalIndicesLongLong()) {
	  return UniFlow2D<long long>(Map, nx, ny, lx, ly, conv, diff, alpha);
  }
  else
#endif
    throw "Galeri::Matrices::UniFlow2D: GlobalIndices type unknown";
}

} // namespace Matrices
} // namespace Galeri
#endif
