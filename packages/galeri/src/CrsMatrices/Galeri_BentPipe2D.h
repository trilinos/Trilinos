// @HEADER
// ************************************************************************
//
//           Galeri: Finite Element and Matrix Generation Package
//                 Copyright (2006) ETHZ/Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions about Galeri? Contact Marzio Sala (marzio.sala _AT_ gmail.com)
//
// ************************************************************************
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

inline
Epetra_CrsMatrix* 
BentPipe2D(const Epetra_Map* Map, const int nx, const int ny,
           const double lx, const double ly,
           const double conv, const double diff)
{
  int NumMyElements = Map->NumMyElements();
  int* MyGlobalElements = Map->MyGlobalElements();

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

} // namespace Matrices
} // namespace Galeri
#endif
