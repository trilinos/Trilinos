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

#ifndef GALERI_STAR2D_H
#define GALERI_STAR2D_H

#include "Galeri_Exception.h"
#include "Galeri_Utils.h"
#include "Epetra_Comm.h"
#include "Epetra_BlockMap.h"
#include "Epetra_CrsMatrix.h"

namespace Galeri {
namespace Matrices {

inline
Epetra_CrsMatrix* 
Star2D(const Epetra_Map* Map, const int nx, const int ny,
       const double a, const double b, const double c,
       const double d, const double e,
       const double z1, const double z2,
       const double z3, const double z4)
{
  Epetra_CrsMatrix* Matrix = new Epetra_CrsMatrix(Copy, *Map,  9);

  int NumMyElements = Map->NumMyElements();
  int* MyGlobalElements = Map->MyGlobalElements();

  int left, right, lower, upper;
  double Values[9];
  int    Indices[9];

  //  z3  e  z4
  //   b  a  c
  //  z1  d  z2
  for (int i = 0 ; i < NumMyElements ; ++i) 
  {
    int NumEntries = 0;
    GetNeighboursCartesian2d(MyGlobalElements[i], nx, ny, 
			     left, right, lower, upper);

    if (left != -1) 
    {
      Values[NumEntries] = b;
      Indices[NumEntries] = left;
      ++NumEntries;
    }
    if (right != -1) 
    {
      Values[NumEntries] = c;
      Indices[NumEntries] = right;
      ++NumEntries;
    }
    if (lower != -1) 
    {
      Values[NumEntries] = d;
      Indices[NumEntries] = lower;
      ++NumEntries;
    }
    if (upper != -1) 
    {
      Values[NumEntries] = e;
      Indices[NumEntries] = upper;
      ++NumEntries;
    }
    if (left != -1 && lower != -1) 
    {
      Values[NumEntries] = z1;
      Indices[NumEntries] = lower - 1;
      ++NumEntries;
    }
    if (right != -1 && lower != -1) 
    {
      Values[NumEntries] = z2;
      Indices[NumEntries] = lower + 1;
      ++NumEntries;
    }
    if (left != -1 && upper != -1) 
    {
      Values[NumEntries] = z3;
      Indices[NumEntries] = upper - 1;
      ++NumEntries;
    }
    if (right != -1 && upper != -1) 
    {
      Values[NumEntries] = z4;
      Indices[NumEntries] = upper + 1;
      ++NumEntries;
    }
    
    Values[NumEntries] = a;
    Indices[NumEntries] = MyGlobalElements[i];
    ++NumEntries;

    Matrix->InsertGlobalValues(MyGlobalElements[i], NumEntries, 
                               &Values[0], &Indices[0]);
  }

  Matrix->FillComplete();
  Matrix->OptimizeStorage();

  return(Matrix);
}

} // namespace Matrices
} // namespace Galeri
#endif
