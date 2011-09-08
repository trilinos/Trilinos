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

#ifndef GALERI_CROSS2D_H
#define GALERI_CROSS2D_H

#include "Galeri_Exception.h"
#include "Galeri_Utils.h"
#include "Epetra_Comm.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"

namespace Galeri {
namespace Matrices {

inline
Epetra_CrsMatrix* 
Cross2D(const Epetra_Map* Map, const int nx, const int ny,
        const double a, const double b, const double c, 
        const double d, const double e)
{
  Epetra_CrsMatrix* Matrix = new Epetra_CrsMatrix(Copy, *Map,  5);

  int NumMyElements = Map->NumMyElements();
  int* MyGlobalElements = Map->MyGlobalElements();

  int left, right, lower, upper;
  vector<double> Values(4);
  vector<int> Indices(4);

  //    e
  //  b a c
  //    d
  for (int i = 0 ; i < NumMyElements ; ++i) 
  {
    int NumEntries = 0;
    GetNeighboursCartesian2d(MyGlobalElements[i], nx, ny, 
			     left, right, lower, upper);

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

Epetra_CrsMatrix* 
Cross2D(const Epetra_Map* Map, const int nx, const int ny,
        const Epetra_Vector& A, const Epetra_Vector& B, const Epetra_Vector& C,
        const Epetra_Vector& D, const Epetra_Vector& E)
{
  Epetra_CrsMatrix* Matrix = new Epetra_CrsMatrix(Copy, *Map,  5);

  int NumMyElements = Map->NumMyElements();
  int* MyGlobalElements = Map->MyGlobalElements();

  int left, right, lower, upper;
  double Values[5];
  int    Indices[5];
    
  for (int i = 0 ; i < NumMyElements ; ++i) 
  {
    int NumEntries = 0;
    GetNeighboursCartesian2d(MyGlobalElements[i], nx, ny, 
			     left, right, lower, upper);

    if (left != -1) 
    {
      Indices[NumEntries] = left;
      Values[NumEntries] = B[i];
      ++NumEntries;
    }
    if (right != -1) 
    {
      Indices[NumEntries] = right;
      Values[NumEntries] = C[i];
      ++NumEntries;
    }
    if (lower != -1) 
    {
      Indices[NumEntries] = lower;
      Values[NumEntries] = D[i];
      ++NumEntries;
    }
    if (upper != -1) 
    {
      Indices[NumEntries] = upper;
      Values[NumEntries] = E[i];
      ++NumEntries;
    }

    Indices[NumEntries] = MyGlobalElements[i];
    Values[NumEntries] = A[i];
    ++NumEntries;

    Matrix->InsertGlobalValues(MyGlobalElements[i], NumEntries, 
                               Values, Indices);
  }

  Matrix->FillComplete();
  Matrix->OptimizeStorage();

  return(Matrix);
}

} // namespace Matrices
} // namespace Galeri
#endif
