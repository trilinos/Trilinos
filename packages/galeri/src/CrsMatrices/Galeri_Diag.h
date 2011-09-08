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

#ifndef GALERI_DIAGMATRIX_H
#define GALERI_DIAGMATRIX_H

#include "Galeri_Exception.h"
#include "Epetra_Comm.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"

namespace Galeri {
namespace Matrices {

inline
Epetra_CrsMatrix* Diag(const Epetra_Map* Map, double Value)
{
  Epetra_CrsMatrix* Matrix = new Epetra_CrsMatrix(Copy, *Map,  1);

  int NumMyElements = Map->NumMyElements();
  int* MyGlobalElements = Map->MyGlobalElements();

  for (int i = 0 ; i < NumMyElements ; ++i) 
  {
    int Indices = MyGlobalElements[i];

    Matrix->InsertGlobalValues(MyGlobalElements[i], 1, &Value, &Indices);
  }

  Matrix->FillComplete();
  Matrix->OptimizeStorage();

  return(Matrix);
}

Epetra_CrsMatrix* Diag(const Epetra_Map* Map, Epetra_Vector& Vector)
{
  Epetra_CrsMatrix* Matrix = new Epetra_CrsMatrix(Copy, *Map,  1);

  if (!Map->SameAs(Vector.Map()))
    throw(Exception(__FILE__, __LINE__,
                    "Vector.Map() is not equivalent to input Map"));

  int NumMyElements = Map->NumMyElements();
  int* MyGlobalElements = Map->MyGlobalElements();

  for (int i = 0 ; i < NumMyElements ; ++i) 
  {
    int Indices = MyGlobalElements[i];
    double Value = Vector[i];

    Matrix->InsertGlobalValues(MyGlobalElements[i], 1, &Value, &Indices);
  }

  Matrix->FillComplete();
  Matrix->OptimizeStorage();

  return(Matrix);
}

} // namespace Matrices
} // namespace Galeri
#endif
