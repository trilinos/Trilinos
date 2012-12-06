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

#ifndef GALERI_UTILS_H
#define GALERI_UTILS_H

#include "Galeri_ConfigDefs.h"
#include "Galeri_Exception.h"

class Epetra_BlockMap;
class Epetra_RowMatrix;
class Epetra_CrsMatrix;
class Epetra_MultiVector;
class Epetra_LinearProblem;
namespace Teuchos {
  class ParameterList;
}

namespace Galeri {

Epetra_MultiVector* 
CreateCartesianCoordinates(const string CoordType,
                           const Epetra_BlockMap* BlockMap,
                           Teuchos::ParameterList& List);

void Solve(const Epetra_LinearProblem Problem);

void Solve(const Epetra_RowMatrix* Matrix, const Epetra_MultiVector* LHS,
             const Epetra_MultiVector* RHS);

double ComputeNorm(const Epetra_MultiVector* LHS,
                   const Epetra_MultiVector* RHS);

double ComputeNorm(const Epetra_RowMatrix* A,
                   const Epetra_MultiVector* LHS, const Epetra_MultiVector* RHS);

string toString(const int& x);

string toString(const unsigned int& x);

string toString(const long int& x);

string toString(const double& x);

void 
GetNeighboursCartesian2d(const int i, const int nx, const int ny,
                         int & left, int & right, int & lower, int & upper);

void 
GetNeighboursCartesian2d(const int i, const int nx, const int ny,
                         int& left, int& right, int& lower, int& upper,
                         int& left2, int& right2, int& lower2, int& upper2);

void 
GetNeighboursCartesian3d(const int i, const int nx, const int ny, const int nz,
                         int& left, int& right, int& lower, int& upper,
                         int& below, int& above);
void
PrintStencil2D(const Epetra_CrsMatrix* Matrix,
               const int nx, const int ny, int GID = -1);

} // namespace Galeri

#endif
