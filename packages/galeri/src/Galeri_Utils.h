// @HEADER
// ************************************************************************
//
//           Galeri: Finite Element and Matrix Generation Package
//                 Copyright (2006) ETHZ/Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
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

void  Solve(const Epetra_LinearProblem Problem);

void Solve(const Epetra_RowMatrix* Matrix, const Epetra_MultiVector* LHS,
             const Epetra_MultiVector* RHS);

double ComputeNorm(const Epetra_MultiVector* LHS,
                   const Epetra_MultiVector* RHS);

double ComputeNorm(const Epetra_RowMatrix* A,
                   const Epetra_MultiVector* LHS, const Epetra_MultiVector* RHS);

string toString(const int& x);

string toString(const unsigned int& x);

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
