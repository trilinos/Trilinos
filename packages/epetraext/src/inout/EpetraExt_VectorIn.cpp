//@HEADER
// ************************************************************************
// 
//          Trilinos: An Object-Oriented Solver Framework
//              Copyright (2001) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
// 
// ************************************************************************
//@HEADER
#include "EpetraExt_VectorIn.h"
#include "EpetraExt_MultiVectorIn.h"
#include "Epetra_Comm.h"
#include "Epetra_Vector.h"
#include "Epetra_BlockMap.h"

using namespace EpetraExt;
namespace EpetraExt {

int MatrixMarketFileToVector( const char *filename, const Epetra_BlockMap & map, Epetra_Vector * & A) {

  Epetra_MultiVector * A1;
  if (MatrixMarketFileToMultiVector(filename, map, A1)) return(-1);
  A = dynamic_cast<Epetra_Vector *>(A1);
  return(0);
}
} // namespace EpetraExt
