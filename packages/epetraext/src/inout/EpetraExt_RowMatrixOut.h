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
#include <stdio.h>
#include "Epetra_RowMatrix.h"
namespace EpetraExt {
 
  int RowMatrixToMatrixMarketFile( const char *filename, const char * matrixName,
				   const char *matrixDescription, const Epetra_RowMatrix & A, 
				   bool writeHeader=true);
  int RowMatrixToMatlabFile( const char *filename, const char * matrixName,
			     const char *matrixDescription, const Epetra_RowMatrix & A);
   
  int RowMatrixToHandle(FILE * handle, const Epetra_RowMatrix & A);
  int writeRowMatrix(FILE * handle, const Epetra_RowMatrix & A);

} // namespace EpetraExt
