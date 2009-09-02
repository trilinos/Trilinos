//@HEADER
// ***********************************************************************
// 
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2001) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER
#include "EpetraExt_VectorOut.h"
#include "EpetraExt_MultiVectorOut.h"
#include "EpetraExt_mmio.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_Import.h"
#include "Epetra_CrsMatrix.h"

using namespace EpetraExt;
namespace EpetraExt {

// Make use of the fact that Vector is a specialization of MultiVector
int VectorToMatlabFile( const char *filename, const Epetra_Vector & A) {
  return(MultiVectorToMatlabFile(filename, A));
}

int VectorToMatrixMarketFile( const char *filename, const Epetra_Vector & A, 
				 const char * matrixName,
				 const char *matrixDescription, 
				 bool writeHeader) {
  return(MultiVectorToMatrixMarketFile( filename, A, matrixName, matrixDescription, writeHeader));
}

int VectorToHandle(FILE * handle, const Epetra_Vector & A) {

  return(MultiVectorToHandle(handle, A, true));
}
int writeVector(FILE * handle, const Epetra_Vector & A) {

  return(writeMultiVector(handle, A, true));
}
} // namespace EpetraExt
