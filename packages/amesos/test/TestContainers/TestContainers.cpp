// @HEADER
// ***********************************************************************
// 
//                Amesos: Direct Sparse Solver Package
//                 Copyright (2004) Sandia Corporation
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
// @HEADER

// Test the two containers:
// - Amesos_ContainerLAPACK
// - Amesos_ContainerEpetraCrs
//
#include "Amesos_config.h"
#include "Amesos.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Trilinos_Util_CrsMatrixGallery.h"
#include "Teuchos_ParameterList.hpp"
#include "Amesos_Container.h"
#include "Amesos_ContainerLAPACK.h"
#include "Amesos_ContainerEpetraCrs.h"
#include "Amesos_Utils.h"
#include "Amesos.h"
#include "Amesos_Preconditioner.h"
#include <stdlib.h>

using namespace Trilinos_Util;

int FillMatrix(Amesos_Container& Container, int NumRows, int NumVectors)
{

  Container.Shape(NumRows,NumVectors);

  double val;

  for (int i = 0 ; i < NumRows ; ++i) {
    // fill matrix
    // add diagonal element
    AMESOS_CHK_ERR(Container.SetMatrixElement(i,i,2.0 * NumRows));

    for (int j = 0 ; j < NumRows ; ++j) {
      if (j != i ) {
	val = 1.0 * rand() / RAND_MAX;
	AMESOS_CHK_ERR(Container.SetMatrixElement(i,j,val));
      }
    }

    // now fill vectors
    for (int j = 0 ; j < NumVectors ; ++j) {
      Container.LHS(i,j) = 0.0;
      Container.RHS(i,j) = 1.0;
    }
  }

  AMESOS_CHK_ERR(Container.Compute());

  return(0);
}
  
//==============================================================================
int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  int NumRows = 100;
  int NumVectors = 10;

  Amesos_ContainerEpetraCrs EpetraCrs;
  Amesos_ContainerLAPACK LAPACK;

  AMESOS_CHK_ERR(FillMatrix(EpetraCrs,NumRows,NumVectors));
  AMESOS_CHK_ERR(FillMatrix(LAPACK,NumRows,NumVectors));

  // for the LAPACK container, I just have to call
  // ApplyInverse().

  LAPACK.ApplyInverse();

  // for the EpetraCrs container, a priori I do not have any ApplyInverse().
  // Here, I build an exact LU factorization, and I use the call
  // Amesos_Preconditioner to define ApplyInverse. This class is stick
  // into the container using SetInversePointer().

  Epetra_CrsMatrix* Matrix;
  EpetraCrs.GetMatrixPointer((void**)&Matrix);

  Teuchos::ParameterList List;
  Amesos_Preconditioner Preconditioner("Amesos_Klu",Matrix,List);

  EpetraCrs.SetInversePointer((void*)&Preconditioner);

  EpetraCrs.ApplyInverse();

  // check that the two approaches are giving approximatively the same
  // results.

  double Norm = 0.0;

  for (int i = 0 ; i < NumRows ; ++i) {
    for (int j = 0 ; j < NumVectors ; ++j) {
      double val =  EpetraCrs.LHS(i,j) - LAPACK.LHS(i,j);
      Norm += val * val;
    }
  }

  cout << Norm << endl;

#ifdef HAVE_MPI
  MPI_Finalize() ; 
#endif

  if (Norm < 1e-2) 
    exit(EXIT_SUCCESS);
  else
    AMESOS_CHK_ERR(-1);

}
