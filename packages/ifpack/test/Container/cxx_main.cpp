// @HEADER
// ***********************************************************************
// 
//                IFPACK
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

#include "Ifpack_ConfigDefs.h"
#if defined(HAVE_IFPACK_AMESOS) && defined(HAVE_IFPACK_TEUCHOS)
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_LinearProblem.h"
#include "Trilinos_Util_CrsMatrixGallery.h"
#include "Teuchos_ParameterList.hpp"
#include "Ifpack_DenseContainer.h"
#include "Ifpack_SparseContainer.h"
#include "Ifpack_Amesos.h"

using namespace Trilinos_Util;

int TestContainer(string Type, CrsMatrixGallery& Gallery)
{

  Epetra_LinearProblem* Problem = Gallery.GetLinearProblem();

  Epetra_MultiVector& RHS = *(Problem->GetRHS());
  Epetra_MultiVector& LHS = *(Problem->GetLHS());
  Epetra_RowMatrix* A = Problem->GetMatrix();

  int NumVectors = RHS.NumVectors();
  int NumMyRows = A->NumMyRows();

  cout << "Container type = " << Type << endl;
  cout << "NumMyRows = " << NumMyRows << ", NumVectors = " << NumVectors << endl;
  LHS.PutScalar(0.0);
  
  Ifpack_Container* Container;

  if (Type == "dense")
    Container = new Ifpack_DenseContainer;
  else
    Container = new Ifpack_SparseContainer<Ifpack_Amesos>;

  assert (Container != 0);

  // shape the container to hold the entire local matrix,
  // and NumVectors vectors.
  IFPACK_CHK_ERR(Container->Shape(A->NumMyRows(), NumVectors));

  // set as ID all the local rows of A
  for (int i = 0 ; i < A->NumMyRows() ; ++i)
    Container->ID(i) = i;

  // extract submatrix (in this case, the entire matrix)
  // and complete setup
  IFPACK_CHK_ERR(Container->Extract(A));
  IFPACK_CHK_ERR(Container->Compute());

  // set the RHS and LHS
  for (int i = 0 ; i < A->NumMyRows() ; ++i)
    for (int j = 0 ; j < NumVectors ; ++j) {
      Container->RHS(i,j) = RHS[j][i];
      Container->LHS(i,j) = LHS[j][i];
    }

  
  // set parameters (empty for dense containers)
  Teuchos::ParameterList List;
  List.set("amesos: solver type", Type);
  IFPACK_CHK_ERR(Container->SetParameters(List));

  // solve the linear system
  IFPACK_CHK_ERR(Container->ApplyInverse());

  // get the computed solution, store it in LHS
  for (int i = 0 ; i < A->NumMyRows() ; ++i)
    for (int j = 0 ; j < NumVectors ; ++j) {
       LHS[j][i] = Container->LHS(i,j);
    }

  // check residual
  double* residual = new double[NumVectors];
  double *diff = new double[NumVectors];

  Gallery.ComputeResidual(residual);
  Gallery.ComputeDiffBetweenStartingAndExactSolutions(diff);

  if( A->Comm().MyPID()==0 ) {
    for (int i = 0 ; i < NumVectors ; ++i) {
      cout << "eq " << i << ", ||b-Ax||_2 = " << residual[i] << endl;
      cout << "eq " << i << ", ||x_exact - x||_2 = " << diff[i] << endl;
    }
  }

  if ((residual[0] < 1e-5) && (diff[0] < 1e-5))
    return(0);
  else
    return(-1);

  delete [] residual;
  delete [] diff;

}

int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  // size of the global matrix. 
  const int NumPoints = 900;

  CrsMatrixGallery Gallery("laplace_2d", Comm);
  Gallery.Set("problem_size", NumPoints);
  Gallery.Set("map_type", "linear");
  Gallery.Set("num_vectors", 5);

  // test the preconditioner
  int TestPassed = true;
  if (TestContainer("dense",Gallery))
    TestPassed = false;
  if (TestContainer("sparse",Gallery))
    TestPassed = false;

#ifdef HAVE_MPI
  MPI_Finalize(); 
#endif

  if (TestPassed)
    exit(EXIT_SUCCESS);
  else
    exit(EXIT_FAILURE);
}

#else

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  puts("please configure IFPACK with --enable-teuchos");
  puts("--enable-amesos to run this test");

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif
  return(EXIT_SUCCESS);
}

#endif
