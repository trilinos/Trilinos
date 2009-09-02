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

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#endif
#include "Epetra_SerialComm.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_LinearProblem.h"
#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"
#include "Galeri_Utils.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Ifpack_DenseContainer.h"
#include "Ifpack_SparseContainer.h"
#include "Ifpack_Amesos.h"
#include "Ifpack_LocalFilter.h"

static bool verbose = false;

// ======================================================================
bool TestContainer(string Type, const Teuchos::RefCountPtr<Epetra_RowMatrix>& A)
{
  int NumVectors = 3;
  int NumMyRows = A->NumMyRows();

  Epetra_MultiVector LHS_exact(A->RowMatrixRowMap(), NumVectors);
  Epetra_MultiVector LHS(A->RowMatrixRowMap(), NumVectors);
  Epetra_MultiVector RHS(A->RowMatrixRowMap(), NumVectors);
  LHS_exact.Random(); LHS.PutScalar(0.0); 
  A->Multiply(false, LHS_exact, RHS);

  Epetra_LinearProblem Problem(&*A, &LHS, &RHS);

  if (verbose) {
    cout << "Container type = " << Type << endl;
    cout << "NumMyRows = " << NumMyRows << ", NumVectors = " << NumVectors << endl;
  }
  LHS.PutScalar(0.0);
  
  Teuchos::RefCountPtr<Ifpack_Container> Container;

  if (Type == "dense")
    Container = Teuchos::rcp( new Ifpack_DenseContainer(A->NumMyRows(), NumVectors) );
  else
    Container = Teuchos::rcp( new Ifpack_SparseContainer<Ifpack_Amesos>(A->NumMyRows(), NumVectors) );

  assert (Container != Teuchos::null);

  IFPACK_CHK_ERR(Container->Initialize());
  // set as ID all the local rows of A
  for (int i = 0 ; i < A->NumMyRows() ; ++i)
    Container->ID(i) = i;

  // extract submatrix (in this case, the entire matrix)
  // and complete setup
  IFPACK_CHK_ERR(Container->Compute(*A));

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

  double residual = Galeri::ComputeNorm(&LHS, &LHS_exact);

  if (A->Comm().MyPID() == 0 && verbose) {
    cout << "||x_exact - x||_2 = " << residual << endl;
    cout << *Container;
  }

  bool passed = false;
  if (residual < 1e-5)
    passed = true;

  return(passed);
}

// ======================================================================
int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif
  Epetra_SerialComm SerialComm;

  verbose = (Comm.MyPID() == 0);

  const int nx = 30;
  Teuchos::ParameterList GaleriList;
  GaleriList.set("n", nx * nx);
  GaleriList.set("nx", nx);
  GaleriList.set("ny", nx);

  Teuchos::RefCountPtr<Epetra_Map> Map = Teuchos::rcp( Galeri::CreateMap("Linear", Comm, GaleriList) );
  Teuchos::RefCountPtr<Epetra_RowMatrix> Matrix = Teuchos::rcp( Galeri::CreateCrsMatrix("Laplace2D", &*Map, GaleriList) );
  
  Teuchos::RefCountPtr<Ifpack_LocalFilter> LocalMatrix = Teuchos::rcp( new Ifpack_LocalFilter(Matrix) );
  int TestPassed = true;

  if (!TestContainer("dense",LocalMatrix))  TestPassed = false;
  if (!TestContainer("sparse",LocalMatrix)) TestPassed = false;

  if (TestPassed)
  {
    if (verbose)
      cout << "Test `TestContainer.exe' passed!" << endl;
  }
  else {
    cout << "Test `TestContainer.exe' FAILED!" << endl;
    exit(EXIT_FAILURE);
  }
  
#ifdef HAVE_MPI
  MPI_Finalize(); 
#endif

  return(EXIT_SUCCESS);
}
