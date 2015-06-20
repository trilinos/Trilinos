/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

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

  // add TriDi

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
