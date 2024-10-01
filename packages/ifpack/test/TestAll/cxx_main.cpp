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
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_LinearProblem.h"
#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Ifpack_AdditiveSchwarz.h"
#include "AztecOO.h"
#include "Ifpack_PointRelaxation.h"
#include "Ifpack_BlockRelaxation.h"
#include "Ifpack_DenseContainer.h"
#include "Ifpack_SparseContainer.h"
#include "Ifpack_Amesos.h"
#include "Ifpack_Utils.h"
#include "Ifpack_Chebyshev.h"
#include "Ifpack_Polynomial.h"
#include "Ifpack_Krylov.h"

template <class CommType>
Teuchos::RefCountPtr<Epetra_RowMatrix>
createTriDiagMat(int NumGlobalElements, CommType Comm, bool str_singular, bool num_singular) {
  // Construct a Map that puts approximatively the same number of 
  // equations on each processor. `0' is the index base (that is,
  // numbering starts from 0.
  Epetra_Map Map(NumGlobalElements, 0, Comm);

  // Create an empty EpetraCrsMatrix 
  Teuchos::RefCountPtr<Epetra_CrsMatrix> A = Teuchos::rcp( new Epetra_CrsMatrix(Copy, Map, 0));

  // Create the structure of the matrix (tridiagonal)
  int NumMyElements = Map.NumMyElements();

  // Add  rows one-at-a-time
  // Need some vectors to help

  double Values[3];
  // Right now, we put zeros only in the matrix.
  int Indices[3];
  int NumEntries;
  /// global ID's of local ID's
  int* MyGlobalElements = Map.MyGlobalElements();

  // At this point we simply set the nonzero structure of A.
  // Actual values will be inserted later (now all zeros)
  for (int i = 0; i < NumMyElements; i++)
  {
    if (MyGlobalElements[i] == 0)
    {
      if (str_singular) {
        NumEntries = 0;
      } else {
        Indices[0] = 0;
        Indices[1] = 1;
        if (num_singular) {
          Values[0] = 0.0;
          Values[1] = 0.0;
        } else {
          Values[0] = 2.0;
          Values[1] = 1.0;
        }
        NumEntries = 2;
      }
    }
    else if (MyGlobalElements[i] == NumGlobalElements-1)
    {
      Indices[0] = NumGlobalElements-1;
      Indices[1] = NumGlobalElements-2;
      Values[0] = 2.0;
      Values[1] = 1.0;
      NumEntries = 2;
    }
    else
    {
      Indices[0] = MyGlobalElements[i]-1;
      Indices[1] = MyGlobalElements[i];
      Indices[2] = MyGlobalElements[i]+1;
      Values[0] = 1.0;
      Values[1] = 2.0;
      Values[2] = 1.0;
      NumEntries = 3;
    }

    if (NumEntries > 0)
    A->InsertGlobalValues(MyGlobalElements[i],
                          NumEntries, Values, Indices);
  }

  // Finish up.
  A->FillComplete();
  //A->Print(std::cout);

  return A;
}

template <class T>
bool Test(const Teuchos::RefCountPtr<Epetra_RowMatrix>& Matrix, Teuchos::ParameterList& List)
{

  int NumVectors = 1;
  bool UseTranspose = false;

  Epetra_MultiVector LHS(Matrix->OperatorDomainMap(),NumVectors);
  Epetra_MultiVector RHS(Matrix->OperatorRangeMap(),NumVectors);
  Epetra_MultiVector LHSexact(Matrix->OperatorDomainMap(),NumVectors);

  LHS.PutScalar(0.0);
  LHSexact.Random();
  Matrix->Multiply(UseTranspose,LHSexact,RHS);

  Epetra_LinearProblem Problem(&*Matrix,&LHS,&RHS);

  Teuchos::RefCountPtr<T> Prec;
  
  Prec = Teuchos::rcp( new T(&*Matrix) );
  assert(Prec != Teuchos::null);

  IFPACK_CHK_ERRB(Prec->SetParameters(List));
  IFPACK_CHK_ERRB(Prec->Initialize());
  IFPACK_CHK_ERRB(Prec->Compute());

  // create the AztecOO solver
  AztecOO AztecOOSolver(Problem);

  // specify solver
  AztecOOSolver.SetAztecOption(AZ_solver,AZ_gmres);
  AztecOOSolver.SetAztecOption(AZ_output,32);

  AztecOOSolver.SetPrecOperator(&*Prec);

  // solver. The solver should converge in one iteration,
  // or maximum two (numerical errors)
  AztecOOSolver.Iterate(1550,1e-8);

  cout << *Prec;
  
  std::vector<double> Norm(NumVectors);
  LHS.Update(1.0,LHSexact,-1.0);
  LHS.Norm2(&Norm[0]);
  for (int i = 0 ; i < NumVectors ; ++i) {
    cout << "Norm[" << i << "] = " << Norm[i] << endl;
    if (Norm[i] > 1e-3)
      return(false);
  }
  return(true);

}

int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  bool verbose = (Comm.MyPID() == 0);

  Teuchos::ParameterList GaleriList;
  const int nx = 30;
  GaleriList.set("n", nx * nx);
  GaleriList.set("nx", nx);
  GaleriList.set("ny", nx);
  Teuchos::RefCountPtr<Epetra_Map> Map = Teuchos::rcp( Galeri::CreateMap("Linear", Comm, GaleriList) );
  Teuchos::RefCountPtr<Epetra_RowMatrix> Matrix = Teuchos::rcp( Galeri::CreateCrsMatrix("Laplace2D", &*Map, GaleriList) );

  Teuchos::ParameterList List, DefaultList;

  // test the preconditioner
  int TestPassed = true;

  if (!Test<Ifpack_Chebyshev>(Matrix,List)) 
  {
    TestPassed = false;
  }

  List.set("polynomial: degree",3);
  if (!Test<Ifpack_Polynomial>(Matrix,List))
  {
    TestPassed = false;
  }

  List = DefaultList;
  List.set("krylov: tolerance", 1e-14);
  List.set("krylov: iterations", 100);
  List.set("krylov: preconditioner", 2);
  List.set("krylov: block size", 9);
  List.set("krylov: number of sweeps", 2);
  if (!Test<Ifpack_Krylov>(Matrix,List))
  {
    TestPassed = false;
  }

  if (!Test< Ifpack_AdditiveSchwarz<Ifpack_Krylov> >(Matrix,List)) {
    TestPassed = false;
  }

  if (!Test<Ifpack_Amesos>(Matrix,List)) 
  {
    TestPassed = false;
  }

  // FIXME
#if 0
  if (!Test<Ifpack_AdditiveSchwarz<Ifpack_BlockRelaxation<Ifpack_DenseContainer> > >(Matrix,List)) {
    TestPassed = false;
  }
#endif

  // this is ok as long as just one sweep is applied
  List = DefaultList;
  List.set("relaxation: type", "Gauss-Seidel");
  if (!Test<Ifpack_PointRelaxation>(Matrix,List)) {
    TestPassed = false;
  }

  // this is ok as long as just one sweep is applied
  List = DefaultList;
  List.set("relaxation: type", "symmetric Gauss-Seidel");
  List.set("relaxation: sweeps", 5);
  List.set("partitioner: local parts", 128);
  List.set("partitioner: type", "linear");
  if (!Test<Ifpack_BlockRelaxation<Ifpack_DenseContainer> >(Matrix,List)) {
    TestPassed = false;
  }
 
  // this is ok as long as just one sweep is applied
  List = DefaultList;
  List.set("relaxation: type", "symmetric Gauss-Seidel");
  List.set("partitioner: local parts", 128);
  List.set("partitioner: type", "linear");
  if (!Test<Ifpack_BlockRelaxation<Ifpack_SparseContainer<Ifpack_Amesos> > >(Matrix,List)) {
    TestPassed = false;
  }

  // this is ok as long as just one sweep is applied
  List = DefaultList;
  List.set("relaxation: type", "symmetric Gauss-Seidel");
  List.set("partitioner: local parts", 128);
  List.set("partitioner: type", "linear");
  if (!Test<Ifpack_AdditiveSchwarz<Ifpack_BlockRelaxation<Ifpack_SparseContainer<Ifpack_Amesos> > > >(Matrix,List)) {
    TestPassed = false;
  }

#ifdef HAVE_IFPACK_AMESOS
  // Additive Schwarz with local Amesos
  // Amesos should fail on MPI-0.
  // So, these tests should fail,
  // but are designed to check that error is propagated to all MPI processes
  // instead of just failing on MPI-0, causing deadlock
  bool check_for_global_error = false;
  if (check_for_global_error) {
    // structurally singular case
    List = DefaultList;
    bool num_singular = false;
    bool str_singular = true;
    Matrix = createTriDiagMat(10, Comm, str_singular, num_singular);
    if (Test<Ifpack_AdditiveSchwarz<Ifpack_Amesos>>(Matrix,List)) {
      TestPassed = false;
    }
    // numerically singular case
    num_singular = true;
    str_singular = false;
    Matrix = createTriDiagMat(10, Comm, str_singular, num_singular);
    if (Test<Ifpack_AdditiveSchwarz<Ifpack_Amesos>>(Matrix,List)) {
      TestPassed = false;
    }
  }
#endif

  if (!TestPassed) {
    cerr << "Test `TestAll.exe' FAILED!" << endl;
    exit(EXIT_FAILURE);
  }

#ifdef HAVE_MPI
  MPI_Finalize(); 
#endif
  if (verbose)
    cout << "Test `TestAll.exe' passed!" << endl;

  return(EXIT_SUCCESS);
}
