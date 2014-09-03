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
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Map.h"
#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"
#include "Galeri_Utils.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Ifpack_PointRelaxation.h"
#include "Ifpack_BlockRelaxation.h"
#include "Ifpack_SparseContainer.h"
#include "Ifpack_TriDiContainer.h"

#include "Ifpack_Amesos.h"
#include "AztecOO.h"

static bool verbose = false;
static bool SymmetricGallery = false;
static bool Solver = AZ_gmres;
const int NumVectors = 3;

// ====================================================================== 
bool ComparePointAndBlock(string PrecType, const Teuchos::RefCountPtr<Epetra_RowMatrix>& A, int sweeps)
{
  Epetra_MultiVector RHS(A->RowMatrixRowMap(), NumVectors);
  Epetra_MultiVector LHS(A->RowMatrixRowMap(), NumVectors);
  LHS.PutScalar(0.0); RHS.Random();

  Epetra_LinearProblem Problem(&*A, &LHS, &RHS);

  // Set up the list
  Teuchos::ParameterList List;
  List.set("relaxation: damping factor", 1.0);
  List.set("relaxation: type", PrecType);
  List.set("relaxation: sweeps",sweeps);
  List.set("partitioner: type", "linear");
  List.set("partitioner: local parts", A->NumMyRows());

  int ItersPoint, ItersBlock;

  // ================================================== //
  // get the number of iterations with point relaxation //
  // ================================================== //
  {
    RHS.PutScalar(1.0);
    LHS.PutScalar(0.0);

    Ifpack_PointRelaxation Point(&*A);
    Point.SetParameters(List);
    Point.Compute();

    // set AztecOO solver object
    AztecOO AztecOOSolver(Problem);
    AztecOOSolver.SetAztecOption(AZ_solver,Solver);
    if (verbose)
      AztecOOSolver.SetAztecOption(AZ_output,32);
    else
      AztecOOSolver.SetAztecOption(AZ_output,AZ_none);
    AztecOOSolver.SetPrecOperator(&Point);

    AztecOOSolver.Iterate(1550,1e-2);

    double TrueResidual = AztecOOSolver.TrueResidual();
    ItersPoint = AztecOOSolver.NumIters();
    // some output
    if (verbose && Problem.GetMatrix()->Comm().MyPID() == 0) {
      cout << "Iterations  = " << ItersPoint << endl;
      cout << "Norm of the true residual = " << TrueResidual << endl;
    }
  }

  // ================================================== //
  // get the number of iterations with block relaxation //
  // ================================================== //
  {

    RHS.PutScalar(1.0);
    LHS.PutScalar(0.0);

    Ifpack_BlockRelaxation<Ifpack_SparseContainer<Ifpack_Amesos> > Block(&*A);
    Block.SetParameters(List);
    Block.Compute();

    // set AztecOO solver object
    AztecOO AztecOOSolver(Problem);
    AztecOOSolver.SetAztecOption(AZ_solver,Solver);
    if (verbose)
      AztecOOSolver.SetAztecOption(AZ_output,32);
    else
      AztecOOSolver.SetAztecOption(AZ_output,AZ_none);
    AztecOOSolver.SetPrecOperator(&Block);

    AztecOOSolver.Iterate(1550,1e-2);

    double TrueResidual = AztecOOSolver.TrueResidual();
    ItersBlock = AztecOOSolver.NumIters();
    // some output
    if (verbose && Problem.GetMatrix()->Comm().MyPID() == 0) {
      cout << "Iterations " << ItersBlock << endl;
      cout << "Norm of the true residual = " << TrueResidual << endl;
    }
  }

  int diff = ItersPoint - ItersBlock;
  if (diff < 0) diff = -diff;
    
  if (diff > 10)
  {
    if (verbose)
      cout << "ComparePointandBlock TEST FAILED!" << endl;
    return(false);
  }
  else {
    if (verbose)
      cout << "ComparePointandBlock TEST PASSED" << endl;
    return(true);
  }
}


// ====================================================================== 
int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  verbose = (Comm.MyPID() == 0);

  int nx = 60;

  for (int i = 1 ; i < argc ; ++i) {
    if (strcmp(argv[i],"-s") == 0) {
      SymmetricGallery = true;
      Solver = AZ_cg;
    }
    if(strcmp(argv[i],"-n") == 0  && i+1 < argc)  {
      i++;
      nx = atoi(argv[i]);
    }

  }

  // size of the global matrix. 
  Teuchos::ParameterList GaleriList;

  GaleriList.set("nx", nx);
  GaleriList.set("ny", nx * Comm.NumProc());
  GaleriList.set("mx", 1);
  GaleriList.set("my", Comm.NumProc());
  Teuchos::RefCountPtr<Epetra_Map> Map = Teuchos::rcp( Galeri::CreateMap("Cartesian2D", Comm, GaleriList) );
  Teuchos::RefCountPtr<Epetra_CrsMatrix> A;
  if (SymmetricGallery)
    A = Teuchos::rcp( Galeri::CreateCrsMatrix("Laplace2D", &*Map, GaleriList) );
  else
    A = Teuchos::rcp( Galeri::CreateCrsMatrix("Recirc2D", &*Map, GaleriList) );

  // coordinates
  Teuchos::RCP<Epetra_MultiVector> coord = Teuchos::rcp( Galeri::CreateCartesianCoordinates("2D",&*Map,GaleriList));

  // test the preconditioner
  int TestPassed = true;


  // ================================== //
  // compare point and block relaxation //
  // ================================== //



  TestPassed = TestPassed && 
    ComparePointAndBlock("Jacobi",A,10);


  TestPassed = TestPassed && 
    ComparePointAndBlock("symmetric Gauss-Seidel",A,10);

  if (!SymmetricGallery) {

    TestPassed = TestPassed && 
      ComparePointAndBlock("Gauss-Seidel",A,10);
  }

  if (!TestPassed) {
    cout << "Test `Performance.exe' failed!" << endl;
    exit(EXIT_FAILURE);
  }
  
#ifdef HAVE_MPI
  MPI_Finalize(); 
#endif

  cout << endl;
  cout << "Test `Performance.exe' passed!" << endl;
  cout << endl;
  return(EXIT_SUCCESS);
}
