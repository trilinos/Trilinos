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
#ifdef HAVE_IFPACK_AZTECOO 
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
#include "Ifpack_PointRelaxation.h"
#include "Ifpack_BlockRelaxation.h"
#include "Ifpack_DenseContainer.h"
#include "Ifpack_AdditiveSchwarz.h"
#include "AztecOO.h"

using namespace Trilinos_Util;
static bool verbose = false;
static bool SymmetricGallery = false;
static bool Solver = AZ_gmres;


// ====================================================================== 
int CompareBlockOverlap(CrsMatrixGallery& Gallery, int Overlap)
{

  Epetra_RowMatrix* A = Gallery.GetMatrix();
  Epetra_LinearProblem* Problem = Gallery.GetLinearProblem();

  Epetra_MultiVector& RHS = *(Problem->GetRHS());
  Epetra_MultiVector& LHS = *(Problem->GetLHS());

  Teuchos::ParameterList List;
  List.set("relaxation: damping factor", 1.0);
  List.set("relaxation: type", "Jacobi");
  List.set("relaxation: sweeps",1);
  List.set("partitioner: overlap", Overlap);
  List.set("partitioner: type", "greedy");
  List.set("partitioner: local parts", 16);

  RHS.PutScalar(1.0);
  LHS.PutScalar(0.0);

  Ifpack_AdditiveSchwarz<Ifpack_BlockRelaxation<Ifpack_DenseContainer> > Prec(A);
  Prec.SetParameters(List);
  Prec.Compute();

  // set AztecOO solver object
  AztecOO AztecOOSolver(*Problem);
  AztecOOSolver.SetAztecOption(AZ_solver,Solver);
  if (verbose)
    AztecOOSolver.SetAztecOption(AZ_output,32);
  else
    AztecOOSolver.SetAztecOption(AZ_output,AZ_none);
  AztecOOSolver.SetPrecOperator(&Prec);

  AztecOOSolver.Iterate(1550,1e-8);

  return(AztecOOSolver.NumIters());
}

// ====================================================================== 
int CompareBlockSizes(string PrecType,
                      CrsMatrixGallery& Gallery, int NumParts)
{

  Epetra_RowMatrix* A = Gallery.GetMatrix();
  Epetra_LinearProblem* Problem = Gallery.GetLinearProblem();

  Epetra_MultiVector& RHS = *(Problem->GetRHS());
  Epetra_MultiVector& LHS = *(Problem->GetLHS());

  Teuchos::ParameterList List;
  List.set("relaxation: damping factor", 1.0);
  List.set("relaxation: type", PrecType);
  List.set("relaxation: sweeps",1);
  List.set("partitioner: type", "greedy");
  List.set("partitioner: local parts", NumParts);

  RHS.PutScalar(1.0);
  LHS.PutScalar(0.0);

  Ifpack_AdditiveSchwarz<Ifpack_BlockRelaxation<Ifpack_DenseContainer> > Prec(A);
  Prec.SetParameters(List);
  Prec.Compute();

  // set AztecOO solver object
  AztecOO AztecOOSolver(*Problem);
  AztecOOSolver.SetAztecOption(AZ_solver,Solver);
  if (verbose)
    AztecOOSolver.SetAztecOption(AZ_output,32);
  else
    AztecOOSolver.SetAztecOption(AZ_output,AZ_none);
  AztecOOSolver.SetPrecOperator(&Prec);

  AztecOOSolver.Iterate(1550,1e-8);

  return(AztecOOSolver.NumIters());
}

// ====================================================================== 
bool ComparePointAndBlock(string PrecType,
                          CrsMatrixGallery& Gallery, int sweeps)
{

  Epetra_RowMatrix* A = Gallery.GetMatrix();
  Epetra_LinearProblem* Problem = Gallery.GetLinearProblem();

  Epetra_MultiVector& RHS = *(Problem->GetRHS());
  Epetra_MultiVector& LHS = *(Problem->GetLHS());

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

    Ifpack_AdditiveSchwarz<Ifpack_PointRelaxation> Point(A);
    Point.SetParameters(List);
    Point.Compute();

    // set AztecOO solver object
    AztecOO AztecOOSolver(*Problem);
    AztecOOSolver.SetAztecOption(AZ_solver,Solver);
    if (verbose)
      AztecOOSolver.SetAztecOption(AZ_output,32);
    else
      AztecOOSolver.SetAztecOption(AZ_output,AZ_none);
    AztecOOSolver.SetPrecOperator(&Point);

    AztecOOSolver.Iterate(1550,1e-8);

    double TrueResidual = AztecOOSolver.TrueResidual();
    ItersPoint = AztecOOSolver.NumIters();
    // some output
    if (verbose && Problem->GetMatrix()->Comm().MyPID() == 0) {
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

    Ifpack_AdditiveSchwarz<Ifpack_BlockRelaxation<Ifpack_DenseContainer> > Block(A);
    Block.SetParameters(List);
    Block.Compute();

    // set AztecOO solver object
    AztecOO AztecOOSolver(*Problem);
    AztecOOSolver.SetAztecOption(AZ_solver,Solver);
    if (verbose)
      AztecOOSolver.SetAztecOption(AZ_output,32);
    else
      AztecOOSolver.SetAztecOption(AZ_output,AZ_none);
    AztecOOSolver.SetPrecOperator(&Block);

    AztecOOSolver.Iterate(1550,1e-8);

    double TrueResidual = AztecOOSolver.TrueResidual();
    ItersBlock = AztecOOSolver.NumIters();
    // some output
    if (verbose && Problem->GetMatrix()->Comm().MyPID() == 0) {
      cout << "Iterations " << ItersBlock << endl;
      cout << "Norm of the true residual = " << TrueResidual << endl;
    }
  }

  if (ItersPoint != ItersBlock) {
    if (verbose)
      cout << "TEST FAILED!" << endl;
    return(false);
  }
  else {
    if (verbose)
      cout << "TEST PASSED" << endl;
    return(true);
  }
}

// ====================================================================== 
bool KrylovTest(string PrecType,
                CrsMatrixGallery& Gallery)
{

  // The following methods of CrsMatrixGallery are used to get pointers
  // to internally stored Epetra_RowMatrix and Epetra_LinearProblem.
  Epetra_RowMatrix* A = Gallery.GetMatrix();
  Epetra_LinearProblem* Problem = Gallery.GetLinearProblem();

  Epetra_MultiVector& LHS = *(Problem->GetLHS());
  LHS.PutScalar(0.0);

  // Set up the list
  Teuchos::ParameterList List;
  List.set("relaxation: damping factor", 1.0);
  List.set("relaxation: type", PrecType);

  int Iters1, Iters10;

  // ============================================== //
  // get the number of iterations with 1 sweep only //
  // ============================================== //
  {

    List.set("relaxation: sweeps",1);
    Ifpack_AdditiveSchwarz<Ifpack_PointRelaxation> Point(A);
    Point.SetParameters(List);
    Point.Compute();

    // set AztecOO solver object
    AztecOO AztecOOSolver(*Problem);
    AztecOOSolver.SetAztecOption(AZ_solver,Solver);
    AztecOOSolver.SetAztecOption(AZ_output,AZ_none);
    AztecOOSolver.SetPrecOperator(&Point);

    AztecOOSolver.Iterate(1550,1e-8);

    double TrueResidual = AztecOOSolver.TrueResidual();
    // some output
    if (verbose && Problem->GetMatrix()->Comm().MyPID() == 0) {
      cout << "Norm of the true residual = " << TrueResidual << endl;
    }
    Iters1 = AztecOOSolver.NumIters();
  }

  // ======================================================== //
  // now re-run with 10 sweeps, solver should converge faster
  // ======================================================== //
  {
    List.set("relaxation: sweeps",10);
    Ifpack_AdditiveSchwarz<Ifpack_PointRelaxation> Point(A);
    Point.SetParameters(List);
    Point.Compute();
    LHS.PutScalar(0.0);

    // set AztecOO solver object
    AztecOO AztecOOSolver(*Problem);
    AztecOOSolver.SetAztecOption(AZ_solver,Solver);
    AztecOOSolver.SetAztecOption(AZ_output,AZ_none);
    AztecOOSolver.SetPrecOperator(&Point);
    AztecOOSolver.Iterate(1550,1e-8);

    double TrueResidual = AztecOOSolver.TrueResidual();
    // some output
    if (verbose && Problem->GetMatrix()->Comm().MyPID() == 0) {
      cout << "Norm of the true residual = " << TrueResidual << endl;
    }
    Iters10 = AztecOOSolver.NumIters();
  }

  if (verbose) {
    cout << "Iters_1 = " << Iters1 << ", Iters_10 = " << Iters10 << endl;
    cout << "(second number should be smaller than first one)" << endl;
  }

  if (Iters10 > Iters1) {
    if (verbose)
      cout << "TEST FAILED!" << endl;
    return(false);
  }
  else {
    if (verbose)
      cout << "TEST PASSED" << endl;
    return(true);
  }
}

// ====================================================================== 
bool BasicTest(string PrecType,
               CrsMatrixGallery& Gallery)
{

  // The following methods of CrsMatrixGallery are used to get pointers
  // to internally stored Epetra_RowMatrix and Epetra_LinearProblem.
  Epetra_RowMatrix* A = Gallery.GetMatrix();
  Epetra_LinearProblem* Problem = Gallery.GetLinearProblem();

  Epetra_MultiVector& RHS = *(Problem->GetRHS());
  Epetra_MultiVector& LHS = *(Problem->GetLHS());

  LHS.PutScalar(0.0);

  // Set up the list
  Teuchos::ParameterList List;
  List.set("relaxation: damping factor", 1.0);
  List.set("relaxation: sweeps",1550);
  List.set("relaxation: type", PrecType);

  Ifpack_PointRelaxation Point(A);

  Point.SetParameters(List);
  Point.Compute();
  // use the preconditioner as solver, with 1550 iterations
  Point.ApplyInverse(RHS,LHS);

  // compute the real residual

  double residual, diff;
  Gallery.ComputeResidual(&residual);
  Gallery.ComputeDiffBetweenStartingAndExactSolutions(&diff);

  if (verbose && A->Comm().MyPID()==0) {
    cout << "||b-Ax||_2 = " << residual << endl;
    cout << "||x_exact - x||_2 = " << diff << endl;
  }
  
  // Jacobi is very slow to converge here
  if (residual < 1e-2) {
    if (verbose)
      cout << "Test passed" << endl;
    return(true);
  }
  else {
    if (verbose)
      cout << "Test failed!" << endl;
    return(false);
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

  bool verbose = (Comm.MyPID() == 0);

  for (int i = 1 ; i < argc ; ++i) {
    if (strcmp(argv[i],"-s") == 0) {
      SymmetricGallery = true;
      Solver = AZ_cg;
    }
  }

  // size of the global matrix. 
  const int NumPoints = 900;

  CrsMatrixGallery Gallery("", Comm);
  Gallery.Set("problem_size", NumPoints);
  Gallery.Set("map_type", "linear");
  if (SymmetricGallery)
    Gallery.Set("problem_type", "laplace_2d");
  else
    Gallery.Set("problem_type", "recirc_2d");

  // test the preconditioner
  int TestPassed = true;

  // ======================================== //
  // first verify that we can get convergence //
  // with all point relaxation methods        //
  // ======================================== //
  if(!BasicTest("Jacobi",Gallery))
    TestPassed = false;

  if(!BasicTest("symmetric Gauss-Seidel",Gallery))
    TestPassed = false;

  if (!SymmetricGallery) {
    if(!BasicTest("Gauss-Seidel",Gallery))
      TestPassed = false;
  }

  // ============================= //
  // check uses as preconditioners //
  // ============================= //
  
  if(!KrylovTest("symmetric Gauss-Seidel",Gallery))
    TestPassed = false;

  if (!SymmetricGallery) {
    if(!KrylovTest("Gauss-Seidel",Gallery))
      TestPassed = false;
  }

  // ================================== //
  // compare point and block relaxation //
  // ================================== //

  TestPassed = TestPassed && 
    ComparePointAndBlock("Jacobi",Gallery,1);

  TestPassed = TestPassed && 
    ComparePointAndBlock("Jacobi",Gallery,10);

  TestPassed = TestPassed && 
    ComparePointAndBlock("symmetric Gauss-Seidel",Gallery,1);

  TestPassed = TestPassed && 
    ComparePointAndBlock("symmetric Gauss-Seidel",Gallery,10);

  if (!SymmetricGallery) {
    TestPassed = TestPassed && 
      ComparePointAndBlock("Gauss-Seidel",Gallery,1);

    TestPassed = TestPassed && 
      ComparePointAndBlock("Gauss-Seidel",Gallery,10);
  }

  // ============================ //
  // verify effect of # of blocks //
  // ============================ //
  
  {
    int Iters4, Iters8, Iters16;
    Iters4 = CompareBlockSizes("Jacobi",Gallery,4);
    Iters8 = CompareBlockSizes("Jacobi",Gallery,8);
    Iters16 = CompareBlockSizes("Jacobi",Gallery,16);
    if ((Iters16 > Iters8) && (Iters8 > Iters4)) {
      if (verbose)
        cout << "Test passed" << endl;
    }
    else {
      if (verbose) 
        cout << "TEST FAILED!" << endl;
      TestPassed = TestPassed && false;
    }
  }

  // ================================== //
  // verify effect of overlap in Jacobi //
  // ================================== //
  
  {
    int Iters0, Iters2, Iters4;
    Iters0 = CompareBlockOverlap(Gallery,0);
    Iters2 = CompareBlockOverlap(Gallery,2);
    Iters4 = CompareBlockOverlap(Gallery,4);
    if ((Iters4 < Iters2) && (Iters2 < Iters0)) {
      if (verbose)
        cout << "Test passed" << endl;
    }
    else {
      if (verbose) 
        cout << "TEST FAILED!" << endl;
      TestPassed = TestPassed && false;
    }
  }

  // ============ //
  // final output //
  // ============ //

  if (!TestPassed) {
    cout << "Test `PointPreconditioner.exe' FAILED!" << endl;
    exit(EXIT_FAILURE);
  }
  
#ifdef HAVE_MPI
  MPI_Finalize(); 
#endif

  cout << "Test `PointPreconditioner.exe' passed!" << endl;
  exit(EXIT_SUCCESS);
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

  puts("please configure IFPACK with --eanble-aztecoo --enable-teuchos");
  puts("to run this test");

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif
  return(EXIT_SUCCESS);
}

#endif
