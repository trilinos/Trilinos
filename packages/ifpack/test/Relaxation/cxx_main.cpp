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
bool TestTriDiVariableBlocking(const Epetra_Comm & Comm) {
  // Basically each processor gets this 5x5 block lower-triangular matrix:
  //
  // [ 2 -1  0  0  0 ;...
  // [-1  2  0  0  0 ;...
  // [ 0 -1  3 -1  0 ;...
  // [ 0  0 -1  3 -1 ;...
  // [ 0  0  0 -1  2  ];
  //

  Epetra_Map RowMap(-1,5,0,Comm); // 5 rows per proc

  Epetra_CrsMatrix A(Copy,RowMap,0);
  
  int num_entries;
  int indices[5];
  double values[5];
  int rb = RowMap.GID(0);

  /*** Fill RHS / LHS ***/
  Epetra_Vector rhs(RowMap), lhs(RowMap), exact_soln(RowMap);
  rhs.PutScalar(2.0);
  lhs.PutScalar(0.0);
  exact_soln.PutScalar(2.0);

  /*** Fill Matrix ****/
  // Row 0 
  num_entries=2;
  indices[0]=rb; indices[1]=rb+1;
  values[0] =2; values[1] =-1;
  A.InsertGlobalValues(rb,num_entries,&values[0],&indices[0]);

  // Row 1
  num_entries=2;
  indices[0]=rb; indices[1]=rb+1; 
  values[0] =-1; values[1] =2;
  A.InsertGlobalValues(rb+1,num_entries,&values[0],&indices[0]);

  // Row 2
  num_entries=3;
  indices[0]=rb+1; indices[1]=rb+2; indices[2]=rb+3; 
  values[0] =-1;   values[1] = 3;   values[2] =-1;   
  A.InsertGlobalValues(rb+2,num_entries,&values[0],&indices[0]);

  // Row 3
  num_entries=3;
  indices[0]=rb+2; indices[1]=rb+3; indices[2]=rb+4;
  values[0] =-1;   values[1] = 3;   values[2] =-1;
  A.InsertGlobalValues(rb+3,num_entries,&values[0],&indices[0]);

  // Row 4
  num_entries=2;
  indices[0]=rb+3; indices[1]=rb+4;
  values[0] =-1;   values[1] = 2;
  A.InsertGlobalValues(rb+4,num_entries,&values[0],&indices[0]); 
  A.FillComplete();

  /* Setup Block Relaxation */
  int PartMap[5]={0,0,1,1,1};

  Teuchos::ParameterList ilist;
  ilist.set("partitioner: type","user");
  ilist.set("partitioner: map",&PartMap[0]);
  ilist.set("partitioner: local parts",2);
  ilist.set("relaxation: sweeps",1);
  ilist.set("relaxation: type","Gauss-Seidel");

  Ifpack_BlockRelaxation<Ifpack_TriDiContainer> TDRelax(&A);

  TDRelax.SetParameters(ilist);
  TDRelax.Initialize();
  TDRelax.Compute();
  TDRelax.ApplyInverse(rhs,lhs);
  
  double norm;
  lhs.Update(1.0,exact_soln,-1.0);
  lhs.Norm2(&norm);

  if(verbose) cout<<"Variable Block Partitioning Test"<<endl;

  if(norm < 1e-14) {
    if(verbose) cout << "Test passed" << endl;
     return true;
  }
  else {
    if(verbose) cout << "Test failed" << endl;
    return false;
  }
}
// ====================================================================== 
bool TestVariableBlocking(const Epetra_Comm & Comm) {
  // Basically each processor gets this 5x5 block lower-triangular matrix:
  //
  // [ 2 -1  0  0  0 ;...
  // [-1  2  0  0  0 ;...
  // [-1 -1  5 -1 -1 ;...
  // [-1 -1 -1  5 -1 ;...
  // [-1 -1 -1 -1  5  ];
  //
  // The nice thing about this matrix is that if the RHS is a constant,the solution is the same constant...

  Epetra_Map RowMap(-1,5,0,Comm); // 5 rows per proc

  Epetra_CrsMatrix A(Copy,RowMap,0);
  
  int num_entries;
  int indices[5];
  double values[5];
  int rb = RowMap.GID(0);

  /*** Fill RHS / LHS ***/
  Epetra_Vector rhs(RowMap), lhs(RowMap), exact_soln(RowMap);
  rhs.PutScalar(2.0);
  lhs.PutScalar(0.0);
  exact_soln.PutScalar(2.0);

  /*** Fill Matrix ****/
  // Row 0 
  num_entries=2;
  indices[0]=rb; indices[1]=rb+1;
  values[0] =2; values[1] =-1;
  A.InsertGlobalValues(rb,num_entries,&values[0],&indices[0]);

  // Row 1
  num_entries=2;
  indices[0]=rb; indices[1]=rb+1;
  values[0] =-1; values[1] =2;
  A.InsertGlobalValues(rb+1,num_entries,&values[0],&indices[0]);

  // Row 2
  num_entries=5;
  indices[0]=rb; indices[1]=rb+1; indices[2]=rb+2; indices[3]=rb+3; indices[4]=rb+4;
  values[0] =-1; values[1] =-1;   values[2] = 5;   values[3] =-1;   values[4] =-1;
  A.InsertGlobalValues(rb+2,num_entries,&values[0],&indices[0]);

  // Row 3
  num_entries=5;
  indices[0]=rb; indices[1]=rb+1; indices[2]=rb+2; indices[3]=rb+3; indices[4]=rb+4;
  values[0] =-1; values[1] =-1;   values[2] =-1;   values[3] = 5;   values[4] =-1;
  A.InsertGlobalValues(rb+3,num_entries,&values[0],&indices[0]);

  // Row 4
  num_entries=5;
  indices[0]=rb; indices[1]=rb+1; indices[2]=rb+2; indices[3]=rb+3; indices[4]=rb+4;
  values[0] =-1; values[1] =-1;   values[2] =-1;   values[3] =-1;   values[4] = 5;
  A.InsertGlobalValues(rb+4,num_entries,&values[0],&indices[0]); 
  A.FillComplete();


  /* Setup Block Relaxation */
  int PartMap[5]={0,0,1,1,1};

  Teuchos::ParameterList ilist;
  ilist.set("partitioner: type","user");
  ilist.set("partitioner: map",&PartMap[0]);
  ilist.set("partitioner: local parts",2);
  ilist.set("relaxation: sweeps",1);
  ilist.set("relaxation: type","Gauss-Seidel");
  Ifpack_BlockRelaxation<Ifpack_DenseContainer> Relax(&A);
  Relax.SetParameters(ilist);
  Relax.Initialize();
  Relax.Compute();

  Relax.ApplyInverse(rhs,lhs);

  
  double norm;
  lhs.Update(1.0,exact_soln,-1.0);
  lhs.Norm2(&norm);

  if(verbose) cout<<"Variable Block Partitioning Test"<<endl;

  if(norm < 1e-14) {
    if(verbose) cout << "Test passed" << endl;
     return true;
  }
  else {
    if(verbose) cout << "Test failed" << endl;
    return false;
  }
}

// ====================================================================== 
int CompareLineSmoother(const Teuchos::RefCountPtr<Epetra_RowMatrix>& A, Teuchos::RCP<Epetra_MultiVector> coord)
{
  Epetra_MultiVector LHS(A->RowMatrixRowMap(), NumVectors);
  Epetra_MultiVector RHS(A->RowMatrixRowMap(), NumVectors);
  LHS.PutScalar(0.0); RHS.Random();

  Epetra_LinearProblem Problem(&*A, &LHS, &RHS);

  Teuchos::ParameterList List;
  List.set("relaxation: damping factor", 1.0);
  List.set("relaxation: type", "symmetric Gauss-Seidel");
  List.set("relaxation: sweeps",1);
  List.set("partitioner: overlap",0);
  List.set("partitioner: type", "line");
  List.set("partitioner: line detection threshold",0.1);
  List.set("partitioner: x-coordinates",&(*coord)[0][0]);
  List.set("partitioner: y-coordinates",&(*coord)[1][0]);
  List.set("partitioner: z-coordinates",(double*) 0);

  RHS.PutScalar(1.0);
  LHS.PutScalar(0.0);

  Ifpack_BlockRelaxation<Ifpack_SparseContainer<Ifpack_Amesos> > Prec(&*A);
  Prec.SetParameters(List);
  Prec.Compute();

  // set AztecOO solver object
  AztecOO AztecOOSolver(Problem);
  AztecOOSolver.SetAztecOption(AZ_solver,Solver);
  if (verbose)
    AztecOOSolver.SetAztecOption(AZ_output,32);
  else
    AztecOOSolver.SetAztecOption(AZ_output,AZ_none);
  AztecOOSolver.SetPrecOperator(&Prec);

  AztecOOSolver.Iterate(1550,1e-5);

  return(AztecOOSolver.NumIters());
}
// ====================================================================== 
int AllSingle(const Teuchos::RefCountPtr<Epetra_RowMatrix>& A, Teuchos::RCP<Epetra_MultiVector> coord)
{
  Epetra_MultiVector LHS(A->RowMatrixRowMap(), NumVectors);
  Epetra_MultiVector RHS(A->RowMatrixRowMap(), NumVectors);
  LHS.PutScalar(0.0); RHS.Random();

  Epetra_LinearProblem Problem(&*A, &LHS, &RHS);

  Teuchos::ParameterList List;
  List.set("relaxation: damping factor", 1.0);
  List.set("relaxation: type", "symmetric Gauss-Seidel");
  List.set("relaxation: sweeps",1);
  List.set("partitioner: overlap",0);
  List.set("partitioner: type", "line");
  List.set("partitioner: line detection threshold",1.0);
  List.set("partitioner: x-coordinates",&(*coord)[0][0]);
  List.set("partitioner: y-coordinates",&(*coord)[1][0]);
  List.set("partitioner: z-coordinates",(double*) 0);

  RHS.PutScalar(1.0);
  LHS.PutScalar(0.0);

  Ifpack_BlockRelaxation<Ifpack_SparseContainer<Ifpack_Amesos> > Prec(&*A);
  Prec.SetParameters(List);
  Prec.Compute();

  // set AztecOO solver object
  AztecOO AztecOOSolver(Problem);
  AztecOOSolver.SetAztecOption(AZ_solver,Solver);
  if (verbose)
    AztecOOSolver.SetAztecOption(AZ_output,32);
  else
    AztecOOSolver.SetAztecOption(AZ_output,AZ_none);
  AztecOOSolver.SetPrecOperator(&Prec);

  AztecOOSolver.Iterate(1550,1e-5);

  printf(" AllSingle  iters %d \n",AztecOOSolver.NumIters());
  return(AztecOOSolver.NumIters());
}

// ====================================================================== 
int CompareBlockOverlap(const Teuchos::RefCountPtr<Epetra_RowMatrix>& A, int Overlap)
{
  Epetra_MultiVector LHS(A->RowMatrixRowMap(), NumVectors);
  Epetra_MultiVector RHS(A->RowMatrixRowMap(), NumVectors);
  LHS.PutScalar(0.0); RHS.Random();

  Epetra_LinearProblem Problem(&*A, &LHS, &RHS);

  Teuchos::ParameterList List;
  List.set("relaxation: damping factor", 1.0);
  List.set("relaxation: type", "Jacobi");
  List.set("relaxation: sweeps",1);
  List.set("partitioner: overlap", Overlap);
  List.set("partitioner: type", "linear");
  List.set("partitioner: local parts", 16);

  RHS.PutScalar(1.0);
  LHS.PutScalar(0.0);

  Ifpack_BlockRelaxation<Ifpack_SparseContainer<Ifpack_Amesos> > Prec(&*A);
  Prec.SetParameters(List);
  Prec.Compute();

  // set AztecOO solver object
  AztecOO AztecOOSolver(Problem);
  AztecOOSolver.SetAztecOption(AZ_solver,Solver);
  if (verbose)
    AztecOOSolver.SetAztecOption(AZ_output,32);
  else
    AztecOOSolver.SetAztecOption(AZ_output,AZ_none);
  AztecOOSolver.SetPrecOperator(&Prec);

  AztecOOSolver.Iterate(1550,1e-5);

  return(AztecOOSolver.NumIters());
}

// ====================================================================== 
int CompareBlockSizes(string PrecType, const Teuchos::RefCountPtr<Epetra_RowMatrix>& A, int NumParts)
{
  Epetra_MultiVector RHS(A->RowMatrixRowMap(), NumVectors);
  Epetra_MultiVector LHS(A->RowMatrixRowMap(), NumVectors);
  LHS.PutScalar(0.0); RHS.Random();

  Epetra_LinearProblem Problem(&*A, &LHS, &RHS);

  Teuchos::ParameterList List;
  List.set("relaxation: damping factor", 1.0);
  List.set("relaxation: type", PrecType);
  List.set("relaxation: sweeps",1);
  List.set("partitioner: type", "linear");
  List.set("partitioner: local parts", NumParts);

  RHS.PutScalar(1.0);
  LHS.PutScalar(0.0);

  Ifpack_BlockRelaxation<Ifpack_SparseContainer<Ifpack_Amesos> > Prec(&*A);
  Prec.SetParameters(List);
  Prec.Compute();

  // set AztecOO solver object
  AztecOO AztecOOSolver(Problem);
  AztecOOSolver.SetAztecOption(AZ_solver,Solver);
  if (verbose)
    AztecOOSolver.SetAztecOption(AZ_output,32);
  else
    AztecOOSolver.SetAztecOption(AZ_output,AZ_none);
  AztecOOSolver.SetPrecOperator(&Prec);

  AztecOOSolver.Iterate(1550,1e-5);

  return(AztecOOSolver.NumIters());
}

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
bool KrylovTest(string PrecType, const Teuchos::RefCountPtr<Epetra_RowMatrix>& A, bool backward, bool reorder=false)
{
  Epetra_MultiVector LHS(A->RowMatrixRowMap(), NumVectors);
  Epetra_MultiVector RHS(A->RowMatrixRowMap(), NumVectors);
  LHS.PutScalar(0.0); RHS.Random();

  Epetra_LinearProblem Problem(&*A, &LHS, &RHS);

  // Set up the list
  Teuchos::ParameterList List;
  List.set("relaxation: damping factor", 1.0);
  List.set("relaxation: type", PrecType);
  if(backward) List.set("relaxation: backward mode",backward);  

  // Reordering if needed
  int NumRows=A->NumMyRows();
  std::vector<int> RowList(NumRows);
  if(reorder) {
    for(int i=0; i<NumRows; i++)
      RowList[i]=i;
    List.set("relaxation: number of local smoothing indices",NumRows);
    List.set("relaxation: local smoothing indices",RowList.size()>0? &RowList[0] : (int*)0);
  }


  int Iters1, Iters10;

  if (verbose) {
    cout << "Krylov test: Using " << PrecType 
         << " with AztecOO" << endl;
  }

  // ============================================== //
  // get the number of iterations with 1 sweep only //
  // ============================================== //
  {

    List.set("relaxation: sweeps",1);
    Ifpack_PointRelaxation Point(&*A);
    Point.SetParameters(List);
    Point.Compute();

    // set AztecOO solver object
    AztecOO AztecOOSolver(Problem);
    AztecOOSolver.SetAztecOption(AZ_solver,Solver);
    AztecOOSolver.SetAztecOption(AZ_output,AZ_none);
    AztecOOSolver.SetPrecOperator(&Point);

    AztecOOSolver.Iterate(1550,1e-5);

    double TrueResidual = AztecOOSolver.TrueResidual();
    // some output
    if (verbose && Problem.GetMatrix()->Comm().MyPID() == 0) {
      cout << "Norm of the true residual = " << TrueResidual << endl;
    }
    Iters1 = AztecOOSolver.NumIters();
  }

  // ======================================================== //
  // now re-run with 10 sweeps, solver should converge faster
  // ======================================================== //
  {
    List.set("relaxation: sweeps",10);
    Ifpack_PointRelaxation Point(&*A);
    Point.SetParameters(List);
    Point.Compute();
    LHS.PutScalar(0.0);

    // set AztecOO solver object
    AztecOO AztecOOSolver(Problem);
    AztecOOSolver.SetAztecOption(AZ_solver,Solver);
    AztecOOSolver.SetAztecOption(AZ_output,AZ_none);
    AztecOOSolver.SetPrecOperator(&Point);
    AztecOOSolver.Iterate(1550,1e-5);

    double TrueResidual = AztecOOSolver.TrueResidual();
    // some output
    if (verbose && Problem.GetMatrix()->Comm().MyPID() == 0) {
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
      cout << "KrylovTest TEST FAILED!" << endl;
    return(false);
  }
  else {
    if (verbose)
      cout << "KrylovTest TEST PASSED" << endl;
    return(true);
  }
}

// ====================================================================== 
bool BasicTest(string PrecType, const Teuchos::RefCountPtr<Epetra_RowMatrix>& A,bool backward, bool reorder=false)
{
  Epetra_MultiVector LHS(A->RowMatrixRowMap(), NumVectors);
  Epetra_MultiVector RHS(A->RowMatrixRowMap(), NumVectors);
  LHS.PutScalar(0.0); RHS.Random();

  double starting_residual = Galeri::ComputeNorm(&*A, &LHS, &RHS);
  Epetra_LinearProblem Problem(&*A, &LHS, &RHS);

  // Set up the list
  Teuchos::ParameterList List;
  List.set("relaxation: damping factor", 1.0);
  List.set("relaxation: sweeps",1550);
  List.set("relaxation: type", PrecType);
  if(backward) List.set("relaxation: backward mode",backward);

  // Reordering if needed
  int NumRows=A->NumMyRows();
  std::vector<int> RowList(NumRows);
  if(reorder) {
    for(int i=0; i<NumRows; i++)
      RowList[i]=i;
    List.set("relaxation: number of local smoothing indices",NumRows);
    List.set("relaxation: local smoothing indices",RowList.size()>0? &RowList[0] : (int*)0);
  }

  Ifpack_PointRelaxation Point(&*A);

  Point.SetParameters(List);
  Point.Compute();
  // use the preconditioner as solver, with 1550 iterations
  Point.ApplyInverse(RHS,LHS);

  // compute the real residual

  double residual = Galeri::ComputeNorm(&*A, &LHS, &RHS);
  
  if (A->Comm().MyPID() == 0 && verbose)
    cout << "||A * x - b||_2 (scaled) = " << residual / starting_residual << endl;
  
  // Jacobi is very slow to converge here
  if (residual / starting_residual < 1e-2) {
    if (verbose)
      cout << "BasicTest Test passed" << endl;
    return(true);
  }
  else {
    if (verbose)
      cout << "BasicTest Test failed!" << endl;
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

  verbose = (Comm.MyPID() == 0);

  for (int i = 1 ; i < argc ; ++i) {
    if (strcmp(argv[i],"-s") == 0) {
      SymmetricGallery = true;
      Solver = AZ_cg;
    }
  }

  // size of the global matrix. 
  Teuchos::ParameterList GaleriList;
  int nx = 3000; 
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

  // ======================================== //
  // first verify that we can get convergence //
  // with all point relaxation methods        //
  // ======================================== //

  if(!BasicTest("Jacobi",A,false))
    TestPassed = false;

  if(!BasicTest("symmetric Gauss-Seidel",A,false))
    TestPassed = false;

  if(!BasicTest("symmetric Gauss-Seidel",A,false,true))
    TestPassed = false;

  if (!SymmetricGallery) {
    if(!BasicTest("Gauss-Seidel",A,false))
      TestPassed = false;
    if(!BasicTest("Gauss-Seidel",A,true))
      TestPassed = false;  

    if(!BasicTest("Gauss-Seidel",A,false,true))
      TestPassed = false;
    if(!BasicTest("Gauss-Seidel",A,true,true))
      TestPassed = false;  

  }

  // ============================= //
  // check uses as preconditioners //
  // ============================= //
  
  if(!KrylovTest("symmetric Gauss-Seidel",A,false))
    TestPassed = false;

  if(!KrylovTest("symmetric Gauss-Seidel",A,false,true))
    TestPassed = false;


  if (!SymmetricGallery) {
    if(!KrylovTest("Gauss-Seidel",A,false))
      TestPassed = false;
    if(!KrylovTest("Gauss-Seidel",A,true))
      TestPassed = false;

    if(!KrylovTest("Gauss-Seidel",A,false,true))
      TestPassed = false;
    if(!KrylovTest("Gauss-Seidel",A,true,true))
      TestPassed = false;

  }

  // ================================== //
  // compare point and block relaxation //
  // ================================== //

  //TestPassed = TestPassed && 
   // ComparePointAndBlock("Jacobi",A,1);

  TestPassed = TestPassed && 
    ComparePointAndBlock("Jacobi",A,10);

  //TestPassed = TestPassed && 
    //ComparePointAndBlock("symmetric Gauss-Seidel",A,1);

  TestPassed = TestPassed && 
    ComparePointAndBlock("symmetric Gauss-Seidel",A,10);

  if (!SymmetricGallery) {
    //TestPassed = TestPassed && 
      //ComparePointAndBlock("Gauss-Seidel",A,1);

    TestPassed = TestPassed && 
      ComparePointAndBlock("Gauss-Seidel",A,10);
  }

  // ============================ //
  // verify effect of # of blocks //
  // ============================ //
  
  {
    int Iters4, Iters8, Iters16;
    Iters4 = CompareBlockSizes("Jacobi",A,4);
    Iters8 = CompareBlockSizes("Jacobi",A,8);
    Iters16 = CompareBlockSizes("Jacobi",A,16);

    if ((Iters16 > Iters8) && (Iters8 > Iters4)) {
      if (verbose)
        cout << "CompareBlockSizes Test passed" << endl;
    }
    else {
      if (verbose) 
        cout << "CompareBlockSizes TEST FAILED!" << endl;
      TestPassed = TestPassed && false;
    }
  }

  // ================================== //
  // verify effect of overlap in Jacobi //
  // ================================== //

  {
    int Iters0, Iters2, Iters4;
    Iters0 = CompareBlockOverlap(A,0);
    Iters2 = CompareBlockOverlap(A,2);
    Iters4 = CompareBlockOverlap(A,4);
    if ((Iters4 < Iters2) && (Iters2 < Iters0)) {
      if (verbose)
        cout << "CompareBlockOverlap Test passed" << endl;
    }
    else {
      if (verbose) 
        cout << "CompareBlockOverlap TEST FAILED!" << endl;
      TestPassed = TestPassed && false;
    }
  }

  // ================================== //
  // check if line smoothing works      //
  // ================================== //
  {
    int Iters1=
    CompareLineSmoother(A,coord);    
    printf(" comparelinesmoother iters %d \n",Iters1);
  }				
 // ================================== //
  // check if All singleton version of CompareLineSmoother    //
  // ================================== //
  {

    AllSingle(A,coord);    

  }				

  // ================================== //
  // test variable blocking             //
  // ================================== //
  {
    TestPassed = TestPassed && TestVariableBlocking(A->Comm());
  }

  // ================================== //
  // test variable blocking             //
  // ================================== //
  {
    TestPassed = TestPassed && TestTriDiVariableBlocking(A->Comm());
  }


  // ============ //
  // final output //
  // ============ //

  if (!TestPassed) {
    cout << "Test `TestRelaxation.exe' failed!" << endl;
    exit(EXIT_FAILURE);
  }
  
#ifdef HAVE_MPI
  MPI_Finalize(); 
#endif

  cout << endl;
  cout << "Test `TestRelaxation.exe' passed!" << endl;
  cout << endl;
  return(EXIT_SUCCESS);
}
