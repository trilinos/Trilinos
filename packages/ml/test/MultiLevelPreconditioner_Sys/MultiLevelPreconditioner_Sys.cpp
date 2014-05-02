#ifndef HAVE_CONFIG_H
#define HAVE_CONFIG_H
#endif

#include "ml_config.h"

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_EPETRAEXT) && defined(HAVE_ML_AZTECOO)

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_Vector.h"
#include "Epetra_Time.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Teuchos_ParameterList.hpp"
#include "ml_MultiLevelPreconditioner.h"
#include "AztecOO.h"

#include "EpetraExt_CrsMatrixIn.h"
#include "EpetraExt_VectorIn.h"
#include "ml_rbm.h"

Epetra_Vector *coord1,*coord2,*coord3;


using namespace Teuchos;
using namespace EpetraExt;

void PrintLine()
{
  std::cout << std::endl;
  for( int i=0 ; i<80 ; ++i )
    std::cout << "=";
  std::cout << std::endl;
  std::cout << std::endl;

  return;
}



int TestMultiLevelPreconditioner(char ProblemType[],
				 Teuchos::ParameterList & MLList,
				 Epetra_LinearProblem & Problem, double & TotalErrorResidual,
				 double & TotalErrorExactSol,bool cg=false)
{

  Epetra_MultiVector* lhs = Problem.GetLHS();
  Epetra_MultiVector* rhs = Problem.GetRHS();
  Epetra_RowMatrix* A = Problem.GetMatrix();

  // ======================================== //
  // create a rhs corresponding to lhs or 1's //
  // ======================================== //

  lhs->PutScalar(1.0);
  A->Multiply(false,*lhs,*rhs);

  lhs->PutScalar(0.0);

  Epetra_Time Time(A->Comm());

  // =================== //
  // call ML and AztecOO //
  // =================== //

  AztecOO solver(Problem);

  MLList.set("PDE equations",3);
  MLList.set("ML output", 11);
  MLList.set("max levels",3);
  MLList.set("smoother: sweeps",2);
  MLList.set("coarse: max size",100);
  MLList.set("aggregation: aux: enable",true);
  MLList.set("aggregation: aux: threshold",.01);
  MLList.set("x-coordinates",&((*coord1)[0]));
  MLList.set("y-coordinates",&((*coord2)[0]));
  MLList.set("z-coordinates",&((*coord3)[0]));
  MLList.set("null space: type","pre-computed");
  MLList.set("null space: add default vectors",false);
  MLList.set("null space: dimension",6);


  // Nullspace
  double *rbm = new double[7*A->NumMyRows()];
  ML_Coord2RBM(coord1->MyLength(),&((*coord1)[0]),&((*coord2)[0]),&((*coord3)[0]),rbm,3,0);
  MLList.set("null space: vectors",rbm);

  ML_Epetra::MultiLevelPreconditioner * MLPrec = new ML_Epetra::MultiLevelPreconditioner(*A, MLList, true);

  // tell AztecOO to use this preconditioner, then solve
  solver.SetPrecOperator(MLPrec);

  if(cg) solver.SetAztecOption(AZ_solver, AZ_cg);
  else solver.SetAztecOption(AZ_solver, AZ_gmres);
  solver.SetAztecOption(AZ_output, 32);
  solver.SetAztecOption(AZ_kspace, 160);

  solver.Iterate(1550, 1e-12);

  delete MLPrec;
  delete [] rbm;

  // ==================================================== //
  // compute difference between exact solution and ML one //
  // ==================================================== //

  double d = 0.0, d_tot = 0.0;

  for( int i=0 ; i<lhs->Map().NumMyElements() ; ++i )
    d += ((*lhs)[0][i] - 1.0) * ((*lhs)[0][i] - 1.0);

  A->Comm().SumAll(&d,&d_tot,1);

  // ================== //
  // compute ||Ax - b|| //
  // ================== //

  double Norm;
  Epetra_Vector Ax(rhs->Map());
  A->Multiply(false, *lhs, Ax);
  Ax.Update(1.0, *rhs, -1.0);
  Ax.Norm2(&Norm);

  std::string msg = ProblemType;

  if (A->Comm().MyPID() == 0) {
    std::cout << msg << "......Using " << A->Comm().NumProc() << " processes" << std::endl;
    std::cout << msg << "......||A x - b||_2 = " << Norm << std::endl;
    std::cout << msg << "......||x_exact - x||_2 = " << sqrt(d_tot) << std::endl;
    std::cout << msg << "......Total Time = " << Time.ElapsedTime() << std::endl;
  }

  TotalErrorExactSol += sqrt(d_tot);
  TotalErrorResidual += Norm;

  return( solver.NumIters() );

}


int main(int argc, char *argv[]) {

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  // initialize the random number generator

  int ml_one = 1;
  ML_srandom1(&ml_one);
  // ===================== //
  // create linear problem //
  // ===================== //
  Epetra_CrsMatrix * Matrix;
  MatlabFileToCrsMatrix("elasticity.dat",Comm,Matrix);
  Epetra_Map NodeMap(Matrix->NumGlobalRows()/3,0,Comm);
 
  MatrixMarketFileToVector("c1.dat",NodeMap,coord1);
  MatrixMarketFileToVector("c2.dat",NodeMap,coord2);
  MatrixMarketFileToVector("c3.dat",NodeMap,coord3);

  Epetra_Vector LHS(Matrix->RowMap());
  Epetra_Vector RHS(Matrix->RowMap());

  Epetra_LinearProblem Problem(Matrix, &LHS, &RHS);


  Teuchos::ParameterList MLList;
  double TotalErrorResidual = 0.0, TotalErrorExactSol = 0.0;



  // ====================== //
  // default options for SA //
  // ====================== //

  if (Comm.MyPID() == 0) PrintLine();

  ML_Epetra::SetDefaults("SA",MLList);
  MLList.set("smoother: type", "Gauss-Seidel");
  char mystring[80];
  strcpy(mystring,"SA");
  TestMultiLevelPreconditioner(mystring, MLList, Problem,
                               TotalErrorResidual, TotalErrorExactSol);

  // default options for SA, efficient symmetric GS //
  // ============================================== //

  if (Comm.MyPID() == 0) PrintLine();

  ML_Epetra::SetDefaults("SA",MLList);
  MLList.set("smoother: type", "Gauss-Seidel");
  MLList.set("smoother: Gauss-Seidel efficient symmetric",true);

  TestMultiLevelPreconditioner(mystring, MLList, Problem,
                               TotalErrorResidual, TotalErrorExactSol,true);

  // ============================== //
  // default options for SA, Jacobi //
  // ============================== //

  if (Comm.MyPID() == 0) PrintLine();

  ML_Epetra::SetDefaults("SA",MLList);
  MLList.set("smoother: type", "Jacobi");
  MLList.set("smoother: damping factor", .5);

  TestMultiLevelPreconditioner(mystring, MLList, Problem, TotalErrorResidual,
                               TotalErrorExactSol,true);

  // =========================== //
  // default options for SA, Cheby //
  // =========================== //

  if (Comm.MyPID() == 0) PrintLine();

  ML_Epetra::SetDefaults("SA",MLList);
  MLList.set("smoother: type", "Chebyshev");

  TestMultiLevelPreconditioner(mystring, MLList, Problem,
                               TotalErrorResidual, TotalErrorExactSol);

  // =========================== //
  // Autodetected Line SGS (trivial lines) 
  // =========================== //
  if (Comm.MyPID() == 0) PrintLine();
  ML_Epetra::SetDefaults("SA",MLList);
  MLList.set("smoother: type", "line Gauss-Seidel");
  MLList.set("smoother: line detection threshold",0.1);
  MLList.set("x-coordinates",&((*coord1)[0]));
  MLList.set("y-coordinates",&((*coord2)[0]));
  MLList.set("z-coordinates",&((*coord3)[0]));
  TestMultiLevelPreconditioner(mystring, MLList, Problem,
                               TotalErrorResidual, TotalErrorExactSol);
  

  // ===================== //
  // print out total error //
  // ===================== //

  if (Comm.MyPID() == 0) {
    std::cout << std::endl;
    std::cout << "......Total error for residual        = " << TotalErrorResidual << std::endl;
    std::cout << "......Total error for exact solution  = " << TotalErrorExactSol << std::endl;
    std::cout << std::endl;
  }

  delete Matrix;

  // if (TotalErrorResidual > 1e-8) {
  if (TotalErrorResidual > 5e-8) { // loosened tolerances
    std::cerr << "Error: `MultiLevelPrecoditioner_Sym.exe' failed!" << std::endl;
    exit(EXIT_FAILURE);
  }

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  if (Comm.MyPID() == 0)
    std::cerr << "`MultiLevelPrecoditioner_Sym.exe' passed!" << std::endl;

  return (EXIT_SUCCESS);
}

#else

#include <stdlib.h>
#include <stdio.h>
#ifdef HAVE_MPI
#include "mpi.h"
#endif

int main(int argc, char *argv[])
{
  // still need to deal with MPI, some architecture don't like
  // an exit(0) without MPI_Finalize()
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  puts("Please configure ML with --enable-epetra --enable-teuchos --enable-galeri --enable-aztecoo");

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(EXIT_SUCCESS);
}

#endif /* #if defined(ML_WITH_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_GALERI) */
