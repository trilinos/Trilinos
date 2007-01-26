#ifndef HAVE_CONFIG_H
#define HAVE_CONFIG_H
#endif

#include "ml_config.h"

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_GALERI) && defined(HAVE_ML_AZTECOO)

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Time.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Teuchos_ParameterList.hpp"
#include "ml_MultiLevelPreconditioner.h"
#include "AztecOO.h"
#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"

using namespace Teuchos;
using namespace Galeri;

void PrintLine() 
{
  cout << endl;
  for( int i=0 ; i<80 ; ++i )
    cout << "=";
  cout << endl;
  cout << endl;
  
  return;
}

int TestMultiLevelPreconditioner(char ProblemType[],
				 Teuchos::ParameterList & MLList,
				 Epetra_LinearProblem & Problem, double & TotalErrorResidual,
				 double & TotalErrorExactSol)
{
  
  Epetra_MultiVector   * lhs     = Problem.GetLHS();
  Epetra_MultiVector   * rhs     = Problem.GetRHS();
  Epetra_RowMatrix     * A       = Problem.GetMatrix();
  
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
  
  ML_Epetra::MultiLevelPreconditioner * MLPrec = new ML_Epetra::MultiLevelPreconditioner(*A, MLList, true);
  
  // tell AztecOO to use this preconditioner, then solve
  solver.SetPrecOperator(MLPrec);
  
  solver.SetAztecOption(AZ_solver, AZ_gmres);
  solver.SetAztecOption(AZ_output, 32);
  solver.SetAztecOption(AZ_kspace, 160);
  
  solver.Iterate(1550, 1e-12);
  
  delete MLPrec;
  
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
  Epetra_MultiVector Ax(rhs->Map(), rhs->NumVectors());
  A->Multiply(false, *lhs, Ax);
  Ax.Update(1.0, *rhs, -1.0);
  Ax.Norm2(&Norm);
  
  string msg = ProblemType;
  
  if (A->Comm().MyPID() == 0) {
    cout << msg << "......Using " << A->Comm().NumProc() << " processes" << endl;
    cout << msg << "......||A x - b||_2 = " << Norm << endl;
    cout << msg << "......||x_exact - x||_2 = " << sqrt(d_tot) << endl;
    cout << msg << "......Total Time = " << Time.ElapsedTime() << endl;
  }
  
  TotalErrorExactSol += sqrt(d_tot);
  TotalErrorResidual += Norm;
  
  return (solver.NumIters());
  
}

using namespace Galeri;

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

  ParameterList GaleriList;
  GaleriList.set("nx", 32);
  GaleriList.set("ny", 32 * Comm.NumProc());
  GaleriList.set("mx", 1);
  GaleriList.set("my", Comm.NumProc());

  Epetra_Map* Map = CreateMap("Cartesian2D", Comm, GaleriList);
  Epetra_CrsMatrix* Matrix = CreateCrsMatrix("Recirc2D", Map, GaleriList);

  Epetra_MultiVector LHS(*Map, 1);
  Epetra_MultiVector RHS(*Map, 1);

  Epetra_LinearProblem Problem(Matrix, &LHS, &RHS);

  Teuchos::ParameterList MLList;

  double TotalErrorResidual = 0.0, TotalErrorExactSol = 0.0;


  // JJH Note: default internal ML options are no longer appropriate for
  // JJH       nonsymmetric problems.  Thus, the user should always populate
  // JJH       the parameter list by calling SetDefaults().

  // ======================== //
  // default options  for NSSA//
  // ======================== //

  if (Comm.MyPID() == 0) PrintLine();

  char mystring[80];
  ML_Epetra::SetDefaults("NSSA",MLList);
  strcpy(mystring,"NSSA");
  TestMultiLevelPreconditioner(mystring, MLList, Problem, 
                               TotalErrorResidual, TotalErrorExactSol );

  // ====================== //
  // default options for DD //
  // ====================== //

  if (Comm.MyPID() == 0) PrintLine();

  ML_Epetra::SetDefaults("DD",MLList);
  strcpy(mystring,"DD");
  TestMultiLevelPreconditioner(mystring, MLList, Problem, 
                               TotalErrorResidual, TotalErrorExactSol );

  // ========================================== //
  // default options for DD -- 16 aggr per proc //
  // ========================================== //

  if (Comm.MyPID() == 0) PrintLine();

  ML_Epetra::SetDefaults("DD",MLList);
  MLList.set("aggregation: local aggregates", 16);
  TestMultiLevelPreconditioner(mystring, MLList, Problem, 
                               TotalErrorResidual, TotalErrorExactSol );

  // ========================= //
  // default options for DD-ML //
  // ========================= //

  if (Comm.MyPID() == 0) PrintLine();

  ML_Epetra::SetDefaults("DD-ML",MLList);
  strcpy(mystring,"DD-ML");
  TestMultiLevelPreconditioner(mystring, MLList, Problem, 
                               TotalErrorResidual, TotalErrorExactSol );

  // ========================= //
  // default options for DD-ML //
  // ========================= //

  if (Comm.MyPID() == 0) PrintLine();

  ML_Epetra::SetDefaults("DD-ML",MLList);
  MLList.set("aggregation: nodes per aggregate (level 0)", 64);
  MLList.set("aggregation: nodes per aggregate (level 1)", 27);
  TestMultiLevelPreconditioner(mystring, MLList, Problem, 
                               TotalErrorResidual, TotalErrorExactSol );

  // ===================== //
  // print out total error //
  // ===================== //

  if (Comm.MyPID() == 0) 
  {
    cout << endl;
    cout << "......Total error for residual        = " << TotalErrorResidual << endl;
    cout << "......Total error for exact solution  = " << TotalErrorExactSol << endl;
    cout << endl;
  }

  if (TotalErrorResidual > 1e-8) {
    if (Comm.MyPID() == 0)
      cout << "ERROR: test `MultiLevelPreconditioner_NonSym.exe failed!" <<endl;
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return( EXIT_FAILURE );
  }

  delete Matrix;
  delete Map;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  if (Comm.MyPID() == 0)
    cout << endl << "Test `MultiLevelPreconditioner_NonSym.exe passed!" << endl;

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

  return(0);
}

#endif /* #if defined(ML_WITH_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_GALERI) HAVE_ML_AZTECOO */
