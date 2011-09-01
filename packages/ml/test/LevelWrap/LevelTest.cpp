#ifndef HAVE_CONFIG_H
#define HAVE_CONFIG_H
#endif

#include "ml_config.h"


#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_GALERI) && defined(HAVE_ML_AZTECOO) && defined(HAVE_ML_IFPACK)

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
#include "ml_LevelWrap.h"
#include "AztecOO.h"
#include "ml_epetra_utils.h"
#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"

#include "EpetraExt_MultiVectorOut.h"

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




void BuildProlongator(Epetra_CrsMatrix &A,Epetra_CrsMatrix *&P) { 
  ML_Comm* ml_comm_;
  ML_Comm_Create(&ml_comm_);
  ML_Operator* A_ML = ML_Operator_Create(ml_comm_);
  
  /* Wrap A in a ML_Operator */
  ML_Operator_WrapEpetraCrsMatrix(&A,A_ML);
  
  
  /* Pull Teuchos Options */
  string CoarsenType = "Uncoupled";
  double Threshold   = 0.0;

  /* Setup the Aggregation */
  ML_Aggregate_Struct * MLAggr;
  ML_Aggregate_Create(&MLAggr);
  ML_Aggregate_Set_MaxLevels(MLAggr, 2);
  ML_Aggregate_Set_StartLevel(MLAggr, 0);
  ML_Aggregate_Set_Threshold(MLAggr, Threshold);
  MLAggr->cur_level = 0;
  ML_Aggregate_Set_Reuse(MLAggr);
  MLAggr->keep_agg_information = 1;  
  ML_Operator *P_ML = ML_Operator_Create(ml_comm_);
  ML_Aggregate_Set_CoarsenScheme_Uncoupled(MLAggr);
  
  /* Aggregate Nodes */
  int NumAggregates = ML_Aggregate_Coarsen(MLAggr, A_ML, &P_ML, ml_comm_);
  if (NumAggregates == 0){
    cerr << "Found 0 aggregates, perhaps the problem is too small." << endl;
    exit(-2);
  }/*end if*/
  
  // Wrap P to Crs
  int nnz;
  double time;
  ML_Operator2EpetraCrsMatrix(P_ML,P,nnz,true,time,0,false);
  P->OptimizeStorage();

  ML_Operator_Destroy(&A_ML);
  ML_Operator_Destroy(&P_ML);
  ML_Comm_Destroy(&ml_comm_);
}





int TestLevelWrapPreconditioner(char ProblemType[],
				Teuchos::ParameterList & MLList,
				Epetra_LinearProblem & Problem, 
				Epetra_CrsMatrix & P0,
				double & TotalErrorResidual,
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
  
  ML_Epetra::LevelWrap *MLPrec= new ML_Epetra::LevelWrap(rcp<Epetra_CrsMatrix>(dynamic_cast<Epetra_CrsMatrix*>(A),false),rcp<Epetra_CrsMatrix>(&P0,false),MLList,true);
  
  // tell AztecOO to use this preconditioner, then solve
  solver.SetPrecOperator(MLPrec);

  solver.SetAztecOption(AZ_solver, AZ_cg);
  solver.SetAztecOption(AZ_output, 10);
  //  solver.SetAztecOption(AZ_kspace, 160);
  
  EpetraExt::MultiVectorToMatrixMarketFile("rhs.dat",(*rhs));

  solver.Iterate(100, 1e-12);
  
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
  Epetra_Vector Ax(rhs->Map());
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
  
  return( solver.NumIters() );
  
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
  GaleriList.set("n", 1200* Comm.NumProc());
  GaleriList.set("m", Comm.NumProc());

  Epetra_Map* Map = CreateMap("Linear", Comm, GaleriList);
  Epetra_CrsMatrix* Matrix = CreateCrsMatrix("Laplace1D", Map, GaleriList);

  Epetra_Vector LHS(*Map);
  Epetra_Vector RHS(*Map);
  
  Epetra_LinearProblem Problem(Matrix, &LHS, &RHS);

  Teuchos::ParameterList MLList;
  double TotalErrorResidual = 0.0, TotalErrorExactSol = 0.0;

  // ===================== //
  // Build faux prolongator - use injection //
  // ===================== //
#ifdef OLD_AND_BUSTED
  Epetra_CrsMatrix P(Copy,Matrix->RowMap(),1);
  int Pmax=Matrix->NumGlobalRows()/3;
  Epetra_Map Pdomain(Pmax,0,Comm);

  for(int i=0;i<P.NumMyRows();i++){
    double one=1.0;
    int grid=P.GRID(i);
    int grid3=P.GRID(i)/3;
    P.InsertGlobalValues(grid,1,&one,&grid3);
  }
  P.FillComplete(Pdomain,Matrix->DomainMap());
#else
  Epetra_CrsMatrix *P;
  BuildProlongator(*Matrix,P);
#endif

  // ====================== //
  // default options for LW+SA
  // ====================== //

  if (Comm.MyPID() == 0) PrintLine();

  ML_Epetra::SetDefaultsLevelWrap(MLList);
  
  char mystring[80];
  strcpy(mystring,"SA");
  TestLevelWrapPreconditioner(mystring, MLList, Problem, *P,
                               TotalErrorResidual, TotalErrorExactSol);



  // ===================== //
  // print out total error //
  // ===================== //

  if (Comm.MyPID() == 0) {
    cout << endl;
    cout << "......Total error for residual        = " << TotalErrorResidual << endl;
    cout << "......Total error for exact solution  = " << TotalErrorExactSol << endl;
    cout << endl;
  }

  delete Matrix;
  delete Map;
  
  if (TotalErrorResidual > 1e-8) {
    cerr << "Error: `MultiLevelPrecoditioner_Sym.exe' failed!" << endl;
    exit(EXIT_FAILURE);
  }

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  if (Comm.MyPID() == 0)
    cerr << "`MultiLevelPrecoditioner_Sym.exe' passed!" << endl;

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
