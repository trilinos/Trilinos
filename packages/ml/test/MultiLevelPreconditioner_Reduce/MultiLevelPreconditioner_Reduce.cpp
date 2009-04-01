#ifndef HAVE_CONFIG_H
#define HAVE_CONFIG_H
#endif

#include "ml_config.h"

#if defined(HAVE_ML_EPETRAEXT) && defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_GALERI) && defined(HAVE_ML_AZTECOO) && defined(HAVE_ML_ZOLTAN) && defined(HAVE_ML_ISORROPIA)

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

#include "Teuchos_RCP.hpp"
#include "Isorropia_Epetra.hpp"
#include "Isorropia_EpetraRedistributor.hpp"

#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_MultiVectorOut.h"

#include "EpetraExt_RestrictedCrsMatrixWrapper.h"
#include "EpetraExt_RestrictedMultiVectorWrapper.h"

using namespace Teuchos;
using namespace Galeri;
using namespace EpetraExt;

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
  
  Epetra_MultiVector* lhs = Problem.GetLHS();
  Epetra_MultiVector* rhs = Problem.GetRHS();
  Epetra_CrsMatrix* A = dynamic_cast<Epetra_CrsMatrix*>(Problem.GetMatrix());
  int PID = A->Comm().MyPID();
  int numProcs = A->Comm().NumProc();
  RCP<const Epetra_RowMatrix> Arcp = Teuchos::rcp(A, false);
  double n1, n2,nf;
  
  // ======================================== //
  // create a rhs corresponding to lhs or 1's //
  // ======================================== //
  
  lhs->PutScalar(1.0);
  A->Multiply(false,*lhs,*rhs);

  lhs->PutScalar(0.0);
  MLList.set("ML output", 0);

  RowMatrixToMatlabFile("mat_f.dat",*A);  
  MultiVectorToMatrixMarketFile("lhs_f.dat",*lhs,0,0,false);
  MultiVectorToMatrixMarketFile("rhs_f.dat",*rhs,0,0,false);

  
  Epetra_Time Time(A->Comm());
  /* Build the Zoltan list - Group #1 */
  ParameterList Zlist1,Sublist1;
  Sublist1.set("DEBUG_LEVEL","0");
  Sublist1.set("NUM_GLOBAL_PARTITIONS","2");
  Zlist1.set("Zoltan",Sublist1);
  
  /* Start Isorropia's Ninja Magic - Group #1 */
  Teuchos::RCP<Isorropia::Epetra::Partitioner> partitioner1 =
    Isorropia::Epetra::create_partitioner(Arcp, Zlist1);
  Isorropia::Epetra::Redistributor rd1(partitioner1);

  Teuchos::RCP<Epetra_CrsMatrix> ResA1=rd1.redistribute(*A);
  Teuchos::RCP<Epetra_MultiVector> ResX1=rd1.redistribute(*lhs);
  Teuchos::RCP<Epetra_MultiVector> ResB1=rd1.redistribute(*rhs);

  RestrictedCrsMatrixWrapper RW1;
  RW1.restrict_comm(ResA1);
  RestrictedMultiVectorWrapper RX1,RB1;
  RX1.restrict_comm(ResX1);
  RB1.restrict_comm(ResB1);

  /* Build the Zoltan list - Group #2 */
  ParameterList Zlist2,Sublist2;
  Sublist2.set("DEBUG_LEVEL","0");
  if(PID > 1) Sublist2.set("NUM_LOCAL_PARTITIONS","1");
  else Sublist2.set("NUM_LOCAL_PARTITIONS","0");
  Zlist2.set("Zoltan",Sublist2);
    
  /* Start Isorropia's Ninja Magic - Group #2 */
  Teuchos::RCP<Isorropia::Epetra::Partitioner> partitioner2 =
    Isorropia::Epetra::create_partitioner(Arcp, Zlist2);
  Isorropia::Epetra::Redistributor rd2(partitioner2);

  Teuchos::RCP<Epetra_CrsMatrix> ResA2=rd2.redistribute(*A);
  Teuchos::RCP<Epetra_MultiVector> ResX2=rd2.redistribute(*lhs);
  Teuchos::RCP<Epetra_MultiVector> ResB2=rd2.redistribute(*rhs);

  RestrictedCrsMatrixWrapper RW2;
  RW2.restrict_comm(ResA2);
  RestrictedMultiVectorWrapper RX2,RB2;
  RX2.restrict_comm(ResX2);
  RB2.restrict_comm(ResB2);

  if(RW1.RestrictedProcIsActive()){
    Teuchos::RCP<Epetra_CrsMatrix> SubA1 = RW1.RestrictedMatrix();
    Teuchos::RCP<Epetra_MultiVector> SubX1 = RX1.RestrictedMultiVector();
    Teuchos::RCP<Epetra_MultiVector> SubB1 = RB1.RestrictedMultiVector();    
    ML_Epetra::MultiLevelPreconditioner * SubPrec1 = new ML_Epetra::MultiLevelPreconditioner(*SubA1, MLList, true);        

    Epetra_LinearProblem Problem1(&*SubA1,&*SubX1,&*SubB1);
    AztecOO solver1(Problem1);
    solver1.SetPrecOperator(SubPrec1);  
    solver1.SetAztecOption(AZ_solver, AZ_gmres);
    solver1.SetAztecOption(AZ_output, 32);
    solver1.SetAztecOption(AZ_kspace, 160);  
    solver1.Iterate(1550, 1e-12);
    delete SubPrec1;

  }
  else{
    Teuchos::RCP<Epetra_CrsMatrix> SubA2 = RW2.RestrictedMatrix();
    Teuchos::RCP<Epetra_MultiVector> SubX2 = RX2.RestrictedMultiVector();
    Teuchos::RCP<Epetra_MultiVector> SubB2 = RB2.RestrictedMultiVector();        
    ML_Epetra::MultiLevelPreconditioner * SubPrec2 = new ML_Epetra::MultiLevelPreconditioner(*SubA2, MLList, true);        
    
    Epetra_LinearProblem Problem2(&*SubA2,&*SubX2,&*SubB2);
    AztecOO solver2(Problem2);
    solver2.SetPrecOperator(SubPrec2);  
    solver2.SetAztecOption(AZ_solver, AZ_gmres);
    solver2.SetAztecOption(AZ_output, 32);
    solver2.SetAztecOption(AZ_kspace, 160);  
    solver2.Iterate(1550, 1e-12);
    delete SubPrec2;

  }

  /* Post-processing exports */
  Epetra_MultiVector ans1(*lhs), ans2(*lhs); 
  rd1.redistribute_reverse(*ResX1,ans1);
  rd2.redistribute_reverse(*ResX2,ans2);
  
  /* Run on Full Problem */
  A->Comm().Barrier();    
  ML_Epetra::MultiLevelPreconditioner * FullPrec = new ML_Epetra::MultiLevelPreconditioner(*A, MLList, true);          
  AztecOO solverF(Problem);
  solverF.SetPrecOperator(FullPrec);  
  solverF.SetAztecOption(AZ_solver, AZ_gmres);
  solverF.SetAztecOption(AZ_output, 32);
  solverF.SetAztecOption(AZ_kspace, 160);  
  solverF.Iterate(1550, 1e-12);
  delete FullPrec;


  /* Solution Comparison */
  ans1.Update(1.0,*lhs,-1.0);
  ans2.Update(1.0,*lhs,-1.0);
  ans1.Norm2(&n1);
  ans2.Norm2(&n2);
  if(!PID) {
    printf("Norm Diff 1 = %6.4e\n",n1);
    printf("Norm Diff 2 = %6.4e\n",n2);
  }

  TotalErrorExactSol += n1 + n2;
    
  
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
  int base=10;
  GaleriList.set("nx", base);
  GaleriList.set("ny", base);
  GaleriList.set("nz", base * Comm.NumProc());
  GaleriList.set("mx", 1);
  GaleriList.set("my", 1);
  GaleriList.set("mz", Comm.NumProc());

  Epetra_Map* Map = CreateMap("Cartesian3D", Comm, GaleriList);
  Epetra_CrsMatrix* Matrix = CreateCrsMatrix("Laplace3D", Map, GaleriList);

  Epetra_Vector LHS(*Map);
  Epetra_Vector RHS(*Map);
  
  Epetra_LinearProblem Problem(Matrix, &LHS, &RHS);

  Teuchos::ParameterList MLList;
  double TotalErrorResidual = 0.0, TotalErrorExactSol = 0.0;

  // ==================n==== //
  // default options for SA //
  // ====================== //

  if (Comm.MyPID() == 0) PrintLine();

  ML_Epetra::SetDefaults("SA",MLList);
  MLList.set("smoother: type", "Gauss-Seidel");
  char mystring[80];
  strcpy(mystring,"SA");
  TestMultiLevelPreconditioner(mystring, MLList, Problem, 
                               TotalErrorResidual, TotalErrorExactSol);

  // ============================== //
  // default options for SA, Jacobi //
  // ============================== //

  if (Comm.MyPID() == 0) PrintLine();

  ML_Epetra::SetDefaults("SA",MLList);
  MLList.set("smoother: type", "Jacobi");

  TestMultiLevelPreconditioner(mystring, MLList, Problem, TotalErrorResidual,
                               TotalErrorExactSol);

  // =========================== //
  // default options for SA, Cheby //
  // =========================== //

  if (Comm.MyPID() == 0) PrintLine();

  ML_Epetra::SetDefaults("SA",MLList);
  MLList.set("smoother: type", "Chebyshev");

  TestMultiLevelPreconditioner(mystring, MLList, Problem, 
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
