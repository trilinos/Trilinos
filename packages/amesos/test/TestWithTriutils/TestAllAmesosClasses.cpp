#ifndef HAVE_CONFIG_H
#define HAVE_CONFIG_H
#endif

#include "Epetra_ConfigDefs.h"

bool TestKLU         = true;
bool TestUMFPACK     = true;
bool TestSuperLU     = true;
bool TestSuperLUdist = true;
bool TestMUMPS       = true;
bool TestScaLAPACK   = true;
const bool verbose         = true;
const int OutputLevel      = 2;

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Vector.h"
#include "Epetra_Time.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"

#include "Amesos_ConfigDefs.h"
#include "Amesos.h"
#include "Teuchos_ParameterList.hpp"

#include "Trilinos_Util_CrsMatrixGallery.h"
#include "Trilinos_Util_CommandLineParser.h"


void PrintLine() 
{
  cout << endl;
  for( int i=0 ; i<80 ; ++i )
    cout << "=";
  cout << endl;
  cout << endl;
  
  return;
}

void TestAmesos(char ProblemType[],
		Teuchos::ParameterList & AmesosList,
		bool UseTranspose,
		Epetra_LinearProblem & Problem, double & TotalErrorResidual,
		double & TotalErrorExactSol)
{

  Epetra_MultiVector   * lhs     = Problem.GetLHS();
  Epetra_MultiVector   * rhs     = Problem.GetRHS();
  Epetra_RowMatrix     * A       = Problem.GetMatrix();

  // create an empty linear problem; Amesos should work with this too
  Epetra_LinearProblem Problem2;
  
  Amesos A_Factory;
  
  Amesos_BaseSolver * Solver = A_Factory.Create(ProblemType, Problem2);

  if( Solver ) {

    // Both sentences should work
    //    Solver->SetUseTranspose(UseTranspose);
    AmesosList.set("UseTranspose",UseTranspose);

    Solver->SetParameters(AmesosList);
    
    // create a rhs corresponding to lhs or 1's
    lhs->PutScalar(1.0);
    A->Multiply(UseTranspose,*lhs,*rhs);
    
    Epetra_Time Time(A->Comm());
    Epetra_Time StartTime(A->Comm());

    // set the problem here
    Problem2.SetOperator(A);
    
    Solver->SymbolicFactorization();
    double TimeForSymbolicFactorization = Time.ElapsedTime();
  
    Time.ResetStartTime();
    Solver->NumericFactorization();
    double TimeForNumericFactorization = Time.ElapsedTime();

    // set sol and rhs here
    Problem2.SetLHS(lhs);
    Problem2.SetRHS(rhs);
    
    Time.ResetStartTime();
    Solver->Solve();
    double TimeForSolve = Time.ElapsedTime();

    // compute difference between exact solution and Amesos
    double d = 0.0, d_tot = 0.0;
    
    for( int i=0 ; i<lhs->Map().NumMyElements() ; ++i )
      d += ((*lhs)[0][i] - 1.0) * ((*lhs)[0][i] - 1.0);

    A->Comm().SumAll(&d,&d_tot,1);

    // compute ||Ax - b||
    double Norm;
    Epetra_Vector Ax(rhs->Map());
    A->Multiply(UseTranspose, *lhs, Ax);
    Ax.Update(1.0, *rhs, -1.0);
    Ax.Norm2(&Norm);
    
    string msg = ProblemType;
    
    if( A->Comm().MyPID() == 0 ) {
      cout << msg << "......Using " << A->Comm().NumProc() << " processes, UseTranspose = " << UseTranspose << endl;
      cout << msg << "......||A x - b||_2 = " << Norm << endl;
      cout << msg << "......||x_exact - x||_2 = " << sqrt(d_tot) << endl;
      cout << msg << "......Time for Symbolic Factorization = " << TimeForSymbolicFactorization << endl;
      cout << msg << "......Time for Numeric Factorization  = " << TimeForNumericFactorization << endl;
      cout << msg << "......Time for Solve                  = " << TimeForSolve << endl;
      cout << msg << "......Total Time = " << StartTime.ElapsedTime() << endl;

      if( Norm > 1e-9 ) {
	cerr << endl << msg << " WARNING : TEST FAILED!" << endl << endl;
      }
    }

    TotalErrorExactSol += sqrt(d_tot);
    TotalErrorResidual += Norm;
    
    delete Solver;
  
  }

  return;
  
}

using namespace Trilinos_Util;

int main(int argc, char *argv[]) {

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  // ===================== //
  // create linear problem //
  // ===================== //
  
  CommandLineParser CLP(argc,argv);
  CrsMatrixGallery Gallery("", Comm);
  
  // default values for problem type and size
  if( CLP.Has("-problem_type") == false ) CLP.Add("-problem_type", "recirc_2d" ); 
  if( CLP.Has("-problem_size") == false ) CLP.Add("-problem_size", "900" ); 

  // disable tests from CLP
  if( CLP.Has("-disable_klu")         == true ) TestKLU       = false;
  if( CLP.Has("-disable_superlu")     == true ) TestUMFPACK   = false;
  if( CLP.Has("-disable_superludist") == true ) TestSuperLU   = false;
  if( CLP.Has("-disable_mumps" )      == true ) TestMUMPS     = false;
  if( CLP.Has("-disable_scalapack" )  == true ) TestScaLAPACK = false;

  Gallery.Set(CLP);
  Epetra_LinearProblem * Problem = Gallery.GetLinearProblem();

  double TotalErrorResidual = 0.0, TotalErrorExactSol = 0.0;
    
  // ====================== //
  // KLU -- default options //
  // ====================== //

  if( Comm.MyPID() == 0 ) PrintLine();

  if( TestKLU ) {
    
    Teuchos::ParameterList AmesosList;
    TestAmesos("Amesos_Klu", AmesosList, false, *Problem, TotalErrorResidual, TotalErrorExactSol );
    
  }

  // ======================================= //
  // KLU -- default options -- use transpose //
  // ======================================= //

  if( Comm.MyPID() == 0 ) PrintLine();
  
  if( TestKLU ) {
    
    Teuchos::ParameterList AmesosList;
    TestAmesos("Amesos_Klu", AmesosList, true, *Problem, TotalErrorResidual, TotalErrorExactSol );
    
  }
  
  // ================ //
  // KLU -- options 1 //
  // ================ //

  if( Comm.MyPID() == 0 ) PrintLine();

  if( TestKLU ) {
    
    Teuchos::ParameterList AmesosList;
    AmesosList.set("ComputeVectorNorms",true);
    AmesosList.set("ComputeTrueResidual",true);
    AmesosList.set("PrintStatus",true);
    AmesosList.set("PrintTiming",true);
    AmesosList.set("OutputLevel",1);
    TestAmesos("Amesos_Klu", AmesosList, false, *Problem, TotalErrorResidual, TotalErrorExactSol );
    
  }

  // ========================== //
  // UMFPACK -- default options //
  // ========================== //

  if( Comm.MyPID() == 0 ) PrintLine();

  if( TestUMFPACK ) {
    
     Teuchos::ParameterList AmesosList;
     TestAmesos("Amesos_Umfpack", AmesosList, false, *Problem, TotalErrorResidual, TotalErrorExactSol );
     
  }
  
  // ==========================================  //
  // UMFPACK -- default options -- use transpose //
  // =========================================== //

  if( Comm.MyPID() == 0 ) PrintLine();

  if( TestUMFPACK ) {
    
     Teuchos::ParameterList AmesosList;
     TestAmesos("Amesos_Umfpack", AmesosList, true, *Problem, TotalErrorResidual, TotalErrorExactSol );
  }
  
  // ==================== //
  // UMFPACK -- options 1 //
  // ===================== //

  if( Comm.MyPID() == 0 ) PrintLine();

  if( TestUMFPACK ) {
    
     Teuchos::ParameterList AmesosList;
     AmesosList.set("ComputeVectorNorms",true);
     AmesosList.set("ComputeTrueResidual",true);
     AmesosList.set("PrintStatus",true);
     AmesosList.set("PrintTiming",true);
     AmesosList.set("OutputLevel",1);
     TestAmesos("Amesos_Umfpack", AmesosList, false, *Problem, TotalErrorResidual, TotalErrorExactSol );
     
  }
  
  // ======= //
  // Superlu //
  // ======= //

  if( Comm.MyPID() == 0 ) PrintLine();

  if( TestSuperLU ) {
  
     Teuchos::ParameterList AmesosList;
     TestAmesos("Amesos_Superlu", AmesosList, false, *Problem, TotalErrorResidual, TotalErrorExactSol );
  }

  // ============================== //
  // Superludist -- default options //
  // ============================== //

  if( Comm.MyPID() == 0 ) PrintLine();

  if( TestSuperLUdist ) { 

    Teuchos::ParameterList AmesosList;
    TestAmesos("Amesos_Superludist", AmesosList, false, *Problem, TotalErrorResidual, TotalErrorExactSol );
    
  }

  // ============================================== //
  // Superludist -- default options -- use tranpose //
  // ============================================== //
  /* this is not supported now
  if( Comm.MyPID() == 0 ) PrintLine();

  if( TestSuperLUdist ) {
  
     Teuchos::ParameterList AmesosList;
     TestAmesos("Amesos_Superludist", AmesosList, true, *Problem, TotalErrorResidual, TotalErrorExactSol );
  }
  */

  // ======================== //
  // Superludist -- options 1 //
  // ======================== //

  if( Comm.MyPID() == 0 ) PrintLine();

  if( TestSuperLUdist ) { 

    Teuchos::ParameterList AmesosList;
    AmesosList.set("PrintTiming",true);
    AmesosList.set("PrintStatus",true);
    AmesosList.set("MaxProcs",7);
    AmesosList.set("Redistribute",false);
    TestAmesos("Amesos_Superludist", AmesosList, false, *Problem, TotalErrorResidual, TotalErrorExactSol );
    
  }

  // ======================== //
  // Superludist -- options 2 //
  // ======================== //

  if( Comm.MyPID() == 0 ) PrintLine();

  if( TestSuperLUdist ) { 

    Teuchos::ParameterList AmesosList;
    AmesosList.set("PrintTiming",true);
    AmesosList.set("PrintStatus",true);
    AmesosList.set("MaxProcs",-2);
    AmesosList.set("Redistribute",true);
    AmesosList.set("ComputeTrueResidual",true);    
    AmesosList.set("ComputeVectorNorms",true);
    TestAmesos("Amesos_Superludist", AmesosList, false, *Problem, TotalErrorResidual, TotalErrorExactSol );
    
  }

  // ======================== //
  // MUMPS -- default options //
  // ======================== //

  if( Comm.MyPID() == 0 ) PrintLine();

  if( TestMUMPS ) {
    
     Teuchos::ParameterList AmesosList;
     TestAmesos("Amesos_Mumps", AmesosList, false, *Problem, TotalErrorResidual, TotalErrorExactSol );
     
  }

  // ======================================== //
  // MUMPS -- default options -- use tranpose //
  // ======================================== //

  if( Comm.MyPID() == 0 ) PrintLine();

  if( TestMUMPS ) {
     Teuchos::ParameterList AmesosList;
     AmesosList.set("OutputLevel",OutputLevel);
     TestAmesos("Amesos_Mumps", AmesosList, true, *Problem, TotalErrorResidual, TotalErrorExactSol );
  }
  
  // =================== //
  // MUMPS --  options 1 //
  // =================== //

  if( Comm.MyPID() == 0 ) PrintLine();
  
  if( TestMUMPS ) {
  
     Teuchos::ParameterList AmesosList;
     AmesosList.set("Threshold",.0);
     AmesosList.set("MaxProcs",-1);
     AmesosList.set("MaxProcsMatrix",-1);
     AmesosList.set("OutputLevel",OutputLevel);
     
     TestAmesos("Amesos_Mumps", AmesosList, false, *Problem, TotalErrorResidual, TotalErrorExactSol );
  }

  // ================== //
  // MUMPS -- options 2 //
  // ================== //

  if( Comm.MyPID() == 0 ) PrintLine();

  if( TestMUMPS ) {
    
     Teuchos::ParameterList AmesosList;
     
     AmesosList.set("OutputLevel",0);
     AmesosList.set("MaxProcs",-2);
     AmesosList.set("MaxProcsMatrix",-2);
     AmesosList.set("OutputLevel",OutputLevel);
          
     TestAmesos("Amesos_Mumps", AmesosList, false, *Problem, TotalErrorResidual, TotalErrorExactSol );
  }
  
  // ================== //
  // MUMPS -- options 2 //
  // ================== //

  if( Comm.MyPID() == 0 ) PrintLine();

  if( TestMUMPS ) {
    
     Teuchos::ParameterList AmesosList;
     
     AmesosList.set("OutputLevel",0);
     AmesosList.set("MaxProcs",-3);
     AmesosList.set("MaxProcsMatrix",-3);
     AmesosList.set("OutputLevel",OutputLevel);
          
     TestAmesos("Amesos_Mumps", AmesosList, false, *Problem, TotalErrorResidual, TotalErrorExactSol );
  }

  // ================== //
  // MUMPS -- options 3 //
  // ================== //

  if( Comm.MyPID() == 0 ) PrintLine();

  if( TestMUMPS ) {
    
     Teuchos::ParameterList AmesosList;
     
     AmesosList.set("OutputLevel",0);
     AmesosList.set("MaxProcs",4);
     AmesosList.set("MaxProcsMatrix",2);
     AmesosList.set("OutputLevel",OutputLevel);
          
     TestAmesos("Amesos_Mumps", AmesosList, false, *Problem, TotalErrorResidual, TotalErrorExactSol );
  }

  // ========= //
  // Scalapack //
  // ========= //
  
  if( Comm.MyPID() == 0 ) PrintLine();

  if( TestScaLAPACK ) {
    
     Teuchos::ParameterList AmesosList;
     TestAmesos("Amesos_Scalapack", AmesosList, false, *Problem, TotalErrorResidual, TotalErrorExactSol );
  }
    
  // ====================== //
  // Scalapack -- transpose //
  // ====================== //
  
  if( Comm.MyPID() == 0 ) PrintLine();

  if( TestScaLAPACK ) {
    
     Teuchos::ParameterList AmesosList;
     AmesosList.set("ComputeVectorNorms",true);
     AmesosList.set("ComputeTrueResidual",true);
     AmesosList.set("PrintStatus",true);
     AmesosList.set("PrintTiming",true);
     AmesosList.set("OutputLevel",1);
     TestAmesos("Amesos_Scalapack", AmesosList, true, *Problem, TotalErrorResidual, TotalErrorExactSol );
  }

  // ====================== //
  // Scalapack -- options 1 //
  // ====================== //
  
  if( Comm.MyPID() == 0 ) PrintLine();

  if( TestScaLAPACK ) {
    
     Teuchos::ParameterList AmesosList;
     TestAmesos("Amesos_Scalapack", AmesosList, false, *Problem, TotalErrorResidual, TotalErrorExactSol );
  }

  // print out total error
  
  if( Comm.MyPID() == 0 ) {
    cout << endl;
    cout << "......Total error for residual = " << TotalErrorResidual << endl;
    cout << "......Total error for exact solution  = " << TotalErrorExactSol << endl;
    cout << endl;
 }
  
#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  if( TotalErrorResidual < 1e-9 ) return( EXIT_SUCCESS );
  else                            return( EXIT_FAILURE );

}
