/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

#include "ml_config.h"
#include "ml_utils.h"

#if defined(HAVE_ML_AMESOS) && defined(HAVE_ML_TEUCHOS)

#include "ml_epetra_utils.h"
#include "ml_xyt.h"

#include "Epetra_Map.h" 
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h" 
#include "Epetra_LinearProblem.h"
#include "Epetra_Time.h"
#include "ml_amesos.h"
#include "ml_amesos_wrap.h"
#include "Amesos_BaseSolver.h"
#include "Amesos.h" 

#include "ml_amesos_wrap.h"

#ifdef EPETRA_MPI
#ifndef ML_MPI
Garbage - ML_MPI and EPETRA_MPI must be the same 
#endif
#include "Epetra_MpiComm.h"
#else
#ifdef ML_MPI
Garbage - ML_MPI and EPETRA_MPI must be the same 
#endif
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Comm.h"
#include "Teuchos_ParameterList.hpp"

static double TimeForSolve__ = 0.0;
static int Level__ = -1;
static int NumSolves__ = 0;

#ifdef DEBUG
static double MaxError__ = 0.0;
#endif

// ================================================ ====== ==== ==== == =

int ML_Amesos_Gen(ML *ml, int curr_level, int choice,
		  int MaxProcs, void **Amesos_Handle)
{

  int MaxNumNonzeros;
  double Time1, Time2;

  ML_Operator *Ke = &(ml->Amat[curr_level]);
  Epetra_CrsMatrix * Amesos_CrsMatrix;

  ML_Operator2EpetraCrsMatrix( Ke, Amesos_CrsMatrix, MaxNumNonzeros,
			       true, Time1);

  if( Amesos_CrsMatrix->Comm().MyPID() == 0 && ML_Get_PrintLevel()>2 ) {
#ifdef TFLOP
    printf("Amesos (level %d) : Time to convert to Epetra_CrsMatrix = %f (s)\n",curr_level,Time1);
#else
    cout << "Amesos (level " << curr_level
	 << ") : Time to convert to Epetra_CrsMatrix  = "
	 << Time1 << " (s)" << endl;
#endif
  }

  
  double NormInf = Amesos_CrsMatrix->NormInf();
  double NormOne = Amesos_CrsMatrix->NormOne();
  int NumGlobalRows = Amesos_CrsMatrix->NumGlobalRows();
  int NumGlobalNonzeros = Amesos_CrsMatrix->NumGlobalNonzeros();

#ifdef TFLOP
  if( Amesos_CrsMatrix->Comm().MyPID() == 0 && ML_Get_PrintLevel() > 2 ) {
    printf("Amesos (level %d) : NumGlobalRows = %d\n",curr_level,NumGlobalRows);
    printf("Amesos (level %d) : NumGlobalNonzeros = %d\n",curr_level,NumGlobalNonzeros);
    printf("Amesos (level %d) : MaxNonzerosPerRow = %d\n",curr_level,MaxNumNonzeros);
    printf("Amesos (level %d) : fill-in = %f %\n",curr_level,100.0*NumGlobalNonzeros/(NumGlobalRows*NumGlobalRows));
    printf("Amesos (level %d) : ||A||_\\infty = %f ||A||_1= %f\n",curr_level,NormInf,NormOne);
  }
#else
  if( Amesos_CrsMatrix->Comm().MyPID() == 0 && ML_Get_PrintLevel() > 2 ) {
    cout << "Amesos (level " << curr_level
	 << ") : NumGlobalRows = "
	 << NumGlobalRows << endl;
    cout << "Amesos (level " << curr_level
	 << ") : NumGlobalNonzeros = "
	 << NumGlobalNonzeros << endl;
    cout << "Amesos (level " << curr_level
	 << ") : MaxNonzerosPerRow = "
	 << MaxNumNonzeros << endl;
    cout << "Amesos (level " << curr_level
	 << ") : Fill-in = "
	 << 100.0*NumGlobalNonzeros/(NumGlobalRows*NumGlobalRows)
	 << " %" << endl;
    cout << "Amesos (level " << curr_level
	 << ") : ||A||_\\infty = "
	 << NormInf << " ||A||_1 = "
	 << NormOne << endl;
  }
  //cout << *Amesos_CrsMatrix;

  // MS // introduce support for Amesos_BaseFactory to
  // MS // allow different Amesos_Solvers
#endif
  
  Epetra_LinearProblem *Amesos_LinearProblem = new Epetra_LinearProblem;
  Amesos_LinearProblem->SetOperator( Amesos_CrsMatrix ) ; 

  Teuchos::ParameterList ParamList ;

  Teuchos::ParameterList & SluParamList=ParamList.sublist("Superludist");

  // this is specific to Superludist-2.0
  if( choice == ML_AMESOS_SUPERLUDIST ) {
    
    if( MaxProcs == -2 ) {
      if( Amesos_CrsMatrix->RowMatrixRowMap().LinearMap() == true ) {
	ParamList.set("Redistribute",false);
      } else {
#ifdef TFLOP
	if( Amesos_CrsMatrix->Comm().MyPID() == 0 ) {
 	  printf("*ML*WRN* in Amesos_Smoother, you can set MaxProcs = -1\n");
	  printf("*ML*WRN* (that is, matrix will not be redistributed)\n");
	  printf("*ML*WRN* ONLY if the matrix map is linear. Now proceeding\n");
	  printf("*ML*WRN* with redistribution of the matrix\n");
	  printf("*ML*WRN* (file %s line %d)\n",__FILE__,__LINE__);
	}
#else
	if( Amesos_CrsMatrix->Comm().MyPID() == 0 ) {
	  cout << "*ML*WRN* in Amesos_Smoother, you can set MaxProcs = -1\n"
	       << "*ML*WRN* (that is, matrix will not be redistributed)\n"
	       << "*ML*WRN* ONLY if the matrix map is linear. Now proceeding\n"
	       << "*ML*WRN* with redistribution of the matrix\n"
	       << "*ML*WRN* (file " << __FILE__ << ", line "
	       << __LINE__ << ")\n";
	}
#endif
      }
    } else {
      ParamList.set("Redistribute",true);
      SluParamList.set("MaxProcesses",MaxProcs);
    }
  }

  Amesos_BaseSolver* A_Base;
  Amesos A_Factory;

  // I got the impression that small problems are "safer"
  // in other hands than superludist ones.

  if( NumGlobalRows < 128 ) choice = ML_AMESOS_KLU;
  
  switch( choice ) {

  case ML_AMESOS_UMFPACK:
#ifdef TFLOP
    if( Amesos_CrsMatrix->Comm().MyPID() == 0 && ML_Get_PrintLevel()>2 )
      printf("Amesos (level %d) : Building UMFPACK\n",curr_level);
#else
    if( Amesos_CrsMatrix->Comm().MyPID() == 0 && ML_Get_PrintLevel()>2 )
      cout << "Amesos (level " << curr_level
	   << ") : Building UMFPACK\n";
#endif
    A_Base = A_Factory.Create("Amesos_Klu", *Amesos_LinearProblem);
    break;

  case ML_AMESOS_SUPERLUDIST:
#ifdef TFLOP
    if( Amesos_CrsMatrix->Comm().MyPID() == 0 && ML_Get_PrintLevel()>2 )
      printf("Amesos (level %d) : Building SUPERLUDIST\n",curr_level);
#else
    if( Amesos_CrsMatrix->Comm().MyPID() == 0 && ML_Get_PrintLevel()>2 )
      cout << "Amesos (level " << curr_level
	   << ") : Building SUPERLUDIST\n";
#endif
    A_Base = A_Factory.Create("Amesos_Superludist", *Amesos_LinearProblem);
    
    break;

  case ML_AMESOS_SCALAPACK:
#ifdef TFLOP
    if( Amesos_CrsMatrix->Comm().MyPID() == 0 && ML_Get_PrintLevel()>2 )
      printf("Amesos (level %d) : Building SCALAPACK\n",curr_level);
#else
    if( Amesos_CrsMatrix->Comm().MyPID() == 0 && ML_Get_PrintLevel()>2 )
      cout << "Amesos (level " << curr_level
	   << ") : Building SCALAPACK\n";
#endif
    A_Base = A_Factory.Create("Amesos_Scalapack", *Amesos_LinearProblem);
    
    break;

  case ML_AMESOS_MUMPS:
#ifdef TFLOP
    if( Amesos_CrsMatrix->Comm().MyPID() == 0 && ML_Get_PrintLevel()>2 )
      printf("Amesos (level %d) : Building MUMPS\n",curr_level);
#else
    if( Amesos_CrsMatrix->Comm().MyPID() == 0 && ML_Get_PrintLevel()>2 )
      cout << "Amesos (level " << curr_level
	   << ") : Building MUMPS\n";
#endif
    A_Base = A_Factory.Create("Amesos_Mumps", *Amesos_LinearProblem);
    break;

  case ML_AMESOS_KLU:
  default:
#ifdef TFLOP
    if( Amesos_CrsMatrix->Comm().MyPID() == 0 && ML_Get_PrintLevel()>2 )
      printf("Amesos (level %d) : Building KLU\n",curr_level);
#else
    if( Amesos_CrsMatrix->Comm().MyPID() == 0 && ML_Get_PrintLevel()>2 )
      cout << "Amesos (level " << curr_level
	   << ") : Building KLU\n";
#endif
    A_Base = A_Factory.Create("Amesos_Klu", *Amesos_LinearProblem);
    break;
  }

  // may happen the desired solver is not available. KLU is almost
  // always compiled, so try this...
  if( A_Base == 0 ) {
#ifdef TFLOP
    if( Amesos_CrsMatrix->Comm().MyPID() == 0 && ML_Get_PrintLevel()>2 )
      printf("Amesos (level %d) : Now re-building with KLU\n",curr_level);
#else
    if( Amesos_CrsMatrix->Comm().MyPID() == 0 && ML_Get_PrintLevel()>2 )
      cout << "Amesos (level " << curr_level
	   << ") : Now re-building with KLU\n";
#endif
    A_Base = A_Factory.Create("Amesos_Klu", *Amesos_LinearProblem);
    if( A_Base == 0 ) {
      if( Amesos_CrsMatrix->Comm().MyPID() == 0 ) {
	cerr << "Amesos:ERROR: not Amesos solver is available!" << endl;
      }
      exit( EXIT_FAILURE );
    }
 }

  A_Base->SetParameters(ParamList);

  Epetra_Time Time(Amesos_CrsMatrix->Comm());

  A_Base->SymbolicFactorization();
  Time1 = Time.ElapsedTime();
  Time.ResetStartTime();
  A_Base->NumericFactorization();
  Time2 = Time.ElapsedTime();

  Level__ = -1;

#ifdef TFLOP
  if( Amesos_CrsMatrix->Comm().MyPID() == 0 && ML_Get_PrintLevel()>2 ) {
    Level__ = curr_level;
    printf("Amesos (level %d) : Time for symbolic fact = %f (s)\n",curr_level,Time1);
    printf("Amesos (level %d) : Time for numerical fact = %f (s)\n",curr_level,Time2);
  }
#else
  if( Amesos_CrsMatrix->Comm().MyPID() == 0 && ML_Get_PrintLevel()>2 ) {
    Level__ = curr_level;
    cout << "Amesos (level " << curr_level
	 << ") : Time for symbolic fact  = "
	 << Time1 << " (s)" << endl;
    cout << "Amesos (level " << curr_level
	 << ") : Time for numerical fact = "
	 << Time2 << " (s)" << endl;
  }
#endif
  
  // those are very simple timing for solution
  TimeForSolve__ = 0.0;
  NumSolves__ = 0;
  
  *Amesos_Handle = (void *) A_Base ;

  return 0;
}

// ================================================ ====== ==== ==== == =

int ML_Amesos_Solve( void *Amesos_Handle, double x[], double rhs[] )
{

  Amesos_BaseSolver *A_Base = (Amesos_BaseSolver *) Amesos_Handle ;

  Epetra_Time Time(A_Base->Comm());  

  Epetra_LinearProblem *Amesos_LinearProblem = (Epetra_LinearProblem *)A_Base->GetProblem() ;
  
  Epetra_BlockMap map = Amesos_LinearProblem->GetOperator()->OperatorDomainMap() ; 

  Epetra_Vector EV_rhs( View, map, rhs ) ;
  Epetra_Vector EV_lhs( View, map, x ) ;

  Amesos_LinearProblem->SetRHS( &EV_rhs ) ; 
  Amesos_LinearProblem->SetLHS( &EV_lhs ) ;

  A_Base->Solve() ; 

  TimeForSolve__ += Time.ElapsedTime();
  NumSolves__++;

#ifdef DEBUG
  // verify that the residual is actually small (and print the max
  // in the destruction phase)
  Epetra_Vector Ax(map);
  
  (Amesos_LinearProblem->GetMatrix())->Multiply(false,EV_lhs,Ax);
  
  assert(Ax.Update(1.0, EV_rhs, -1.0)==0);
  
  double residual;
  assert(Ax.Norm2(&residual)==0);
  if( residual > MaxError__ ) MaxError__ = residual;
#endif
  
  return 0;
}

#ifdef DEBUG
#include <iomanip>
#endif

// ================================================ ====== ==== ==== == =

void ML_Amesos_Destroy(void *Amesos_Handle)
{

#ifdef TFLOP
  if( Level__ != -1 ) {
    printf("Amesos (level %d) : Time for solve = %f (s)\n",Level__,TimeForSolve__);
    if( NumSolves__ ) 
      printf("Amesos (level %d) : avg time for solve = %f (s) ( # solve = %d)\n",Level__,TimeForSolve__/NumSolves__,NumSolves__);
    else
      printf("Amesos (level %d) : no solve\n",Level__);
  }
#else
  if( Level__ != -1 ) {
    cout << endl;
    cout << "Amesos (level " << Level__
	 << ") : Time for solve = "
	 << TimeForSolve__ << " (s)" << endl;
    if( NumSolves__ ) 
      cout << "Amesos (level " << Level__
	   << ") : avg time for solve = "
	   << TimeForSolve__/NumSolves__ << " (s) ( # solves = "
	   << NumSolves__ << ")" << endl;
    else
      cout << "Amesos (level " << Level__
	   << ") : no solve" << endl;

#ifdef DEBUG
    cout << "Amesos (level " << Level__
	 << ") : max (over solves) ||Ax - b|| = " << setiosflags(ios::scientific) << MaxError__ << endl;
#endif
    cout << endl;

  }
#endif
  
  Amesos_BaseSolver *A_Base = (Amesos_BaseSolver *) Amesos_Handle ;
  const Epetra_LinearProblem *Amesos_LinearProblem;
  Amesos_LinearProblem = A_Base->GetProblem() ; 

  delete Amesos_LinearProblem->GetOperator() ; 

  delete Amesos_LinearProblem ;
  delete A_Base ;

}


#else

#include "ml_include.h"
#include "ml_amesos_wrap.h"
#include <stdio.h>

int ML_Amesos_Gen(ML *ml, int curr_level, int choice,
		  int MaxProcs, void **Amesos_Handle)
{
  puts("You must configure with --with-ml_amesos.");
  exit( EXIT_FAILURE );
  return EXIT_FAILURE;
}

int ML_Amesos_Solve( void *Amesos_Handle, double x[], double rhs[] )
{
  puts("You must configure with --with-ml_amesos.");
  exit( EXIT_FAILURE );
  return EXIT_FAILURE;
}

void ML_Amesos_Destroy(void *Amesos_Handle)
{
  puts("You must configure with --with-ml_amesos.");
  exit( EXIT_FAILURE );
}


#endif
