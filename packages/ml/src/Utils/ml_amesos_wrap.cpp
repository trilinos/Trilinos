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
#include "Epetra_LinearProblem.h"
#include "Epetra_Time.h"
#include "ml_amesos.h"
#include "ml_amesos_wrap.h"
#include "ml_RowMatrix.h"
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

#ifdef ML_AMESOS_DEBUG
static double MaxError__ = 0.0;
#endif

// ================================================ ====== ==== ==== == =

int ML_Amesos_Gen(ML *ml, int curr_level, int choice,
		  int MaxProcs, void **Amesos_Handle)
{

#ifdef HAVE_MPI
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  ML_Operator *Ke = &(ml->Amat[curr_level]);
  ML_Epetra::RowMatrix* Amesos_Matrix = new ML_Epetra::RowMatrix(Ke,Comm);
  assert (Amesos_Matrix != 0);
  
  int NumGlobalRows = Amesos_Matrix->NumGlobalRows();
  int NumGlobalNonzeros = Amesos_Matrix->NumGlobalNonzeros();

  // sanity check, coarse matrix should not be empty
  if( NumGlobalRows == 0 && Amesos_Matrix->Comm().MyPID() == 0 ) {
    cerr << endl;
    cerr << "ERROR : Coarse matrix has no rows!" << endl;
    cerr << endl;
  }
  if( NumGlobalNonzeros == 0 && Amesos_Matrix->Comm().MyPID() == 0 ) {
    cerr << endl;
    cerr << "ERROR : Coarse matrix has no nonzero elements!" << endl;
    cerr << endl;
  }

#ifdef TFLOP
  if( Amesos_Matrix->Comm().MyPID() == 0 && ML_Get_PrintLevel() > 2 ) {
    printf("Amesos (level %d) : NumGlobalRows = %d\n",curr_level,NumGlobalRows);
    printf("Amesos (level %d) : NumGlobalNonzeros = %d\n",curr_level,NumGlobalNonzeros);
    printf("Amesos (level %d) : fill-in = %f %\n",curr_level,100.0*NumGlobalNonzeros/(NumGlobalRows*NumGlobalRows));
  }
#else
  if( Amesos_Matrix->Comm().MyPID() == 0 && ML_Get_PrintLevel() > 2 ) {
    cout << "Amesos (level " << curr_level
	 << ") : NumGlobalRows = "
	 << NumGlobalRows << endl;
    cout << "Amesos (level " << curr_level
	 << ") : NumGlobalNonzeros = "
	 << NumGlobalNonzeros << endl;
    cout << "Amesos (level " << curr_level
	 << ") : Fill-in = "
	 << 100.0*NumGlobalNonzeros/(NumGlobalRows*NumGlobalRows)
	 << " %" << endl;
  }

#endif
  
  Epetra_LinearProblem *Amesos_LinearProblem = new Epetra_LinearProblem;
  Amesos_LinearProblem->SetOperator(Amesos_Matrix); 

  Teuchos::ParameterList AmesosList;

  if( ML_Get_PrintLevel() > 8 ) {
    AmesosList.set("PrintTiming",true);
    AmesosList.set("PrintStatus",true);
  }
  AmesosList.set("MaxProcs",MaxProcs);

  // don't use iterative refinement for Superludist only
  Teuchos::ParameterList & SuperludistList = AmesosList.sublist("Superludist");
  SuperludistList.set("IterRefine","NO");

  Amesos_BaseSolver* A_Base;
  Amesos A_Factory;

  switch( choice ) {

  case ML_AMESOS_UMFPACK:
#ifdef TFLOP
    if( Amesos_Matrix->Comm().MyPID() == 0 && ML_Get_PrintLevel()>2 )
      printf("Amesos (level %d) : Building UMFPACK\n",curr_level);
#else
    if( Amesos_Matrix->Comm().MyPID() == 0 && ML_Get_PrintLevel()>2 )
      cout << "Amesos (level " << curr_level
	   << ") : Building UMFPACK\n";
#endif
    A_Base = A_Factory.Create("Amesos_Klu", *Amesos_LinearProblem);
    break;

  case ML_AMESOS_SUPERLUDIST:
#ifdef TFLOP
    if( Amesos_Matrix->Comm().MyPID() == 0 && ML_Get_PrintLevel()>2 )
      printf("Amesos (level %d) : Building SUPERLUDIST\n",curr_level);
#else
    if( Amesos_Matrix->Comm().MyPID() == 0 && ML_Get_PrintLevel()>2 )
      cout << "Amesos (level " << curr_level
	   << ") : Building SUPERLUDIST\n";
#endif
    A_Base = A_Factory.Create("Amesos_Superludist", *Amesos_LinearProblem);
    
    break;

  case ML_AMESOS_SCALAPACK:
#ifdef TFLOP
    if( Amesos_Matrix->Comm().MyPID() == 0 && ML_Get_PrintLevel()>2 )
      printf("Amesos (level %d) : Building SCALAPACK\n",curr_level);
#else
    if( Amesos_Matrix->Comm().MyPID() == 0 && ML_Get_PrintLevel()>2 )
      cout << "Amesos (level " << curr_level
	   << ") : Building SCALAPACK\n";
#endif
    A_Base = A_Factory.Create("Amesos_Scalapack", *Amesos_LinearProblem);
    
    break;

  case ML_AMESOS_MUMPS:
#ifdef TFLOP
    if( Amesos_Matrix->Comm().MyPID() == 0 && ML_Get_PrintLevel()>2 )
      printf("Amesos (level %d) : Building MUMPS\n",curr_level);
#else
    if( Amesos_Matrix->Comm().MyPID() == 0 && ML_Get_PrintLevel()>2 )
      cout << "Amesos (level " << curr_level
	   << ") : Building MUMPS\n";
#endif
    A_Base = A_Factory.Create("Amesos_Mumps", *Amesos_LinearProblem);
    break;

  case ML_AMESOS_KLU:
  default:
#ifdef TFLOP
    if( Amesos_Matrix->Comm().MyPID() == 0 && ML_Get_PrintLevel()>2 )
      printf("Amesos (level %d) : Building KLU\n",curr_level);
#else
    if( Amesos_Matrix->Comm().MyPID() == 0 && ML_Get_PrintLevel()>2 )
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
    if( Amesos_Matrix->Comm().MyPID() == 0 && ML_Get_PrintLevel()>2 )
      printf("Amesos (level %d) : Now re-building with KLU\n",curr_level);
#else
    if( Amesos_Matrix->Comm().MyPID() == 0 && ML_Get_PrintLevel()>2 ) {
      cout << "Amesos (level " << curr_level
	   << ") : This coarse solver is not available." << endl;
      cout << "Amesos (level " << curr_level
	   << ") : Now re-building with KLU\n";
    }
#endif
    A_Base = A_Factory.Create("Amesos_Klu", *Amesos_LinearProblem);
    if( A_Base == 0 ) {
      if( Amesos_Matrix->Comm().MyPID() == 0 ) {
	cerr << "*ML*ERR* no Amesos solver is available!" << endl;
      }
      exit( EXIT_FAILURE );
    }
 }

  A_Base->SetParameters(AmesosList);

  Epetra_Time Time(Amesos_Matrix->Comm());

  A_Base->SymbolicFactorization();
  double Time1 = Time.ElapsedTime();
  Time.ResetStartTime();
  A_Base->NumericFactorization();
  double Time2 = Time.ElapsedTime();

  Level__ = -1;

#ifdef TFLOP
  if( Amesos_Matrix->Comm().MyPID() == 0 && ML_Get_PrintLevel()>2 ) {
    Level__ = curr_level;
    printf("Amesos (level %d) : Time for symbolic fact = %f (s)\n",curr_level,Time1);
    printf("Amesos (level %d) : Time for numerical fact = %f (s)\n",curr_level,Time2);
  }
#else
  if( Amesos_Matrix->Comm().MyPID() == 0 && ML_Get_PrintLevel()>2 ) {
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
  
  const Epetra_BlockMap & map = Amesos_LinearProblem->GetOperator()->OperatorDomainMap() ; 

  Epetra_Vector EV_rhs( View, map, rhs ) ;
  Epetra_Vector EV_lhs( View, map, x ) ;

  Amesos_LinearProblem->SetRHS( &EV_rhs ) ; 
  Amesos_LinearProblem->SetLHS( &EV_lhs ) ;

  A_Base->Solve() ; 

  TimeForSolve__ += Time.ElapsedTime();
  NumSolves__++;

#ifdef ML_AMESOS_DEBUG
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

#ifdef ML_AMESOS_DEBUG
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

#ifdef ML_AMESOS_DEBUG
    cout << "Amesos (level " << Level__
	 << ") : max (over solves) ||Ax - b|| = " << setiosflags(ios::scientific) << MaxError__ << endl;
#endif
    cout << endl;

  }
#endif
  
  Amesos_BaseSolver *A_Base = (Amesos_BaseSolver *) Amesos_Handle;
  const Epetra_LinearProblem *Amesos_LinearProblem;
  Amesos_LinearProblem = A_Base->GetProblem(); 

  delete Amesos_LinearProblem->GetOperator(); 

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
