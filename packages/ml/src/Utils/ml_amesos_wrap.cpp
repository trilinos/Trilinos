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
#include "Teuchos_ParameterList.hpp"

static double TimeForSolve__ = 0.0;
static int Level__ = -1;
static int NumSolves__ = 0;

#ifdef ML_AMESOS_DEBUG
static double MaxError__ = 0.0;
#endif

int ML_isKLUAvailable()
{
  Amesos Factory;
  if (Factory.Query("Amesos_Klu"))
    return(1); // this is true
  else
    return(0); // this is false
}
  
static void print_out(const Epetra_Comm& Comm, const int level, const char* what)
{
  if (Comm.MyPID() == 0 && ML_Get_PrintLevel() > 2)
#ifdef TFLOP
    printf("Amesos (level %d) : Building %s\n", level, what);
#else
    cout << "Amesos (level " << level << ") : Building " << what << "\n";
#endif
}

// ================================================ ====== ==== ==== == =

int ML_Amesos_Gen(ML *ml, int curr_level, int choice, int MaxProcs, 
                  double AddToDiag, void **Amesos_Handle)
{

  ML_Operator *Ke = &(ml->Amat[curr_level]);
  ML_Epetra::RowMatrix* Amesos_Matrix = 
    new ML_Epetra::RowMatrix(Ke, 0, false, ml->comm->USR_comm);
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

#if 0
  // MS // I don't like this output any more 
  if( ML_Get_PrintLevel() > 8 ) {
    AmesosList.set("PrintTiming",true);
    AmesosList.set("PrintStatus",true);
  }
#endif
  AmesosList.set("MaxProcs",MaxProcs);
  AmesosList.set("AddToDiag", AddToDiag);

  // don't use iterative refinement for Superludist only
  Teuchos::ParameterList & SuperludistList = AmesosList.sublist("Superludist");
  SuperludistList.set("IterRefine","NO");

  Amesos_BaseSolver* A_Base;
  Amesos A_Factory;
  const Epetra_Comm& Comm = Amesos_Matrix->Comm();

  switch (choice) {

  case ML_AMESOS_LAPACK:
    print_out(Comm, curr_level, "LAPACK");
    A_Base = A_Factory.Create("Amesos_Lapack", *Amesos_LinearProblem);
    break;

  case ML_AMESOS_UMFPACK:
    print_out(Comm, curr_level, "UMFPACK");
    A_Base = A_Factory.Create("Amesos_Klu", *Amesos_LinearProblem);
    break;

  case ML_AMESOS_SUPERLUDIST:
    print_out(Comm, curr_level, "SuperLU_DIST");
    A_Base = A_Factory.Create("Amesos_Superludist", *Amesos_LinearProblem);
    
    break;

  case ML_AMESOS_SUPERLU:
    print_out(Comm, curr_level, "SuperLU");
    A_Base = A_Factory.Create("Amesos_Superlu", *Amesos_LinearProblem);
    
    break;

  case ML_AMESOS_SCALAPACK:
    print_out(Comm, curr_level, "ScaLAPACK");
    A_Base = A_Factory.Create("Amesos_Scalapack", *Amesos_LinearProblem);
    
    break;

  case ML_AMESOS_MUMPS:
    print_out(Comm, curr_level, "MUMPS");
    A_Base = A_Factory.Create("Amesos_Mumps", *Amesos_LinearProblem);
    break;

  case ML_AMESOS_KLU:
  default:
    print_out(Comm, curr_level, "KLU");
    A_Base = A_Factory.Create("Amesos_Klu", *Amesos_LinearProblem);
    break;
  }

  // may happen the desired solver is not available. KLU is almost
  // always compiled, so try this first. If not, then LAPACK is
  // the last choice before quitting
  if (A_Base == 0) 
  {
    if (choice != ML_AMESOS_KLU)
    {
      if (Amesos_Matrix->Comm().MyPID() == 0 && ML_Get_PrintLevel() > 2)
      {
        cerr << "Amesos (level " << curr_level
             << ") : This coarse solver is not available." << endl;
        cerr << "Amesos (level " << curr_level
             << ") : Now re-building with KLU" << endl;
      }
      A_Base = A_Factory.Create("Amesos_Klu", *Amesos_LinearProblem);
    }
    if (A_Base == 0) 
    {
      if (Amesos_Matrix->Comm().MyPID() == 0) 
      {
        cerr << "Amesos (level " << curr_level
             << ") : This coarse solver is not available." << endl;
        cerr << "Amesos (level " << curr_level
             << ") : Now re-building with LAPACK" << endl;
      }
      A_Base = A_Factory.Create("Amesos_Lapack", *Amesos_LinearProblem);
      if (A_Base == 0) 
      {
        if (Amesos_Matrix->Comm().MyPID() == 0) 
        {
          cerr << "*ML*ERR* no Amesos solver is available!" << endl;
        }
        exit( EXIT_FAILURE );
      }
    }
  }

  A_Base->SetParameters(AmesosList);

  Epetra_Time Time(Amesos_Matrix->Comm());

  // Changed on 27-Nov-05, MS
  // It is faster to just call NumericFactorization(), otherwise the
  // code might have to ship the matrix twice, first to gather the
  // structure, then to gather the numerical values.
  //A_Base->SymbolicFactorization();
  //double Time1 = Time.ElapsedTime();
  Time.ResetStartTime();
  A_Base->NumericFactorization();
  double Time2 = Time.ElapsedTime();

  Level__ = -1;

#ifdef TFLOP
  if( Amesos_Matrix->Comm().MyPID() == 0 && ML_Get_PrintLevel()>2 ) {
    Level__ = curr_level;
    printf("Amesos (level %d) : Time for factorization = %f (s)\n",curr_level,Time2);
  }
#else
  if( Amesos_Matrix->Comm().MyPID() == 0 && ML_Get_PrintLevel()>2 ) {
    Level__ = curr_level;
    cout << "Amesos (level " << curr_level << ") : Time for factorization  = "
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
  
  ML_CHK_ERR(Ax.Update(1.0, EV_rhs, -1.0));
  
  double residual;
  ML_CHK_ERR(Ax.Norm2(&residual));
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

  delete A_Base ;
  delete Amesos_LinearProblem->GetOperator(); 
  delete Amesos_LinearProblem ;

}

#else

#include "ml_include.h"
#include "ml_amesos_wrap.h"
#include <stdio.h>

int ML_isKLUAvailable()
{
  return(0); // this is false
}

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
