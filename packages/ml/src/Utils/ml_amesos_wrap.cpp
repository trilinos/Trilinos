/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */
#include "ml_utils.h"
#include "ml_epetra_utils.h"
#include "ml_xyt.h"

#include "Epetra_Map.h" 
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h" 
#include "Epetra_LinearProblem.h"
#include "Epetra_Time.h"
#ifdef HAVE_ML_AMESOS
#include "ml_amesos.h"
#include "ml_amesos_wrap.h"
#include "Amesos_BaseSolver.h"
#include "Amesos_Factory.h" 
#include "AmesosClassType.h"

//  Jonathan - I need to convert an ml to an ML_Operator 
//  Did I pick off the right ML_Operator?

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
#include "Amesos_Parameter_List.h"

  extern "C" { 
  double *SendThisToEMV ; 
  }

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
    cout << "Amesos (level " << curr_level
	 << "): Time to convert to Epetra_CrsMatrix  = "
	 << Time1 << " (s)" << endl;
  }
  
  double NormInf = Amesos_CrsMatrix->NormInf();
  double NormOne = Amesos_CrsMatrix->NormOne();
  int NumGlobalRows = Amesos_CrsMatrix->NumGlobalRows();
  int NumGlobalNonzeros = Amesos_CrsMatrix->NumGlobalNonzeros();

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
	 << ") : fill-in = "
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
  
  Epetra_LinearProblem *Amesos_LinearProblem = new Epetra_LinearProblem;
  Amesos_LinearProblem->SetOperator( Amesos_CrsMatrix ) ; 

  AMESOS::Parameter::List ParamList ;

  AMESOS::Parameter::List SluParamList=ParamList.sublist("Superludist");

  // this is specific to Superludist-2.0
  if( choice == ML_AMESOS_SUPERLUDIST ) {
    
    if( MaxProcs == -2 ) {
      if( Amesos_CrsMatrix->RowMatrixRowMap().LinearMap() == true ) {
	ParamList.setParameter("Redistribute",false);
      } else {
	if( Amesos_CrsMatrix->Comm().MyPID() == 0 ) {
	  cout << "*ML*WRN* in Amesos_Smoother, you can set MaxProcs = -1\n"
	       << "*ML*WRN* (that is, matrix will not be redistributed)\n"
	       << "*ML*WRN* ONLY if the matrix map is linear. Now proceeding\n"
	       << "*ML*WRN* with redistribution of the matrix\n"
	       << "*ML*WRN* (file " << __FILE__ << ", line "
	       << __LINE__ << ")\n";
	}
      }
    } else {
      ParamList.setParameter("Redistribute",true);
      SluParamList.setParameter("MaxProcesses",MaxProcs);
    }
  }

  Amesos_BaseSolver* A_Base;
  Amesos_Factory A_Factory;

  // I got the impression that small problems are "safer"
  // in other hands than superludist ones.
  // Certo che 'stp superludist e` proprio 'na schifezza ;)
  // ????? brrrrr, what the hell is this ???????
  if( NumGlobalRows < 4*MaxProcs || NumGlobalRows < 16 ) choice = ML_AMESOS_KLU;
  
  switch( choice ) {

  case ML_AMESOS_UMFPACK:
    if( Amesos_CrsMatrix->Comm().MyPID() == 0 && ML_Get_PrintLevel()>2 )
      cout << "Amesos (level " << curr_level
	   << ") : building UMFPACK\n";
    A_Base = A_Factory.Create(AMESOS_UMFPACK, *Amesos_LinearProblem, ParamList );
    assert(A_Base!=0);
    break;

  case ML_AMESOS_SUPERLUDIST:
    if( Amesos_CrsMatrix->Comm().MyPID() == 0 && ML_Get_PrintLevel()>2 )
      cout << "Amesos (level " << curr_level
	   << ") : building SUPERLUDIST\n";
    A_Base = A_Factory.Create( AMESOS_SUPERLUDIST, *Amesos_LinearProblem, ParamList );
    
    assert(A_Base!=0);
    break;

  case ML_AMESOS_KLU:
  default:
    if( Amesos_CrsMatrix->Comm().MyPID() == 0 && ML_Get_PrintLevel()>2 )
      cout << "Amesos (level " << curr_level
	   << ") : building KLU\n";
    A_Base = A_Factory.Create(AMESOS_KLU, *Amesos_LinearProblem, ParamList );
    assert(A_Base!=0);
    break;
  }


  Epetra_Time Time(Amesos_CrsMatrix->Comm());

  A_Base->SymbolicFactorization();
  Time1 = Time.ElapsedTime();
  Time.ResetStartTime();
  A_Base->NumericFactorization();
  Time2 = Time.ElapsedTime();
  
  if( Amesos_CrsMatrix->Comm().MyPID() == 0 && ML_Get_PrintLevel()>2 ) {
    cout << "Amesos (level " << curr_level
	 << ") : Time for symbolic fact  = "
	 << Time1 << " (s)" << endl;
    cout << "Amesos (level " << curr_level
	 << ") : Time for numerical fact = "
	 << Time2 << " (s)" << endl;
  }

  *Amesos_Handle = (void *) A_Base ;

  return 0;
}

int ML_Amesos_Solve( void *Amesos_Handle, double x[], double rhs[] )
{

  Amesos_BaseSolver *A_Base = (Amesos_BaseSolver *) Amesos_Handle ;
  Epetra_LinearProblem *Amesos_LinearProblem = (Epetra_LinearProblem *) A_Base->GetProblem() ; 

  Epetra_BlockMap map = Amesos_LinearProblem->GetOperator()->OperatorDomainMap() ; 

  Epetra_Vector EV_rhs( View, map, rhs ) ;
  Epetra_Vector EV_lhs( View, map, x ) ;

  Amesos_LinearProblem->SetRHS( &EV_rhs ) ; 
  Amesos_LinearProblem->SetLHS( &EV_lhs ) ;

  A_Base->Solve() ; 

  return 0;
}

void ML_Amesos_Destroy(void *Amesos_Handle)
{

  Amesos_BaseSolver *A_Base = (Amesos_BaseSolver *) Amesos_Handle ;
  const Epetra_LinearProblem *Amesos_LinearProblem;
  Amesos_LinearProblem = A_Base->GetProblem() ; 
  const Epetra_Operator *EO = Amesos_LinearProblem->GetOperator() ; 

  delete Amesos_LinearProblem->GetOperator() ; 

  delete Amesos_LinearProblem ;
  delete A_Base ;
}


#else

int ciao=0;

#endif
