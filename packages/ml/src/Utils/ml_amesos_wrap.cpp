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
#include "Epetra_Operator.h"

#include "Ifpack_CrsIct.h"
#include "Ifpack_CrsRiluk.h"
//#include "Ifpack_CrsRick.h"
#include "Ifpack_IlukGraph.h"

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

int ML_Ifpack_Gen(ML *ml, int curr_level, int choice, int * options,
		  double * params, void ** Ifpack_Handle)
{

  char msg[80];
  sprintf( msg, "Ifpack (level %d) : ", curr_level);
  string Msg(msg);
  
  int MaxNumNonzeros;
  double Time1, Time2;
  
  ML_Operator *Ke = &(ml->Amat[curr_level]);
  Epetra_CrsMatrix * Ifpack_CrsMatrix;

  ML_Operator2EpetraCrsMatrix( Ke, Ifpack_CrsMatrix, MaxNumNonzeros,
			       true, Time1);
  
  if( Ifpack_CrsMatrix->Comm().MyPID() == 0 && ML_Get_PrintLevel()>2 ) {
    cout << Msg << "Time to convert to Epetra_CrsMatrix  = "
	 << Time1 << " (s)" << endl;
  }
  
  double NormInf = Ifpack_CrsMatrix->NormInf();
  double NormOne = Ifpack_CrsMatrix->NormOne();
  int NumGlobalRows = Ifpack_CrsMatrix->NumGlobalRows();
  int NumGlobalNonzeros = Ifpack_CrsMatrix->NumGlobalNonzeros();

  if( Ifpack_CrsMatrix->Comm().MyPID() == 0 && ML_Get_PrintLevel() > 2 ) {
    cout << Msg << "NumGlobalRows = " << NumGlobalRows << endl;
    cout << Msg << "NumGlobalNonzeros = " << NumGlobalNonzeros << endl;
    cout << Msg << "MaxNonzerosPerRow = " << MaxNumNonzeros << endl;
    cout << Msg << "fill-in = " << 100.0*NumGlobalNonzeros/(NumGlobalRows*NumGlobalRows)
	 << " %" << endl;
    cout << Msg << "||A||_\\infty = " << NormInf << " ||A||_1 = " << NormOne << endl;
  }
   
  // ============================ //
  // Construct ILU preconditioner //
  // ---------------------------- //

  double RelaxValue;
  int Overlap;
  int LevelOfFill;
  
  if( params != NULL ) {
    RelaxValue  = params[ML_IFPACK_RELAX_VALUE];
  } else {
    RelaxValue = 1.0;
  }

  if( options != NULL ) {
    Overlap     = options[ML_IFPACK_OVERLAP];
    LevelOfFill = options[ML_IFPACK_LEVEL_OF_FILL];
  } else {
    Overlap = 1;
    LevelOfFill = 1;
  }

  double Condest;
  int    ierr;
  
  switch( choice ) {

  case ML_IFPACK_RILUK:
    {
      
    Ifpack_IlukGraph * IlukGraph = new Ifpack_IlukGraph(Ifpack_CrsMatrix->Graph(),
							LevelOfFill,Overlap);
    
    assert(IlukGraph->ConstructFilledGraph()==0);
    
    Ifpack_CrsRiluk * RILU = NULL;
    RILU = new Ifpack_CrsRiluk(*IlukGraph);
    
    RILU->SetRelaxValue(RelaxValue);
    
    ierr = RILU->InitValues(*Ifpack_CrsMatrix);
    if (ierr!=0) cout << "*ERR* InitValues = " << ierr;

    assert(RILU->Factor()==0);

    // and now estimate the condition number
    RILU->Condest(false,Condest);

    *Ifpack_Handle = (void *) RILU ;
    }
    
    break;

  case ML_IFPACK_ICT:


    break;
    
  case ML_IFPACK_RICK:

    break;

  default:
    cerr << "*ML*ERR* choice in input to Ifpack not valid" << endl;
    exit( EXIT_FAILURE );
  }
  
  // qui metti la media dei condnum
  if( Ifpack_CrsMatrix->Comm().MyPID() == 0 && ML_Get_PrintLevel()>5 ) {
    cout << Msg << "overlap = "  << Overlap <<  endl;
    cout << Msg << "LevelOfFill = "  << LevelOfFill <<  endl;
    cout << Msg << "Condnum estimate " << Condest << endl;
  }

  delete Ifpack_CrsMatrix; Ifpack_CrsMatrix = NULL;

  return 0;
  
} /* ML_Ifpack_Gen */


int ML_Ifpack_Solve( void * Ifpack_Handle, double * x, double * rhs )
{

  Epetra_Operator * Handle = (Epetra_Operator *)Ifpack_Handle;
  Ifpack_CrsIct   * ICT   = dynamic_cast<Ifpack_CrsIct *>(Handle);
  Ifpack_CrsRiluk * RILUK = dynamic_cast<Ifpack_CrsRiluk *>(Handle);
  //  Ifpack_CrsRick  * RICK  = (Ifpack_CrsRick *)  Ifpack_Handle ;
  
  if( ICT != NULL ) {
    Epetra_Vector Erhs( View, (ICT->OperatorRangeMap()), rhs ) ;
    Epetra_Vector Ex( View, (ICT->OperatorDomainMap()), x ) ;
    ICT->Solve(false,Erhs,Ex); 
  } else if( RILUK != NULL ) {
    Epetra_Vector Erhs( View, (RILUK->OperatorRangeMap()), rhs ) ;
    Epetra_Vector Ex( View, (RILUK->OperatorDomainMap()), x ) ;
    RILUK->Solve(false,Erhs,Ex); 
    //  } else if( RICK != NULL ) {
    //    Epetra_Vector Erhs( View, (RICK->OperatorRangeMap()), rhs ) ;
    //    Epetra_Vector Ex( View, (RICK->OperatorDomainMap()), x ) ;
    //    RICK->Solve(false,Erhs,Ex); 
  } else {
    cerr << "*ML*ERR* Something wrong in `ML_Ifpack_Solve'" << endl;
    exit( EXIT_FAILURE );
  }
    
  return 0;

} /* ML_Ifpack_Solve */

void ML_Ifpack_Destroy(void * Ifpack_Handle)
{

  cout << "to deallocate  Ifpack_OverlapGraph * = */ ???? " << endl;

  Epetra_Operator * Handle = (Epetra_Operator *)Ifpack_Handle;
  
  Ifpack_CrsIct * ICT   = dynamic_cast<Ifpack_CrsIct *>(Handle);
  Ifpack_CrsRiluk * RILUK = dynamic_cast<Ifpack_CrsRiluk *>(Handle);
  //  Ifpack_CrsRick  * RICK  = (Ifpack_CrsRick *)  Ifpack_Handle ;

  if( ICT != NULL ) {
    delete ICT; ICT = NULL;
  } else if( RILUK != NULL ) {
    delete RILUK; RILUK = NULL;
    //  } else if( RICK != NULL ) {
    //    delete RICK; RICK = NULL;
  } else {
    cerr << "*ML*ERR* Something wrong in `ML_Ifpack_Destroy'" << endl;
    exit( EXIT_FAILURE );
  }
  
} /* ML_Ifpack_Destroy */

#else

int ciao=0;

#endif
