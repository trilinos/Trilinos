/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

#include "ml_config.h"
#include "ml_include.h"
#ifdef HAVE_ML_IFPACK
#include "ml_utils.h"
#include "ml_epetra_utils.h"

#include "Epetra_Map.h" 
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h" 
#include "Epetra_LinearProblem.h"
#include "Epetra_Time.h"
#include "ml_ifpack.h"
#include "ml_ifpack_wrap.h"

#include "Epetra_Operator.h"

#include "Ifpack_CrsIct.h"
#include "Ifpack_CrsRiluk.h"
//#include "Ifpack_CrsRick.h"
#include "Ifpack_IlukGraph.h"

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

// ================================================ ====== ==== ==== == =

int ML_Ifpack_Gen(ML *ml, int curr_level, int choice, int * options,
		  double * params, void ** Ifpack_Handle)
{

  char msg[80];
  sprintf( msg, "Ifpack (level %d) : ", curr_level);
  string Msg(msg);
  
  int MaxNumNonzeros;
  double Time1;
  
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
    Overlap = 0;
    LevelOfFill = 0;
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

// ================================================ ====== ==== ==== == =

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

// ================================================ ====== ==== ==== == =

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

#ifndef DOXYGEN_SHOULD_SKIP_THIS
static int ciao=0;
#endif


#endif /* #ifdef HAVE_ML_IFPACK */
