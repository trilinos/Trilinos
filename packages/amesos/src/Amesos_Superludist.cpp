  /* Copyright (2003) Sandia Corportation. Under the terms of Contract 
   * DE-AC04-94AL85000, there is a non-exclusive license for use of this 
   * work by or on behalf of the U.S. Government.  Export of this program
   * may require a license from the United States Government. */


  /* NOTICE:  The United States Government is granted for itself and others
   * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
   * license in ths data to reproduce, prepare derivative works, and
   * perform publicly and display publicly.  Beginning five (5) years from
   * July 25, 2001, the United States Government is granted for itself and
   * others acting on its behalf a paid-up, nonexclusive, irrevocable
   * worldwide license in this data to reproduce, prepare derivative works,
   * distribute copies to the public, perform publicly and display
   * publicly, and to permit others to do so.
   * 
   * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
   * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
   * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
   * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
   * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
   * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */

#include "Amesos_Superludist.h"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_Util.h"
#include "CrsMatrixTranspose.h"

int SLU_NumProcRows( int NumProcs ) {
#ifdef TFLOP
  //  Else, parameter 6 of DTRSV CTNLU is incorrect 
  return 1;
#else
  int i;
  int NumProcRows ;
  for ( i = 1; i*i <= NumProcs; i++ ) 
    ;
  bool done = false ;
  for ( NumProcRows = i-1 ; done == false ; ) {
    int NumCols = NumProcs / NumProcRows ; 
    if ( NumProcRows * NumCols == NumProcs ) 
      done = true; 
    else 
      NumProcRows-- ; 
  }
  return NumProcRows;
#endif
}

  //=============================================================================
  Amesos_Superludist::Amesos_Superludist(const Epetra_LinearProblem &prob, 
				 const AMESOS::Parameter::List &ParameterList ) :  

    UseTranspose_(false), 
    GridCreated_(0), 
    FactorizationDone_(0), 
    UniformMap_(0) 
{

  Problem_ = &prob ; 
  ParameterList_ = &ParameterList ; 
}

//=============================================================================
Amesos_Superludist::~Amesos_Superludist(void) {

  if ( FactorizationDone_ ) {
    SUPERLU_FREE( superluA_.Store );
    ScalePermstructFree(&ScalePermstruct_);
    Destroy_LU(numrows_, &grid_, &LUstruct_);
    LUstructFree(&LUstruct_);
    if ( options_.SolveInitialized ) {
      dSolveFinalize(&options_, &SOLVEstruct_ ) ; 
    }
  }
  if ( GridCreated_ ) {
    superlu_gridexit(&grid_);
  }
  if (UniformMap_ ) delete UniformMap_ ; 
}



int Amesos_Superludist::PerformSymbolicFactorization() {

  //  Superludist does not offer SymbolicFactorization

  return 0;
}

int Amesos_Superludist::PerformNumericFactorization( ) {

  //
  //  Cast input matrix to a CrsMatrix 
  //
  Epetra_RowMatrix *RowMatrixA = 
    dynamic_cast<Epetra_RowMatrix *>(Problem_->GetOperator());
  EPETRA_CHK_ERR( RowMatrixA == 0 ) ; 

  Epetra_CrsMatrix *CastCrsMatrixA = dynamic_cast<Epetra_CrsMatrix*>(RowMatrixA) ; 
  EPETRA_CHK_ERR( CastCrsMatrixA == 0 ) ; 

  // Ken note:  If the following does not work, move OrignalMap_ to Amesos_Superludist.h
  // 
  const Epetra_Map OriginalMap = CastCrsMatrixA->RowMap() ; 
  const Epetra_Comm &Comm = RowMatrixA->Comm();
  
  int iam = Comm.MyPID() ;

  //
  //  Compute a cannonical uniform distribution:
  //    MyFirstElement - The first element which this processor would have
  //    NumExpectedElemetns - The number of elements which this processor would have
  //
  numrows_ = OriginalMap.NumGlobalElements() ; 
  assert( numrows_ == CastCrsMatrixA->NumGlobalRows() ) ; 
  int numcols = CastCrsMatrixA->NumGlobalCols() ; 
  assert( numrows_ == numcols ) ; 

  int m_per_p = numrows_ / Comm.NumProc() ;
  int remainder = numrows_ - ( m_per_p * Comm.NumProc() );
  int MyFirstElement = iam * m_per_p + EPETRA_MIN( iam, remainder );
  int MyFirstNonElement = (iam+1) * m_per_p + EPETRA_MIN( iam+1, remainder );
  int NumExpectedElements = MyFirstNonElement - MyFirstElement ; 

  //
  //  Convert the matrix to the form needed by SuperLU
  //  If the matrix is already in the right form, 
  //    SuperLUmat points to CastCrsMatrixA 
  //  Else
  //    The matrix is redistributed to UniformMatrix 
  //    SuperLUmat points to UniformMatrix
  //
  //  KEN NOTE:  debugxx work - redistribute should be set automatically
  //
  Epetra_CrsMatrix *SuperLUmat = 0 ;
  UniformMap_ = new Epetra_Map(  OriginalMap.NumGlobalElements(), NumExpectedElements, 0, Comm );
  Epetra_CrsMatrix UniformMatrix(Copy, *UniformMap_, 0);

  redistribute_ = true ;
  if ( redistribute_ ) {

    Epetra_Export export_to_dist( OriginalMap, *UniformMap_);

    UniformMatrix.Export( *CastCrsMatrixA, export_to_dist, Add ); 
    UniformMatrix.TransformToLocal() ; 
    SuperLUmat = &UniformMatrix ;

    
  } else {
    {  
      EPETRA_CHK_ERR( ! ( OriginalMap.LinearMap())  ) ; // Map must be contiguously divided
      //
      //  This is another way to check that the distribution is as pdgssvx expects it
      //  (i.e. a linear map)
      //
      int Phase2NumElements = OriginalMap.NumMyElements() ; 
      vector <int> MyRows( Phase2NumElements ) ; 
      OriginalMap.MyGlobalElements( &MyRows[0] ) ; 
      for (int row = 0 ; row < Phase2NumElements ; row++ ) {
	EPETRA_CHK_ERR( MyFirstElement+row != MyRows[row] ) ;
      }
    }
    SuperLUmat = CastCrsMatrixA ;
  }

  int MyActualFirstElement = SuperLUmat->RowMap().MinMyGID() ; 
  int NumMyElements = SuperLUmat->NumMyRows() ; 
  vector <int> MyRowPtr( NumMyElements+1 ) ;  
  //
  //  This is actually redundant with the Ap below
  //
  MyRowPtr[0] = 0 ; 
  int CurrentRowPtr = 0 ;
  for ( int i = 0; i < NumMyElements ; i++ ) { 
    CurrentRowPtr += SuperLUmat->NumMyEntries( i ) ; 
    MyRowPtr[i+1] = CurrentRowPtr ; 
  }

  //
  //  Ap, Ai, Aval form the compressed row storage used by Klu
  //
  vector <int> Ap;
  vector <int> Ai;
  vector <double> Aval;

  //
  //  Extract Ai, Ap and Aval from SuperLUmat
  //

  int nnz_loc = SuperLUmat->NumMyNonzeros() ;
  //
  //  Step 6) Convert the matrix to Ap, Ai, Aval
  //
  Ap.resize( NumMyElements+1 );
  Ai.resize( EPETRA_MAX( NumMyElements, nnz_loc) ) ; 
  Aval.resize( EPETRA_MAX( NumMyElements, nnz_loc) ) ; 
  
  int NzThisRow ;
  double *RowValues;
  int *ColIndices;
  int Ai_index = 0 ; 
  int MyRow;
  int num_my_cols = SuperLUmat->NumMyCols() ; 
  vector <int>Global_Columns( num_my_cols ) ; 
  for ( int i = 0 ; i < num_my_cols ; i ++ ) { 
    Global_Columns[i] = SuperLUmat->GCID( i ) ; 
  }
  
  for ( MyRow = 0; MyRow < NumMyElements ; MyRow++ ) {
    int status = SuperLUmat->ExtractMyRowView( MyRow, NzThisRow, RowValues, ColIndices ) ;
    assert( status == 0 ) ; 
    Ap[MyRow] = Ai_index ; 
    assert( Ap[MyRow] == MyRowPtr[MyRow] ) ; 
    for ( int j = 0; j < NzThisRow; j++ ) { 
      Ai[Ai_index] = Global_Columns[ColIndices[j]] ; 
      Aval[Ai_index] = RowValues[j] ; 
      Ai_index++;
    }
  }
  assert( NumMyElements == MyRow );
  Ap[ NumMyElements ] = Ai_index ; 

  //
  //  Call SuperludistOO
  //
  //
  //
  const Epetra_MpiComm & comm1 = dynamic_cast<const Epetra_MpiComm &> (Comm);
  MPI_Comm MPIC = comm1.Comm() ;

  if ( ! GridCreated_ ) {
    GridCreated_ = true; 
    int numprocs = Comm.NumProc() ;                 
    nprow_ = SLU_NumProcRows( numprocs ) ; 
    npcol_ = numprocs / nprow_ ;
    assert ( nprow_ * npcol_ == numprocs ) ; 
    superlu_gridinit( MPIC, nprow_, npcol_, &grid_);
  } else {
    assert( nprow_ * npcol_ == Comm.NumProc() ) ; 
  }

  /* Bail out if I do not belong in the grid. */
  if ( iam < nprow_ * npcol_ ) {
	
    if ( FactorizationDone_ ) { 
      assert( false ) ; 
      //      free a:  ken note debugxx work - free stuff here
      //  Copy from the destructor - or not - do we really need to?
    } 
    set_default_options(&options_);
    
    dCreate_CompRowLoc_Matrix_dist( &superluA_, numrows_, numcols, 
				    nnz_loc, NumMyElements, MyActualFirstElement,
				    &Aval[0], &Ai[0], &Ap[0], 
				    SLU_NR_loc, SLU_D, SLU_GE );
    
    FactorizationDone_ = true; 

    /* Initialize ScalePermstruct and LUstruct. */
    ScalePermstructInit(numrows_, numcols, &ScalePermstruct_);
    LUstructInit(numrows_, numcols, &LUstruct_);
    
    assert( options_.Fact == DOFACT );  
    options_.Fact = DOFACT ;       
    
    SuperLUStat_t stat;
    PStatInit(&stat);    /* Initialize the statistics variables. */

    int info ;
    double berr ;    //  Should be untouched
    double xValues;  //  Should be untouched
    int nrhs = 0 ;   //  Prevents forward and back solves
    int ldx = 1;     //  Should be untouched
    pdgssvx(&options_, &superluA_, &ScalePermstruct_, &xValues, ldx, nrhs, &grid_,
	    &LUstruct_, &SOLVEstruct_, &berr, &stat, &info);
    EPETRA_CHK_ERR( info ) ; 

    PStatFree(&stat);
  }

  return 0;
}




bool Amesos_Superludist::MatrixShapeOK() const { 
  bool OK ;

  if ( GetProblem()->GetOperator()->OperatorRangeMap().NumGlobalPoints() != 
       GetProblem()->GetOperator()->OperatorDomainMap().NumGlobalPoints() ) OK = false;
  return OK; 
}


int Amesos_Superludist::SymbolicFactorization() {

  PerformSymbolicFactorization();

  return 0;
}

int Amesos_Superludist::NumericFactorization() {
  
  PerformNumericFactorization( );
  return 0;
}


int Amesos_Superludist::Solve() { 


  Epetra_MultiVector   *vecX = Problem_->GetLHS() ; 
  Epetra_MultiVector   *vecB = Problem_->GetRHS() ; 

  //
  //  Compute the number of right hands sides (and check that X and B have the same shape) 
  //
  int nrhs; 
  if ( vecX == 0 ) { 
    nrhs = 0 ;
    EPETRA_CHK_ERR( vecB != 0 ) ; 
  } else { 
    nrhs = vecX->NumVectors() ; 
    EPETRA_CHK_ERR( vecB->NumVectors() != nrhs ) ; 
  }

  //
  //  Set up vecXptr and vecBptr 
  //
  double *bValues ;
  double *xValues ;
  int ldb, ldx ; 

  Epetra_MultiVector vecXdistributed( *UniformMap_, nrhs ) ; 
  Epetra_MultiVector vecBdistributed( *UniformMap_, nrhs ) ; 
  Epetra_MultiVector* vecXptr; 
  Epetra_MultiVector* vecBptr; 

  
  if ( redistribute_ ) {
    Epetra_RowMatrix *RowMatrixA = dynamic_cast<Epetra_RowMatrix *>(Problem_->GetOperator());
    Epetra_CrsMatrix *CastCrsMatrixA = dynamic_cast<Epetra_CrsMatrix*>(RowMatrixA) ; 
    const Epetra_Map OriginalMap = CastCrsMatrixA->RowMap() ; 
    Epetra_Import ImportToDistributed( *UniformMap_, OriginalMap);

    vecXdistributed.Import( *vecX, ImportToDistributed, Insert ) ;
    vecBdistributed.Import( *vecB, ImportToDistributed, Insert ) ;

    vecXptr = &vecXdistributed ; 
    vecBptr = &vecBdistributed ; 
  } else {
    vecXptr = vecX ; 
    vecBptr = vecB ; 
  }


  int NumMyElements = vecBptr->MyLength(); 
  EPETRA_CHK_ERR( vecBptr->ExtractView( &bValues, &ldb ) )  ; 
  EPETRA_CHK_ERR( vecXptr->ExtractView( &xValues, &ldx ) ) ; 
  EPETRA_CHK_ERR( ! ( ldx == ldb ) ) ; 

  //
  //  pdgssvx returns x in b, so we copy b into x.  
  //
  for ( int j = 0 ; j < nrhs; j++ )
    for ( int i = 0 ; i < NumMyElements; i++ ) xValues[i+j*ldx] = bValues[i+j*ldb]; 
  
  Epetra_RowMatrix *RowMatrixA = dynamic_cast<Epetra_RowMatrix *>(Problem_->GetOperator());
  const Epetra_Comm &Comm = RowMatrixA->Comm();
  int iam = Comm.MyPID() ;
  
  /* Bail out if I do not belong in the grid. */
  if ( iam < nprow_ * npcol_ ) {
    int info ;
    vector<double>berr(nrhs);
    SuperLUStat_t stat;
    PStatInit(&stat);    /* Initialize the statistics variables. */
    
    assert( GridCreated_ ) ; 
    assert( FactorizationDone_ ) ; 
    options_.Fact = FACTORED ;       
    pdgssvx(&options_, &superluA_, &ScalePermstruct_, &xValues[0], ldx, nrhs, &grid_,
	    &LUstruct_, &SOLVEstruct_, &berr[0], &stat, &info);
    EPETRA_CHK_ERR( info ) ; 
    
    PStatFree(&stat);
  }
  
  if ( redistribute_ ) { 
    Epetra_RowMatrix *RowMatrixA = dynamic_cast<Epetra_RowMatrix *>(Problem_->GetOperator());
    Epetra_CrsMatrix *CastCrsMatrixA = dynamic_cast<Epetra_CrsMatrix*>(RowMatrixA) ; 
    const Epetra_Map OriginalMap = CastCrsMatrixA->RowMap() ; 
    Epetra_Import ImportBackToOriginal( OriginalMap,*UniformMap_);
    
    vecX->Import( *vecXptr, ImportBackToOriginal, Insert ) ;
  }

  return(0) ; 
}
