#include "SuperludistOO.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
    This code cannot be compiled without mpi.h.
#include "Epetra_Comm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "Epetra_Vector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_Operator.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_CrsMatrix.h"
#include "CrsMatrixTranspose.h"
#include <vector>

#ifdef DEBUG
#include "Comm_assert_equal.h"
    // #include "CrsMatricesAreIdentical.h"
#endif
#ifdef EPETRA_CRSMATRIX_CONSTRUCT_FROM_ROWMATRIX
#include "ExtractCrsFromRowMatrix.h"
#endif
/*
  Returns the largest number of rows that allows NumProcs to be used within a 
  rectangular grid with no more rows than columns.

  i.e. max i such that there exists j > i such that i*j = NumProcs
*/
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
SuperludistOO::SuperludistOO(const Epetra_LinearProblem &prob ) {
  //  AllocAzArrays();

  Problem_ = &prob ; 
  Transpose_ = false; 
  A_and_LU_built = false ; 
  Factored_ = false ; 
  FirstCallToSolve_ = true ; 
  //
  //  The following are initialized just on general principle
  //
  numprocs = -13 ; 
  nprow = -13 ; 
  npcol = -13 ; 
  numrows = -13 ; 

}

//=============================================================================
SuperludistOO::~SuperludistOO(void) {
  //  DeleteMemory();
  //  DeleteAzArrays();

  if ( false & A_and_LU_built ) { 
    // Destroy_CompCol_Matrix_dist(&A);
    SUPERLU_FREE(A.Store);
    Destroy_LU(numrows, &grid, &LUstruct);
    ScalePermstructFree(&ScalePermstruct);
    LUstructFree(&LUstruct);
    SUPERLU_FREE(berr);
  }

}

//=============================================================================

//
//  Solve() uses several intermediate matrices to convert the input matrix
//  to one that we can pass to the Sparse Direct Solver
//
//  Epetra_RowMatrix *RowMatrixA - The input matrix
//  Epetra_CrsMatrix *CastCrsMatrixA - The input matrix casted to a crs matrix
//  Epetra_CrsMatrix ExtractCrsMatrixA - Converted to a Crs matrix 
//                                 (Unused if RowMatrix is an Epetra_CrsMatrix)
//  Epetra_CrsMatrix *Phase2Mat - Guaranteed to be a CrsMatrix
//  Epetra_CrsMatrix SerialCrsMatrixA - Phase2Mat coalesced to one process
//  Epetra_CrsMatrix *Phase3Mat - A pseudo-distributed CrsMatrix
//    (i.e. stored exclusively on one process)
//  Epetra_CrsMatrix Phase3MatTrans - The transpose of Phase3Mat
//  Epetra_CrsMatrix *Phase4Mat - A pseudo-serial CrsMatrix with the 
//    proper transposition
//  Epetra_CrsMatrix Phase5Mat - A replicated CrsMatrix with the 
//    proper transposition 
//
//  This is what the code does:
//  Step 1)  Convert the matrix to an Epetra_CrsMatrix
//  Step 2)  Coalesce the matrix onto process 0
//  Step 3)  Transpose the matrix 
//  Step 4)  Replicate the matrix
//  Step 5)  Convert vector b to a replicated vector
//  Step 6)  Convert the matrix to Ap, Ai, Aval
//  Step 7)  Call SuperLUdist
//  Step 8)  Convert vector x back to a distributed vector
//
//  Remaining tasks:
//  1)  I still need to make it possible for SuperludistOO to accept a 
//  replicated matrix without using any additional memory.
//  I believe that all that I need is to move the definition
//  of ExtractCrsMatrixA,  SerialCrsMatrixA and Phase3MatTrans up 
//  to the top of the code and add a conditional which tests to 
//  see if RowMatrixA is actually in the exact format that we need,
//  and if so, skip all the matrix transformations.  Also, Phase5Mat 
//  needs to be changed to Phase4Replicated with an added pointer, named
//  *Phase5Mat.
//  2)  Can we handle a transposed matrix any cleaner?  We build Ap, 
//  Ai and Avalues - can we build that as a transpose more efficiently
//  than doing a full CrsMatrix to CrsMatrix transpose?
//  
//  Memory usage:
//    ExtractCrsMatrixA - 1 if RowMAtrixA is not a CrsMatrix
//    SerialCrsMatrixA - 1 if RowMatrixA is not a serial matrix
//    Phase3MatTrans - 1 if RowMatrixA unless a transpose solve is requested
//    Phase5Mat - 1 
//  If we need SerialCrsMAttrixA, then ExtractCrsMatrixA will not be a
//  serial matrix, hence three extra serial copies is the maximum.
//

#undef EPETRA_CHK_ERR
#define EPETRA_CHK_ERR(xxx) assert( (xxx) == 0 ) 

int SuperludistOO::Solve(bool factor) { 
  //
  //  I am going to put these here until I determine that I need them in 
  //  SuperludistOO.h 
  //


  SOLVEstruct_t SOLVEstruct;    // This - and many other variables will 
                                // need to move to SuperludistOO.h before we can 
                                // make multiple solves work.
  bool CheckExtraction = false;    //  Set to true to force extraction for unit test

  Epetra_RowMatrix *RowMatrixA = 
    dynamic_cast<Epetra_RowMatrix *>(Problem_->GetOperator());
  
  EPETRA_CHK_ERR( RowMatrixA == 0 ) ; 

  Epetra_CrsMatrix *CastCrsMatrixA = dynamic_cast<Epetra_CrsMatrix*>(RowMatrixA) ; 
#ifdef EPETRA_CRSMATRIX_CONSTRUCT_FROM_ROWMATRIX
  Epetra_CrsMatrix *ExtractCrsMatrixA = 0;
#endif
  Epetra_CrsMatrix *Phase2Mat = 0 ;
  const Epetra_Comm &Comm = RowMatrixA->Comm();
  
  int iam = Comm.MyPID() ;
#if 0
  //
  //  The following lines allow us time to attach the debugger
  //
  int hatever;
  if ( iam == 0 )  cin >> hatever ; 
  Comm.Barrier();
#endif

  //  return 0; // WORK GXX BOGUS KEN 

  //
  //  Step 1)  Convert the matrix to an Epetra_CrsMatrix
  //
  //  If RowMatrixA is not a CrsMatrix, i.e. the dynamic cast fails, 
  //  extract a CrsMatrix from the RowMatrix.
  //
  if ( CastCrsMatrixA != 0 && ! CheckExtraction ) { 
    Phase2Mat = CastCrsMatrixA ; 
  } else {
#ifdef EPETRA_CRSMATRIX_CONSTRUCT_FROM_ROWMATRIX
    ExtractCrsMatrixA = new Epetra_CrsMatrix( *RowMatrixA ) ; 

    Phase2Mat = ExtractCrsMatrixA ; 
#ifdef DEBUG
    if ( CheckExtraction ) 
      assert( CrsMatricesAreIdentical( CastCrsMatrixA, ExtractCrsMatrixA ) ) ; 
#endif
#else
    assert( false ) ;
#endif
  }
  assert( Phase2Mat != NULL ) ; 
  const Epetra_Map &Phase2Matmap = Phase2Mat->RowMap() ; 

  //  return 0; // WORK GXX BOGUS KEN -  OK, if we put the return we get a big error, but at least we don't 
            // do that bizarre double return thing.

  //
  //  Glossary:
  //    numrows, numcols = m,n, i.e. the size of the full, global, matrix
  //    m_loc            = the number of rows owned by this process
  //    nnz_loc          = the number of non zeros in the rows owned by this process
  //    Ap, Aval, Ai,    = rowptr, nzval, colind, = sparse row storage 
  //    MyRowPtr         = a redundant computation of Ap  (rowptr) 

  //
  //  Compute nnz_loc, m_loc, and MyRowPtr from the map
  //  

#if 0
  cout << " A  Comm.NumProc() = " <<  Comm.NumProc() << endl ; 

  cout << " true = " << true << " LInMap = " <<  Phase2Matmap.LinearMap() << endl  ; // Map must be contiguously divided
  cout << " SuperludistOO.cpp::   traceback mode = " << Epetra_Object::GetTracebackMode() << endl ; 
  cerr << " Send this to cerr cerr cerr   traceback mode = " << Epetra_Object::GetTracebackMode() << endl ; 
#endif
  EPETRA_CHK_ERR( ! ( Phase2Matmap.LinearMap())  ) ; // Map must be contiguously divided

  int numrows = Phase2Matmap.NumGlobalElements() ; 
  assert( numrows == Phase2Mat->NumGlobalRows() ) ; 
  int numcols = Phase2Mat->NumGlobalCols() ; 
  assert( numrows == numcols ) ; 

  int m_loc = Phase2Matmap.NumMyElements() ; 
  int nnz_loc = Phase2Mat->NumMyNonzeros() ;
  vector <int> MyRowPtr( m_loc+1 ) ;  

  //  Here is another attempted hack
  int MyFirstElement ;
  //  assert( Comm.NumProc() <= 2 ) ; 
  //  if ( iam == 0 ) MyFirstElement = 0 ;
  //  if ( iam == 1 ) MyFirstElement = numrows - m_loc ; // Proc 0 has the rest
  

#if 1
  //
  //  Here I compute what a uniform distribution should look like
  //  so that I can compare what I think MyFirstElement should be 
  //  against what it really is.
  //
  int m_per_p = numrows / Comm.NumProc() ;
  cout << " m_per_p = " << m_per_p << endl ; 
  int remainder = numrows - ( m_per_p * Comm.NumProc() );
  MyFirstElement = iam * m_per_p + EPETRA_MIN( iam, remainder );
#endif
  cout << " iam = " << iam << " MyFirstElement = " << MyFirstElement << endl ; 
  if ( ( numrows == 5 ) && ( Comm.NumProc() == 2)  ) {
    assert( iam ==0 || MyFirstElement == 3 ) ; 
    assert( iam ==1 || MyFirstElement == 0 ) ; 
  }

  //
  //  Check to make sure that we have exactly the rows (elements) that an uniform
  //  distribution would give us.
  //
  vector <int> MyRows( m_loc ) ; 
  Phase2Matmap.MyGlobalElements( &MyRows[0] ) ; 
  cout << " iam = " << iam << "MyRows = " ;
  for ( int i = 0; i < m_loc ; i++) cout << MyRows[i] << " " ;
  cout << endl ; 

  cout << "NumMyElements = " <<  Phase2Matmap.NumMyElements() << endl ; 
  cout << "iam = " << iam << " My GIDs = " <<  Phase2Matmap.MinMyGID() << " ..  " 
       << Phase2Matmap.MaxMyGID() << endl ; 
  assert( MyFirstElement ==  Phase2Matmap.MinMyGID() ) ; 
  assert( MyFirstElement+m_loc-1 ==  Phase2Matmap.MaxMyGID() ) ; 


  //
  //  This is actually redundant with the Ap below
  //
  MyRowPtr[0] = 0 ; 
  int CurrentRowPtr = 0 ;
  for ( int i = 0; i < m_loc ; i++ ) { 
    CurrentRowPtr += Phase2Mat->NumMyEntries( i ) ; 
    MyRowPtr[i+1] = CurrentRowPtr ; 
  }

  //  return 0; // WORK GXX BOGUS KEN - Fails before here 

  //
  //  Extract nzval(Aval) and colind(Ai) from the CrsMatrix (Phase2Mat) 
  //

  if ( factor ) { 
  //
  //  Step 6) Convert the matrix to Ap, Ai, Aval
  //
    Ap.resize( m_loc+1 );
    Ai.resize( EPETRA_MAX( m_loc, nnz_loc) ) ; 
    Aval.resize( EPETRA_MAX( m_loc, nnz_loc) ) ; 
    
    int NzThisRow ;
    double *RowValues;
    int *ColIndices;
    int Ai_index = 0 ; 
    int MyRow;
    int num_my_cols = Phase2Mat->NumMyCols() ; 
    vector <int>Global_Columns( num_my_cols ) ; 
    for ( int i = 0 ; i < num_my_cols ; i ++ ) { 
      Global_Columns[i] = Phase2Mat->GCID( i ) ; 
    }

    for ( MyRow = 0; MyRow <m_loc; MyRow++ ) {
      int status = Phase2Mat->ExtractMyRowView( MyRow, NzThisRow, RowValues, ColIndices ) ;
      assert( status == 0 ) ; 
      Ap[MyRow] = Ai_index ; 
      assert( Ap[MyRow] == MyRowPtr[MyRow] ) ; 
      for ( int j = 0; j < NzThisRow; j++ ) { 
	Ai[Ai_index] = Global_Columns[ColIndices[j]] ; 
	Aval[Ai_index] = RowValues[j] ; 
	Ai_index++;
      }
    }
    assert( m_loc == MyRow );
    Ap[ m_loc ] = Ai_index ; 
  }
  

  //
  //  WORK GXX - DONE - trivial actually 
  //  Pull B out of the Epetra_vector 
  //

  Epetra_MultiVector   *vecX = Problem_->GetLHS() ; 
  Epetra_MultiVector   *vecB = Problem_->GetRHS() ; 

  double *bValues ;
  double *xValues ;
  int ldb, ldx ; 

  EPETRA_CHK_ERR( vecB->ExtractView( &bValues, &ldb ) )  ; 
  EPETRA_CHK_ERR( vecX->ExtractView( &xValues, &ldx ) ) ; 
  EPETRA_CHK_ERR( ! ( ldx == ldb ) ) ; 
  EPETRA_CHK_ERR( ! ( ldx == m_loc ) ) ; 

  int nrhs; 
  if ( vecX == 0 ) { 
    nrhs = 0 ;
    EPETRA_CHK_ERR( vecB != 0 ) ; 
  } else { 
    nrhs = vecX->NumVectors() ; 
    EPETRA_CHK_ERR( vecB->NumVectors() != nrhs ) ; 
  }


  //  return 0; // WORK GXX BOGUS KEN 

  //
  //  Step 7)  Call SuperLUdist
  //  

  //
  //  This really belongs somewhere else - perhaps in the constructor
  //
  const Epetra_MpiComm & comm1 = dynamic_cast<const Epetra_MpiComm &> (Comm);
  MPI_Comm MPIC = comm1.Comm() ;

  if ( factor ) { 
    numprocs = Comm.NumProc() ;                 
    nprow = SLU_NumProcRows( numprocs ) ; 
    npcol = numprocs / nprow ;
    assert ( nprow * npcol == numprocs ) ; 
    superlu_gridinit( MPIC, nprow, npcol, &grid);
  
#ifdef DEBUG
    assert( Comm_assert_equal( &Comm, numrows ) );
    assert( Comm_assert_equal( &Comm, numcols ) );
#endif
    
  } else {
    assert( numprocs == Comm.NumProc() ) ; 
  }

  /* Bail out if I do not belong in the grid. */
  if ( iam < nprow * npcol ) {
	
    if ( factor ) { 
      set_default_options(&options);
      
      if ( !(berr = doubleMalloc_dist(nrhs)) )
	EPETRA_CHK_ERR( -1 ) ; 
      
#if 0
    /* Set up the local A in NR_loc format */
    cout << " iam = " << iam << " nnz_loc = " << nnz_loc 
	   << " m_loc = " << m_loc 
	   << " MyFirstElement = " << MyFirstElement 
	   << " numrows = " << numrows 
	   << " numcols = " << numcols  << endl ; 
    cout << " Ap = " ; 
    for (int i = 0; i < m_loc+1 ; i++ ) { 
      cout << Ap[i] << " " ; 
    } ; 
    cout << endl ; 

    cout << " Ai = " ;
    for (int i = 0; i < nnz_loc ; i++ ) { 
      cout << Ai[i] << " "  ; 
    } ; 
    cout << endl ; 

    cout << " Aval = "; 
    for (int i = 0; i < nnz_loc ; i++ ) { 
      cout << Aval[i]  << " " ; 
    } ; 
    cout << endl ; 
    cout << " END OF MIDDLE PRINTOUT of A " << endl ; 
#endif
    dCreate_CompRowLoc_Matrix_dist( &A, numrows, numcols, 
				    nnz_loc, m_loc, MyFirstElement,
				    &Aval[0], &Ai[0], &Ap[0], 
				    SLU_NR_loc, SLU_D, SLU_GE );
    
#if 0
  Comm.Barrier();
    cout << " iam = " << iam << " nnz_loc = " << nnz_loc 
	   << " m_loc = " << m_loc 
	   << " MyFirstElement = " << MyFirstElement 
	   << " numrows = " << numrows 
	   << " numcols = " << numcols  << endl ; 
    cout << " iam = " << iam << " Ap = " ; 
    for (int i = 0; i < m_loc+1 ; i++ ) { 
      cout << Ap[i] << " " ; 
    } ; 
    cout << endl ; 

    cout << " iam = " << iam << " Ai = " ;
    for (int i = 0; i < nnz_loc ; i++ ) { 
      cout << Ai[i] << " "  ; 
    } ; 
    cout << endl ; 

    cout << " iam = " << iam << " Aval = " ;
    for (int i = 0; i < nnz_loc ; i++ ) { 
      cout << Aval[i]  << " " ; 
    } ; 
    cout << endl ; 
    cout << " END OF BOTTOM PRINTOUT of A " << endl ; 
#endif
#if 1
  Comm.Barrier();
  //  Print for matlab 
  //
  for (int i = 0; i < m_loc ; i++ ) { 
    for ( int j = Ap[i]; j < Ap[i+1] ; j++ ) { 
      cout << "A(" << i + MyFirstElement +1  << "," << Ai[j]+1 << " ) = " << Aval[j] << "; % iam = " << iam <<endl ; 
    } 
  }
#endif


    //      /* Create compressed column matrix for A. */
    //      dCreate_CompCol_Matrix_dist(&A, numrows, numcols, 
    //				  nnz_loc, &Aval[0], &Ai[0], 
    //				  &Ap[0], SLU_NC, SLU_D, SLU_GE);
      A_and_LU_built = true; 
#if 0
      cout << " Here is A " << endl ; 
      dPrint_CompCol_Matrix_dist( &A ); 
      cout << " That was A " << "numrows = " << numrows <<  endl ; 
      cout << "numcols = " << numcols <<  endl ; 
#endif
      
      /* Initialize ScalePermstruct and LUstruct. */
      ScalePermstructInit(numrows, numcols, &ScalePermstruct);
      LUstructInit(numrows, numcols, &LUstruct);
      
      assert( options.Fact == DOFACT );  
      options.Fact = DOFACT ;       

      Factored_ = true; 
    } else {
      assert( Factored_ == true ) ; 
      EPETRA_CHK_ERR( Factored_ == false ) ; 
      options.Fact = FACTORED ; 
    }

    //
    //  WORK GXX - We may not need to do anything over than copy the 
    //  data from b into x  - This looks good to me. 
    //  pdgssvx_ABglobal returns x in b, so we copy b into x and pass x as b.
    //
    for ( int j = 0 ; j < nrhs; j++ )
      for ( int i = 0 ; i < m_loc; i++ ) xValues[i+j*ldx] = bValues[i+j*ldx]; 

#if 0
    cout << " B = [ " ;
    for ( int i = 0 ; i < m_loc; i++ ) cout << xValues[i+0*ldx] << " " ; 
    cout << "];"  << endl ; 
#endif
    
    /* Initialize the statistics variables. */
    PStatInit(&stat);

    /* Call the linear equation solver. */
    int info ;
    //    dPrint_CompCol_Matrix_dist(&A);
    pdgssvx(&options, &A, &ScalePermstruct, &xValues[0], ldx, nrhs, &grid,
	    &LUstruct, &SOLVEstruct, berr, &stat, &info);

    //    pdgssvx_ABglobal(&options, &A, &ScalePermstruct, &xValues[0], 
    //		     ldb, nrhs, &grid, &LUstruct, berr, &stat, &info);
    EPETRA_CHK_ERR( info ) ; 

#if 0
    cout << " X = [ " ;
    for ( int i = 0 ; i < m_loc; i++ ) cout << xValues[i+0*ldx] << " " ; 
    cout << "];"  << endl ; 
#endif
    
    PStatFree(&stat);

#if 0
    cout << " Here is X: " ;
    vecX->Print(cout ) ; 

    cout << endl << " Here is B: " ;
    vecB->Print(cout) ; 
    cout << endl ; 
#endif

  }

  //
  //  NO WORK GXX - I don't think that we need to do anything to convert X 
  //

  return(0) ; 
}
