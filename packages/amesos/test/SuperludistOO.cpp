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
#include "Epetra_Operator.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_CrsMatrix.h"
#ifdef SPARSE_DIRECT_TIMINGS
#include "SparseDirectTimingVars.h"
#endif
#include "supermatrix.h"      
#include "superlu_ddefs.h"
#include "CrsMatrixTranspose.h"
#include <vector>
#include "Epetra_LinearProblemRedistor.h"


#ifdef DEBUG
#include "Comm_assert_equal.h"
#endif
/*
  Returns the largest number of rows that allows NumProcs to be used within a 
  rectangular grid with no more rows than columns.

  i.e. max i such that there exists j > i such that i*j = NumProcs
*/
int SLU_NumRows( int NumProcs ) {
#ifdef TFLOP
  //  Else, parameter 6 of DTRSV CTNLU is incorrect 
  return 1;
#else
  int i;
  int numrows ;
  for ( i = 1; i*i <= NumProcs; i++ ) 
    ;
  bool done = false ;
  for ( numrows = i-1 ; done == false ; ) {
    int NumCols = NumProcs / numrows ; 
    if ( numrows * NumCols == NumProcs ) 
      done = true; 
    else 
      numrows-- ; 
  }
  return numrows;
#endif
}

//=============================================================================
SuperludistOO::SuperludistOO(const Epetra_LinearProblem &prob ) {
  //  AllocAzArrays();

  Problem_ = &prob ; 
  A_and_LU_built = false ; 
  Transpose_ = false ; 
  Factored_ = false ; 
  FirstCallToSolve_ = true ; 
  //
  //  The following are initialized just on general principle
  //
  numprocs = -13 ; 
  nprow = -13 ; 
  npcol = -13 ; 
  M = -13 ; 
  N = -13 ; 
  nz = -13;
  ptr = 0; 
  ind = 0;
  val = 0;
  rhs = 0 ; 
  lhs = 0 ;
  Nrhs = -13;
  ldrhs = -13;
  ldlhs = -13; 
  

}

//=============================================================================
SuperludistOO::~SuperludistOO(void) {
  //  DeleteMemory();
  //  DeleteAzArrays();

  if ( A_and_LU_built ) { 
    // Destry_CompCol_Matrix seg faults because it tries to free memory 
    // that it does not own: ind, ptr and val. 
    // Destroy_CompCol_Matrix_dist(&A);
    assert( Factored_ ) ; 
    //    The following failure is belived to be in ~Epetra_LinearProblemRedistor()
    //    if (Factored_) delete redistor ;  THIS FAILS TODAY (26 Nov 2002)
    SUPERLU_FREE(A.Store);
    Destroy_LU( M, &grid, &LUstruct);
    ScalePermstructFree(&ScalePermstruct);
    LUstructFree(&LUstruct);
    SUPERLU_FREE(berr);
  }

}

//=============================================================================

//
//
int SuperludistOO::Solve(bool factor) { 

  const Epetra_Comm &Comm = Problem_->GetOperator()->Comm();

  //
  //  This really belongs somewhere else - perhaps in the constructor
  //
  const Epetra_MpiComm & comm1 = dynamic_cast<const Epetra_MpiComm &> (Comm);
  MPI_Comm MPIC = comm1.Comm() ;

  if ( FirstCallToSolve_ ) { 
    numprocs = Comm.NumProc() ;                 
    nprow = SLU_NumRows( numprocs ) ; 
    npcol = numprocs / nprow ;
    assert ( nprow * npcol == numprocs ) ; 
    superlu_gridinit( MPIC, nprow, npcol, &grid);
    FirstCallToSolve_ = false ; 
  } else {
    assert( numprocs == Comm.NumProc() ) ; 
  }

  if ( factor ) { 
    //    if (Factored_) delete redistor ;  THIS FAILS TODAY 26 Nov 2002 
    redistor = new Epetra_LinearProblemRedistor( (Epetra_LinearProblem *)Problem_, 
						 Comm.NumProc(), true);
    bool ConstructTranspose = true; 
    bool MakeDataContiguous = true;
  redistor->CreateRedistProblem(ConstructTranspose, MakeDataContiguous, redistProblem);
  } else {
    redistor->UpdateRedistRHS(Problem_->GetRHS());
  }
  redistor->ExtractHbData( M, N, nz, ptr, ind, 
			   val, Nrhs, rhs, ldrhs, 
			   lhs, ldlhs);
  int iam = Comm.MyPID() ;

  //
  //  Step 7)  Call SuperLUdist
  //  

  /* Bail out if I do not belong in the grid. */
  if ( iam < nprow * npcol ) {
    //
    //  All processes need to have identical values of:
    //    numrows(m), numcols(n), nnz(NumEntries), 
    //    Aval(a), Ap(xa), Ai(asub)
    //    double(numentries), int(n+1), int(numentries) 

#ifdef DEBUG
    for ( int ii = 0; ii < min( nz, 10 ) ; ii++ ) { 
      assert( Comm_assert_equal( &Comm, val[ii] ) ) ; 
    }
    for ( int ii = 0; ii < min( N+1, 10 ) ; ii++ ) { 
      assert( Comm_assert_equal( &Comm, ind[ii] ) ); 
      assert( Comm_assert_equal( &Comm, ptr[ii] ) ); 
    }
    for ( int ii = 0; ii < min( M, 10 ) ; ii++ ) { 
      assert( Comm_assert_equal( &Comm, val[ii] ) ); 
    }
#endif
	
    if ( factor ) { 
      set_default_options(&options);
      
      if ( !(berr = doubleMalloc_dist(Nrhs)) )
	EPETRA_CHK_ERR( -1 ) ; 
      
      /* Create compressed column matrix for A. */
      dCreate_CompCol_Matrix_dist(&A, M, N, 
				  nz, &val[0], &ind[0], 
				  &ptr[0], SLU_NC, SLU_D, SLU_GE);
      A_and_LU_built = true; 
      
#if 0
      cout << " Here is A " << endl ; 
      dPrint_CompCol_Matrix_dist( &A ); 
      cout << " That was A " << "M = " << M <<  endl ; 
      cout << "N = " << N <<  endl ; 
#endif
      
      /* Initialize ScalePermstruct and LUstruct. */
      ScalePermstructInit(M, N, &ScalePermstruct);
      LUstructInit(M, N, &LUstruct);
      
      assert( options.Fact == DOFACT );  
      Factored_ = true; 
    } else {
      assert( Factored_ == true ) ; 
      EPETRA_CHK_ERR( Factored_ == false ) ; 
      options.Fact = FACTORED ; 
    }
    //
    //  pdgssvx_ABglobal returns x in b, so we copy b into x and pass x as b.
    //
    for ( int j = 0 ; j < Nrhs; j++ )
      for ( int i = 0 ; i < M; i++ ) lhs[i+j*ldlhs] = rhs[i+j*ldrhs]; 
    
    /* Initialize the statistics variables. */
    PStatInit(&stat);

    /* Call the linear equation solver. */
    int info ;
    pdgssvx_ABglobal(&options, &A, &ScalePermstruct, &lhs[0], 
		     ldlhs, Nrhs, &grid, &LUstruct, berr, &stat, &info);
    EPETRA_CHK_ERR( info ) ; 

    PStatFree(&stat);
  }

  redistor->UpdateOriginalLHS( Problem_->GetLHS() ) ; 

  return(0) ; 
}
