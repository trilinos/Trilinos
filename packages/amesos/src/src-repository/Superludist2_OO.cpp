// @HEADER
// ***********************************************************************
// 
//                Amesos: Direct Sparse Solver Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#define KLOPTER
#include "Superludist2_OO.h"
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
Superludist2_OO::Superludist2_OO(const Epetra_LinearProblem &prob ) {
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
Superludist2_OO::~Superludist2_OO(void) {
  //  DeleteMemory();
  //  DeleteAzArrays();

  if ( A_and_LU_built ) { 
    //    Destroy_CompRowLoc_Matrix_dist(&A);
    SUPERLU_FREE( A.Store );
    ScalePermstructFree(&ScalePermstruct);
    Destroy_LU(numrows, &grid, &LUstruct);
    LUstructFree(&LUstruct);
    if ( options.SolveInitialized ) {
      dSolveFinalize(&options, &SOLVEstruct ) ; 
    }
    superlu_gridexit(&grid);
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
//  1)  I still need to make it possible for Superludist2_OO to accept a 
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

int Superludist2_OO::Solve(bool factor) { 
  //
  //  I am going to put these here until I determine that I need them in 
  //  Superludist2_OO.h 
  //


  bool CheckExtraction = false;    //  Set to true to force extraction for unit test

  assert( GetTrans() == false ) ; 

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

  //
  //  Old Glossary:
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
  cout << " Superludist2_OO.cpp::   traceback mode = " << Epetra_Object::GetTracebackMode() << endl ; 
  cerr << " Send this to cerr cerr cerr   traceback mode = " << Epetra_Object::GetTracebackMode() << endl ; 
#endif

  numrows = Phase2Matmap.NumGlobalElements() ; 
  assert( numrows == Phase2Mat->NumGlobalRows() ) ; 
  int numcols = Phase2Mat->NumGlobalCols() ; 
  assert( numrows == numcols ) ; 

  //
  //  Here I compute what a uniform distribution should look like
  //  so that I can compare what I think MyFirstElement should be 
  //  against what it really is.
  //
  int m_per_p = numrows / Comm.NumProc() ;
  int remainder = numrows - ( m_per_p * Comm.NumProc() );
  int MyFirstElement = iam * m_per_p + EPETRA_MIN( iam, remainder );
  int MyFirstNonElement = (iam+1) * m_per_p + EPETRA_MIN( iam+1, remainder );
  int NumExpectedElements = MyFirstNonElement - MyFirstElement ; 


  int IsLocal = ( Phase2Matmap.NumMyElements() == 
		  Phase2Matmap.NumGlobalElements() )?1:0;
  Comm.Broadcast( &IsLocal, 1, 0 ) ; 
  //
  //  Step ?)  Convert to a distributed matrix, if appropriate
  //
  Epetra_CrsMatrix *Phase3Mat = 0 ;
  Epetra_Map DistMap(  Phase2Matmap.NumGlobalElements(), NumExpectedElements, 0, Comm );
  Epetra_CrsMatrix DistCrsMatrixA(Copy, DistMap, 0);

  bool redistribute = true ;
  if ( redistribute ) {

    Epetra_Export export_to_dist( Phase2Matmap, DistMap);

    DistCrsMatrixA.Export( *Phase2Mat, export_to_dist, Add ); 
    
    DistCrsMatrixA.TransformToLocal() ; 
    Phase3Mat = &DistCrsMatrixA ;

    
  } else {
    {  
      EPETRA_CHK_ERR( ! ( Phase2Matmap.LinearMap())  ) ; // Map must be contiguously divided
      //
      //  This is another way to check that the distribution is as pdgssvx expects it
      //  (i.e. a linear map)
      //
      int Phase2NumElements = Phase2Matmap.NumMyElements() ; 
      vector <int> MyRows( Phase2NumElements ) ; 
      Phase2Matmap.MyGlobalElements( &MyRows[0] ) ; 
      for (int row = 0 ; row < Phase2NumElements ; row++ ) {
	EPETRA_CHK_ERR( MyFirstElement+row != MyRows[row] ) ;
      }
    }
    Phase3Mat = Phase2Mat ;
  }

#if 0
  assert( MyFirstElement ==  Phase3Mat->RowMap().MinMyGID() ) ; 
  assert( NumExpectedElements == Phase3Mat->RowMap().NumMyElements() ) ; 
  assert( MyFirstElement+NumExpectedElements-1 ==  Phase3Mat->RowMap().MaxMyGID() ) ; 
#endif
  //  Comm.Barrier(); 
  int MyActualFirstElement = Phase3Mat->RowMap().MinMyGID() ; 
  int NumMyElements = Phase3Mat->NumMyRows() ; 
  vector <int> MyRowPtr( NumMyElements+1 ) ;  
  //
  //  This is actually redundant with the Ap below
  //
  MyRowPtr[0] = 0 ; 
  int CurrentRowPtr = 0 ;
  for ( int i = 0; i < NumMyElements ; i++ ) { 
    CurrentRowPtr += Phase3Mat->NumMyEntries( i ) ; 
    MyRowPtr[i+1] = CurrentRowPtr ; 
  }

  //
  //  Extract nzval(Aval) and colind(Ai) from the CrsMatrix (Phase3Mat) 
  //

  int nnz_loc = Phase3Mat->NumMyNonzeros() ;
  if ( factor ) { 
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
    int num_my_cols = Phase3Mat->NumMyCols() ; 
    vector <int>Global_Columns( num_my_cols ) ; 
    for ( int i = 0 ; i < num_my_cols ; i ++ ) { 
      Global_Columns[i] = Phase3Mat->GCID( i ) ; 
    }

    for ( MyRow = 0; MyRow < NumMyElements ; MyRow++ ) {
      int status = Phase3Mat->ExtractMyRowView( MyRow, NzThisRow, RowValues, ColIndices ) ;
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
  }

  //
  //  Pull B out of the Epetra_vector 
  //

  Epetra_MultiVector   *vecX = Problem_->GetLHS() ; 
  Epetra_MultiVector   *vecB = Problem_->GetRHS() ; 

  int nrhs; 
  if ( vecX == 0 ) { 
    nrhs = 0 ;
    EPETRA_CHK_ERR( vecB != 0 ) ; 
  } else { 
    nrhs = vecX->NumVectors() ; 
    EPETRA_CHK_ERR( vecB->NumVectors() != nrhs ) ; 
  }

  Epetra_MultiVector vecXdistributed( DistMap, nrhs ) ; 
  Epetra_MultiVector vecBdistributed( DistMap, nrhs ) ; 


  double *bValues ;
  double *xValues ;
  int ldb, ldx ; 

  Epetra_MultiVector* vecXptr; 
  Epetra_MultiVector* vecBptr; 

  if ( redistribute ) {
    Epetra_Import ImportToDistributed( DistMap, Phase2Matmap);

    vecXdistributed.Import( *vecX, ImportToDistributed, Insert ) ;
    vecBdistributed.Import( *vecB, ImportToDistributed, Insert ) ;

    vecXptr = &vecXdistributed ; 
    vecBptr = &vecBdistributed ; 
  } else {
    vecXptr = vecX ; 
    vecBptr = vecB ; 
  }


  EPETRA_CHK_ERR( vecBptr->ExtractView( &bValues, &ldb ) )  ; 
  EPETRA_CHK_ERR( vecXptr->ExtractView( &xValues, &ldx ) ) ; 
  EPETRA_CHK_ERR( ! ( ldx == ldb ) ) ; 
  EPETRA_CHK_ERR( ! ( ldx == NumMyElements ) ) ; 

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
  

  } else {
    assert( numprocs == Comm.NumProc() ) ; 
  }

  /* Bail out if I do not belong in the grid. */
  if ( iam < nprow * npcol ) {
	
    if ( factor ) { 
      set_default_options(&options);
      
      dCreate_CompRowLoc_Matrix_dist( &A, numrows, numcols, 
				      nnz_loc, NumMyElements, MyActualFirstElement,
				      &Aval[0], &Ai[0], &Ap[0], 
				      SLU_NR_loc, SLU_D, SLU_GE );

      A_and_LU_built = true; 
      
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
    //  pdgssvx returns x in b, so we copy b into x.  
    //
    for ( int j = 0 ; j < nrhs; j++ )
      for ( int i = 0 ; i < NumMyElements; i++ ) xValues[i+j*ldx] = bValues[i+j*ldb]; 

    PStatInit(&stat);    /* Initialize the statistics variables. */

    int info ;
    vector<double>berr(nrhs);
    pdgssvx(&options, &A, &ScalePermstruct, &xValues[0], ldx, nrhs, &grid,
	    &LUstruct, &SOLVEstruct, &berr[0], &stat, &info);
    EPETRA_CHK_ERR( info ) ; 

    PStatFree(&stat);

  }

  if ( redistribute ) { 
    Epetra_Import ImportBackToOriginal( Phase2Matmap,DistMap);
    
    vecX->Import( *vecXptr, ImportBackToOriginal, Insert ) ;
  }


  return(0) ; 
}
