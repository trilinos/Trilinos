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
#include "supermatrix.h"      
#include "superlu_ddefs.h"
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

  if ( A_and_LU_built ) { 
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
int SuperludistOO::Solve(bool factor) { 
  //
  //  I am going to put these here until I determine that I need them in 
  //  SuperludistOO.h 
  //

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
  //  Step 2)  Coalesce the matrix onto process 0
  //
  const Epetra_MpiComm & comm1 = dynamic_cast<const Epetra_MpiComm &> (Comm);
  MPI_Comm MPIC = comm1.Comm() ;

  int IsLocal = ( Phase2Matmap.NumMyElements() == 
		  Phase2Matmap.NumGlobalElements() )?1:0;
  Comm.Broadcast( &IsLocal, 1, 0 ) ; 
#ifdef DEBUG
  assert( Comm_assert_equal( &Comm, IsLocal ) );
#endif

  Epetra_CrsMatrix *Phase3Mat = 0 ;

  int NumGlobalElements_ = Phase2Matmap.NumGlobalElements() ;
  //  Create a serial map in case we end up needing it 
  //  If it is created inside the else block below it would have to
  //  be with a call to new().
  int NumMyElements_ = 0 ;
  if (iam==0) NumMyElements_ = NumGlobalElements_;
  Epetra_Map SerialMap( NumGlobalElements_, NumMyElements_, 0, Comm );
  Epetra_CrsMatrix SerialCrsMatrixA(Copy, SerialMap, 0);


  if ( IsLocal==1 ) {
     Phase3Mat = Phase2Mat ;
  } else {

    Epetra_Export export_to_serial( Phase2Matmap, SerialMap);

    SerialCrsMatrixA.Export( *Phase2Mat, export_to_serial, Add ); 
    
    SerialCrsMatrixA.TransformToLocal() ; 
    Phase3Mat = &SerialCrsMatrixA ;

  }
  Comm.Barrier() ; 



  //
  //  Step 3)  Transpose the matrix 
  //
  const Epetra_Map &Phase3Matmap = Phase3Mat->RowMap() ; 

  numrows = Phase3Mat->NumGlobalRows();
  int numcols = Phase3Mat->NumGlobalCols();
  int numentries = Phase3Mat->NumGlobalNonzeros();

  Epetra_CrsMatrix Phase3MatTrans(Copy, Phase3Matmap, 0);
  Epetra_CrsMatrix *Phase4Mat;

  if ( GetTrans() ) { 
    Phase4Mat = Phase3Mat ; 
  } else {
    assert( CrsMatrixTranspose( Phase3Mat, &Phase3MatTrans ) == 0 ) ; 
    Phase4Mat = &Phase3MatTrans ;
  }


  //
  //  Step 4)  Replicate the matrix
  //
  int * AllIDs = new int[numrows];
  for ( int i = 0; i < numrows ; i++ ) AllIDs[i] = i ; 

  // Create a replicated map and matrix
  Epetra_Map ReplicatedMap( -1, numrows, AllIDs, 0, Comm);
  //  Epetra_LocalMap ReplicatedMap( numrows, 0, Comm);   // Compiles,runs, fails to replicate

  delete [] AllIDs;

  Epetra_Import importer( ReplicatedMap, Phase3Matmap );

  Epetra_CrsMatrix Phase5Mat(Copy, ReplicatedMap, 0);
  EPETRA_CHK_ERR( Phase5Mat.Import( *Phase4Mat, importer, Insert) );
  EPETRA_CHK_ERR( Phase5Mat.TransformToLocal() ) ; 

  assert( Phase5Mat.NumMyRows() == Phase4Mat->NumGlobalRows() ) ;

  //
  //  Step 5)  Convert vector b to a replicated vector
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

  int nArows = Phase3Mat->NumGlobalRows() ; 
  int nAcols = Phase3Mat->NumGlobalCols() ; 

#ifdef ONE_VECTOR_ONLY
  assert( vecX->NumVectors() == 1 ) ; 
  assert( vecB->NumVectors() == 1 ) ; 

  Epetra_Vector *vecXvector = dynamic_cast<Epetra_Vector*>(vecX) ; 
  Epetra_Vector *vecBvector = dynamic_cast<Epetra_Vector*>(vecB) ; 

  assert( vecXvector != 0 ) ; 
  assert( vecBvector != 0 ) ; 

  Epetra_Vector vecXreplicated( ReplicatedMap ) ; 
  Epetra_Vector vecBreplicated( ReplicatedMap ) ; 
#else
  Epetra_MultiVector *vecXvector = (vecX) ; 
  Epetra_MultiVector *vecBvector = (vecB) ; 

  Epetra_MultiVector vecXreplicated( ReplicatedMap, nrhs ) ; 
  Epetra_MultiVector vecBreplicated( ReplicatedMap, nrhs ) ; 
#endif

  Epetra_Import ImportToReplicated( ReplicatedMap, Phase2Matmap);

  vecXreplicated.Import( *vecXvector, ImportToReplicated, Insert ) ;
  vecBreplicated.Import( *vecBvector, ImportToReplicated, Insert ) ;

  assert( nArows == vecXreplicated.MyLength() ) ; 
  assert( nAcols == vecBreplicated.MyLength() ) ;

  double *bValues ;
  double *xValues ;
  int bLda, xLda ; 

  assert( vecBreplicated.ExtractView( &bValues, &bLda ) == 0 )  ; 
  assert( vecXreplicated.ExtractView( &xValues, &xLda ) == 0 ) ; 

  if ( factor ) { 
  //
  //  Step 6) Convert the matrix to Ap, Ai, Aval
  //
  Ap.resize( numrows+1 );
  Ai.resize( EPETRA_MAX( numrows, numentries) ) ; 
  Aval.resize( EPETRA_MAX( numrows, numentries) ) ; 

  int NumEntries ;
  double *RowValues;
  int *ColIndices;
  int Ai_index = 0 ; 
  int MyRow;
  for ( MyRow = 0; MyRow <numrows; MyRow++ ) {
    int status = Phase5Mat.ExtractMyRowView( MyRow, NumEntries, RowValues, ColIndices ) ;
    assert( status == 0 ) ; 
    Ap[MyRow] = Ai_index ; 
    for ( int j = 0; j < NumEntries; j++ ) { 
      Ai[Ai_index] = ColIndices[j] ; 
      Aval[Ai_index] = RowValues[j] ; 
      Ai_index++;
    }
  }
  assert( numrows == MyRow );
  Ap[ numrows ] = Ai_index ; 
  }

  //
  //  Step 7)  Call SuperLUdist
  //  

  //
  //  This really belongs somewhere else - perhaps in the constructor
  //
  if ( factor ) { 
    numprocs = Comm.NumProc() ;                 
    nprow = SLU_NumRows( numprocs ) ; 
    npcol = numprocs / nprow ;
    assert ( nprow * npcol == numprocs ) ; 
    superlu_gridinit( MPIC, nprow, npcol, &grid);
  
#ifdef DEBUG
    assert( Comm_assert_equal( &Comm, numentries ) );
    assert( Comm_assert_equal( &Comm, numrows ) );
    assert( Comm_assert_equal( &Comm, numcols ) );
#endif
    
  } else {
    assert( numprocs == Comm.NumProc() ) ; 
  }
  int ldb = numrows ; 
  /* Bail out if I do not belong in the grid. */
  if ( iam < nprow * npcol ) {
    //
    //  All processes need to have identical values of:
    //    numrows(m), numcols(n), nnz(NumEntries), 
    //    Aval(a), Ap(xa), Ai(asub)
    //    double(numentries), int(n+1), int(numentries) 

#ifdef DEBUG
    for ( int ii = 0; ii < min( numentries, 10 ) ; ii++ ) { 
      assert( Comm_assert_equal( &Comm, Aval[ii] ) ) ; 
    }
    for ( int ii = 0; ii < min( numcols+1, 10 ) ; ii++ ) { 
      assert( Comm_assert_equal( &Comm, Ai[ii] ) ); 
      assert( Comm_assert_equal( &Comm, Ap[ii] ) ); 
    }
    for ( int ii = 0; ii < min( numrows, 10 ) ; ii++ ) { 
      assert( Comm_assert_equal( &Comm, bValues[ii] ) ); 
    }
#endif
	
    if ( factor ) { 
      set_default_options(&options);
      
      if ( !(berr = doubleMalloc_dist(nrhs)) )
	EPETRA_CHK_ERR( -1 ) ; 
      
      /* Create compressed column matrix for A. */
      dCreate_CompCol_Matrix_dist(&A, numrows, numcols, 
				  numentries, &Aval[0], &Ai[0], 
				  &Ap[0], SLU_NC, SLU_D, SLU_GE);
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
    //  pdgssvx_ABglobal returns x in b, so we copy b into x and pass x as b.
    //
    for ( int j = 0 ; j < nrhs; j++ )
      for ( int i = 0 ; i < numrows; i++ ) xValues[i+j*xLda] = bValues[i+j*xLda]; 
    
    /* Initialize the statistics variables. */
    PStatInit(&stat);

    /* Call the linear equation solver. */
    int info ;
    pdgssvx_ABglobal(&options, &A, &ScalePermstruct, &xValues[0], 
		     ldb, nrhs, &grid, &LUstruct, berr, &stat, &info);
    EPETRA_CHK_ERR( info ) ; 

    PStatFree(&stat);

  }


  //
  //  Step 8)  Convert vector x back to a distributed vector
  //
  //  This is an ugly hack - it should be cleaned up someday
  //
  for (int i = 0 ; i < numrows; i++ ) { 
    int lid[1] ; 
    lid[0] = Phase2Matmap.LID( i ) ; 
    if ( lid[0] >= 0 ) { 
      for (int j = 0 ; j < nrhs; j++ ) { 
	vecXvector->ReplaceMyValue(  lid[0], j, xValues[i + ldb * j]  ) ; 
      }
    }
  }

  return(0) ; 
}
