//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright 2011 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER



//
//
//      valgrind usage:
//         valgrind --suppressions=Suppressions --leak-check=yes --show-reachable=yes ./VbrMatrix_test.exe
//
//      mpirun -np 2 valgrind --suppressions=Suppressions --logfile=valg.out --leak-check=yes --show-reachable=yes ./VbrMatrix_test.exe
//      The file Suppressions can be found in packages/epetra/test/VbrMatrix/Suppressions.in
//
//
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_Time.h"
#include "Epetra_Vector.h"
#include "Epetra_Flops.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_VbrRowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include <vector>
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "../epetra_test_err.h"
#include "../src/Epetra_matrix_data.h"
#include "Epetra_Version.h"

// prototypes

int checkValues( double x, double y, string message = "", bool verbose = false) {
  if (fabs((x-y)/x) > 0.01) {
    return(1);
    if (verbose) cout << "********** " << message << " check failed.********** " << endl;
  }
  else {
    if (verbose) cout << message << " check OK." << endl;
    return(0);
  }
}

int checkMultiVectors( Epetra_MultiVector & X, Epetra_MultiVector & Y, string message = "", bool verbose = false) {
  int numVectors = X.NumVectors();
  int length = Y.MyLength();
  int badvalue = 0;
  int globalbadvalue = 0;
  for (int j=0; j<numVectors; j++)
    for (int i=0; i< length; i++)
      if (checkValues(X[j][i], Y[j][i])==1) badvalue = 1;
  X.Map().Comm().MaxAll(&badvalue, &globalbadvalue, 1);

  if (verbose) {
    if (globalbadvalue==0) cout << message << " check OK." << endl;
    else cout << "********* " << message << " check failed.********** " << endl;
  }
  return(globalbadvalue);
}

int checkVbrMatrixOptimizedGraph(Epetra_Comm& comm, bool verbose);

int checkVbrRowMatrix(Epetra_RowMatrix& A, Epetra_RowMatrix & B, bool verbose);

int CompareValues(double * A, int LDA, int NumRowsA, int NumColsA, 
		  double * B, int LDB, int NumRowsB, int NumColsB);

int check(Epetra_VbrMatrix& A, 
	  int NumMyRows1, int NumGlobalRows1, int NumMyNonzeros1, int NumGlobalNonzeros1, 
	  int NumMyBlockRows1, int NumGlobalBlockRows1, int NumMyBlockNonzeros1, int NumGlobalBlockNonzeros1, 
	  int * MyGlobalElements, bool verbose);

int power_method(bool TransA, Epetra_VbrMatrix& A, 
		 Epetra_MultiVector& q,
		 Epetra_MultiVector& z, 
		 Epetra_MultiVector& resid, 
		 double * lambda, int niters, double tolerance,
		 bool verbose);

int checkMergeRedundantEntries(Epetra_Comm& comm, bool verbose);

int checkExtractMyRowCopy(Epetra_Comm& comm, bool verbose);

int checkMatvecSameVectors(Epetra_Comm& comm, bool verbose);

int checkEarlyDelete(Epetra_Comm& comm, bool verbose);

//
//  ConvertVbrToCrs is a crude but effective way to convert a Vbr matrix to a Crs matrix
//  Caller is responsible for deleting the CrsMatrix CrsOut
//
void ConvertVbrToCrs( Epetra_VbrMatrix* VbrIn, Epetra_CrsMatrix*& CrsOut ) {

  const Epetra_Map &E_Vbr_RowMap = VbrIn->RowMatrixRowMap() ; 
  const Epetra_Map &E_Vbr_ColMap = VbrIn->RowMatrixColMap() ; 
  //    const Epetra_Map &E_Vbr_RowMap = VbrIn->OperatorRangeMap() ; 
  //    const Epetra_Map &E_Vbr_ColMap = VbrIn->OperatorDomainMap() ; 
    
    CrsOut = new Epetra_CrsMatrix( Copy, E_Vbr_RowMap, E_Vbr_ColMap, 0 ); 
    //  CrsOut = new Epetra_CrsMatrix( Copy, E_Vbr_RowMap, 0 ); 
    int NumMyElements = VbrIn->RowMatrixRowMap().NumMyElements() ;
    vector<int> MyGlobalElements( NumMyElements );
    VbrIn->RowMatrixRowMap().MyGlobalElements( &MyGlobalElements[0] ) ;

    int NumMyColumns = VbrIn->RowMatrixColMap().NumMyElements() ;
    vector<int> MyGlobalColumns( NumMyColumns );
    VbrIn->RowMatrixColMap().MyGlobalElements( &MyGlobalColumns[0] ) ;

    int MaxNumIndices = VbrIn->MaxNumEntries();
    int NumIndices;
    vector<int> LocalColumnIndices(MaxNumIndices);
    vector<int> GlobalColumnIndices(MaxNumIndices);
    vector<double> MatrixValues(MaxNumIndices); 

    for( int LocalRow=0; LocalRow<NumMyElements; ++LocalRow ) {

      VbrIn->ExtractMyRowCopy( LocalRow, 
					MaxNumIndices,
					NumIndices, 
					&MatrixValues[0],
					&LocalColumnIndices[0] );
      
      for (int j = 0 ; j < NumIndices ; j++ ) 
	{ 
	  GlobalColumnIndices[j] = MyGlobalColumns[ LocalColumnIndices[j] ]  ;
	}
      
#if 0
      if ( CrsOut->InsertGlobalValues( MyGlobalElements[LocalRow], 
				      NumIndices, 
				      &MatrixValues[0],
				      &GlobalColumnIndices[0] )!=0)abort();
#else
      if ( CrsOut->InsertMyValues( LocalRow, 
				      NumIndices, 
				      &MatrixValues[0],
				      &LocalColumnIndices[0] )!= 0) abort();
#endif
      
      
    }
    CrsOut->FillComplete();
}

//
//  checkmultiply checks to make sure that AX=Y using both multiply 
//  and both Crs Matrices, multiply1
//

int checkmultiply( bool transpose, Epetra_VbrMatrix& A, Epetra_MultiVector& X, Epetra_MultiVector& Check_Y ) {

  int numerrors = 0 ; 

#if 1
  //
  //  If X and Y are Epetra_Vectors, we first test Multiply1
  //
  Epetra_Vector *vecY =  dynamic_cast<Epetra_Vector *>( &Check_Y );
  Epetra_Vector *vecX =  dynamic_cast<Epetra_Vector *>( &X );
  assert( (vecX && vecY) || (!vecX && !vecY) ) ;

  if ( vecX && vecY ) {
    double normY, NormError;
    Check_Y.NormInf( &normY ) ; 
    Epetra_Vector Y = *vecX ; 
    Y.PutScalar( -13.0 ) ; 
    Epetra_Vector Error = *vecX ; 
    A.Multiply1( transpose, *vecX, Y ) ;  

    Error.Update( 1.0, Y, -1.0, *vecY, 0.0 ) ; 
      
    Error.NormInf( &NormError ) ; 
    if ( NormError / normY > 1e-13 ) {
       numerrors++; 
       //cout << "Y = " << Y << endl;
       //cout << "vecY " << *vecY << endl;
       //cout << "Error " << Error << endl;
       //abort();
    }
    //
    //  Check x = Ax
    //
    Epetra_Vector Z = *vecX ; 

    A.Multiply1( transpose, Z, Z ) ;  
    Error.Update( 1.0, Z, -1.0, *vecY, 0.0 ) ; 
    //    Error.Update( 1.0, Y, -1.0, *vecY, 0.0 ) ; 
      
    Error.NormInf( &NormError ) ; 
    if ( NormError / normY > 1e-13 ) numerrors++; 
  }
  //
  //  Here we test Multiply 
  //
  Epetra_MultiVector Y = X ; 
  Epetra_MultiVector Error = X ; 
  A.Multiply( transpose, X, Y ) ; 

  int NumVecs = X.NumVectors() ; 
  vector<double> NormError(NumVecs);
  vector<double> Normy(NumVecs);

  Check_Y.NormInf( &Normy[0] ) ;

  Error.Update( 1.0, Y, -1.0, Check_Y, 0.0 ) ; 
  Error.NormInf( &NormError[0] ) ; 

  bool LoopError = false ; 
  for ( int ii = 0 ; ii < NumVecs ; ii++ ) {
    if( NormError[ii] / Normy[ii] > 1e-13 ) {
      LoopError = true ;
    }
  }
  if ( LoopError ) {
    numerrors++ ; 
  }

  //
  //  Check X = AX
  //

#endif
  
  return numerrors;
}
//
//  TestMatrix contains the bulk of the testing.  
//    MinSize and MaxSize control the range of the block sizes - which are chosen randomly
//    ConstructWithNumNz controls whether a Column Map is passed to the VbrMatrix contructor
//      if false, the underlying graph will not be optimized
//    ExtraBlocks, if true, causes extra blocks to be added to the matrix, further 
//      guaranteeing that the matrix and graph will not have optimal storage.  
//      ExtraBlocks is only supported for fixed block sizes (i.e. when minsize=maxsize)
//
//  If TestMatrix() is called with any of the ConsTypes that use PreviousA, the values used to 
//  build A, i.e. MinSize, MaxSize must be the same as they were on the previous call to 
//  TestMatrix (i.e. the call that built PreviousA).  Furthermore, the ConsType must be a 
//  matching ConsType.  
//  The ConsType values that cause TestMatrix to
//  use PreviousA are:  
//		RowMapColMap_VEPR, RowMapColMap_FEPR, RowMapColMap_NEPR, 
//		WithGraph, CopyConstructor
//  The matching ConsTypes are:  
//               VariableEntriesPerRow, RowMapColMap_VEPR, NoEntriesPerRow, RowMapColMap_NEPR, WithGraph
//               FixedEntriesPerRow, RowMapColMap_FEPR
//
//  TestMatrix() when called with ConsType=WithGraph, returns with PreviousA = 0 ;
//              
//
enum ConsType { VariableEntriesPerRow, FixedEntriesPerRow, NoEntriesPerRow,
		RowMapColMap_VEPR, RowMapColMap_FEPR, RowMapColMap_NEPR, 
		WithGraph, CopyConstructor } ;

int TestMatrix( Epetra_Comm& Comm, bool verbose, bool debug, 
		int NumMyElements, int MinSize, int MaxSize, 
		ConsType ConstructorType, bool ExtraBlocks, 
		bool insertlocal, 
		bool symmetric,
		Epetra_VbrMatrix** PreviousA ) {


  if (verbose) 
    cout << "MinSize         = " << MinSize << endl
	 << "MaxSize         = " << MaxSize << endl
	 << "ConstructorType = " << ConstructorType << endl
	 << "ExtraBlocks     = " << ExtraBlocks << endl
	 << "insertlocal     = " << insertlocal << endl
	 << "symmetric       = " << symmetric << endl
	 << "PreviousA       = " << PreviousA << endl;

  int ierr = 0, forierr = 0;
  int MyPID = Comm.MyPID();
  if (MyPID < 3) NumMyElements++;
  if (NumMyElements<2) NumMyElements = 2; // This value must be greater than one on each processor

  // Define pseudo-random block sizes using an Epetra Vector of random numbers
  Epetra_Map randmap(-1, NumMyElements, 0, Comm);
  Epetra_Vector randvec(randmap);
  randvec.Random(); // Fill with random numbers
  int * ElementSizeList = new int[NumMyElements];
  int SizeRange = MaxSize - MinSize + 1;
  double DSizeRange = SizeRange;
  
  const Epetra_BlockMap* rowmap = 0 ; 
  const Epetra_BlockMap* colmap = 0 ; 
  Epetra_CrsGraph* graph = 0 ; 
  if ( *PreviousA != 0 ) {
    
    rowmap = &((*PreviousA)->RowMap());
    colmap = &((*PreviousA)->ColMap());
  }

  ElementSizeList[0] = MaxSize;
  for (int i=1; i<NumMyElements-1; i++) {
    int curSize = MinSize + (int) (DSizeRange * fabs(randvec[i]) + .99);
    ElementSizeList[i] = EPETRA_MAX(MinSize, EPETRA_MIN(MaxSize, curSize));
  }
  ElementSizeList[NumMyElements-1] = MaxSize;

  

  // Construct a Map

  int *randMyGlobalElements = randmap.MyGlobalElements();

  if ( ConstructorType == RowMapColMap_VEPR || 
       ConstructorType == RowMapColMap_FEPR ||
       ConstructorType == RowMapColMap_NEPR || 
       ConstructorType == WithGraph ) {
    rowmap->ElementSizeList( ElementSizeList ) ; 
  }
    
  Epetra_BlockMap Map (-1, NumMyElements, randMyGlobalElements, ElementSizeList, 0, Comm);
  
  // Get update list and number of local elements from newly created Map
  int NumGlobalElements = Map.NumGlobalElements();
  int * MyGlobalElements = Map.MyGlobalElements();

  // Create an integer vector NumNz that is used to build the Petra Matrix.
  // NumNz[i] is the Number of OFF-DIAGONAL term for the ith global equation on this processor

  int * NumNz = new int[NumMyElements];

  // We are building a block tridiagonal matrix

  for (int i=0; i<NumMyElements; i++)
    if (MyGlobalElements[i]==0 || MyGlobalElements[i] == NumGlobalElements-1)
      NumNz[i] = 2;
    else
      NumNz[i] = 3;
  // Create a Epetra_Matrix

  bool FixedNumEntries = false ;    // If FixedNumEntries == true, we add the upper left and lower right hand corners to 
                                    // our tri-diagonal matrix so that each row has exactly three entries
  bool HaveColMap = false ;         // Matrices constructed without a column map cannot use insertmy to create the matrix
  bool FixedBlockSize = ( MinSize == MaxSize ) ; 
  bool HaveGraph = false ; 


  Epetra_VbrMatrix* A ; 

  switch( ConstructorType ) { 
  case VariableEntriesPerRow:
    A = new Epetra_VbrMatrix( Copy, Map, NumNz ) ; 
    break; 
  case FixedEntriesPerRow:
    A = new Epetra_VbrMatrix( Copy, Map, 3 ) ; 
    FixedNumEntries = true; 
    break; 
  case NoEntriesPerRow:
    A = new Epetra_VbrMatrix( Copy, Map, 0 ) ; 
    break; 
  case RowMapColMap_VEPR: 
    A = new Epetra_VbrMatrix( Copy, *rowmap, *colmap, NumNz ) ; 
    HaveColMap = true; 
    break; 
  case RowMapColMap_FEPR: 
    A = new Epetra_VbrMatrix( Copy, *rowmap, *colmap, 3 ) ; 
    FixedNumEntries = true; 
    HaveColMap = true; 
    break; 
  case RowMapColMap_NEPR: 
    A = new Epetra_VbrMatrix( Copy, *rowmap, *colmap, 0 ) ; 
    HaveColMap = true; 
    break; 
  case WithGraph:
    graph =  new Epetra_CrsGraph( (*PreviousA)->Graph() );
    A = new Epetra_VbrMatrix( Copy, *graph );
    HaveGraph = true; 
    HaveColMap = true; 
    break; 
  case CopyConstructor:
    A = new Epetra_VbrMatrix( (**PreviousA) );
    HaveColMap = true; 
    break; 
  default:
    assert(false); 
 } 

  if ( insertlocal ) assert( HaveColMap );            // you can't insert local without a column map
  if ( ExtraBlocks ) assert( FixedBlockSize );        // ExtraBlocks is only supported for fixed block sizes
  if ( ExtraBlocks ) assert( ! HaveColMap );          // Matrices constructed with a column map cannot be extended
  if ( FixedNumEntries) assert( FixedBlockSize ) ;   // Can't handle a Fixed Number of Entries and a variable block size
  if ( insertlocal && HaveGraph ) assert( ! FixedNumEntries ) ;   // HaveGraph assumes the standard matrix shape
  if ( insertlocal && HaveGraph ) assert( ! ExtraBlocks ) ;   // HaveGraph assumes the standard matrix shape


  // WORK    Insert/Replace/Suminto  MY should fail here when there is no colmap 


  EPETRA_TEST_ERR(A->IndicesAreGlobal(),ierr);
  if ( ! HaveGraph ) EPETRA_TEST_ERR(A->IndicesAreLocal(),ierr);
  
  // Use an array of Epetra_SerialDenseMatrix objects to build VBR matrix

  Epetra_SerialDenseMatrix ** BlockEntries = new Epetra_SerialDenseMatrix*[SizeRange];

  // The array of dense matrices will increase in size from MinSize to MaxSize (defined above)
  for (int kr=0; kr<SizeRange; kr++) {
    BlockEntries[kr] = new Epetra_SerialDenseMatrix[SizeRange];
    int RowDim = MinSize+kr;
    for (int kc = 0; kc<SizeRange; kc++) {
      int ColDim = MinSize+kc;
      Epetra_SerialDenseMatrix * curmat = &(BlockEntries[kr][kc]);
      curmat->Shape(RowDim,ColDim);
      for (int j=0; j < ColDim; j++)
	for (int i=0; i < RowDim; i++) {
	  BlockEntries[kr][kc][j][i] = -1.0;
	  if (i==j && kr==kc) BlockEntries[kr][kc][j][i] = 9.0;
	  else BlockEntries[kr][kc][j][i] = -1.0;

	  if ( ! symmetric )  BlockEntries[kr][kc][j][i] += ((double) j)/10000.0;
	}
    }
  }

  // Add  rows one-at-a-time

  int *Indices = new int[3];
  int *MyIndices = new int[3];
  int *ColDims = new int[3];
  int NumEntries;
  int NumMyNonzeros = 0, NumMyEquations = 0;



  for (int i=0; i<NumMyElements; i++) {
    int CurRow = MyGlobalElements[i];
    if ( HaveColMap ) { 
      assert ( i == rowmap->LID( CurRow ) ) ; 
      assert ( CurRow == rowmap->GID( i ) ) ; 
    }
    int RowDim = ElementSizeList[i]-MinSize;
    NumMyEquations += BlockEntries[RowDim][RowDim].M();

    if (CurRow==0)
      {
	Indices[0] = CurRow;
	Indices[1] = CurRow+1;
	if ( FixedNumEntries ) {
	  Indices[2] = NumGlobalElements-1;
	  ColDims[2] = 0 ; 
	  assert( ElementSizeList[i] == MinSize );   // This is actually enforced above as well 
	  NumEntries = 3;
	} else {
	  NumEntries = 2;
	}
	ColDims[0] = ElementSizeList[i] - MinSize;
	ColDims[1] = ElementSizeList[i+1] - MinSize; // Assumes linear global ordering and > 1 row/proc.
      }
    else if (CurRow == NumGlobalElements-1)
      {
	Indices[0] = CurRow-1;
	Indices[1] = CurRow;
	if ( FixedNumEntries ) {
	  Indices[2] = 0;
	  ColDims[2] = 0 ; 
	  assert( ElementSizeList[i] == MinSize );   // This is actually enforced above as well 
	  NumEntries = 3;
	} else {
	  NumEntries = 2;
	}
	ColDims[0] = ElementSizeList[i-1] - MinSize;
	ColDims[1] = ElementSizeList[i] - MinSize; // Assumes linear global ordering and > 1 row/proc.
      }
      else {
	Indices[0] = CurRow-1;
	Indices[1] = CurRow;
	Indices[2] = CurRow+1;
	NumEntries = 3;
	if (i==0) ColDims[0] = MaxSize - MinSize; // ElementSize on MyPID-1
	else ColDims[0] = ElementSizeList[i-1] - MinSize;
	ColDims[1] = ElementSizeList[i] - MinSize;
	// ElementSize on MyPID+1
	if (i==NumMyElements-1) ColDims[2] = MaxSize - MinSize;
	else ColDims[2] = ElementSizeList[i+1] - MinSize;
      }
    if ( insertlocal ) { 
      for ( int ii=0; ii < NumEntries; ii++ ) 
	MyIndices[ii] = colmap->LID( Indices[ii] ) ; 
      //      Epetra_MpiComm* MComm = dynamic_cast<Epetra_MpiComm*>( &Comm ) ;
      //  MComm->SetTracebackMode(1); // This should enable error traceback reporting
      if ( HaveGraph ) {
	EPETRA_TEST_ERR(!(A->BeginReplaceMyValues(rowmap->LID(CurRow), NumEntries, MyIndices)==0),ierr);
      } else { 
	EPETRA_TEST_ERR(!(A->BeginInsertMyValues(rowmap->LID(CurRow), NumEntries, MyIndices)==0),ierr);
      }
      //  MComm->SetTracebackMode(0); // This should shut down any error traceback reporting
    } else { 
      if ( HaveGraph ) { 
	EPETRA_TEST_ERR(!(A->BeginReplaceGlobalValues(CurRow, NumEntries, Indices)==0),ierr);
      } else { 
	//
	//  I, Ken Stanley, think the following call should return an error since it 
	//  makes no sense to insert a value with a local index in the absence of a 
	//  map indicating what that index means.  Instead, this call returns with an 
	//  error code of 0, but problems appear later.  
	//
	//  EPETRA_TEST_ERR((A->BeginInsertMyValues(CurRow, NumEntries, Indices)==0),ierr);  // Should fail
	EPETRA_TEST_ERR(!(A->BeginInsertGlobalValues(CurRow, NumEntries, Indices)==0),ierr);
      }
    }
    forierr = 0;
    for (int j=0; j < NumEntries; j++) {
      Epetra_SerialDenseMatrix * AD = &(BlockEntries[RowDim][ColDims[j]]);
      NumMyNonzeros += AD->M() * AD->N();	  
      forierr += !(A->SubmitBlockEntry(AD->A(), AD->LDA(), AD->M(), AD->N())==0);
    }
    EPETRA_TEST_ERR(forierr,ierr);

    A->EndSubmitEntries();
  }

  int NumMyBlockEntries = 3*NumMyElements;
  if ( ! FixedNumEntries ) { 
    if (A->LRID(0)>=0) NumMyBlockEntries--; // If I own first global row, then there is one less nonzero
    if (A->LRID(NumGlobalElements-1)>=0) NumMyBlockEntries--; // If I own last global row, then there is one less nonzero
  }

  if ( ExtraBlocks ) { 


    //
    //  Add a block to the matrix on each process.  
    //  The i index is chosen from among the block rows that this process owns (i.e. MyGlobalElements[0..NumMyElements-1])
    //  The j index is chosen from among all block columns 
    //  
    //  Bugs - does not support non-contiguous matrices
    //
    //  Adding more than one off diagonal block could have resulted in adding an off diagonal block
    //  twice, resulting in errors in NumMyNonzeros and NumMyBlockEntries
    //
    const int NumTries = 100; 
    Epetra_SerialDenseMatrix BlockIindex = Epetra_SerialDenseMatrix( NumTries, 1 );
    Epetra_SerialDenseMatrix BlockJindex = Epetra_SerialDenseMatrix( NumTries, 1 );

    BlockIindex.Random(); 
    BlockJindex.Random();

    BlockIindex.Scale( 1.0 * NumMyElements );
    BlockJindex.Scale( 1.0 * A->NumGlobalBlockCols() );
    bool OffDiagonalBlockAdded = false ; 
    for ( int ii=0; ii < NumTries && ! OffDiagonalBlockAdded; ii++ ) {
      int i = (int) BlockIindex[0][ii]; 
      int j = (int) BlockJindex[0][ii]; 
      if ( i < 0 ) i = - i ; 
      if ( j < 0 ) j = - j ; 
      assert( i >= 0 ) ;
      assert( j >= 0 ) ;
      assert( i < NumMyElements ) ;
      assert( j < A->NumGlobalBlockCols() ) ;
      int CurRow = MyGlobalElements[i];
      int IndicesL[1] ; 
      IndicesL[0] = j ; 
      Epetra_SerialDenseMatrix * AD = &(BlockEntries[0][0]);
      EPETRA_TEST_ERR(!(A->BeginInsertGlobalValues( CurRow, 1, IndicesL)==0),ierr);

      if ( CurRow < j-1 || CurRow > j+1 ) {
	OffDiagonalBlockAdded = true ; 
	NumMyNonzeros += AD->M() * AD->N();	  
	NumMyBlockEntries++ ; 
      }

      //  EPETRA_TEST_ERR(!(A->SubmitBlockEntry(AD->A(), AD->LDA(), AD->M(), AD->N())==0), ierr);
      EPETRA_TEST_ERR(!(A->SubmitBlockEntry(*AD)==0), ierr);
      A->EndSubmitEntries();
    }
  }

  // Finish up
  if ( ! HaveGraph && ! insertlocal ) 
    EPETRA_TEST_ERR(!(A->IndicesAreGlobal()),ierr);
  EPETRA_TEST_ERR(!(A->FillComplete()==0),ierr);
  EPETRA_TEST_ERR(!(A->IndicesAreLocal()),ierr);
  EPETRA_TEST_ERR((A->IndicesAreGlobal()),ierr);
  EPETRA_TEST_ERR(A->StorageOptimized(),ierr);

  A->OptimizeStorage();
  if ( FixedBlockSize ) {
    EPETRA_TEST_ERR(!(A->StorageOptimized()),ierr);
  } else { 
    //    EPETRA_TEST_ERR(A->StorageOptimized(),ierr);    //  Commented out until I figure out why it occasionally fails on one process
  }
  EPETRA_TEST_ERR(A->UpperTriangular(),ierr);
  EPETRA_TEST_ERR(A->LowerTriangular(),ierr);







  int NumGlobalBlockEntries ;
  int NumGlobalNonzeros, NumGlobalEquations;
  Comm.SumAll(&NumMyBlockEntries, &NumGlobalBlockEntries, 1);
  Comm.SumAll(&NumMyNonzeros, &NumGlobalNonzeros, 1);

  Comm.SumAll(&NumMyEquations, &NumGlobalEquations, 1);

  if (! ExtraBlocks ) {
    if ( FixedNumEntries ) {
      EPETRA_TEST_ERR( !( NumGlobalBlockEntries == (3*NumGlobalElements)), ierr );
    } else {
      EPETRA_TEST_ERR( !( NumGlobalBlockEntries == (3*NumGlobalElements-2)), ierr ); 
    }
  }


  EPETRA_TEST_ERR(!(check(*A, NumMyEquations, NumGlobalEquations, NumMyNonzeros, NumGlobalNonzeros, 
	       NumMyElements, NumGlobalElements, NumMyBlockEntries, NumGlobalBlockEntries, 
	       MyGlobalElements, verbose)==0),ierr);

  forierr = 0;
  if ( ! ExtraBlocks ) {
    if ( FixedNumEntries ) 
      for (int i=0; i<NumMyElements; i++) forierr += !(A->NumGlobalBlockEntries(MyGlobalElements[i])==3);
    else
      for (int i=0; i<NumMyElements; i++) forierr += !(A->NumGlobalBlockEntries(MyGlobalElements[i])==NumNz[i]);
  }
  EPETRA_TEST_ERR(forierr,ierr);
  forierr = 0;
  if ( ! ExtraBlocks ) {
    if ( FixedNumEntries ) 
      for (int i=0; i<NumMyElements; i++) forierr += !(A->NumMyBlockEntries(i)==3);
    else
      for (int i=0; i<NumMyElements; i++) forierr += !(A->NumMyBlockEntries(i)==NumNz[i]);
  }
  EPETRA_TEST_ERR(forierr,ierr);

  if (verbose) cout << "\n\nNumEntries function check OK" << endl<< endl;

  delete [] NumNz;


  // Create vectors for Power method

  Epetra_Vector q(Map);
  Epetra_Vector z(Map);
  Epetra_Vector z_initial(Map);
  Epetra_Vector resid(Map);

  
  // Fill z with random Numbers 
  z_initial.Random();

  // variable needed for iteration
  double lambda = 0.0;
  int niters = 100;
  // int niters = 200;
  double tolerance = 1.0e-3;

  /////////////////////////////////////////////////////////////////////////////////////////////////

  // Iterate
  Epetra_Flops flopcounter;
  A->SetFlopCounter(flopcounter);
  q.SetFlopCounter(*A);
  z.SetFlopCounter(*A);
  resid.SetFlopCounter(*A);
  z = z_initial;  // Start with common initial guess
  Epetra_Time timer(Comm);
  int ierr1 = power_method(false, *A, q, z, resid, &lambda, niters, tolerance, verbose);
  double elapsed_time = timer.ElapsedTime();
  double total_flops = flopcounter.Flops();
  double MFLOPs = total_flops/elapsed_time/1000000.0;

  if (verbose) cout << "\n\nTotal MFLOPs for first solve = " << MFLOPs << endl<< endl;
  if (verbose && ierr1==1) cout << "***** Power Method did not converge. *****" << endl << endl;

  /////////////////////////////////////////////////////////////////////////////////////////////////

  // Solve transpose problem

  if (verbose) cout << "\n\nUsing transpose of matrix and solving again (should give same result).\n\n"
		    << endl;
  // Iterate
  lambda = 0.0;
  z = z_initial;  // Start with common initial guess
  flopcounter.ResetFlops();
  timer.ResetStartTime();
  ierr1 = power_method(true, *A, q, z, resid, &lambda, niters, tolerance, verbose);
  elapsed_time = timer.ElapsedTime();
  total_flops = flopcounter.Flops();
  MFLOPs = total_flops/elapsed_time/1000000.0;

  if (verbose) cout << "\n\nTotal MFLOPs for transpose solve = " << MFLOPs << endl<< endl;
  if (verbose && ierr1==1) cout << "***** Power Method did not converge. *****" << endl << endl;

  /////////////////////////////////////////////////////////////////////////////////////////////////

  Epetra_CrsMatrix* OrigCrsA;
  ConvertVbrToCrs( A, OrigCrsA ) ; 

  // Increase diagonal dominance

  if (verbose) cout << "\n\nIncreasing the magnitude of first diagonal term and solving again\n\n"
		    << endl;

  double AnormInf = -13 ;
  double AnormOne = -13 ;
  AnormInf = A->NormInf( ) ; 
  AnormOne = A->NormOne( ) ; 

  EPETRA_TEST_ERR(!(A->IndicesAreLocal()),ierr);
  EPETRA_TEST_ERR((A->IndicesAreGlobal()),ierr);

  if (A->MyGlobalBlockRow(0)) {
    int numvals = A->NumGlobalBlockEntries(0);
    Epetra_SerialDenseMatrix ** Rowvals;
    int* Rowinds = new int[numvals];
    int  RowDim;
    A->ExtractGlobalBlockRowPointers(0, numvals, RowDim, numvals, Rowinds, 
				    Rowvals); // Get A[0,:]

    for (int i=0; i<numvals; i++) {
      if (Rowinds[i] == 0) {
	//	Rowvals[i]->A()[0] *= 10.0; // Multiply first diag value by 10.0
	Rowvals[i]->A()[0] += 1000.0; // Add 1000 to first diag value
      }
    }
    delete [] Rowinds;
  }
  //
  //  NormOne() and NormInf() will NOT return cached values
  //  See bug #1151
  //

  EPETRA_TEST_ERR( ! (AnormOne != A->NormOne( )), ierr ); 
  EPETRA_TEST_ERR( ! (AnormInf != A->NormInf( )), ierr );
  //
  //  On Process 0, let the class know that NormInf_ and NormOne_ are
  //  out of date.  
  // 
  if ( MyPID == 0 ) {
    EPETRA_TEST_ERR(!(A->BeginSumIntoGlobalValues( 0, 0, 0 )==0),ierr);
    EPETRA_TEST_ERR(  A->EndSubmitEntries(), ierr );
  }
  EPETRA_TEST_ERR( ! (AnormOne != A->NormOne( )), ierr ); 
  EPETRA_TEST_ERR( ! (AnormInf != A->NormInf( )), ierr ); 
  if ( MyPID == 0 ) 
  // Iterate (again)
  lambda = 0.0;
  z = z_initial;  // Start with common initial guess
  flopcounter.ResetFlops();
  timer.ResetStartTime();
  ierr1 = power_method(false, *A, q, z, resid, &lambda, niters, tolerance, verbose);
  elapsed_time = timer.ElapsedTime();
  total_flops = flopcounter.Flops();
  MFLOPs = total_flops/elapsed_time/1000000.0;

  if (verbose) cout << "\n\nTotal MFLOPs for second solve = " << MFLOPs << endl<< endl;
  if (verbose && ierr1==1) cout << "***** Power Method did not converge. *****" << endl << endl;



  /////////////////////////////////////////////////////////////////////////////////////////////////

  EPETRA_TEST_ERR(!(A->IndicesAreLocal()),ierr);
  EPETRA_TEST_ERR((A->IndicesAreGlobal()),ierr);
  // Solve transpose problem

  if (verbose) cout << "\n\nUsing transpose of matrix and solving again (should give same result).\n\n"
		    << endl;

  // Iterate (again)
  lambda = 0.0;
  z = z_initial;  // Start with common initial guess
  flopcounter.ResetFlops();
  timer.ResetStartTime();
  ierr1 = power_method(true, *A, q, z, resid, &lambda, niters, tolerance, verbose);
  elapsed_time = timer.ElapsedTime();
  total_flops = flopcounter.Flops();
  MFLOPs = total_flops/elapsed_time/1000000.0;

  if (verbose) cout << "\n\nTotal MFLOPs for tranpose of second solve = " << MFLOPs << endl<< endl;
  if (verbose && ierr1==1) cout << "***** Power Method did not converge. *****" << endl << endl;


  if (debug) Comm.Barrier();

  if (verbose) cout << "\n\n*****Comparing against CrsMatrix " << endl<< endl;

  Epetra_CrsMatrix* CrsA;
  ConvertVbrToCrs( A, CrsA ) ; 

  Epetra_Vector CrsX = Epetra_Vector( A->OperatorDomainMap(), false ) ; 
  Epetra_Vector CrsY = Epetra_Vector( A->OperatorRangeMap(), false ) ; 
  Epetra_Vector OrigCrsY = Epetra_Vector( A->OperatorRangeMap(), false ) ; 
  Epetra_Vector Y_Apply = Epetra_Vector( A->OperatorRangeMap(), false ) ; 
  Epetra_Vector x(Map);
  Epetra_Vector y(Map);
  Epetra_Vector orig_check_y(Map);
  Epetra_Vector Apply_check_y(Map);
  Epetra_Vector check_y(Map);
  Epetra_Vector check_ytranspose(Map);

  x.Random() ; 
  CrsX = x; 
  CrsY = CrsX;
  CrsA->Multiply( false, CrsX, CrsY ) ; 
  OrigCrsA->Multiply( false, CrsX, OrigCrsY ) ; 


  check_y = CrsY ; 
  orig_check_y = OrigCrsY ; 
  CrsA->Multiply( true, CrsX, CrsY ) ; 
  check_ytranspose = CrsY ; 

  EPETRA_TEST_ERR( checkmultiply( false, *A, x, check_y ), ierr ); 

  EPETRA_TEST_ERR( checkmultiply( true, *A, x, check_ytranspose ), ierr ); 

  if (! symmetric ) EPETRA_TEST_ERR( !checkmultiply( false, *A, x, check_ytranspose ), ierr );   // Just to confirm that A is not symmetric

  EPETRA_TEST_ERR(!(A->IndicesAreLocal()),ierr);
  EPETRA_TEST_ERR((A->IndicesAreGlobal()),ierr);

  if (verbose) cout << "\n\n*****Try the Apply method " << endl<< endl;


  A->Apply( CrsX, Y_Apply ) ; 
  Apply_check_y = Y_Apply ; 
  EPETRA_TEST_ERR( checkmultiply( false, *A, x, Apply_check_y ), ierr ); 

  if (verbose) cout << "\n\n*****Multiply multivectors " << endl<< endl;

  const int NumVecs = 4 ; 
  
  Epetra_Map Amap = A->OperatorDomainMap() ; 

  //  Epetra_MultiVector CrsMX = Epetra_MultiVector( A->OperatorDomainMap(), false ) ; 
  Epetra_MultiVector CrsMX( Amap, NumVecs, false ) ; 
  Epetra_MultiVector CrsMY( A->OperatorRangeMap(), NumVecs, false ) ; 
  Epetra_MultiVector mx(Map, NumVecs);
  Epetra_MultiVector my(Map, NumVecs);
  Epetra_MultiVector check_my(Map, NumVecs);
  Epetra_MultiVector check_mytranspose(Map, NumVecs);
  mx.Random(); // Fill mx with random numbers
#if 0 
  CrsMX = mx; 
  CrsA->Multiply( false, CrsMX, CrsMY ) ; 
#else
  CrsMY = mx; 
  CrsA->Multiply( false, CrsMY, CrsMY ) ; 
#endif
  check_my = CrsMY ; 
  CrsMY = mx; 
  CrsA->Multiply( true, CrsMY, CrsMY ) ; 
  check_mytranspose = CrsMY ; 


  EPETRA_TEST_ERR( checkmultiply( false, *A, mx, check_my ), ierr ); 

  EPETRA_TEST_ERR( checkmultiply( true, *A, mx, check_mytranspose ), ierr ); 

  
  EPETRA_TEST_ERR(!(A->IndicesAreLocal()),ierr);
  EPETRA_TEST_ERR((A->IndicesAreGlobal()),ierr);
  if (verbose) cout << "\n\n*****Testing copy constructor" << endl<< endl;

  //  B was changed to a pointer so that we could delete it before the
  //  underlying graph is deleted.  This was necessary before Bug #1116 was
  //  fixed, to avoid seg-faults. Bug #1116 is now fixed so this is no longer
  //  an issue.
  Epetra_VbrMatrix* B = new Epetra_VbrMatrix(*A);

  EPETRA_TEST_ERR(!(check(*B, NumMyEquations, NumGlobalEquations, NumMyNonzeros, NumGlobalNonzeros, 
	       NumMyElements, NumGlobalElements, NumMyBlockEntries, NumGlobalBlockEntries, 
	       MyGlobalElements, verbose)==0),ierr);

  EPETRA_TEST_ERR( ! ( A->StorageOptimized() == B->StorageOptimized() ), ierr ) ; 

  EPETRA_TEST_ERR( checkmultiply( false, *B, mx, check_my ), ierr ); 


  *B = *B;    // Should be harmless - check to make sure that it is

  EPETRA_TEST_ERR(!(check(*B, NumMyEquations, NumGlobalEquations, NumMyNonzeros, NumGlobalNonzeros, 
	       NumMyElements, NumGlobalElements, NumMyBlockEntries, NumGlobalBlockEntries, 
	       MyGlobalElements, verbose)==0),ierr);

  EPETRA_TEST_ERR( ! ( A->StorageOptimized() == B->StorageOptimized() ), ierr ) ; 

  EPETRA_TEST_ERR( checkmultiply( false, *B, mx, check_my ), ierr ); 

  AnormInf =  A->NormInf( );
  AnormOne =  A->NormOne( );
  EPETRA_TEST_ERR( ! (AnormOne == B->NormOne( )), ierr ); 
  EPETRA_TEST_ERR( ! (AnormInf == B->NormInf( )), ierr ); 


  Epetra_CrsMatrix* CrsB;
  ConvertVbrToCrs( B, CrsB ) ; 

  if (verbose) cout << "\n\n*****Testing PutScalar, LeftScale, RightScale, and ReplaceDiagonalValues" << endl<< endl;
  //
  //  Check PutScalar, 
  //
  B->PutScalar( 1.0 ) ; 
  CrsB->PutScalar( 1.0 ) ; 
  EPETRA_TEST_ERR(! ( B->NormOne() == CrsB->NormOne() ), ierr ) ; 


  check_y = CrsY ; 
  //
  EPETRA_TEST_ERR( B->ReplaceDiagonalValues( check_y ), ierr ) ; 
  EPETRA_TEST_ERR( CrsB->ReplaceDiagonalValues( CrsY ), ierr ) ; 
  EPETRA_TEST_ERR(! ( B->NormOne() == CrsB->NormOne() ), ierr ) ; 

  EPETRA_TEST_ERR( B->LeftScale( check_y ), ierr ) ; 
  EPETRA_TEST_ERR( CrsB->LeftScale( CrsY ), ierr ) ; 
  EPETRA_TEST_ERR(! ( B->NormOne() == CrsB->NormOne() ), ierr ) ; 

  EPETRA_TEST_ERR( B->RightScale( check_y ), ierr ) ; 
  EPETRA_TEST_ERR( CrsB->RightScale( CrsY ), ierr ) ; 
  EPETRA_TEST_ERR(! ( B->NormOne() == CrsB->NormOne() ), ierr ) ; 

  double B_norm_frob = B->NormFrobenius();
  double CrsB_norm_frob = CrsB->NormFrobenius();
  //need to use a fairly large tolerance when comparing the frobenius
  //norm from a vbr-matrix and a crs-matrix, because the point-entries
  //are visited in different orders during the norm calculation, causing
  //round-off differences to accumulate. That's why we choose 5.e-5
  //instead of a smaller number like 1.e-13 or so.
  if (fabs(B_norm_frob-CrsB_norm_frob) > 5.e-5) {
    std::cout << "fabs(B_norm_frob-CrsB_norm_frob): "
      << fabs(B_norm_frob-CrsB_norm_frob) << std::endl;
    std::cout << "VbrMatrix::NormFrobenius test FAILED."<<std::endl;
    EPETRA_TEST_ERR(-99, ierr);
  }
  if (verbose) std::cout << "\n\nVbrMatrix::NormFrobenius OK"<<std::endl<<std::endl;

  if (debug) Comm.Barrier();

  if (verbose) cout << "\n\n*****Testing post construction modifications" << endl<< endl;
  if (verbose) cout << "\n\n*****Replace methods should be OK" << endl<< endl;

  //  Check to see if we can restore the matrix to its original value
  // Does not work if ExtraBlocks is true

  EPETRA_TEST_ERR(!(A->IndicesAreLocal()),ierr);
  EPETRA_TEST_ERR((A->IndicesAreGlobal()),ierr);

  if ( ! ExtraBlocks ) 
  {
    // Add  rows one-at-a-time

    NumMyNonzeros = NumMyEquations = 0;
    
    for (int i=0; i<NumMyElements; i++) {
      int CurRow = MyGlobalElements[i];
      int RowDim = ElementSizeList[i]-MinSize;
      NumMyEquations += BlockEntries[RowDim][RowDim].M();
      
      if (CurRow==0)
	{
	  Indices[0] = CurRow;
	  Indices[1] = CurRow+1;
	  NumEntries = 2;
	  ColDims[0] = ElementSizeList[i] - MinSize;
	  ColDims[1] = ElementSizeList[i+1] - MinSize; // Assumes linear global ordering and > 1 row/proc.
	}
      else if (CurRow == NumGlobalElements-1)
	{
	  Indices[0] = CurRow-1;
	  Indices[1] = CurRow;
	  NumEntries = 2;
	  ColDims[0] = ElementSizeList[i-1] - MinSize;
	  ColDims[1] = ElementSizeList[i] - MinSize; // Assumes linear global ordering and > 1 row/proc.
	}
      else {
	Indices[0] = CurRow-1;
	Indices[1] = CurRow;
	Indices[2] = CurRow+1;
	NumEntries = 3;
	if (i==0) ColDims[0] = MaxSize - MinSize; // ElementSize on MyPID-1
	else ColDims[0] = ElementSizeList[i-1] - MinSize;
	ColDims[1] = ElementSizeList[i] - MinSize;
	// ElementSize on MyPID+1
	if (i==NumMyElements-1) ColDims[2] = MaxSize - MinSize;
	else ColDims[2] = ElementSizeList[i+1] - MinSize;
      }
      EPETRA_TEST_ERR(!(A->IndicesAreLocal()),ierr);
      EPETRA_TEST_ERR((A->IndicesAreGlobal()),ierr);
      EPETRA_TEST_ERR(!(A->BeginReplaceGlobalValues(CurRow, NumEntries, Indices)==0),ierr);
      forierr = 0;
      for (int j=0; j < NumEntries; j++) {
	Epetra_SerialDenseMatrix * AD = &(BlockEntries[RowDim][ColDims[j]]);
	NumMyNonzeros += AD->M() * AD->N();	  
	forierr += !(A->SubmitBlockEntry(AD->A(), AD->LDA(), AD->M(), AD->N())==0);
      }
      EPETRA_TEST_ERR(forierr,ierr);

      A->EndSubmitEntries();
    }
    EPETRA_TEST_ERR( checkmultiply( false, *A, x, orig_check_y ), ierr );     
  }

  //
  //  Suminto should cause the matrix to be doubled 
  //
  if ( ! ExtraBlocks )   {
    // Add  rows one-at-a-time

    NumMyNonzeros = NumMyEquations = 0;
    
    for (int i=0; i<NumMyElements; i++) {
      int CurRow = MyGlobalElements[i];
      int RowDim = ElementSizeList[i]-MinSize;
      NumMyEquations += BlockEntries[RowDim][RowDim].M();
      
      if (CurRow==0)
	{
	  Indices[0] = CurRow;
	  Indices[1] = CurRow+1;
	  NumEntries = 2;
	  ColDims[0] = ElementSizeList[i] - MinSize;
	  ColDims[1] = ElementSizeList[i+1] - MinSize; // Assumes linear global ordering and > 1 row/proc.
	}
      else if (CurRow == NumGlobalElements-1)
	{
	  Indices[0] = CurRow-1;
	  Indices[1] = CurRow;
	  NumEntries = 2;
	  ColDims[0] = ElementSizeList[i-1] - MinSize;
	  ColDims[1] = ElementSizeList[i] - MinSize; // Assumes linear global ordering and > 1 row/proc.
	}
      else {
	Indices[0] = CurRow-1;
	Indices[1] = CurRow;
	Indices[2] = CurRow+1;
	NumEntries = 3;
	if (i==0) ColDims[0] = MaxSize - MinSize; // ElementSize on MyPID-1
	else ColDims[0] = ElementSizeList[i-1] - MinSize;
	ColDims[1] = ElementSizeList[i] - MinSize;
	// ElementSize on MyPID+1
	if (i==NumMyElements-1) ColDims[2] = MaxSize - MinSize;
	else ColDims[2] = ElementSizeList[i+1] - MinSize;
      }
      if ( insertlocal ) {
	for ( int ii=0; ii < NumEntries; ii++ ) 
	  MyIndices[ii] = colmap->LID( Indices[ii] ) ; 
	EPETRA_TEST_ERR(!(A->BeginSumIntoMyValues( rowmap->LID(CurRow), NumEntries, MyIndices)==0),ierr);
      } else { 
	EPETRA_TEST_ERR(!(A->BeginSumIntoGlobalValues(CurRow, NumEntries, Indices)==0),ierr);
      }
      forierr = 0;
      for (int j=0; j < NumEntries; j++) {
	Epetra_SerialDenseMatrix * AD = &(BlockEntries[RowDim][ColDims[j]]);
	NumMyNonzeros += AD->M() * AD->N();	  
	//  This has nothing to do with insertlocal, but that is a convenient bool to key on 
	if ( insertlocal ) {
	  forierr += !(A->SubmitBlockEntry( *AD )==0);
	} else { 
	  forierr += !(A->SubmitBlockEntry(AD->A(), AD->LDA(), AD->M(), AD->N())==0);
	}
      }
      EPETRA_TEST_ERR(forierr,ierr);

      A->EndSubmitEntries();
    }
    
    orig_check_y.Scale(2.0) ;

    //  This will not work with FixedNumEntries unless we fix the fix the above loop to add the sorner elements to the tridi matrix
    if ( ! FixedNumEntries ) 
      EPETRA_TEST_ERR( checkmultiply( false, *A, x, orig_check_y ), ierr ); 
  }


  {for (int kr=0; kr<SizeRange; kr++) delete [] BlockEntries[kr];}
  delete [] BlockEntries;
  delete [] ColDims;
  delete [] MyIndices ;
  delete [] Indices;
  delete [] ElementSizeList;




  if (verbose) cout << "\n\n*****Insert methods should not be accepted" << endl<< endl;

  int One = 1;
  if (B->MyGRID(0)) EPETRA_TEST_ERR(!(B->BeginInsertGlobalValues(0, 1, &One)==-2),ierr);

  Epetra_Vector checkDiag(B->RowMap());
  forierr = 0;
  int NumMyEquations1 = B->NumMyRows();
  double two1 = 2.0;

  // Test diagonal replacement and extraction methods

  forierr = 0;
  for (int i=0; i<NumMyEquations1; i++) checkDiag[i]=two1;
  EPETRA_TEST_ERR(forierr,ierr);
  
  EPETRA_TEST_ERR(!(B->ReplaceDiagonalValues(checkDiag)==0),ierr);
  
  Epetra_Vector checkDiag1(B->RowMap());
  EPETRA_TEST_ERR(!(B->ExtractDiagonalCopy(checkDiag1)==0),ierr);
  
  forierr = 0;
  for (int i=0; i<NumMyEquations1; i++) forierr += !(checkDiag[i]==checkDiag1[i]);
  EPETRA_TEST_ERR(forierr,ierr);

  if (verbose) cout << "\n\nDiagonal extraction and replacement OK.\n\n" << endl;

  double originfnorm = B->NormInf();
  double origonenorm = B->NormOne();
  EPETRA_TEST_ERR(!(B->Scale(4.0)==0),ierr);
  //  EPETRA_TEST_ERR((B->NormOne()!=origonenorm),ierr);
  EPETRA_TEST_ERR(!(B->NormInf()==4.0 * originfnorm),ierr);
  EPETRA_TEST_ERR(!(B->NormOne()==4.0 * origonenorm),ierr);

  if (verbose) cout << "\n\nMatrix scale OK.\n\n" << endl;
  
  // Testing VbrRowMatrix adapter
  Epetra_VbrRowMatrix ARowMatrix(A);

  EPETRA_TEST_ERR(checkVbrRowMatrix(*A, ARowMatrix, verbose),ierr);

  if ( PreviousA ) 
    delete *PreviousA; 
  
  //
  //  The following code deals with the fact that A has to be delete before graph is, when 
  //  A is built with a contructor that takes a graph as input.
  //
  delete B;
  if ( graph ) {
    delete A;
    delete graph ; 
    *PreviousA = 0 ; 
  } else { 
    *PreviousA = A ; 
  }
  
  delete CrsA;
  delete CrsB;
  delete OrigCrsA ;

  return ierr; 

}


int main(int argc, char *argv[])
{
  int ierr = 0;
  bool debug = false;

#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();
  bool verbose = false;

  // Check if we should print results to standard out
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;

  // char tmp;
  //  int rank = Comm.MyPID();
  //  if (rank==0) cout << "Press any key to continue..."<< endl;
  //  if (rank==0) cin >> tmp;
  //  Comm.Barrier();

  //  if ( ! verbose )  
  Comm.SetTracebackMode(0); // This should shut down any error traceback reporting

  if (verbose && MyPID==0)
    cout << Epetra_Version() << endl << endl;

  if (verbose) cout << "Processor "<<MyPID<<" of "<< NumProc
              << " is alive."<<endl;

  // Redefine verbose to only print on PE 0
  if (verbose && Comm.MyPID()!=0) verbose = false;


  int NumMyElements = 400;
  //  int NumMyElements = 3; 
  int MinSize = 2;
  int MaxSize = 8;
  bool NoExtraBlocks = false; 
  bool symmetric = true; 
  bool NonSymmetric = false;
  bool NoInsertLocal = false ; 
  bool InsertLocal = true ; 

  Epetra_VbrMatrix* PreviousA = 0 ; 

  ierr +=TestMatrix( Comm, verbose, debug, NumMyElements, 1, 1, VariableEntriesPerRow, NoExtraBlocks, NoInsertLocal, symmetric, &PreviousA );  

  //
  //  Check the various constructors
  //  
  
  ierr +=TestMatrix( Comm, verbose, debug, NumMyElements, MinSize, MinSize, VariableEntriesPerRow, NoExtraBlocks, NoInsertLocal, symmetric, &PreviousA );

  ierr +=TestMatrix( Comm, verbose, debug, NumMyElements, MinSize, MaxSize, VariableEntriesPerRow, NoExtraBlocks, NoInsertLocal, symmetric, &PreviousA );

  ierr +=TestMatrix( Comm, verbose, debug, NumMyElements, MinSize, MinSize, FixedEntriesPerRow, NoExtraBlocks, NoInsertLocal, symmetric, &PreviousA );   

  ierr +=TestMatrix( Comm, verbose, debug, NumMyElements, MinSize, MaxSize, NoEntriesPerRow, NoExtraBlocks, NoInsertLocal, symmetric, &PreviousA ); 

  ierr +=TestMatrix( Comm, verbose, debug, NumMyElements, MinSize, MaxSize, RowMapColMap_VEPR, NoExtraBlocks, NoInsertLocal, symmetric, &PreviousA );

  ierr +=TestMatrix( Comm, verbose, debug, NumMyElements, MinSize, MaxSize, RowMapColMap_NEPR, NoExtraBlocks, NoInsertLocal, NonSymmetric, &PreviousA ); 

  ierr +=TestMatrix( Comm, verbose, debug, NumMyElements, MinSize, MaxSize, RowMapColMap_VEPR, NoExtraBlocks, InsertLocal, NonSymmetric, &PreviousA );

  ierr +=TestMatrix( Comm, verbose, debug, NumMyElements, MinSize, MaxSize, RowMapColMap_NEPR, NoExtraBlocks, InsertLocal, NonSymmetric, &PreviousA ); 

  ierr +=TestMatrix( Comm, verbose, debug, NumMyElements, MinSize, MaxSize, WithGraph, NoExtraBlocks, NoInsertLocal, NonSymmetric, &PreviousA ); 
  assert ( PreviousA == 0 ) ; 


  //
  //  Check some various options
  //
  ierr +=TestMatrix( Comm, verbose, debug, NumMyElements, MinSize, MaxSize, NoEntriesPerRow, NoExtraBlocks, NoInsertLocal, NonSymmetric, &PreviousA ); 

  ierr +=TestMatrix( Comm, verbose, debug, NumMyElements, MinSize, MaxSize, RowMapColMap_NEPR, NoExtraBlocks, InsertLocal, NonSymmetric, &PreviousA ); 

  ierr +=TestMatrix( Comm, verbose, debug, NumMyElements, MinSize, MaxSize, WithGraph, NoExtraBlocks, InsertLocal, NonSymmetric, &PreviousA ); 
  assert ( PreviousA == 0 ) ; 

  ierr +=TestMatrix( Comm, verbose, debug, NumMyElements, 4, 4, NoEntriesPerRow, NoExtraBlocks, NoInsertLocal, NonSymmetric, &PreviousA ); 

  ierr +=TestMatrix( Comm, verbose, debug, NumMyElements, 4, 4, RowMapColMap_NEPR, NoExtraBlocks, InsertLocal, NonSymmetric, &PreviousA ); 

  ierr +=TestMatrix( Comm, verbose, debug, NumMyElements, 4, 4, WithGraph, NoExtraBlocks, InsertLocal, NonSymmetric, &PreviousA ); 
  assert ( PreviousA == 0 ) ; 



  ierr +=TestMatrix( Comm, verbose, debug, NumMyElements, MinSize, MinSize, FixedEntriesPerRow, NoExtraBlocks, NoInsertLocal, NonSymmetric, &PreviousA );   

  ierr +=TestMatrix( Comm, verbose, debug, NumMyElements, MinSize, MinSize, RowMapColMap_FEPR, NoExtraBlocks, NoInsertLocal, NonSymmetric, &PreviousA );   



  ierr +=TestMatrix( Comm, verbose, debug, NumMyElements, 2, 2, NoEntriesPerRow, NoExtraBlocks, NoInsertLocal, NonSymmetric, &PreviousA ); 

  ierr +=TestMatrix( Comm, verbose, debug, NumMyElements, 3, 3, NoEntriesPerRow, NoExtraBlocks, NoInsertLocal, NonSymmetric, &PreviousA ); 

  ierr +=TestMatrix( Comm, verbose, debug, NumMyElements, 4, 4, NoEntriesPerRow, NoExtraBlocks, NoInsertLocal, NonSymmetric, &PreviousA ); 

  ierr +=TestMatrix( Comm, verbose, debug, NumMyElements, 5, 5, NoEntriesPerRow, NoExtraBlocks, NoInsertLocal, NonSymmetric, &PreviousA ); 

  ierr +=TestMatrix( Comm, verbose, debug, NumMyElements, 6, 6, NoEntriesPerRow, NoExtraBlocks, NoInsertLocal, NonSymmetric, &PreviousA ); 

  ierr +=TestMatrix( Comm, verbose, debug, NumMyElements, 7, 7, NoEntriesPerRow, NoExtraBlocks, NoInsertLocal, NonSymmetric, &PreviousA ); 

  ierr +=TestMatrix( Comm, verbose, debug, NumMyElements, 8, 8, NoEntriesPerRow, NoExtraBlocks, NoInsertLocal, NonSymmetric, &PreviousA ); 

  //  ierr +=TestMatrix( Comm, verbose, debug, NumMyElements, 2, 2, NoEntriesPerRow, NoExtraBlocks, NoInsertLocal, NonSymmetric, &PreviousA ); 


  delete PreviousA;

  /*
  if (verbose) {
    // Test ostream << operator (if verbose)
    // Construct a Map that puts 2 equations on each PE
    
    int NumMyElements1 = 2;
    int NumMyElements1 = NumMyElements1;
    int NumGlobalElements1 = NumMyElements1*NumProc;

    Epetra_Map Map1(-1, NumMyElements1, 0, Comm);
    
    // Get update list and number of local equations from newly created Map
    int * MyGlobalElements1 = new int[Map1.NumMyElements()];
    Map1.MyGlobalElements(MyGlobalElements1);
    
    // Create an integer vector NumNz that is used to build the Petra Matrix.
    // NumNz[i] is the Number of OFF-DIAGONAL term for the ith global equation on this processor
    
    int * NumNz1 = new int[NumMyElements1];
    
    // We are building a tridiagonal matrix where each row has (-1 2 -1)
    // So we need 2 off-diagonal terms (except for the first and last equation)
    
    for (int i=0; i<NumMyElements1; i++)
      if (MyGlobalElements1[i]==0 || MyGlobalElements1[i] == NumGlobalElements1-1)
	NumNz1[i] = 1;
      else
	NumNz1[i] = 2;
    
    // Create a Epetra_Matrix
    
    Epetra_VbrMatrix A1(Copy, Map1, NumNz1);
    
    // Add  rows one-at-a-time
    // Need some vectors to help
    // Off diagonal Values will always be -1
    
    
    int *Indices1 = new int[2];
    double two1 = 2.0;
    int NumEntries1;

    forierr = 0;
    for (int i=0; i<NumMyElements1; i++)
      {
	if (MyGlobalElements1[i]==0)
	  {
	    Indices1[0] = 1;
	    NumEntries1 = 1;
	  }
	else if (MyGlobalElements1[i] == NumGlobalElements1-1)
	  {
	    Indices1[0] = NumGlobalElements1-2;
	    NumEntries1 = 1;
	  }
	else
	  {
	    Indices1[0] = MyGlobalElements1[i]-1;
	    Indices1[1] = MyGlobalElements1[i]+1;
	    NumEntries1 = 2;
	  }
        forierr += !(A1.InsertGlobalValues(MyGlobalElements1[i], NumEntries1, Values1, Indices1)==0);
	forierr += !(A1.InsertGlobalValues(MyGlobalElements1[i], 1, &two1, MyGlobalElements1+i)>0); // Put in the diagonal entry
      }
    EPETRA_TEST_ERR(forierr,ierr);
    // Finish up
    EPETRA_TEST_ERR(!(A1.FillComplete()==0),ierr);
    
    if (verbose) cout << "\n\nPrint out tridiagonal matrix, each part on each processor.\n\n" << endl;
    cout << A1 << endl;
    
    // Release all objects
    delete [] NumNz1;
    delete [] Values1;
    delete [] Indices1;
    delete [] MyGlobalElements1;

  }
  */

  /* Active Issue: 5744: EPETRA_TEST_ERR( checkVbrMatrixOptimizedGraph(Comm, verbose), ierr); */

  EPETRA_TEST_ERR( checkMergeRedundantEntries(Comm, verbose), ierr);

  EPETRA_TEST_ERR( checkExtractMyRowCopy(Comm, verbose), ierr);

  EPETRA_TEST_ERR( checkMatvecSameVectors(Comm, verbose), ierr);

  EPETRA_TEST_ERR( checkEarlyDelete(Comm, verbose), ierr);

  if (verbose) {
    if (ierr==0) cout << "All VbrMatrix tests OK" << endl;
    else cout << ierr << " errors encountered." << endl;
  }

#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

/* end main
*/
return ierr ;
}

int power_method(bool TransA, Epetra_VbrMatrix& A, 
		 Epetra_MultiVector& q,
		 Epetra_MultiVector& z, 
		 Epetra_MultiVector& resid, 
		 double * lambda, int niters, double tolerance,
		 bool verbose) {  

  // variable needed for iteration
  double normz, residual;

  int ierr = 1;

  for (int iter = 0; iter < niters; iter++)
    {
      z.Norm2(&normz); // Compute 2-norm of z
      q.Scale(1.0/normz, z);
      A.Multiply(TransA, q, z); // Compute z = A*q
      q.Dot(z, lambda); // Approximate maximum eigenvaluE
      if (iter%100==0 || iter+1==niters)
	{
	  resid.Update(1.0, z, -(*lambda), q, 0.0); // Compute A*q - lambda*q
	  resid.Norm2(&residual);
	  if (verbose) cout << "Iter = " << iter << "  Lambda = " << *lambda 
			     << "  Residual of A*q - lambda*q = " << residual << endl;
	} 
      if (residual < tolerance) {
	ierr = 0;
	break;
      }
    }
  return(ierr);
}
int check(Epetra_VbrMatrix& A, 
	  int NumMyRows1, int NumGlobalRows1, int NumMyNonzeros1, int NumGlobalNonzeros1, 
	  int NumMyBlockRows1, int NumGlobalBlockRows1, int NumMyBlockNonzeros1, int NumGlobalBlockNonzeros1, 
	  int * MyGlobalElements, bool verbose) {

  int ierr = 0, forierr = 0;
  // Test query functions

  int NumMyRows = A.NumMyRows();
  if (verbose) cout << "\n\nNumber of local Rows = " << NumMyRows << endl<< endl;
  // TEMP
  //if (verbose) cout << "\n\nNumber of local Rows should = " << NumMyRows1 << endl<< endl;
  if (verbose) cout << "\n\nNumber of local Rows is = " << NumMyRows << endl<< endl;
  if (verbose) cout << "\n\nNumber of local Rows should = " << NumMyRows1 << endl<< endl;

  EPETRA_TEST_ERR(!(NumMyRows==NumMyRows1),ierr);

  int NumMyNonzeros = A.NumMyNonzeros();
  if (verbose) cout << "\n\nNumber of local Nonzero entries = " << NumMyNonzeros << endl<< endl;
  //if (verbose) cout << "                            Should  = " << NumMyNonzeros1 << endl<< endl;


  if ( NumMyNonzeros != NumMyNonzeros1 ) {
    cout << " MyPID = " << A.Comm().MyPID() 
	 << " NumMyNonzeros = " << NumMyNonzeros 
	 << " NumMyNonzeros1 = " << NumMyNonzeros1 << endl ; 
  }

  EPETRA_TEST_ERR(!(NumMyNonzeros==NumMyNonzeros1),ierr);

  int NumGlobalRows = A.NumGlobalRows();
  if (verbose) cout << "\n\nNumber of global Rows = " << NumGlobalRows << endl<< endl;

  EPETRA_TEST_ERR(!(NumGlobalRows==NumGlobalRows1),ierr);

  int NumGlobalNonzeros = A.NumGlobalNonzeros();
  if (verbose) cout << "\n\nNumber of global Nonzero entries = " << NumGlobalNonzeros << endl<< endl;

  if ( NumGlobalNonzeros != NumGlobalNonzeros1 ) {
    cout << " MyPID = " << A.Comm().MyPID() 
	 << " NumGlobalNonzeros = " << NumGlobalNonzeros 
	 << " NumGlobalNonzeros1 = " << NumGlobalNonzeros1 << endl ; 
  }
  EPETRA_TEST_ERR(!(NumGlobalNonzeros==NumGlobalNonzeros1),ierr);

  int NumMyBlockRows = A.NumMyBlockRows();
  if (verbose) cout << "\n\nNumber of local Block Rows = " << NumMyBlockRows << endl<< endl;

  EPETRA_TEST_ERR(!(NumMyBlockRows==NumMyBlockRows1),ierr);

  int NumMyBlockNonzeros = A.NumMyBlockEntries();

  if (verbose) cout << "\n\nNumber of local Nonzero Block entries = " << NumMyBlockNonzeros << endl<< endl;
  if (verbose) cout << "\n\nNumber of local Nonzero Block entries 1 = " << NumMyBlockNonzeros1 << endl<< endl;

  EPETRA_TEST_ERR(!(NumMyBlockNonzeros==NumMyBlockNonzeros1),ierr);

  int NumGlobalBlockRows = A.NumGlobalBlockRows();
  if (verbose) cout << "\n\nNumber of global Block Rows = " << NumGlobalBlockRows << endl<< endl;

  EPETRA_TEST_ERR(!(NumGlobalBlockRows==NumGlobalBlockRows1),ierr);

  int NumGlobalBlockNonzeros = A.NumGlobalBlockEntries();
  if (verbose) cout << "\n\nNumber of global Nonzero Block entries = " << NumGlobalBlockNonzeros << endl<< endl;

  EPETRA_TEST_ERR(!(NumGlobalBlockNonzeros==NumGlobalBlockNonzeros1),ierr);

  
  // Test RowMatrix interface implementations
  int RowDim, NumBlockEntries, * BlockIndices;
  Epetra_SerialDenseMatrix ** Values;
  // Get View of last block row
  A.ExtractMyBlockRowView(NumMyBlockRows-1, RowDim, NumBlockEntries,
			  BlockIndices, Values);
  int NumMyEntries1 = 0;
  {for (int i=0; i < NumBlockEntries; i++) NumMyEntries1 += Values[i]->N();}
  int NumMyEntries;
  A.NumMyRowEntries(NumMyRows-1, NumMyEntries);
  if (verbose) {
    cout << "\n\nNumber of nonzero values in last row = "
	 << NumMyEntries << endl<< endl;
  }

  EPETRA_TEST_ERR(!(NumMyEntries==NumMyEntries1),ierr);
  
  // Other binary tests

  EPETRA_TEST_ERR(A.NoDiagonal(),ierr);
  EPETRA_TEST_ERR(!(A.Filled()),ierr);
  EPETRA_TEST_ERR(!(A.MyGRID(A.RowMap().MaxMyGID())),ierr);
  EPETRA_TEST_ERR(!(A.MyGRID(A.RowMap().MinMyGID())),ierr);
  EPETRA_TEST_ERR(A.MyGRID(1+A.RowMap().MaxMyGID()),ierr);
  EPETRA_TEST_ERR(A.MyGRID(-1+A.RowMap().MinMyGID()),ierr);
  EPETRA_TEST_ERR(!(A.MyLRID(0)),ierr);
  EPETRA_TEST_ERR(!(A.MyLRID(NumMyBlockRows-1)),ierr);
  EPETRA_TEST_ERR(A.MyLRID(-1),ierr);
  EPETRA_TEST_ERR(A.MyLRID(NumMyBlockRows),ierr);

    
  int MaxNumBlockEntries = A.MaxNumBlockEntries();

  // Pointer Extraction approach

  //   Local index
  int MyPointersRowDim, MyPointersNumBlockEntries;
  int * MyPointersBlockIndices = new int[MaxNumBlockEntries];
  Epetra_SerialDenseMatrix **MyPointersValuesPointers;
  //   Global Index
  int GlobalPointersRowDim, GlobalPointersNumBlockEntries;
  int * GlobalPointersBlockIndices = new int[MaxNumBlockEntries];
  Epetra_SerialDenseMatrix **GlobalPointersValuesPointers;

  // Copy Extraction approach

  //   Local index
  int MyCopyRowDim, MyCopyNumBlockEntries;
  int * MyCopyBlockIndices = new int[MaxNumBlockEntries];
  int * MyCopyColDims = new int[MaxNumBlockEntries];
  int * MyCopyLDAs = new int[MaxNumBlockEntries];
  int MaxRowDim = A.MaxRowDim();
  int MaxColDim = A.MaxColDim();
  int MyCopySizeOfValues = MaxRowDim*MaxColDim;
  double ** MyCopyValuesPointers = new double*[MaxNumBlockEntries];
  for (int i=0; i<MaxNumBlockEntries; i++)
    MyCopyValuesPointers[i] = new double[MaxRowDim*MaxColDim];

  //   Global Index
  int GlobalCopyRowDim, GlobalCopyNumBlockEntries;
  int * GlobalCopyBlockIndices = new int[MaxNumBlockEntries];
  int * GlobalCopyColDims = new int[MaxNumBlockEntries];
  int * GlobalCopyLDAs = new int[MaxNumBlockEntries];
  
  int GlobalMaxRowDim = A.GlobalMaxRowDim();
  int GlobalMaxColDim = A.GlobalMaxColDim();
  int GlobalCopySizeOfValues = GlobalMaxRowDim*GlobalMaxColDim;
  double ** GlobalCopyValuesPointers = new double*[MaxNumBlockEntries];
  for (int i=0; i<MaxNumBlockEntries; i++)
    GlobalCopyValuesPointers[i] = new double[GlobalMaxRowDim*GlobalMaxColDim];

  // View Extraction approaches

  //   Local index (There is no global view available)
  int MyView1RowDim, MyView1NumBlockEntries;
  int * MyView1BlockIndices;
  Epetra_SerialDenseMatrix **MyView1ValuesPointers = new Epetra_SerialDenseMatrix*[MaxNumBlockEntries];

  //   Local index version 2 (There is no global view available)
  int MyView2RowDim, MyView2NumBlockEntries;
  int * MyView2BlockIndices;
  Epetra_SerialDenseMatrix **MyView2ValuesPointers;


  // For each row, test six approaches to extracting data from a given local index matrix
  forierr = 0;
  for (int i=0; i<NumMyBlockRows; i++) {
    int MyRow = i;
    int GlobalRow = A.GRID(i);
    // Get a copy of block indices in local index space, pointers to everything else
    EPETRA_TEST_ERR( A.ExtractMyBlockRowPointers(MyRow, MaxNumBlockEntries, MyPointersRowDim, 
				MyPointersNumBlockEntries, MyPointersBlockIndices,
				MyPointersValuesPointers), ierr );
    // Get a copy of block indices in local index space, pointers to everything else
    EPETRA_TEST_ERR( A.ExtractGlobalBlockRowPointers(GlobalRow, MaxNumBlockEntries, GlobalPointersRowDim, 
				    GlobalPointersNumBlockEntries, GlobalPointersBlockIndices,
				    GlobalPointersValuesPointers), ierr ) ;

    // Initiate a copy of block row in local index space.
    EPETRA_TEST_ERR( A.BeginExtractMyBlockRowCopy(MyRow, MaxNumBlockEntries, MyCopyRowDim, 
				 MyCopyNumBlockEntries, MyCopyBlockIndices,
				 MyCopyColDims), ierr);
    // Copy Values
    for (int j=0; j<MyCopyNumBlockEntries; j++) {
      EPETRA_TEST_ERR( A.ExtractEntryCopy(MyCopySizeOfValues, MyCopyValuesPointers[j], MaxRowDim, false), ierr) ;
      MyCopyLDAs[j] = MaxRowDim;
    }

    // Initiate a copy of block row in global index space.
    EPETRA_TEST_ERR( A.BeginExtractGlobalBlockRowCopy(GlobalRow, MaxNumBlockEntries, GlobalCopyRowDim, 
				    GlobalCopyNumBlockEntries, GlobalCopyBlockIndices,
				    GlobalCopyColDims), ierr ) ;
    // Copy Values
    for (int j=0; j<GlobalCopyNumBlockEntries; j++) {
      EPETRA_TEST_ERR( A.ExtractEntryCopy(GlobalCopySizeOfValues, GlobalCopyValuesPointers[j], GlobalMaxRowDim, false), ierr );
      GlobalCopyLDAs[j] = GlobalMaxRowDim;
    }

    // Initiate a view of block row in local index space (Version 1)
    EPETRA_TEST_ERR( A.BeginExtractMyBlockRowView(MyRow, MyView1RowDim, 
				 MyView1NumBlockEntries, MyView1BlockIndices), ierr ) ;
    // Set pointers to values
    for (int j=0; j<MyView1NumBlockEntries; j++) 
      EPETRA_TEST_ERR ( A.ExtractEntryView(MyView1ValuesPointers[j]), ierr ) ;


    // Extract a view of block row in local index space (version 2)
    EPETRA_TEST_ERR( A.ExtractMyBlockRowView(MyRow, MyView2RowDim, 
			    MyView2NumBlockEntries, MyView2BlockIndices,
			    MyView2ValuesPointers), ierr );

    forierr += !(MyPointersNumBlockEntries==GlobalPointersNumBlockEntries);
    forierr += !(MyPointersNumBlockEntries==MyCopyNumBlockEntries);
    forierr += !(MyPointersNumBlockEntries==GlobalCopyNumBlockEntries);
    forierr += !(MyPointersNumBlockEntries==MyView1NumBlockEntries);
    forierr += !(MyPointersNumBlockEntries==MyView2NumBlockEntries);
    for (int j=1; j<MyPointersNumBlockEntries; j++) {
      forierr += !(MyCopyBlockIndices[j-1]<MyCopyBlockIndices[j]);
      forierr += !(MyView1BlockIndices[j-1]<MyView1BlockIndices[j]);
      forierr += !(MyView2BlockIndices[j-1]<MyView2BlockIndices[j]);

      forierr += !(GlobalPointersBlockIndices[j]==A.GCID(MyPointersBlockIndices[j]));
      forierr += !(A.LCID(GlobalPointersBlockIndices[j])==MyPointersBlockIndices[j]);
      forierr += !(GlobalPointersBlockIndices[j]==GlobalCopyBlockIndices[j]);
      
      Epetra_SerialDenseMatrix* my = MyPointersValuesPointers[j];
      Epetra_SerialDenseMatrix* global = GlobalPointersValuesPointers[j];

      Epetra_SerialDenseMatrix* myview1 = MyView1ValuesPointers[j];
      Epetra_SerialDenseMatrix* myview2 = MyView2ValuesPointers[j];

      forierr += !(CompareValues(my->A(), my->LDA(), 
			   MyPointersRowDim, my->N(), 
			   global->A(), global->LDA(), 
			   GlobalPointersRowDim, global->N())==0);
      forierr += !(CompareValues(my->A(), my->LDA(), 
			   MyPointersRowDim, my->N(), 
			   MyCopyValuesPointers[j], MyCopyLDAs[j], 
			   MyCopyRowDim, MyCopyColDims[j])==0);
      forierr += !(CompareValues(my->A(), my->LDA(), 
			   MyPointersRowDim, my->N(), 
			   GlobalCopyValuesPointers[j], GlobalCopyLDAs[j], 
			   GlobalCopyRowDim, GlobalCopyColDims[j])==0);
      forierr += !(CompareValues(my->A(), my->LDA(), 
			   MyPointersRowDim, my->N(), 
			   myview1->A(), myview1->LDA(), 
			   MyView1RowDim, myview1->N())==0);
      forierr += !(CompareValues(my->A(), my->LDA(),
			   MyPointersRowDim, my->N(),
				 myview2->A(), myview2->LDA(),
			   MyView2RowDim, myview2->N())==0);
    }
  }
  EPETRA_TEST_ERR(forierr,ierr);

  // GlobalRowView should be illegal (since we have local indices)
  EPETRA_TEST_ERR(!(A.BeginExtractGlobalBlockRowView(A.GRID(0), MyView1RowDim, 
						     MyView1NumBlockEntries,
						     MyView1BlockIndices)==-2),ierr);
  
  // Extract a view of block row in local index space (version 2)
  EPETRA_TEST_ERR(!(A.ExtractGlobalBlockRowView(A.GRID(0), MyView2RowDim, 
				     MyView2NumBlockEntries, MyView2BlockIndices,
				     MyView2ValuesPointers)==-2),ierr);
  
  delete [] MyPointersBlockIndices;
  delete [] GlobalPointersBlockIndices;
  delete [] MyCopyBlockIndices;
  delete [] MyCopyColDims;
  delete [] MyCopyLDAs;
  for (int i=0; i<MaxNumBlockEntries; i++) delete [] MyCopyValuesPointers[i];
  delete [] MyCopyValuesPointers;
  delete [] GlobalCopyBlockIndices;
  delete [] GlobalCopyColDims;
  delete [] GlobalCopyLDAs;
  for (int i=0; i<MaxNumBlockEntries; i++) delete [] GlobalCopyValuesPointers[i];
  delete [] GlobalCopyValuesPointers;
  delete [] MyView1ValuesPointers;
  if (verbose) cout << "\n\nRows sorted check OK" << endl<< endl;
  
  return ierr;
}

//=============================================================================
int CompareValues(double * A, int LDA, int NumRowsA, int NumColsA, 
		  double * B, int LDB, int NumRowsB, int NumColsB) {
  
  int ierr = 0, forierr = 0;
  double * ptr1 = B;
  double * ptr2;
  
  if (NumRowsA!=NumRowsB) EPETRA_TEST_ERR(-2,ierr);
  if (NumColsA!=NumColsB) EPETRA_TEST_ERR(-3,ierr);
 

  forierr = 0;
  for (int j=0; j<NumColsA; j++) {
    ptr1 = B + j*LDB;
    ptr2 = A + j*LDA;
    for (int i=0; i<NumRowsA; i++) forierr += (*ptr1++ != *ptr2++);
  }
  EPETRA_TEST_ERR(forierr,ierr);
  return ierr;
}

int checkMergeRedundantEntries(Epetra_Comm& comm, bool verbose)
{
  int numProcs = comm.NumProc();
  int localProc = comm.MyPID();

  int myFirstRow = localProc*3;
  int myLastRow = myFirstRow+2;
  int numMyRows = myLastRow - myFirstRow + 1;
  int numGlobalRows = numProcs*numMyRows;
  int ierr;

  //We'll set up a matrix which is globally block-diagonal, i.e., on each
  //processor the list of columns == list of rows.
  //Set up a list of column-indices which is twice as long as it should be,
  //and its contents will be the list of local rows, repeated twice.
  int numCols = 2*numMyRows;
  int* myCols = new int[numCols];

  int col = myFirstRow;
  for(int i=0; i<numCols; ++i) {
    myCols[i] = col++;
    if (col > myLastRow) col = myFirstRow;
  }

  int elemSize = 2;
  int indexBase = 0;

  Epetra_BlockMap map(numGlobalRows, numMyRows,
                      elemSize, indexBase, comm);

  Epetra_VbrMatrix A(Copy, map, numCols);

  Epetra_MultiVector x(map, 1), y(map, 1);
  x.PutScalar(1.0);

  Epetra_MultiVector x3(map, 3), y3(map, 3);
  x.PutScalar(1.0);

  double* coef = new double[elemSize*elemSize];
  for(int i=0; i<elemSize*elemSize; ++i) {
    coef[i] = 0.5;
  }

  //we're going to insert each row twice, with coef values of 0.5. So after
  //FillComplete, which internally calls MergeRedundantEntries, the
  //matrix should contain 1.0 in each entry.

  for(int i=myFirstRow; i<=myLastRow; ++i) {
    EPETRA_TEST_ERR( A.BeginInsertGlobalValues(i, numCols, myCols), ierr);

    for(int j=0; j<numCols; ++j) {
      EPETRA_TEST_ERR( A.SubmitBlockEntry(coef, elemSize,
                                          elemSize, elemSize), ierr);
    }

    EPETRA_TEST_ERR( A.EndSubmitEntries(), ierr);
  }

  EPETRA_TEST_ERR( A.FillComplete(), ierr);

  delete [] coef;

  if (verbose) cout << "Multiply x"<<endl;
  EPETRA_TEST_ERR( A.Multiply(false, x, y), ierr );


  //Next we're going to extract pointers-to-block-rows and check values to make
  //sure that the internal method Epetra_VbrMatrix::mergeRedundantEntries()
  //worked correctly. 
  //At the same time, we're also going to create another VbrMatrix which will
  //be a View of the matrix we've already created. This will serve to provide
  //more test coverage of the VbrMatrix code.

  int numBlockEntries = 0;
  int RowDim;
  int** BlockIndices = new int*[numMyRows];
  Epetra_SerialDenseMatrix** Values;
  Epetra_VbrMatrix Aview(View, map, numMyRows);

  for(int i=myFirstRow; i<=myLastRow; ++i) {
    BlockIndices[i-myFirstRow] = new int[numCols];
    EPETRA_TEST_ERR( A.ExtractGlobalBlockRowPointers(i, numCols,
                                                     RowDim, numBlockEntries,
                                                     BlockIndices[i-myFirstRow],
                                                     Values), ierr);

    EPETRA_TEST_ERR( Aview.BeginInsertGlobalValues(i, numBlockEntries,
                                                   BlockIndices[i-myFirstRow]), ierr);

    if (numMyRows != numBlockEntries) return(-1);
    if (RowDim != elemSize) return(-2);
    for(int j=0; j<numBlockEntries; ++j) {
      if (Values[j]->A()[0] != 1.0) {
        cout << "Row " << i << " Values["<<j<<"][0]: "<< Values[j][0]
             << " should be 1.0" << endl;
        return(-3); //comment-out this return to de-activate this test
      }

      EPETRA_TEST_ERR( Aview.SubmitBlockEntry(Values[j]->A(),
                                              Values[j]->LDA(),
                                              Values[j]->M(),
                                              Values[j]->N()), ierr);
    }

    EPETRA_TEST_ERR( Aview.EndSubmitEntries(), ierr);
  }

  EPETRA_TEST_ERR( Aview.FillComplete(), ierr);

  //So the test appears to have passed for the original matrix A. Now check the
  //values of our second "view" of the matrix, 'Aview'.

  for(int i=myFirstRow; i<=myLastRow; ++i) {
    EPETRA_TEST_ERR( Aview.ExtractGlobalBlockRowPointers(i, numMyRows,
                                                         RowDim, numBlockEntries,
                                                         BlockIndices[i-myFirstRow],
                                                         Values), ierr);

    if (numMyRows != numBlockEntries) return(-1);
    if (RowDim != elemSize) return(-2);
    for(int j=0; j<numBlockEntries; ++j) {
      if (Values[j]->A()[0] != 1.0) {
        cout << "Aview: Row " << i << " Values["<<j<<"][0]: "<< Values[j][0]
             << " should be 1.0" << endl;
        return(-3); //comment-out this return to de-activate this test
      }
    }

    delete [] BlockIndices[i-myFirstRow];
  }


  if (verbose&&localProc==0) {
    cout << "checkMergeRedundantEntries, A:" << endl;
  }


  delete [] BlockIndices;
  delete [] myCols;

  return(0);
}

int checkExtractMyRowCopy(Epetra_Comm& comm, bool verbose)
{
  int numProcs = comm.NumProc();
  int localProc = comm.MyPID();

  int myFirstRow = localProc*3;
  int myLastRow = myFirstRow+2;
  int numMyRows = myLastRow - myFirstRow + 1;
  int numGlobalRows = numProcs*numMyRows;
  int ierr;

  int numCols = numMyRows;
  int* myCols = new int[numCols];

  int col = myFirstRow;
  for(int i=0; i<numCols; ++i) {
    myCols[i] = col++;
    if (col > myLastRow) col = myFirstRow;
  }

  int elemSize = 2;
  int indexBase = 0;

  Epetra_BlockMap map(numGlobalRows, numMyRows,
		      elemSize, indexBase, comm);

  Epetra_VbrMatrix A(Copy, map, numCols);

  double* coef = new double[elemSize*elemSize];

  for(int i=myFirstRow; i<=myLastRow; ++i) {
    int myPointRow = i*elemSize;

    //The coefficients need to be laid out in column-major order. i.e., the
    //coefficients in a column occur contiguously.
    for(int ii=0; ii<elemSize; ++ii) {
      for(int jj=0; jj<elemSize; ++jj) {
	double val = (myPointRow+ii)*1.0;
	coef[ii+elemSize*jj] = val;
      }
    }

    EPETRA_TEST_ERR( A.BeginInsertGlobalValues(i, numCols, myCols), ierr);

    for(int j=0; j<numCols; ++j) {
      EPETRA_TEST_ERR( A.SubmitBlockEntry(coef, elemSize,
					  elemSize, elemSize), ierr);
    }

    EPETRA_TEST_ERR( A.EndSubmitEntries(), ierr);
  }

  EPETRA_TEST_ERR( A.FillComplete(), ierr);

  delete [] coef;
  delete [] myCols;

  Epetra_SerialDenseMatrix** blockEntries;
  int len = elemSize*numCols, checkLen;
  double* values = new double[len];
  int* indices = new int[len];
  int RowDim, numBlockEntries;

  for(int i=myFirstRow; i<=myLastRow; ++i) {
    EPETRA_TEST_ERR( A.ExtractGlobalBlockRowPointers(i, numMyRows,
						     RowDim, numBlockEntries,
						     indices,
						     blockEntries), ierr);
    if (numMyRows != numBlockEntries) return(-1);
    if (RowDim != elemSize) return(-2);

    int myPointRow = i*elemSize - myFirstRow*elemSize;
    for(int ii=0; ii<elemSize; ++ii) {
      EPETRA_TEST_ERR( A.ExtractMyRowCopy(myPointRow+ii, len,
					  checkLen, values, indices), ierr);
      if (len != checkLen) return(-3);

      double val = (i*elemSize+ii)*1.0;
      double blockvalue = blockEntries[0]->A()[ii];

      for(int jj=0; jj<len; ++jj) {
	if (values[jj] != val) return(-4);
	if (values[jj] != blockvalue) return(-5);
      }
    }
  }

  delete [] values;
  delete [] indices;

  return(0);
}

int checkMatvecSameVectors(Epetra_Comm& comm, bool verbose)
{
  int numProcs = comm.NumProc();
  int localProc = comm.MyPID();

  int myFirstRow = localProc*3;
  int myLastRow = myFirstRow+2;
  int numMyRows = myLastRow - myFirstRow + 1;
  int numGlobalRows = numProcs*numMyRows;
  int ierr;

  int elemSize = 2;
  int num_off_diagonals = 1;

  epetra_test::matrix_data matdata(numGlobalRows, numGlobalRows,
				   num_off_diagonals, elemSize);

  Epetra_BlockMap map(numGlobalRows, numMyRows, elemSize, 0, comm);

  Epetra_VbrMatrix A(Copy, map, num_off_diagonals*2+1);

  int* rowlengths = matdata.rowlengths();
  int** colindices = matdata.colindices();

  for(int i=myFirstRow; i<=myLastRow; ++i) {

    EPETRA_TEST_ERR( A.BeginInsertGlobalValues(i, rowlengths[i],
                                               colindices[i]), ierr);

    for(int j=0; j<rowlengths[i]; ++j) {
      EPETRA_TEST_ERR( A.SubmitBlockEntry(matdata.coefs(i,colindices[i][j]), elemSize,
                                          elemSize, elemSize), ierr);
    }

    EPETRA_TEST_ERR( A.EndSubmitEntries(), ierr);
  }

  EPETRA_TEST_ERR( A.FillComplete(), ierr);

  Epetra_Vector x(map), y(map);

  x.PutScalar(1.0);

  A.Multiply(false, x, y);
  A.Multiply(false, x, x);

  double* xptr = x.Values();
  double* yptr = y.Values();

  for(int i=0; i<numMyRows; ++i) {
    if (xptr[i] != yptr[i]) {
      return(-1);
    }
  }

  return(0);
}

int checkEarlyDelete(Epetra_Comm& comm, bool verbose)
{
  int localProc = comm.MyPID();
  int numProcs = comm.NumProc();
  int myFirstRow = localProc*3;
  int myLastRow = myFirstRow+2;
  int numMyRows = myLastRow - myFirstRow + 1;
  int numGlobalRows = numProcs*numMyRows;
  int ierr;

  int elemSize = 2;
  int num_off_diagonals = 1;

  epetra_test::matrix_data matdata(numGlobalRows, numGlobalRows,
                                   num_off_diagonals, elemSize);

  Epetra_BlockMap map(numGlobalRows, numMyRows, elemSize, 0, comm);

  Epetra_VbrMatrix* A = new Epetra_VbrMatrix(Copy, map, num_off_diagonals*2+1);

  int* rowlengths = matdata.rowlengths();
  int** colindices = matdata.colindices();

  for(int i=myFirstRow; i<=myLastRow; ++i) {

    EPETRA_TEST_ERR( A->BeginInsertGlobalValues(i, rowlengths[i],
                                               colindices[i]), ierr);

    for(int j=0; j<rowlengths[i]; ++j) {
      EPETRA_TEST_ERR( A->SubmitBlockEntry(matdata.coefs(i,colindices[i][j]),
                                           elemSize, elemSize, elemSize), ierr);
    }

    EPETRA_TEST_ERR( A->EndSubmitEntries(), ierr);
  }

  //A call to BeginReplaceMyValues should produce an error at this
  //point, since IndicesAreLocal should be false.
  int errcode = A->BeginReplaceMyValues(0, 0, 0);
  if (errcode == 0) EPETRA_TEST_ERR(-1, ierr);

  EPETRA_TEST_ERR( A->FillComplete(), ierr);

  Epetra_VbrMatrix B(Copy, A->Graph());

  delete A;

  for(int i=myFirstRow; i<=myLastRow; ++i) {

    EPETRA_TEST_ERR( B.BeginReplaceGlobalValues(i, rowlengths[i],
                                               colindices[i]), ierr);

    for(int j=0; j<rowlengths[i]; ++j) {
      EPETRA_TEST_ERR( B.SubmitBlockEntry(matdata.coefs(i,colindices[i][j]),
                                           elemSize, elemSize, elemSize), ierr);
    }

    EPETRA_TEST_ERR( B.EndSubmitEntries(), ierr);
  }

  EPETRA_TEST_ERR( B.FillComplete(), ierr);

  return(0);
}

int checkVbrMatrixOptimizedGraph(Epetra_Comm& comm, bool verbose)
{
  using std::vector;

  int ierr = 0;

  const int node = comm.MyPID();
  const int nodes = comm.NumProc();

  int Ni    = 4;
  int Nj    = 4;
  int Gi    = 4;
  // int Gj    = 4;
  int i_off = 0;
  int j_off = 0;

  int first = 0;
  int last  = 3;

  Epetra_BlockMap* map;
  if (nodes == 1)
  {
    map = new Epetra_BlockMap(-1, 16, 3, 0, comm);
  }
  else if (nodes == 2)
  {
    Ni = 2;
    vector<int> l2g(8);
    if (node == 1) 
    {
      i_off = 2;
    }
    for (int j = 0; j < 4; ++j)
    {
      for (int i = 0; i < 2; ++i)
      {
        l2g[i + 2 * j] = (i + i_off) + (j + j_off) * Gi;
      }
    }
    map = new Epetra_BlockMap(-1, Ni*Nj, &l2g[0], 3, 0, comm); 
  }
  else if (nodes == 4)
  {
    Ni = 2;
    Nj = 2;
    vector<int> l2g(4);
    if (node == 1) 
    {
      i_off = 2;
    }
    else if (node == 2)
    {
      j_off = 2;
    }
    else if (node == 3)
    {
      i_off = 2;
      j_off = 2;
    }
    for (int j = 0; j < 2; ++j)
    {
      for (int i = 0; i < 2; ++i)
      {
        l2g[i + 2 * j] = (i + i_off) + (j + j_off) * Gi;
      }
    }
    map = new Epetra_BlockMap(-1, Ni*Nj, &l2g[0], 3, 0, comm);
  }
  else {
    map = NULL;
    return 0;
  }

  // graph
  Epetra_CrsGraph *graph = new Epetra_CrsGraph(Copy, *map, 5);
  int indices[5];

  for (int j = 0; j < Nj; ++j)
  {
    for (int i = 0; i < Ni; ++i)
    {
      int ctr = 0;
      int gi  = i + i_off;
      int gj  = j + j_off;
      indices[ctr++] = gi + gj * Gi;
      if (gi > first)
        indices[ctr++] = (gi - 1) + gj * Gi;
      if (gi < last)
        indices[ctr++] = (gi + 1) + gj * Gi;
      if (gj > first)
        indices[ctr++] = gi + (gj - 1) * Gi;
      if (gj < last)
        indices[ctr++] = gi + (gj + 1) * Gi;
      EPETRA_TEST_ERR( ! (ctr <= 5), ierr );
      // assign the indices to the graph
      graph->InsertGlobalIndices(indices[0], ctr, &indices[0]);
    }
  }

  // complete the graph
  int result = graph->FillComplete();
  EPETRA_TEST_ERR( ! result == 0,               ierr );
  EPETRA_TEST_ERR( ! graph->Filled(),           ierr );
  EPETRA_TEST_ERR(   graph->StorageOptimized(), ierr );
  EPETRA_TEST_ERR( ! graph->IndicesAreLocal(),  ierr );
  result = graph->OptimizeStorage(); 
  EPETRA_TEST_ERR( ! result == 0,               ierr );
  EPETRA_TEST_ERR( ! graph->StorageOptimized(), ierr );

  EPETRA_TEST_ERR( ! Ni*Nj   == graph->NumMyBlockRows(), ierr );
  EPETRA_TEST_ERR( ! Ni*Nj*3 == graph->NumMyRows(),      ierr );

  Epetra_VbrMatrix *matrix = new Epetra_VbrMatrix(Copy, *graph);
  EPETRA_TEST_ERR( ! matrix->IndicesAreLocal(),          ierr );

  Epetra_SerialDenseMatrix C(3, 3); 
  Epetra_SerialDenseMatrix L(3, 3); 
  Epetra_SerialDenseMatrix R(3, 3); 
  EPETRA_TEST_ERR( ! 3 == C.LDA(), ierr );
  EPETRA_TEST_ERR( ! 3 == L.LDA(), ierr );
  EPETRA_TEST_ERR( ! 3 == R.LDA(), ierr );
  std::fill(C.A(), C.A()+9, -4.0); 
  std::fill(L.A(), L.A()+9,  2.0); 
  std::fill(R.A(), R.A()+9,  2.0); 

  // fill matrix
  {
    for (int j = 0; j < Nj; ++j)
    {
      for (int i = 0; i < Ni; ++i)
      {
        int ctr = 0;
        int gi  = i + i_off;
        int gj  = j + j_off;

        int local  = i + j * Ni;
        int global = gi + gj * Gi;

        int left   = (gi - 1) + gj * Gi;
        int right  = (gi + 1) + gj * Gi;
        int bottom = gi + (gj - 1) * Gi;
        int top    = gi + (gj + 1) * Gi;

        EPETRA_TEST_ERR( ! local == matrix->LCID(global), ierr );

        indices[ctr++] = local;
        if (gi > first)
          indices[ctr++] = left;
        if (gi < last)
          indices[ctr++] = right;;
        if (gj > first)
          indices[ctr++] = bottom;
        if (gj < last)
          indices[ctr++] = top;

        matrix->BeginReplaceMyValues(local, ctr, &indices[0]);
        matrix->SubmitBlockEntry(C);
        if (gi > first) matrix->SubmitBlockEntry(L);
        if (gi < last)  matrix->SubmitBlockEntry(R);
        if (gj > first) matrix->SubmitBlockEntry(L);
        if (gj < last)  matrix->SubmitBlockEntry(R);
        matrix->EndSubmitEntries();
      }
    }
  }
  matrix->FillComplete();
  EPETRA_TEST_ERR( ! matrix->Filled(),               ierr );
  EPETRA_TEST_ERR( ! matrix->StorageOptimized(),     ierr );
  EPETRA_TEST_ERR( ! matrix->IndicesAreContiguous(), ierr );
  // EPETRA_TEST_ERR( matrix->StorageOptimized(),     ierr );
  // EPETRA_TEST_ERR( matrix->IndicesAreContiguous(), ierr );
  
  delete matrix;
  delete graph;
  delete map;

  return(0);
}


int checkVbrRowMatrix(Epetra_RowMatrix& A, Epetra_RowMatrix & B, bool verbose)  {  

  if (verbose) cout << "Checking VbrRowMatrix Adapter..." << endl;
  int ierr = 0;
  EPETRA_TEST_ERR(!A.Comm().NumProc()==B.Comm().NumProc(),ierr);
  EPETRA_TEST_ERR(!A.Comm().MyPID()==B.Comm().MyPID(),ierr);
  EPETRA_TEST_ERR(!A.Filled()==B.Filled(),ierr);
  EPETRA_TEST_ERR(!A.HasNormInf()==B.HasNormInf(),ierr);
  // EPETRA_TEST_ERR(!A.LowerTriangular()==B.LowerTriangular(),ierr);
  // EPETRA_TEST_ERR(!A.Map().SameAs(B.Map()),ierr);
  EPETRA_TEST_ERR(!A.MaxNumEntries()==B.MaxNumEntries(),ierr);
  EPETRA_TEST_ERR(!A.NumGlobalCols()==B.NumGlobalCols(),ierr);
  EPETRA_TEST_ERR(!A.NumGlobalDiagonals()==B.NumGlobalDiagonals(),ierr);
  EPETRA_TEST_ERR(!A.NumGlobalNonzeros()==B.NumGlobalNonzeros(),ierr);
  EPETRA_TEST_ERR(!A.NumGlobalRows()==B.NumGlobalRows(),ierr);
  EPETRA_TEST_ERR(!A.NumMyCols()==B.NumMyCols(),ierr);
  EPETRA_TEST_ERR(!A.NumMyDiagonals()==B.NumMyDiagonals(),ierr);
  EPETRA_TEST_ERR(!A.NumMyNonzeros()==B.NumMyNonzeros(),ierr);
  for (int i=0; i<A.NumMyRows(); i++) {
    int nA, nB;
    A.NumMyRowEntries(i,nA); B.NumMyRowEntries(i,nB);
    EPETRA_TEST_ERR(!nA==nB,ierr);
  }
  EPETRA_TEST_ERR(!A.NumMyRows()==B.NumMyRows(),ierr);
  EPETRA_TEST_ERR(!A.OperatorDomainMap().SameAs(B.OperatorDomainMap()),ierr);
  EPETRA_TEST_ERR(!A.OperatorRangeMap().SameAs(B.OperatorRangeMap()),ierr);
  EPETRA_TEST_ERR(!A.RowMatrixColMap().SameAs(B.RowMatrixColMap()),ierr);
  EPETRA_TEST_ERR(!A.RowMatrixRowMap().SameAs(B.RowMatrixRowMap()),ierr);
  // EPETRA_TEST_ERR(!A.UpperTriangular()==B.UpperTriangular(),ierr);
  EPETRA_TEST_ERR(!A.UseTranspose()==B.UseTranspose(),ierr);

  int NumVectors = 5;
  { // No transpose case
    Epetra_MultiVector X(A.OperatorDomainMap(), NumVectors);
    Epetra_MultiVector YA1(A.OperatorRangeMap(), NumVectors);
    Epetra_MultiVector YA2(YA1);
    Epetra_MultiVector YB1(YA1);
    Epetra_MultiVector YB2(YA1);
    X.Random();
    
    bool transA = false;
    A.SetUseTranspose(transA);
    B.SetUseTranspose(transA);
    A.Apply(X,YA1);
    A.Multiply(transA, X, YA2);
    EPETRA_TEST_ERR(checkMultiVectors(YA1,YA2,"A Multiply and A Apply", verbose),ierr);
    B.Apply(X,YB1);
    EPETRA_TEST_ERR(checkMultiVectors(YA1,YB1,"A Multiply and B Multiply", verbose),ierr);
    B.Multiply(transA, X, YB2);
    EPETRA_TEST_ERR(checkMultiVectors(YA1,YB2,"A Multiply and B Apply", verbose), ierr);
    
  }
  {// transpose case
    Epetra_MultiVector X(A.OperatorRangeMap(), NumVectors);
    Epetra_MultiVector YA1(A.OperatorDomainMap(), NumVectors);
    Epetra_MultiVector YA2(YA1);
    Epetra_MultiVector YB1(YA1);
    Epetra_MultiVector YB2(YA1);
    X.Random();
    
    bool transA = true;
    A.SetUseTranspose(transA);
    B.SetUseTranspose(transA);
    A.Apply(X,YA1);
    A.Multiply(transA, X, YA2);
    EPETRA_TEST_ERR(checkMultiVectors(YA1,YA2, "A Multiply and A Apply (transpose)", verbose),ierr);
    B.Apply(X,YB1);
    EPETRA_TEST_ERR(checkMultiVectors(YA1,YB1, "A Multiply and B Multiply (transpose)", verbose),ierr);
    B.Multiply(transA, X,YB2);
    EPETRA_TEST_ERR(checkMultiVectors(YA1,YB2, "A Multiply and B Apply (transpose)", verbose),ierr);
    
  }

  /*
    Epetra_Vector diagA(A.RowMatrixRowMap());
    EPETRA_TEST_ERR(A.ExtractDiagonalCopy(diagA),ierr);
    Epetra_Vector diagB(B.RowMatrixRowMap());
    EPETRA_TEST_ERR(B.ExtractDiagonalCopy(diagB),ierr);
    EPETRA_TEST_ERR(checkMultiVectors(diagA,diagB, "ExtractDiagonalCopy", verbose),ierr);
  */
  Epetra_Vector rowA(A.RowMatrixRowMap());
  EPETRA_TEST_ERR(A.InvRowSums(rowA),ierr);
  Epetra_Vector rowB(B.RowMatrixRowMap());
  EPETRA_TEST_ERR(B.InvRowSums(rowB),ierr)
  EPETRA_TEST_ERR(checkMultiVectors(rowA,rowB, "InvRowSums", verbose),ierr);

  Epetra_Vector colA(A.OperatorDomainMap());
  EPETRA_TEST_ERR(A.InvColSums(colA),ierr);  // -2 error code
  Epetra_Vector colB(B.OperatorDomainMap());
  EPETRA_TEST_ERR(B.InvColSums(colB),ierr);
  EPETRA_TEST_ERR(checkMultiVectors(colA,colB, "InvColSums", verbose),ierr); // 1 error code

  EPETRA_TEST_ERR(checkValues(A.NormInf(), B.NormInf(), "NormInf before scaling", verbose), ierr);
  EPETRA_TEST_ERR(checkValues(A.NormOne(), B.NormOne(), "NormOne before scaling", verbose),ierr);

  EPETRA_TEST_ERR(A.RightScale(colA),ierr);  // -3 error code
  EPETRA_TEST_ERR(B.RightScale(colB),ierr);  // -3 error code


  EPETRA_TEST_ERR(A.LeftScale(rowA),ierr);
  EPETRA_TEST_ERR(B.LeftScale(rowB),ierr);


  EPETRA_TEST_ERR(checkValues(A.NormInf(), B.NormInf(), "NormInf after scaling", verbose), ierr);
  EPETRA_TEST_ERR(checkValues(A.NormOne(), B.NormOne(), "NormOne after scaling", verbose),ierr);

  vector<double> valuesA(A.MaxNumEntries());
  vector<int> indicesA(A.MaxNumEntries());  
  vector<double> valuesB(B.MaxNumEntries());
  vector<int> indicesB(B.MaxNumEntries());
  //return(0);
  for (int i=0; i<A.NumMyRows(); i++) {
    int nA, nB;
    EPETRA_TEST_ERR(A.ExtractMyRowCopy(i, A.MaxNumEntries(), nA, &valuesA[0], &indicesA[0]),ierr); 
    EPETRA_TEST_ERR(B.ExtractMyRowCopy(i, B.MaxNumEntries(), nB, &valuesB[0], &indicesB[0]),ierr);
    EPETRA_TEST_ERR(!nA==nB,ierr);
    for (int j=0; j<nA; j++) {
      double curVal = valuesA[j];
      int curIndex = indicesA[j];
      bool notfound = true;
      int jj = 0;
      while (notfound && jj< nB) {
	if (!checkValues(curVal, valuesB[jj])) notfound = false;
	jj++;
      }
      EPETRA_TEST_ERR(notfound, ierr);
      vector<int>::iterator p = find(indicesB.begin(),indicesB.end(),curIndex);  // find curIndex in indicesB
      EPETRA_TEST_ERR(p==indicesB.end(), ierr);
    }

  }

  if (verbose) {
    if (ierr==0)
      cout << "RowMatrix Methods check OK" << endl;
    else
      cout << "ierr = " << ierr << ". RowMatrix Methods check failed." << endl;
  }
  return (ierr);
}

