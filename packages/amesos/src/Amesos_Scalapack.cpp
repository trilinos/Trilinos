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

#include "Amesos_Scalapack.h"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_Util.h"
#include "CrsMatrixTranspose.h"

//
//  pcolnum computes the process column which a given column belongs to
//  in the ScaLAPACK 2D grid.
//
inline int pcolnum( int j, int nb, int npcol ) {
  return ((j/nb)%npcol) ; }


//=============================================================================
Amesos_Scalapack::Amesos_Scalapack(const Epetra_LinearProblem &prob ):
  ictxt_(-1313),
  ScaLAPACK1DMap_(0), 
  ScaLAPACK1DMatrix_(0), 
  VectorMap_(0),
  UseTranspose_(false),   // Overwritten by call to SetParameters below
  Problem_(&prob), 
  PrintTiming_(false),   // Overwritten by call to SetParameters below
  PrintStatus_(false),   // Overwritten by call to SetParameters below
  ComputeVectorNorms_(false),   // Overwritten by call to SetParameters below
  ComputeTrueResidual_(false),   // Overwritten by call to SetParameters below
  verbose_(1),   // Overwritten by call to SetParameters below
  debug_(0),   // Overwritten by call to SetParameters below
  ConTime_(0.0),
  SymTime_(0.0),
  NumTime_(0.0),
  SolTime_(0.0),
  VecTime_(0.0),
  MatTime_(0.0),
  Time_(0),
  NumSymbolicFact_(0),
  NumNumericFact_(0),
  NumSolve_(0)
{
    Teuchos::ParameterList ParamList ;
    SetParameters( ParamList ) ; 
}

//=============================================================================
Amesos_Scalapack::~Amesos_Scalapack(void) {

  if ( ScaLAPACK1DMap_ ) delete ScaLAPACK1DMap_ ; 
  if ( ScaLAPACK1DMatrix_ ) delete ScaLAPACK1DMatrix_ ; 
  // print out some information if required by the user
  if( (verbose_ && PrintTiming_) || verbose_ == 2 ) PrintTiming();
  if( (verbose_ && PrintStatus_) || verbose_ == 2 ) PrintStatus();

  if( Time_ ) delete Time_;
  
}
//  See  pre and post conditions in Amesos_Scalapack.h

//
//
//  Distribution of the matrix:
//    Amesos_Scalapack uses one of two data distributions:  
//    1)  A 1D blocked data distribution in which each row is assigned
//        to a single process.  The first (n/p) rows are assigned to 
//        the first process and the i^th (n/p) rows are assigned to 
//        process i.   A^T X = B is computed.
//    2)  A 2D data distribution (A X = B is computed)
//
//    The 1D data distribution should be reasonably efficient for the matrices
//    that we expect to deal with, i.e. n on the order of 4,000 
//    with up to 16 processes.  For n >= 10,000 we should switch to 
//    a 2D block-cyclis data distribution.  
//
//  Distribution of the vector(s)
//    1)  1D data distribution
//      ScaLAPACK requires that the vectors be distributed in the same manner
//      as the matrices.  Since PDGETRF factors the transpose of the matrix, 
//      using NPROW=1, the vectors must also be distributed with NPROW=1, i.e.
//      each vector fully owned by a particular process.  And, the first nb
//      vectors belong on the first process. 
//
//      The easiest way to deal with this is to limit ourselves to nb right hand 
//      sides at a time (this is not a significant limitation as nb >= n/p ) 
//      and using our basic heuristic for the number of processes to use 
//      nb >= min(200,n) as well. 
//   
//      Limiting the number of right hand sides to <= nb means that all right hand 
//      sides are stored on process 0.  Hence, they can be treated as a serial 
//      matrix of vectors.
//
//    2)  2D data distribution
//      If we restrict ourselves to nb (typically 32) right hand sides,
//      we can use a simple epetra exporter to export the vector a single
//      ScaLAPACK process column (i.e. to nprow processes)
//

//
int Amesos_Scalapack::RedistributeA( ) {

  if( debug_ == 1 ) cout << "Entering `RedistributeA()'" << endl;

  Time_->ResetStartTime();
  
  Epetra_RowMatrix *RowMatrixA = dynamic_cast<Epetra_RowMatrix *>(Problem_->GetOperator());
  EPETRA_CHK_ERR( RowMatrixA == 0 ) ; 

  const Epetra_Map &OriginalMap = RowMatrixA->RowMatrixRowMap() ; 
  int NumberOfProcesses = Comm().NumProc() ; 

  //
  //  Compute a uniform distribution as ScaLAPACK would want it
  //    MyFirstElement - The first element which this processor would have
  //    NumExpectedElemetns - The number of elements which this processor would have
  //

  int NumRows_ = RowMatrixA->NumGlobalRows() ; 
  int NumColumns_  = RowMatrixA->NumGlobalCols() ; 
  if ( MaxProcesses_ > 0 ) {
    NumberOfProcesses = EPETRA_MIN( NumberOfProcesses, MaxProcesses_ ) ; 
  }
  else {
    int ProcessNumHeuristic = (1+NumRows_/200)*(1+NumRows_/200);
    NumberOfProcesses = EPETRA_MIN( NumberOfProcesses,  ProcessNumHeuristic );
  }

  if( debug_ == 1 ) cout << "Amesos_Scalaapack.cpp:171" << endl;
  //
  // Create the ScaLAPACK data distribution.
  // The TwoD data distribution is created in a completely different
  // manner and is not transposed (whereas the SaLAPACK 1D data
  // distribution was transposed) 
  //
  if ( TwoD_distribution_ ) { 
    npcol_ = EPETRA_MIN( NumberOfProcesses, 
			 EPETRA_MAX ( 2, (int) sqrt( NumberOfProcesses * 0.5 ) ) ) ; 
    nprow_ = NumberOfProcesses / npcol_ ;

    //
    //  Create the map for FatA - our first intermediate matrix
    //
    int NumMyElements = RowMatrixA->RowMatrixRowMap().NumMyElements() ;
    vector<int> MyGlobalElements( NumMyElements );
    RowMatrixA->RowMatrixRowMap().MyGlobalElements( &MyGlobalElements[0] ) ;

    int NumMyColumns = RowMatrixA->RowMatrixColMap().NumMyElements() ;
    vector<int> MyGlobalColumns( NumMyColumns );
    RowMatrixA->RowMatrixColMap().MyGlobalElements( &MyGlobalColumns[0] ) ;

  if( debug_ == 1 ) cout << "Amesos_Scalaapack.cpp:194" << endl;

    vector<int> MyFatElements( NumMyElements * npcol_ ); 

    for( int LocalRow=0; LocalRow<NumMyElements; LocalRow++ ) {
      for (int i = 0 ; i < npcol_; i++ ){
	MyFatElements[LocalRow*npcol_+i] = MyGlobalElements[LocalRow]*npcol_+i;
      } 
    }

    Epetra_Map FatInMap( npcol_*NumRows_, NumMyElements*npcol_, 
			 &MyFatElements[0], 0, Comm() ); 

    //
    //  Create FatIn, our first intermediate matrix
    //
    Epetra_CrsMatrix FatIn( Copy, FatInMap, 0 );


    const int INITIAL_SIZE = 1; 
    vector<vector<int> > FatColumnIndices(npcol_,INITIAL_SIZE);
    vector<vector<double> > FatMatrixValues(npcol_,INITIAL_SIZE);
    vector<int> FatRowPtrs(npcol_);  // A FatRowPtrs[i] = the number 
    // of entries in local row LocalRow*npcol_ + i 

  if( debug_ == 1 ) cout << "Amesos_Scalaapack.cpp:219" << endl;
    //
    mypcol_ = iam_%npcol_;
    myprow_ = (iam_/npcol_)%nprow_;
    //  Each row is split into npcol_ rows, with each of the 
    //  new rows containing only those elements belonging to 
    //  its process column (in the ScaLAPACK 2D process grid)
    //
    int MaxNumIndices = RowMatrixA->MaxNumEntries();
    int NumIndices;
    vector<int> ColumnIndices(MaxNumIndices);
    vector<double> MatrixValues(MaxNumIndices); 

  if( debug_ == 1 ) cout << "Amesos_Scalaapack.cpp:232" << endl;

    for( int LocalRow=0; LocalRow<NumMyElements; ++LocalRow ) {

      RowMatrixA->ExtractMyRowCopy( LocalRow, 
				    MaxNumIndices,
				    NumIndices, 
				    &MatrixValues[0],
				    &ColumnIndices[0] );

      for (int i=0; i<npcol_; i++ )  FatRowPtrs[i] = 0 ; 

      //
      //  Deal the individual matrix entries out to the row owned by
      //  the process to which this matrix entry will belong.
      //
      for( int i=0 ; i<NumIndices ; ++i ) {
	mb_ = grid_mb_;
	nb_ = grid_nb_;
	int GlobalCol = MyGlobalColumns[ ColumnIndices[i] ];
	int pcol_i = pcolnum( GlobalCol, nb_, npcol_ ) ;
	if ( FatRowPtrs[ pcol_i ]+1 >= FatColumnIndices[ pcol_i ].size() ) {
	  FatColumnIndices[ pcol_i ]. resize( 2 * FatRowPtrs[ pcol_i ]+1 );
	  FatMatrixValues[ pcol_i ]. resize( 2 * FatRowPtrs[ pcol_i ]+1 );
	}
	FatColumnIndices[pcol_i][FatRowPtrs[pcol_i]] = GlobalCol ;
	FatMatrixValues[pcol_i][FatRowPtrs[pcol_i]] = MatrixValues[i];

	FatRowPtrs[ pcol_i ]++;
      }

      //
      //  Insert each of the npcol_ rows individually
      //
      for ( int pcol_i = 0 ; pcol_i < npcol_ ; pcol_i++ ) { 
	FatIn.InsertGlobalValues( MyGlobalElements[LocalRow]*npcol_ + pcol_i, 
				  FatRowPtrs[ pcol_i ],
				  &FatMatrixValues[ pcol_i ][0], 
				  &FatColumnIndices[ pcol_i ][0] );
      }
    }
    FatIn.FillComplete();

  if( debug_ == 1 ) cout << "Amesos_Scalaapack.cpp:276" << endl;
    //
    //  Compute the map for our second intermediate matrix, FatOut
    //
    //  Compute directly
    int UniformRows =  ( NumRows_ / ( nprow_ * mb_ ) ) * mb_ ; 
    int AllExcessRows = NumRows_ - UniformRows * nprow_ ; 
    int OurExcessRows = EPETRA_MIN( mb_, AllExcessRows - ( myprow_ * mb_ ) ) ; 
    OurExcessRows = EPETRA_MAX( 0, OurExcessRows );
    NumOurRows_ = UniformRows + OurExcessRows ; 

    int UniformColumns =  ( NumColumns_ / ( npcol_ * mb_ ) ) * mb_ ; 
    int AllExcessColumns = NumColumns_ - UniformColumns * npcol_ ; 
    int OurExcessColumns = EPETRA_MIN( mb_, AllExcessColumns - ( mypcol_ * mb_ ) ) ; 
    OurExcessColumns = EPETRA_MAX( 0, OurExcessColumns );
    NumOurColumns_ = UniformColumns + OurExcessColumns ; 

  if( debug_ == 1 ) cout << "Amesos_Scalaapack.cpp:295" << endl;
#if 0
    //  Compute using ScaLAPACK's numroc routine, assert agreement  
    int izero = 0; // All matrices start at process 0
    int NumRocSays = numroc_( &NumRows_, &mb_, &myprow_, &izero, &nprow_ );
    assert( NumOurRows_ == NumRocSays );
#endif
    //
    //  Compute the rows which this process row owns in the ScaLAPACK 2D
    //  process grid.
    //
    vector<int> AllOurRows(NumOurRows_);

    int RowIndex = 0 ; 
    int BlockRow = 0 ;
    for ( ; BlockRow < UniformRows / mb_ ; BlockRow++ ) {
      for ( int RowOffset = 0; RowOffset < mb_ ; RowOffset++ ) {
	AllOurRows[RowIndex++] = BlockRow*mb_*nprow_  + myprow_*mb_ + RowOffset ;
      } 
    }
  if( debug_ == 1 ) cout << "Amesos_Scalaapack.cpp:315" << endl;
    assert ( BlockRow == UniformRows / mb_ ) ; 
    for ( int RowOffset = 0; RowOffset < OurExcessRows ; RowOffset++ ) {
      AllOurRows[RowIndex++] = BlockRow*mb_*nprow_ + myprow_*mb_ + RowOffset ;
    } 
    assert( RowIndex == NumOurRows_ );
    //
    //  Distribute those rows amongst all the processes in that process row
    //  This is an artificial distribution with the following properties:
    //  1)  It is a 1D data distribution (each row belogs entirely to 
    //      a single process
    //  2)  All data which will eventually belong to a given process row, 
    //      is entirely contained within the processes in that row.
    //

  if( debug_ == 1 ) cout << "Amesos_Scalaapack.cpp:330" << endl;
    //
    //  Compute MyRows directly
    //
    vector<int>MyRows(NumOurRows_);
    RowIndex = 0 ; 
    BlockRow = 0 ;
    for ( ; BlockRow < UniformRows / mb_ ; BlockRow++ ) {
      for ( int RowOffset = 0; RowOffset < mb_ ; RowOffset++ ) {
	MyRows[RowIndex++] = BlockRow*mb_*nprow_*npcol_  + 
	  myprow_*mb_*npcol_ + RowOffset*npcol_  + mypcol_ ;
      } 
    }
    assert ( BlockRow == UniformRows / mb_ ) ; 
    for ( int RowOffset = 0; RowOffset < OurExcessRows ; RowOffset++ ) {
      MyRows[RowIndex++] = BlockRow*mb_*nprow_*npcol_  + 
	myprow_*mb_*npcol_ + RowOffset*npcol_  + mypcol_ ;
    } 

    for (int i=0; i < NumOurRows_; i++ ) { 
      assert( MyRows[i] == AllOurRows[i]*npcol_+mypcol_ );
    } 

    Epetra_Map FatOutMap( npcol_*NumRows_, NumOurRows_, &MyRows[0], 0, Comm() ); 


    FatOut_ = new Epetra_CrsMatrix( Copy, FatOutMap, 0 ) ;
  
    Epetra_Export ExportToFatOut( FatInMap, FatOutMap ) ;

  if( debug_ == 1 ) cout << "Amesos_Scalaapack.cpp:360" << endl;

    FatOut_->Export( FatIn, ExportToFatOut, Add );
    FatOut_->FillComplete();

    //
    //  Create a map to allow us to redistribute the vectors X and B 
    //
    Epetra_RowMatrix *RowMatrixA = dynamic_cast<Epetra_RowMatrix *>(Problem_->GetOperator());
    const Epetra_Map &OriginalMap = RowMatrixA->RowMatrixRowMap() ; 
    assert( NumGlobalElements_ == OriginalMap.NumGlobalElements() ) ;
    int NumMyVecElements = 0 ;
    if ( mypcol_ == 0 ) { 
      NumMyVecElements = NumOurRows_;
    }
    if (VectorMap_) { delete VectorMap_ ; VectorMap_ = 0 ; } 
    VectorMap_ = new Epetra_Map( NumGlobalElements_, 
				 NumMyVecElements, 
				 &AllOurRows[0], 
				 0, 
				 Comm() );
  } else {
    nprow_ = 1 ;
    npcol_ = NumberOfProcesses / nprow_ ;
    assert ( nprow_ * npcol_ == NumberOfProcesses ) ; 

    m_per_p_ = ( NumRows_ + NumberOfProcesses - 1 ) / NumberOfProcesses ;
    int MyFirstElement = EPETRA_MIN( iam_ * m_per_p_, NumRows_ ) ;
    int MyFirstNonElement = EPETRA_MIN( (iam_+1) * m_per_p_, NumRows_ ) ;
    int NumExpectedElements = MyFirstNonElement - MyFirstElement ; 
    
    assert( NumRows_ ==  RowMatrixA->NumGlobalRows() ) ; 
    if ( ScaLAPACK1DMap_ ) delete( ScaLAPACK1DMap_ ) ; 
    ScaLAPACK1DMap_ = new Epetra_Map( NumRows_, NumExpectedElements, 0, Comm() );
    if ( ScaLAPACK1DMatrix_ ) delete( ScaLAPACK1DMatrix_ ) ; 
    ScaLAPACK1DMatrix_ = new Epetra_CrsMatrix(Copy, *ScaLAPACK1DMap_, 0);
    Epetra_Export ExportToScaLAPACK1D_( OriginalMap, *ScaLAPACK1DMap_);
    
    ScaLAPACK1DMatrix_->Export( *RowMatrixA, ExportToScaLAPACK1D_, Add ); 
    
    ScaLAPACK1DMatrix_->FillComplete() ; 
  }
  if( debug_ == 1 ) cout << "Amesos_Scalaapack.cpp:402" << endl;
  cout << " nprow_ = " << nprow_ << endl ; 
  cout << " npcol_ = " << npcol_ << endl ; 
  int info; 
  const int zero = 0 ; 
  if ( ictxt_ == -1313 ) {
    ictxt_ = 0 ; 
  if( debug_ == 1 ) cout << "Amesos_Scalaapack.cpp:408" << endl;
    SL_INIT_F77(&ictxt_, &nprow_, &npcol_) ; 
  }
  if( debug_ == 1 ) cout << "Amesos_Scalaapack.cpp:410" << endl;

  int nprow;
  int npcol;
  int myrow;
  int mycol;
  BLACS_GRIDINFO_F77(&ictxt_, &nprow, &npcol, &myrow, &mycol) ; 
  if ( iam_ < nprow_ * npcol_ ) { 
    assert( nprow == nprow_ ) ; 
    assert( npcol == npcol_ ) ; 
    if ( TwoD_distribution_ ) {
      assert( myrow == myprow_ ) ; 
      assert( mycol == mypcol_ ) ; 
      lda_ = NumOurRows_ ;
    } else { 
      assert( myrow == 0 ) ; 
      assert( mycol == iam_ ) ; 
      mb_ = m_per_p_;            //  Irrelevant as nprow_ = 1, but may affect blocking
      nb_ = m_per_p_;
      lda_ = NumGlobalElements_;
    }
  if( debug_ == 1 ) cout << "Amesos_Scalaapack.cpp:430" << endl;
    DESCINIT_F77(DescA_, 
		 &NumGlobalElements_, 
		 &NumGlobalElements_, 
		 &mb_,
		 &nb_,
		 &zero,
		 &zero,
		 &ictxt_,
		 &lda_,
		 &info) ;
  if( debug_ == 1 ) cout << "Amesos_Scalaapack.cpp:441" << endl;
    assert( info == 0 ) ; 
  } else {
    DescA_[0] = -13;
    assert( nprow == -1 ) ; 
  }

  if( debug_ == 1 ) cout << "Amesos_Scalaapack.cpp:446" << endl;
  MatTime_ += Time_->ElapsedTime();
  
  return 0;
}


int Amesos_Scalapack::ConvertToScalapack(){
  
  //
  //  Convert matrix and vector to the form that Scalapack expects
  //  ScaLAPACK accepts the matrix to be in any 2D block-cyclic form
  //
  //  Amesos_ScaLAPACK uses one of two 2D data distributions: 
  //  a simple 1D non-cyclic data distribution with npcol= 1, or a 
  //  full  2D block-cyclic data distribution.
  //
  //  2D data distribvution:
  //    Because the Epetra export operation is oriented toward a 1D 
  //    data distribution in which each row is entirely stored on 
  //    a single process, we create two intermediate matrices: FatIn and
  //    FatOut, both of which have dimension:  
  //      NumGlobalElements * nprow by NumGlobalElements
  //    This allows each row of FatOut to be owned by a single process.
  //    The larger dimension does not significantly increase the 
  //    storage requirements and allows the export operation to be 
  //    efficient.  
  //
  //  1D data distribution:
  //  We have chosen the simplest 2D block-cyclic form, a 1D blocked (not-cyclic)
  //  data distribution, for the matrix A.
  //  We use the same distribution for the multivectors X and B.  However, 
  //  except for very large numbers of right hand sides, this places all of X and B
  //  on process 0, making it effectively a serial matrix.  
  //  
  //  For now, we simply treat X and B as serial matrices (as viewed from epetra)
  //  though ScaLAPACK treats them as distributed matrices. 
  //

  if( debug_ == 1 ) cout << "Entering `ConvertToScalapack()'" << endl;

  Time_->ResetStartTime();
  
  if ( iam_ < nprow_ * npcol_ ) { 
    if ( TwoD_distribution_ ) { 

      DenseA_.resize( NumOurRows_ * NumOurColumns_ ); 
      for ( int i = 0 ; i < (int)DenseA_.size() ; i++ ) DenseA_[i] = 0 ; 
      int lda = NumOurRows_ ;
      assert( DescA_[8] == lda ) ;
      
      int NzThisRow ;
      int MyRow;
    
      double *RowValues;
      int *ColIndices;
      int MaxNumEntries = FatOut_->MaxNumEntries();
     
      assert( DescA_[8] == NumGlobalElements_ ) ; //  Double check Lda
      vector<int>ColIndicesV(MaxNumEntries);
      vector<double>RowValuesV(MaxNumEntries);
    
      int NumMyElements = FatOut_->NumMyRows() ; 
      for ( MyRow = 0; MyRow < NumMyElements ; MyRow++ ) {
	EPETRA_CHK_ERR( FatOut_->
			ExtractMyRowView( MyRow, NzThisRow, RowValues, 
					  ColIndices ) != 0 ) ;
	//
	//  The following eight lines are just a sanity check on MyRow:
	//
	int MyGlobalRow =  FatOut_->GRID( MyRow );
	assert( (MyGlobalRow/mb_)%nprow_ == myprow_ ) ; 
	int UniformRows =  ( MyGlobalRow / ( nprow_ * mb_ ) ) * mb_ ; 
	int AllExcessRows = MyGlobalRow - UniformRows * nprow_ ; 
	int OurExcessRows =  AllExcessRows - ( myprow_ * mb_ ) ; 
	assert( OurExcessRows >= 0 &&  OurExcessRows < mb_ );
	assert( MyRow == UniformRows + OurExcessRows ) ; 

	for ( int j = 0; j < NzThisRow; j++ ) { 
	  assert(  FatOut_->RowMatrixColMap().GID( ColIndices[j] ) ==
		   FatOut_->GCID( ColIndices[j] ) );

	  int MyGlobalCol =  FatOut_->GCID( ColIndices[j] );
	  assert( (MyGlobalCol/nb_)%npcol_ == mypcol_ ) ; 
	  int UniformCols =  ( MyGlobalCol / ( npcol_ * nb_ ) ) * nb_ ; 
	  int AllExcessCols = MyGlobalCol - UniformCols * npcol_ ; 
	  int OurExcessCols =  AllExcessCols - ( mypcol_ * nb_ ) ; 
	  assert( OurExcessCols >= 0 &&  OurExcessCols < nb_ );
	  int MyCol = UniformCols + OurExcessCols ; 

	  DenseA_[ MyCol * lda + MyRow ] = RowValues[j] ; 
	}
      }
      


    } else { 
      
      int NumMyElements = ScaLAPACK1DMatrix_->NumMyRows() ; 
      
      assert( NumGlobalElements_ ==ScaLAPACK1DMatrix_->NumGlobalRows());
      assert( NumGlobalElements_ ==ScaLAPACK1DMatrix_->NumGlobalCols());
      DenseA_.resize( NumGlobalElements_ * NumMyElements ) ;
      for ( int i = 0 ; i < (int)DenseA_.size() ; i++ ) DenseA_[i] = 0 ; 
    
      int NzThisRow ;
      int MyRow;
    
      double *RowValues;
      int *ColIndices;
      int MaxNumEntries = ScaLAPACK1DMatrix_->MaxNumEntries();
     
      assert( DescA_[8] == NumGlobalElements_ ) ; //  Double check Lda
      vector<int>ColIndicesV(MaxNumEntries);
      vector<double>RowValuesV(MaxNumEntries);
    
      for ( MyRow = 0; MyRow < NumMyElements ; MyRow++ ) {
	EPETRA_CHK_ERR( ScaLAPACK1DMatrix_->
			ExtractMyRowView( MyRow, NzThisRow, RowValues, 
					  ColIndices ) != 0 ) ;
	
	for ( int j = 0; j < NzThisRow; j++ ) { 
	  DenseA_[ ( ScaLAPACK1DMatrix_->RowMatrixColMap().GID( ColIndices[j] ) ) 
		   + MyRow * NumGlobalElements_ ] = RowValues[j] ; 
	}
      }
      //
      //  Create a map to allow us to redistribute the vectors X and B 
      //
      Epetra_RowMatrix *RowMatrixA = dynamic_cast<Epetra_RowMatrix *>(Problem_->GetOperator());
      const Epetra_Map &OriginalMap = RowMatrixA->RowMatrixRowMap() ; 
      assert( NumGlobalElements_ == OriginalMap.NumGlobalElements() ) ;
      int NumMyElements_ = 0 ;
      if (iam_==0) NumMyElements_ = NumGlobalElements_;
      
      if (VectorMap_) { delete VectorMap_ ; VectorMap_ = 0 ; } 
      VectorMap_ = new Epetra_Map( NumGlobalElements_, NumMyElements_, 0, Comm() );
    }
  }
  ConTime_ += Time_->ElapsedTime();
  
return 0;
}   


int Amesos_Scalapack::SetParameters( Teuchos::ParameterList &ParameterList ) {

  if( debug_ == 1 ) cout << "Entering `SetParameters()'" << endl;

  //
  //  We have to set these to their defaults here because user codes 
  //  are not guaranteed to have a "Scalapack" parameter list.
  //
  MaxProcesses_ = - 1; 
  bool UseTrans = false ; 
  PrintTiming_ = false ; 
  PrintStatus_ = false ; 
  ComputeVectorNorms_ = false ; 
  ComputeTrueResidual_ = false ; 
  verbose_ = 1; 
  //  debug_ = 0;    // This is set below because debug is used in this routine
  TwoD_distribution_ = true; 
  grid_mb_ = 32; 
  grid_nb_ = 32; 

  //  Some compilers reject the following cast:
  //  if( &ParameterList == 0 ) return 0;

  // define how many processes to use in the ScaLAPACK factor and solve
  // if (-1), a heuristic is used to determine the number of processes to use 
  if( ParameterList.isParameter("MaxProcs") )
    MaxProcesses_ = ParameterList.get("MaxProcs",MaxProcesses_);
  
  // ========================================= //
  // retrive ScaLAPACK's parameters from list. //
  // ========================================= //
  
  // retrive general parameters

  // solve problem with transpose
  if( ParameterList.isParameter("UseTranspose") )
    UseTrans = ParameterList.get("UseTranspose",UseTrans);
  SetUseTranspose(UseTrans);

  // print some timing information (on process 0)
  if( ParameterList.isParameter("PrintTiming") )
    PrintTiming_ = ParameterList.get("PrintTiming", PrintTiming_ );

  // print some statistics (on process 0). Do not include timing
  if( ParameterList.isParameter("PrintStatus") )
    PrintStatus_ = ParameterList.get("PrintStatus", PrintStatus_);

  // compute norms of some vectors
  if( ParameterList.isParameter("ComputeVectorNorms") )
    ComputeVectorNorms_ = ParameterList.get("ComputeVectorNorms",false);

  // compute the true residual Ax-b after solution
  if( ParameterList.isParameter("ComputeTrueResidual") )
    ComputeTrueResidual_ = ParameterList.get("ComputeTrueResidual",false);

  // some verbose output:
  // 0 - no output at all
  // 1 - output as specified by other parameters
  // 2 - all possible output
  if( ParameterList.isParameter("OutputLevel") )
    verbose_ = ParameterList.get("OutputLevel",verbose_);

  // possible debug statements
  // 0 - no debug
  // 1 - debug
  debug_ = 0 ;
  if( ParameterList.isParameter("DebugLevel") )
    debug_ = ParameterList.get("DebugLevel",debug_);


  if (ParameterList.isSublist("Scalapack") ) {
    Teuchos::ParameterList ScalapackParams = ParameterList.sublist("Scalapack") ;
    TwoD_distribution_ = ScalapackParams.get("2D distribution",TwoD_distribution_);
    grid_mb_ = ScalapackParams.get("grid_mb",grid_mb_);
    grid_nb_ = ScalapackParams.get("grid_nb",grid_nb_);
  }  
  
  return 0;
}

int Amesos_Scalapack::PerformNumericFactorization( ) {

  if( debug_ == 1 ) cout << "Entering `PerformNumericFactorization()'" << endl;
  
  Time_->ResetStartTime();  

  Ipiv_.resize(NumGlobalElements_) ;

  int Ierr[1] ; 
  Ierr[0] = 0 ; 
  const int one = 1 ; 
  if ( iam_ < nprow_ * npcol_ ) {
    PDGETRF_F77(&NumGlobalElements_,  
		&NumGlobalElements_, 
		&DenseA_[0],
		&one,
		&one, 
		DescA_,
		&Ipiv_[0],
		Ierr) ;
  }

  //  All processes should return the same error code
  if ( nprow_ * npcol_ < Comm().NumProc() ) 
    Comm().Broadcast( Ierr, 1, 0 ) ; 

  NumTime_ += Time_->ElapsedTime();

  return Ierr[0];
}




bool Amesos_Scalapack::MatrixShapeOK() const { 
  bool OK ;

  if ( GetProblem()->GetOperator()->OperatorRangeMap().NumGlobalPoints() != 
       GetProblem()->GetOperator()->OperatorDomainMap().NumGlobalPoints() ) OK = false;
  return OK; 
}


int Amesos_Scalapack::SymbolicFactorization() {

  if( debug_ == 1 ) cout << "Entering `PerformSymbolicFactorization()'" << endl;

  if( Time_ == 0 ) Time_ = new Epetra_Time( Comm() );

  NumSymbolicFact_++;
  
  return 0;
}

int Amesos_Scalapack::NumericFactorization() {

  if( debug_ == 1 ) cout << "Entering `NumericFactorization()'" << endl;

  NumNumericFact_++;

  iam_ = Comm().MyPID();
  
  Epetra_RowMatrix *RowMatrixA = dynamic_cast<Epetra_RowMatrix *>(Problem_->GetOperator());
  const Epetra_Map &OriginalMap = RowMatrixA->RowMatrixRowMap() ; 
  NumGlobalElements_ = OriginalMap.NumGlobalElements();

  NumGlobalNonzeros_ = RowMatrixA->NumGlobalNonzeros();

  RedistributeA();
  ConvertToScalapack();

  return PerformNumericFactorization( );
}


int Amesos_Scalapack::Solve() { 

  if( debug_ == 1 ) cout << "Entering `Solve()'" << endl;

  NumSolve_++;

  Epetra_MultiVector   *vecX = Problem_->GetLHS() ; 
  Epetra_MultiVector   *vecB = Problem_->GetRHS() ; 

  //
  //  Compute the number of right hands sides 
  //  (and check that X and B have the same shape) 
  //
  int nrhs; 
  if ( vecX == 0 ) { 
    nrhs = 0 ;
    EPETRA_CHK_ERR( vecB != 0 ) ; 
  } else { 
    nrhs = vecX->NumVectors() ; 
    EPETRA_CHK_ERR( vecB->NumVectors() != nrhs ) ; 
  }

  Epetra_MultiVector *ScalapackB =0;
  Epetra_MultiVector *ScalapackX =0;
  //
  //  Extract Scalapack versions of X and B 
  //
  double *ScalapackXvalues ;

  Epetra_RowMatrix *RowMatrixA = dynamic_cast<Epetra_RowMatrix *>(Problem_->GetOperator());
  Epetra_CrsMatrix *CastCrsMatrixA = dynamic_cast<Epetra_CrsMatrix*>(RowMatrixA) ; 
  Time_->ResetStartTime(); // track time to broadcast vectors
  //
  //  Copy B to the scalapack version of B
  //
  const Epetra_Map &OriginalMap = CastCrsMatrixA->RowMap();
  Epetra_MultiVector *ScalapackXextract = new Epetra_MultiVector( *VectorMap_, nrhs ) ; 
  Epetra_MultiVector *ScalapackBextract = new Epetra_MultiVector( *VectorMap_, nrhs ) ; 
  
  Epetra_Import ImportToScalapack( *VectorMap_, OriginalMap );
  ScalapackBextract->Import( *vecB, ImportToScalapack, Insert ) ;
  ScalapackB = ScalapackBextract ; 
  ScalapackX = ScalapackXextract ; 

  VecTime_ += Time_->ElapsedTime();

  //
  //  Call SCALAPACKs PDGETRS to perform the solve
  //

  int DescX[10];  
  
  ScalapackX->Scale(1.0, *ScalapackB) ;  

  int ScalapackXlda ; 

  Time_->ResetStartTime(); // tract time to solve

  //
  //  Setup DescX 
  //
  cout << " nrhs = " << nrhs << endl;
  cout << " nb_ = " << nb_ << endl;

  if( nrhs > nb_ ) {
    EPETRA_CHK_ERR( -2 );  
  }

  int Ierr[1] ; 
  Ierr[0] = 0 ; 
  const int zero = 0 ; 
  const int one = 1 ; 
  if ( iam_ < nprow_ * npcol_ ) {
    assert( ScalapackX->ExtractView( &ScalapackXvalues, &ScalapackXlda ) == 0 ) ; 

    assert( iam_ >0 || ScalapackXlda == NumGlobalElements_ ) ; 
    
    DESCINIT_F77(DescX, 
		 &NumGlobalElements_, 
		 &nrhs, 
		 &mb_,
		 &nb_,
		 &zero,
		 &zero,
		 &ictxt_,
		 &lda_,
		 Ierr ) ;
    assert( Ierr[0] == 0 ) ; 
		
    //
    //  For the 1D data distribution, we factor the transposed 
    //  matrix, hence we must invert the sense of the transposition
    //
    char trans = 'N';
    if ( TwoD_distribution_ ) {
      if ( UseTranspose() ) trans = 'Y' ;
    } else {
      if ( ! UseTranspose() ) trans = 'Y' ;
    }

    PDGETRS_F77(&trans,
		&NumGlobalElements_,  
		&nrhs, 
		&DenseA_[0],
		&one,
		&one, 
		DescA_,
		&Ipiv_[0],
		ScalapackXvalues,
		&one,
		&one, 
		DescX,
		Ierr ) ;
  }

  SolTime_ += Time_->ElapsedTime();

  Time_->ResetStartTime();  // track time to broadcast vectors
  //
  //  Copy X back to the original vector
  // 
  Epetra_Import ImportFromScalapack( OriginalMap, *VectorMap_ );
  vecX->Import( *ScalapackX, ImportFromScalapack, Insert ) ;
  delete ScalapackBextract ;
  delete ScalapackXextract ;

  VecTime_ += Time_->ElapsedTime();

  //  All processes should return the same error code
  if ( nprow_ * npcol_ < Comm().NumProc() ) 
    Comm().Broadcast( Ierr, 1, 0 ) ; 

  // MS // compute vector norms
  if( ComputeVectorNorms_ == true || verbose_ == 2 ) {
    double NormLHS, NormRHS;
    for( int i=0 ; i<nrhs ; ++i ) {
      assert((*vecX)(i)->Norm2(&NormLHS)==0);
      assert((*vecB)(i)->Norm2(&NormRHS)==0);
      if( verbose_ && Comm().MyPID() == 0 ) {
	cout << "Amesos_Scalapack : vector " << i << ", ||x|| = " << NormLHS
	     << ", ||b|| = " << NormRHS << endl;
      }
    }
  }
  
  // MS // compute true residual
  if( ComputeTrueResidual_ == true || verbose_ == 2  ) {
    double Norm;
    Epetra_MultiVector Ax(vecB->Map(),nrhs);
    for( int i=0 ; i<nrhs ; ++i ) {
      (Problem_->GetMatrix()->Multiply(UseTranspose(), *((*vecX)(i)), Ax));
      (Ax.Update(1.0, *((*vecB)(i)), -1.0));
      (Ax.Norm2(&Norm));
      
      if( verbose_ && Comm().MyPID() == 0 ) {
	cout << "Amesos_Scalapack : vector " << i << ", ||Ax - b|| = " << Norm << endl;
      }
    }
  }
  
  return Ierr[0];

}


// ================================================ ====== ==== ==== == =

void Amesos_Scalapack::PrintStatus() 
{

  if( iam_ != 0  ) return;

  cout << "----------------------------------------------------------------------------" << endl;
  cout << "Amesos_Scalapack : Matrix has " << NumGlobalElements_ << " rows"
       << " and " << NumGlobalNonzeros_ << " nonzeros" << endl;
  cout << "Amesos_Scalapack : Nonzero elements per row = "
       << 1.0*NumGlobalNonzeros_/NumGlobalElements_ << endl;
  cout << "Amesos_Scalapack : Percentage of nonzero elements = "
       << 100.0*NumGlobalNonzeros_/(pow(NumGlobalElements_,2.0)) << endl;
  cout << "Amesos_Scalapack : Use transpose = " << UseTranspose_ << endl;
  cout << "----------------------------------------------------------------------------" << endl;

  return;
  
}

// ================================================ ====== ==== ==== == =

void Amesos_Scalapack::PrintTiming()
{
  if( iam_ ) return;
  
  cout << "----------------------------------------------------------------------------" << endl;
  cout << "Amesos_Scalapack : Time to convert matrix to ScaLAPACK format = "
       << ConTime_ << " (s)" << endl;
  cout << "Amesos_Scalapack : Time to redistribute matrix = "
       << MatTime_ << " (s)" << endl;
  cout << "Amesos_Scalapack : Time to redistribute vectors = "
       << VecTime_ << " (s)" << endl;
  cout << "Amesos_Scalapack : Number of symbolic factorizations = "
       << NumSymbolicFact_ << endl;
  cout << "Amesos_Scalapack : Time for sym fact = "
       << SymTime_ << " (s), avg = " << SymTime_/NumSymbolicFact_
       << " (s)" << endl;
  cout << "Amesos_Scalapack : Number of numeric factorizations = "
       << NumNumericFact_ << endl;
  cout << "Amesos_Scalapack : Time for num fact = "
       << NumTime_ << " (s), avg = " << NumTime_/NumNumericFact_
       << " (s)" << endl;
  cout << "Amesos_Scalapack : Number of solve phases = "
       << NumSolve_ << endl;
  cout << "Amesos_Scalapack : Time for solve = "
       << SolTime_ << " (s), avg = " << SolTime_/NumSolve_
       << " (s)" << endl;
  cout << "----------------------------------------------------------------------------" << endl;
   
  return;
}
