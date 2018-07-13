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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
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
#include "Amesos_SCALAPACK_wrappers.h"
//  #include "Epetra_LAPACK_wrappers.h"


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
  ConTime_(0.0),
  SymTime_(0.0),
  NumTime_(0.0),
  SolTime_(0.0),
  VecTime_(0.0),
  MatTime_(0.0),
  FatOut_(0),
  Time_(0)
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
  if( FatOut_ ) delete FatOut_;
  if( VectorMap_ ) delete VectorMap_;
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

  if( debug_ == 1 ) std::cout << "Entering `RedistributeA()'" << std::endl;

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
  
  if ( debug_ == 1) std::cout  << "iam_ = " << iam_  << "Amesos_Scalapack.cpp:171" << std::endl;
  //
  // Create the ScaLAPACK data distribution.
  // The TwoD data distribution is created in a completely different
  // manner and is not transposed (whereas the SaLAPACK 1D data
  // distribution was transposed) 
  //
  if ( debug_ == 1) std::cout  << "iam_ = " << iam_  << "Amesos_Scalapack.cpp:163" << std::endl;
  Comm().Barrier(); 
  if ( TwoD_distribution_ ) { 
    if ( debug_ == 1) std::cout  << "iam_ = " << iam_  << "Amesos_Scalapack.cpp:166" << std::endl;
    Comm().Barrier(); 
    npcol_ = EPETRA_MIN( NumberOfProcesses, 
			 EPETRA_MAX ( 2, (int) sqrt( NumberOfProcesses * 0.5 ) ) ) ; 
    nprow_ = NumberOfProcesses / npcol_ ;

    //
    //  Create the map for FatA - our first intermediate matrix
    //
    int NumMyElements = RowMatrixA->RowMatrixRowMap().NumMyElements() ;
    std::vector<int> MyGlobalElements( NumMyElements );
    RowMatrixA->RowMatrixRowMap().MyGlobalElements( &MyGlobalElements[0] ) ;

    int NumMyColumns = RowMatrixA->RowMatrixColMap().NumMyElements() ;
    std::vector<int> MyGlobalColumns( NumMyColumns );
    RowMatrixA->RowMatrixColMap().MyGlobalElements( &MyGlobalColumns[0] ) ;

    if ( debug_ == 1) std::cout  << "iam_ = " << iam_  << "Amesos_Scalapack.cpp:194" << std::endl;

    std::vector<int> MyFatElements( NumMyElements * npcol_ ); 
    
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
    
    
    std::vector<std::vector<int> > FatColumnIndices(npcol_,std::vector<int>(1));
    std::vector<std::vector<double> > FatMatrixValues(npcol_,std::vector<double>(1));
    std::vector<int> FatRowPtrs(npcol_);  // A FatRowPtrs[i] = the number 
    // of entries in local row LocalRow*npcol_ + i 
    
    if ( debug_ == 1) std::cout  << "iam_ = " << iam_  << "Amesos_Scalapack.cpp:219" << std::endl;
    //
    mypcol_ = iam_%npcol_;
    myprow_ = (iam_/npcol_)%nprow_;
    if ( iam_ >= nprow_ * npcol_ ) {
      myprow_ = nprow_; 
      mypcol_ = npcol_; 
    }
    //  Each row is split into npcol_ rows, with each of the 
    //  new rows containing only those elements belonging to 
    //  its process column (in the ScaLAPACK 2D process grid)
    //
    int MaxNumIndices = RowMatrixA->MaxNumEntries();
    int NumIndices;
    std::vector<int> ColumnIndices(MaxNumIndices);
    std::vector<double> MatrixValues(MaxNumIndices); 
    
    if ( debug_ == 1) std::cout  << "iam_ = " << iam_  << "Amesos_Scalapack.cpp:232 NumMyElements = " 
			    << NumMyElements 
			    << std::endl;
    
    nb_ = grid_nb_;
    
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
    FatIn.FillComplete( false );
    
    if (  debug_ == 1) std::cout  << "iam_ = " << iam_  << "Amesos_Scalapack.cpp:260" << std::endl;
    if (  debug_ == 1) std::cout  << "Amesos_Scalapack.cpp:265B" 
			     << " iam_ = " << iam_ 
			     << " nb_ = " << nb_ 
			     << " nprow_ = " << nprow_ 
			     << " npcol_ = " << npcol_ 
			     << std::endl;
    
    //
    //  Compute the map for our second intermediate matrix, FatOut
    //
    //  Compute directly
    int UniformRows =  ( NumRows_ / ( nprow_ * nb_ ) ) * nb_ ; 
    int AllExcessRows = NumRows_ - UniformRows * nprow_ ; 
    int OurExcessRows = EPETRA_MIN( nb_, AllExcessRows - ( myprow_ * nb_ ) ) ; 
    OurExcessRows = EPETRA_MAX( 0, OurExcessRows );
    NumOurRows_ = UniformRows + OurExcessRows ; 
    
    if (  debug_ == 1) std::cout  << "iam_ = " << iam_  << "Amesos_Scalapack.cpp:277" << std::endl;
    int UniformColumns =  ( NumColumns_ / ( npcol_ * nb_ ) ) * nb_ ; 
    int AllExcessColumns = NumColumns_ - UniformColumns * npcol_ ; 
    int OurExcessColumns = EPETRA_MIN( nb_, AllExcessColumns - ( mypcol_ * nb_ ) ) ; 
    if (  debug_ == 1) std::cout  << "iam_ = " << iam_  << "Amesos_Scalapack.cpp:281" << std::endl;
    OurExcessColumns = EPETRA_MAX( 0, OurExcessColumns );
    NumOurColumns_ = UniformColumns + OurExcessColumns ; 
    
    if ( iam_ >= nprow_ * npcol_ ) {
      UniformRows = 0;
      NumOurRows_ = 0;
      NumOurColumns_ = 0;
    }
    
    Comm().Barrier(); 
    
    if (  debug_ == 1) std::cout  << "iam_ = " << iam_  << "Amesos_Scalapack.cpp:295" << std::endl;
#if 0
    //  Compute using ScaLAPACK's numroc routine, assert agreement  
    int izero = 0; // All matrices start at process 0
    int NumRocSays = numroc_( &NumRows_, &nb_, &myprow_, &izero, &nprow_ );
    assert( NumOurRows_ == NumRocSays );
#endif
    //
    //  Compute the rows which this process row owns in the ScaLAPACK 2D
    //  process grid.
    //
    std::vector<int> AllOurRows(NumOurRows_);
    
    int RowIndex = 0 ; 
    int BlockRow = 0 ;
    for ( ; BlockRow < UniformRows / nb_ ; BlockRow++ ) {
      for ( int RowOffset = 0; RowOffset < nb_ ; RowOffset++ ) {
	AllOurRows[RowIndex++] = BlockRow*nb_*nprow_  + myprow_*nb_ + RowOffset ;
      } 
    }
    Comm().Barrier(); 
    if (  debug_ == 1) std::cout  << "iam_ = " << iam_  << "Amesos_Scalapack.cpp:315" << std::endl;
    assert ( BlockRow == UniformRows / nb_ ) ; 
    for ( int RowOffset = 0; RowOffset < OurExcessRows ; RowOffset++ ) {
      AllOurRows[RowIndex++] = BlockRow*nb_*nprow_ + myprow_*nb_ + RowOffset ;
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
    
    Comm().Barrier(); 
    if (  debug_ == 1) std::cout  << "iam_ = " << iam_  << "Amesos_Scalapack.cpp:312" << std::endl;
    //
    //  Compute MyRows directly
    //
    std::vector<int>MyRows(NumOurRows_);
    RowIndex = 0 ; 
    BlockRow = 0 ;
    for ( ; BlockRow < UniformRows / nb_ ; BlockRow++ ) {
      for ( int RowOffset = 0; RowOffset < nb_ ; RowOffset++ ) {
	MyRows[RowIndex++] = BlockRow*nb_*nprow_*npcol_  + 
	  myprow_*nb_*npcol_ + RowOffset*npcol_  + mypcol_ ;
      } 
    }
    
    Comm().Barrier(); 
    if (  debug_ == 1) std::cout  << "iam_ = " << iam_  << "Amesos_Scalapack.cpp:326" << std::endl;
    
    assert ( BlockRow == UniformRows / nb_ ) ; 
    for ( int RowOffset = 0; RowOffset < OurExcessRows ; RowOffset++ ) {
      MyRows[RowIndex++] = BlockRow*nb_*nprow_*npcol_  + 
	myprow_*nb_*npcol_ + RowOffset*npcol_  + mypcol_ ;
    } 
    
    Comm().Barrier(); 
    if (  debug_ == 1) std::cout  << "iam_ = " << iam_  << "Amesos_Scalapack.cpp:334" << std::endl;
    Comm().Barrier(); 
    
    for (int i=0; i < NumOurRows_; i++ ) { 
      assert( MyRows[i] == AllOurRows[i]*npcol_+mypcol_ );
    } 
    
    Comm().Barrier(); 
    if (  debug_ == 1) std::cout  << "Amesos_Scalapack.cpp:340" 
			     << " iam_ = " << iam_ 
			     << " myprow_ = " << myprow_ 
			     << " mypcol_ = " << mypcol_ 
			     << " NumRows_ = " << NumRows_ 
			     << " NumOurRows_ = " << NumOurRows_ 
			     << std::endl;
    
    Comm().Barrier(); 
    Epetra_Map FatOutMap( npcol_*NumRows_, NumOurRows_, &MyRows[0], 0, Comm() ); 
    
    if (  debug_ == 1) std::cout  << "iam_ = " << iam_  << "Amesos_Scalapack.cpp:344" << std::endl;
    Comm().Barrier(); 
    
    if ( FatOut_ ) delete FatOut_ ; 
    FatOut_ = new Epetra_CrsMatrix( Copy, FatOutMap, 0 ) ;
    
    if (  debug_ == 1) std::cout  << "iam_ = " << iam_  << "Amesos_Scalapack.cpp:348" << std::endl;
    
    Epetra_Export ExportToFatOut( FatInMap, FatOutMap ) ;
    
    if (  debug_ == 1) std::cout  << "iam_ = " << iam_  << "Amesos_Scalapack.cpp:360" << std::endl;
    
    FatOut_->Export( FatIn, ExportToFatOut, Add );
    FatOut_->FillComplete( false );
    
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
    
    if (  debug_ == 1) std::cout  << "iam_ = " << iam_  << "Amesos_Scalapack.cpp:385" << std::endl;
    
    if (VectorMap_) { delete VectorMap_ ; VectorMap_ = 0 ; } 
    VectorMap_ = new Epetra_Map( NumGlobalElements_, 
				 NumMyVecElements, 
				 &AllOurRows[0], 
				 0, 
				 Comm() );
    if (  debug_ == 1) std::cout  << "iam_ = " << iam_  << " Amesos_Scalapack.cpp:393 debug_ = "
			     << debug_ << std::endl;
    
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
    
    ScaLAPACK1DMatrix_->FillComplete( false ) ; 
  }
  if (  debug_ == 1) std::cout  << "iam_ = " << iam_  << " Amesos_Scalapack.cpp:417 debug_ = "
			   << debug_ << std::endl;
  if (  debug_ == 1) std::cout  << "iam_ = " << iam_  << "Amesos_Scalapack.cpp:402"
			   << " nprow_ = " << nprow_
			   << " npcol_ = " << npcol_ << std::endl ; 
  int info; 
  const int zero = 0 ; 
  if ( ictxt_ == -1313 ) {
    ictxt_ = 0 ; 
    if (  debug_ == 1) std::cout  << "iam_ = " << iam_  << "Amesos_Scalapack.cpp:408" << std::endl;
    SL_INIT_F77(&ictxt_, &nprow_, &npcol_) ; 
  }
  if (  debug_ == 1) std::cout  << "iam_ = " << iam_  << "Amesos_Scalapack.cpp:410A" << std::endl;
  
  int nprow;
  int npcol;
  int myrow;
  int mycol;
  BLACS_GRIDINFO_F77(&ictxt_, &nprow, &npcol, &myrow, &mycol) ; 
  if (  debug_ == 1) std::cout  << "iam_ = " << iam_  << "iam_ = " << iam_ << " Amesos_Scalapack.cpp:410" << std::endl;
  if ( iam_ < nprow_ * npcol_ ) { 
    assert( nprow == nprow_ ) ; 
    if ( npcol != npcol_ ) std::cout << "Amesos_Scalapack.cpp:430 npcol = " << 
      npcol << " npcol_ = " << npcol_ << std::endl ; 
    assert( npcol == npcol_ ) ; 
    if ( TwoD_distribution_ ) {
      assert( myrow == myprow_ ) ; 
      assert( mycol == mypcol_ ) ; 
      lda_ = EPETRA_MAX(1,NumOurRows_) ;
    } else { 
      assert( myrow == 0 ) ; 
      assert( mycol == iam_ ) ; 
      nb_ = m_per_p_;
      lda_ = EPETRA_MAX(1,NumGlobalElements_);
    }
    if (  debug_ == 1) std::cout  << "iam_ = " << iam_  
			     << "Amesos_Scalapack.cpp: " << __LINE__ 
			     << " TwoD_distribution_ = "  << TwoD_distribution_ 
			     << " NumGlobalElements_ = "  << NumGlobalElements_ 
			     << " debug_ = "  << debug_ 
			     << " nb_ = "  << nb_ 
			     << " lda_ = "  << lda_ 
			     << " nprow_ = "  << nprow_ 
			     << " npcol_ = "  << npcol_ 
			     << " myprow_ = "  << myprow_ 
			     << " mypcol_ = "  << mypcol_ 
			     << " iam_ = "  << iam_ << std::endl ;
    AMESOS_PRINT( myprow_ );
    DESCINIT_F77(DescA_, 
		 &NumGlobalElements_, 
		 &NumGlobalElements_, 
		 &nb_,
		 &nb_,
		 &zero,
		 &zero,
		 &ictxt_,
		 &lda_,
		 &info) ;
    if (  debug_ == 1) std::cout  << "iam_ = " << iam_  << "Amesos_Scalapack.cpp:441" << std::endl;
    assert( info == 0 ) ; 
  } else {
    DescA_[0] = -13;
    if (  debug_ == 1) std::cout  << "iam_ = " << iam_  << "Amesos_Scalapack.cpp:458 nprow = " << nprow << std::endl;
    assert( nprow == -1 ) ; 
  }
  
  if (  debug_ == 1) std::cout  << "Amesos_Scalapack.cpp:446" << std::endl;
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
  
  if( debug_ == 1 ) std::cout << "Entering `ConvertToScalapack()'" << std::endl;
  
  Time_->ResetStartTime();
  
  if ( iam_ < nprow_ * npcol_ ) { 
    if ( TwoD_distribution_ ) { 
      
      DenseA_.resize( NumOurRows_ * NumOurColumns_ ); 
      for ( int i = 0 ; i < (int)DenseA_.size() ; i++ ) DenseA_[i] = 0 ; 
      assert( lda_ == EPETRA_MAX(1,NumOurRows_) ) ;
      assert( DescA_[8] == lda_ ) ;
      
      int NzThisRow ;
      int MyRow;
      
      double *RowValues;
      int *ColIndices;
      int MaxNumEntries = FatOut_->MaxNumEntries();
      
      std::vector<int>ColIndicesV(MaxNumEntries);
      std::vector<double>RowValuesV(MaxNumEntries);
      
      int NumMyElements = FatOut_->NumMyRows() ; 
      for ( MyRow = 0; MyRow < NumMyElements ; MyRow++ ) {
	EPETRA_CHK_ERR( FatOut_->
			ExtractMyRowView( MyRow, NzThisRow, RowValues, 
					  ColIndices ) != 0 ) ;
	//
	//  The following eight lines are just a sanity check on MyRow:
	//
	int MyGlobalRow =  FatOut_->GRID( MyRow );
	assert( MyGlobalRow%npcol_ == mypcol_ ) ;   // I should only own rows belonging to my processor column
	int MyTrueRow = MyGlobalRow/npcol_ ;  // This is the original row
	int UniformRows =  ( MyTrueRow / ( nprow_ * nb_ ) ) * nb_ ; 
	int AllExcessRows = MyTrueRow - UniformRows * nprow_ ; 
	int OurExcessRows =  AllExcessRows - ( myprow_ * nb_ ) ; 
	
	if (  MyRow != UniformRows + OurExcessRows ) {
	  std::cout << " iam _ = " << iam_ 
	       << " MyGlobalRow = " << MyGlobalRow 
	       << " MyTrueRow = " << MyTrueRow 
	       << " UniformRows = " << UniformRows 
	       << " AllExcessRows = " << AllExcessRows 
	       << " OurExcessRows = " << OurExcessRows 
	       << " MyRow = " << MyRow  << std::endl ;  
	}
	
	assert( OurExcessRows >= 0 &&  OurExcessRows < nb_ );
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
	  
	  DenseA_[ MyCol * lda_ + MyRow ] = RowValues[j] ; 
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
      
      assert( DescA_[8] == lda_ ) ; //  Double check Lda
      std::vector<int>ColIndicesV(MaxNumEntries);
      std::vector<double>RowValuesV(MaxNumEntries);
      
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
  
  if( debug_ == 1 ) std::cout << "Entering `SetParameters()'" << std::endl;
  
  // retrive general parameters
  SetStatusParameters( ParameterList );
  
  SetControlParameters( ParameterList );
  
  //
  //  We have to set these to their defaults here because user codes 
  //  are not guaranteed to have a "Scalapack" parameter list.
  //
  TwoD_distribution_ = true; 
  grid_nb_ = 32; 
  
  //  Some compilers reject the following cast:
  //  if( &ParameterList == 0 ) return 0;
  
  // ========================================= //
  // retrive ScaLAPACK's parameters from list. //
  // ========================================= //
  
  // retrive general parameters
  //  check to see if they exist before trying to retrieve them
  if (ParameterList.isSublist("Scalapack") ) {
    const Teuchos::ParameterList ScalapackParams = ParameterList.sublist("Scalapack") ;
    //  Fix Bug #3251 
    if ( ScalapackParams.isParameter("2D distribution") )
      TwoD_distribution_ = ScalapackParams.get<bool>("2D distribution");
    if ( ScalapackParams.isParameter("grid_nb") )
      grid_nb_ = ScalapackParams.get<int>("grid_nb");
  }  
  
  return 0;
}

int Amesos_Scalapack::PerformNumericFactorization( ) {
  
  if( debug_ == 1 ) std::cout << "Entering `PerformNumericFactorization()'" << std::endl;
  
  Time_->ResetStartTime();  
  
  Ipiv_.resize(NumGlobalElements_) ;
  for (int i=0; i <NumGlobalElements_ ; i++) 
    Ipiv_[i] = 0 ; // kludge - just to see if this makes valgrind happy
  
  if ( false) std::cout  << " Amesos_Scalapack.cpp: 711 iam_ = " << iam_ << " DescA = " 
		    << DescA_[0] << " " 
		    << DescA_[1] << " " 
		    << DescA_[2] << " " 
		    << DescA_[3] << " " 
		    << DescA_[4] << " " 
		    << DescA_[5] << " " 
		    << DescA_[6] << " " 
		    << DescA_[7] << " " 
		    << DescA_[8] << " " 
		    << std::endl ; 
  
#if 1
  if( NumGlobalElements_ < 10 && nprow_ == 1 && npcol_ == 1 && debug_ == 1 ) {
    assert( lda_ == NumGlobalElements_ ) ; 
    std::cout << " DenseA = " << std::endl ; 
    for (int i=0 ; i < NumGlobalElements_; i++ ) {
      for (int j=0 ; j < NumGlobalElements_; j++ ) {
	if ( DenseA_[ i+j*lda_ ] < 0 ) {
	  DenseA_[ i+j*lda_ ] *= (1+1e-5) ;   // kludge fixme debugxx - just to let vaglrind check to be sure that DenseA is initialized
	}
	//	std::cout << DenseA_[ i+j*lda_ ] << "\t"; 
      }
      //      std::cout << std::endl ; 
    }
  }
#endif
  
  int Ierr[1] ; 
  Ierr[0] = 0 ; 
  const int one = 1 ; 
  if ( iam_ < nprow_ * npcol_ ) {
    if ( nprow_ * npcol_ == 1 ) { 
      DGETRF_F77(&NumGlobalElements_,  
		 &NumGlobalElements_, 
		 &DenseA_[0],
		 &lda_,
		 &Ipiv_[0],
		 Ierr) ;
    } else { 
      PDGETRF_F77(&NumGlobalElements_,  
		  &NumGlobalElements_, 
		  &DenseA_[0],
		  &one,
		  &one, 
		  DescA_,
		  &Ipiv_[0],
		  Ierr) ;
    }
  }
  
  if ( debug_ == 1  && Ierr[0] != 0 ) 
    std::cout << " Amesos_Scalapack.cpp:738 iam_ = " << iam_ 
	 << " Ierr = " << Ierr[0]  << std::endl ; 
  
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
  
  if( debug_ == 1 ) std::cout << "Entering `PerformSymbolicFactorization()'" << std::endl;
  
  if( Time_ == 0 ) Time_ = new Epetra_Time( Comm() );
  
  NumSymbolicFact_++;
  
  return 0;
}

int Amesos_Scalapack::NumericFactorization() {
  
  if( debug_ == 1 ) std::cout << "Entering `NumericFactorization()'" << std::endl;
  
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
  
  if( debug_ == 1 ) std::cout << "Entering `Solve()'" << std::endl;
  
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
  Time_->ResetStartTime(); // track time to broadcast vectors
  //
  //  Copy B to the scalapack version of B
  //
  const Epetra_Map &OriginalMap = RowMatrixA->RowMatrixRowMap();
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
  
  if( nrhs > nb_ ) {
    EPETRA_CHK_ERR( -2 );  
  }
  
  int Ierr[1] ; 
  Ierr[0] = 0 ; 
  const int zero = 0 ; 
  const int one = 1 ; 
  if ( iam_ < nprow_ * npcol_ ) {
    EPETRA_CHK_ERR( ScalapackX->ExtractView( &ScalapackXvalues, &ScalapackXlda ) ) ;
    
    if ( false ) std::cout << "Amesos_Scalapack.cpp: " << __LINE__ << " ScalapackXlda = "  <<  ScalapackXlda 
		      << " lda_ = "  << lda_ 
		      << " nprow_ = "  << nprow_ 
		      << " npcol_ = "  << npcol_ 
		      << " myprow_ = "  << myprow_ 
		      << " mypcol_ = "  << mypcol_ 
		      << " iam_ = "  << iam_ << std::endl ;
    if (  TwoD_distribution_ )    assert( mypcol_ >0 || EPETRA_MAX(ScalapackXlda,1) == lda_ ) ; 
    
    DESCINIT_F77(DescX, 
		 &NumGlobalElements_, 
		 &nrhs, 
		 &nb_,
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
      if ( UseTranspose() ) trans = 'T' ;
    } else {
      
      if ( ! UseTranspose() ) trans = 'T' ;
    }
    
    if ( nprow_ * npcol_ == 1 ) { 
      DGETRS_F77(&trans,
		 &NumGlobalElements_,  
		 &nrhs, 
		 &DenseA_[0],
		 &lda_,
		 &Ipiv_[0],
		 ScalapackXvalues,
		 &lda_,
		 Ierr ) ;
    } else { 
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
      EPETRA_CHK_ERR((*vecX)(i)->Norm2(&NormLHS));
      EPETRA_CHK_ERR((*vecB)(i)->Norm2(&NormRHS));
      if( verbose_ && Comm().MyPID() == 0 ) {
	std::cout << "Amesos_Scalapack : vector " << i << ", ||x|| = " << NormLHS
	     << ", ||b|| = " << NormRHS << std::endl;
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
	std::cout << "Amesos_Scalapack : vector " << i << ", ||Ax - b|| = " << Norm << std::endl;
      }
    }
  }
  
  return Ierr[0];
  
}


// ================================================ ====== ==== ==== == =

void Amesos_Scalapack::PrintStatus() const
{
  
  if( iam_ != 0  ) return;
  
  std::cout << "----------------------------------------------------------------------------" << std::endl;
  std::cout << "Amesos_Scalapack : Matrix has " << NumGlobalElements_ << " rows"
       << " and " << NumGlobalNonzeros_ << " nonzeros" << std::endl;
  std::cout << "Amesos_Scalapack : Nonzero elements per row = "
       << 1.0*NumGlobalNonzeros_/NumGlobalElements_ << std::endl;
  std::cout << "Amesos_Scalapack : Percentage of nonzero elements = "
       << 100.0*NumGlobalNonzeros_/(pow(NumGlobalElements_,2.0)) << std::endl;
  std::cout << "Amesos_Scalapack : Use transpose = " << UseTranspose_ << std::endl;
  std::cout << "----------------------------------------------------------------------------" << std::endl;
  
  return;
  
}

// ================================================ ====== ==== ==== == =

void Amesos_Scalapack::PrintTiming() const
{
  if( iam_ ) return;
  
  std::cout << "----------------------------------------------------------------------------" << std::endl;
  std::cout << "Amesos_Scalapack : Time to convert matrix to ScaLAPACK format = "
       << ConTime_ << " (s)" << std::endl;
  std::cout << "Amesos_Scalapack : Time to redistribute matrix = "
       << MatTime_ << " (s)" << std::endl;
  std::cout << "Amesos_Scalapack : Time to redistribute vectors = "
       << VecTime_ << " (s)" << std::endl;
  std::cout << "Amesos_Scalapack : Number of symbolic factorizations = "
       << NumSymbolicFact_ << std::endl;
  std::cout << "Amesos_Scalapack : Time for sym fact = "
       << SymTime_ << " (s), avg = " << SymTime_/NumSymbolicFact_
       << " (s)" << std::endl;
  std::cout << "Amesos_Scalapack : Number of numeric factorizations = "
       << NumNumericFact_ << std::endl;
  std::cout << "Amesos_Scalapack : Time for num fact = "
       << NumTime_ << " (s), avg = " << NumTime_/NumNumericFact_
       << " (s)" << std::endl;
  std::cout << "Amesos_Scalapack : Number of solve phases = "
       << NumSolve_ << std::endl;
  std::cout << "Amesos_Scalapack : Time for solve = "
       << SolTime_ << " (s), avg = " << SolTime_/NumSolve_
       << " (s)" << std::endl;
  std::cout << "----------------------------------------------------------------------------" << std::endl;
  
  return;
}
