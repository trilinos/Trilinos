b;
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


  //=============================================================================
  Amesos_Scalapack::Amesos_Scalapack(const Epetra_LinearProblem &prob, 
				 const Teuchos::ParameterList &ParameterList ) :  
    ScaLAPACK1DMap_(0), 
    ScaLAPACK1DMatrix_(0), 
    ictxt_(-1313),
    UseTranspose_(false) {

  Problem_ = &prob ; 
  ParameterList_ = &ParameterList ; 
}

//=============================================================================
Amesos_Scalapack::~Amesos_Scalapack(void) {

  if ( ScaLAPACK1DMap_ ) delete ScaLAPACK1DMap_ ; 
  if ( ScaLAPACK1DMatrix_ ) delete ScaLAPACK1DMatrix_ ; 

}
//  See  pre and post conditions in Amesos_Scalapack.h




//
//  Note:
//    The current code uses a data distribution which both ScaLAPACK
//    and Trilinos support:  A 1D data distribution with NPROW=1 
//    and NPCOL = NumberOfProcesses.  
//    
//    With NPROW=1, each column of the matrix is owned by a single 
//    process.  And, the vectors of the matrix are 
//
//    Question:  Is there any way to convert a MultiVector to a CrsMatrix?
//    Assume:  No
//
//    Question:  Is DESCB allowed to differ from DESCA in the call to PSGETRS?
//    Assume:  No 
//
//    Two options:  
//      1)  Create two Epetra_Maps.  One for converting the Epetra_CrsMatrix
//      to a ScaLAPACK matrix, the other for converting the Epetra_Multivector 
//      to a ScaLAPACK matrix.  
//      I can't yet find a clean way to distibute the vectors to the
//      ScaLAPACK matrix that we need from within Trilinos.  
//      I suspect that for the right hand side and left hand side, we
//      have no choice but to convert first to a ScaLAPACK matrix and then 
//      call a ScaLAPACK redistribution routine.  
//      2)  
//      3)  Accept only Epetra_Vectors for now.  
//
//  Goals:
//    1)  Single process, single vector 
//        A)  Write code to convert the Epetra_CrsMatrix to a 1D 
//        multi-process ScaLAPACK matrix.  
//        B)  Write code to convert the Epetra_Vector to a ScaLAPACK 
//        matrix in 1D.  (This code will have no value for the 
//        general multi-vector problem.  Actually, it will have some 
//        value, albeit limited.)  
//    2)  Multiple process, single vector, 1D 
//        No change to the code.  (Barring bugs)
//    3)  Multiple process, 2D 
//        Use the data redistribution routines to convert the data to 
//        a 2D data distribution.  
//        The vectors can be serialized and then redistributed.  
//        We can leave as a bug the serial mature of this particular step.
//    4)  Test with transpose
//    5)  Send the vectors to a different intermediate format:
//        (i.e. remove the serial bottleneck)
//    
//    Coding tasks:
//                                                           .h DOC   .cpp CODE
//        i)  SymbolicFactorization() is a NOP                DONE     NONE
//        ii)  NumericFactorization() - A call to PDGETRF     DONE     DONE 
//        iii) Solve() - a call to PDGETRS                    DONE     DONE
//        iv)  RedistributeA -                                         DONE
//        v)  ConvertToScaLAPACK                                       DONE
//        vi)  Constructor                                             DONE
//        vii)  Destructor                                             DONE
//        viii)  Amesos_SCALAPACK_wrappers.h -                         DONE 
//               what is this CHAR_MACRO thing? For passing characters - such as UPLO.
//        ix)  Where do we want to put the construction of the serialmap and do we 
//             want to call it a serial map?  
//             Let's put it in ConvertToScalapack1D
//
//    Since Trilinos stores its matrices in compressed column 
//    storage with nprow=1, they are fairly close to being in 
//    scalapack 1D form.  
//
//  Note:
//    An efficient ScaLAPACK solution would differ from this in two ways:
//    1)  We would first redistribute the data in a ScaLAPACK form, but 
//    on all processors.
//    2)  We would then use the ScaLAPACK redistribution routine 
//    to redistribute the data to a 2D block cyclic data distribution
//    3)  We would move the vector(s) to 
//


//
int Amesos_Scalapack::RedistributeA( ) {

  cout << "Amesos_Scalapack.cpp::140" << endl ;

  Epetra_RowMatrix *RowMatrixA = dynamic_cast<Epetra_RowMatrix *>(Problem_->GetOperator());
  EPETRA_CHK_ERR( RowMatrixA == 0 ) ; 

  cout << "Amesos_Scalapack.cpp::145" << endl ;

  const Epetra_Map &OriginalMap = RowMatrixA->RowMatrixRowMap() ; 
  int NumberOfProcesses = Comm().NumProc() ; 

  cout << "Amesos_Scalapack.cpp::150" << endl ;
  //
  //  Compute a uniform distribution as ScaLAPACK would want it
  //    MyFirstElement - The first element which this processor would have
  //    NumExpectedElemetns - The number of elements which this processor would have
  //

  int NumRows_ = RowMatrixA->NumGlobalRows() ; 
  if ( MaxProcesses_ > 0 ) {
    NumberOfProcesses = EPETRA_MIN( NumberOfProcesses, MaxProcesses_ ) ; 
  }
  else {
    int ProcessNumHeuristic = (1+NumRows_/200)^2 ;
    NumberOfProcesses = EPETRA_MIN( NumberOfProcesses,  ProcessNumHeuristic );
  }
  nprow_ = 1 ;
  npcol_ = NumberOfProcesses / nprow_ ;
  assert ( nprow_ * npcol_ == NumberOfProcesses ) ; 

  cout << "Amesos_Scalapack.cpp::170" << endl ;

  int m_per_p_ = ( NumRows_ + NumberOfProcesses - 1 ) / NumberOfProcesses ;
  int MyFirstElement = EPETRA_MIN( iam_ * m_per_p_, NumRows_ ) ;
  int MyFirstNonElement = EPETRA_MIN( (iam_+1) * m_per_p_, NumRows_ ) ;
  int NumExpectedElements = MyFirstNonElement - MyFirstElement ; 

  cout << " iam_ = " << iam_ << endl ;
  cout << " m_per_p_ = " << m_per_p_ << endl ;
  cout << " NumberOfPRoceses = " << NumberOfProcesses << endl ;
  cout << " MyFirstElement = " << MyFirstElement << endl ;
  cout << " MyFirstNonElement = " << MyFirstNonElement << endl ;


  assert( NumRows_ ==  RowMatrixA->NumGlobalRows() ) ; 
  if ( ScaLAPACK1DMap_ ) delete( ScaLAPACK1DMap_ ) ; 
  cout << "Amesos_Scalapack.cpp::178" 
       << " NumRows_= " << NumRows_ 
       << " NumExpectedElements= " << NumExpectedElements 
       << endl ;
  ScaLAPACK1DMap_ = new Epetra_Map( NumRows_, NumExpectedElements, 0, Comm() );
  if ( ScaLAPACK1DMatrix_ ) delete( ScaLAPACK1DMatrix_ ) ; 
  cout << "Amesos_Scalapack.cpp::181" << endl ;
  ScaLAPACK1DMatrix_ = new Epetra_CrsMatrix(Copy, *ScaLAPACK1DMap_, 0);
  //  if ( ExportToScaLAPACK1D_ ) delete ExportToScaLAPACK1D_;
  Epetra_Export ExportToScaLAPACK1D_( OriginalMap, *ScaLAPACK1DMap_);
  cout << "Amesos_Scalapack.cpp::187" << endl ;

  ScaLAPACK1DMatrix_->Export( *RowMatrixA, ExportToScaLAPACK1D_, Add ); 
  
  cout << "Amesos_Scalapack.cpp::194" << endl ;

  ScaLAPACK1DMatrix_->FillComplete() ; 

  cout << "Amesos_Scalapack.cpp::203" << endl ;
  cout << " nprow_ = " << nprow_ << endl ;
  cout << " npcol_ = " << npcol_ << endl ;

  int info; 
  const int zero = 0 ; 
  if ( iam_ < nprow_ * npcol_ ) { 
    cout << "Amesos_Scalapack.cpp::208" << endl ;
    if ( ictxt_ == -1313 ) {
      ictxt_ = 0 ; 
      SL_INIT_F77(&ictxt_, &nprow_, &npcol_) ; 
    }
    int nprow;
    int npcol;
    int myrow;
    int mycol;
  cout << "Amesos_Scalapack.cpp::210" << endl ;
    BLACS_GRIDINFO_F77(&ictxt_, &nprow, &npcol, &myrow, &mycol) ; 
    cout << " ictxt_ = " << ictxt_ << endl ; 
    cout << " nprow = " << nprow << endl ; 
    cout << " nprow_ = " << nprow_ << endl ; 
    assert( nprow == nprow_ ) ; 
    assert( npcol == npcol_ ) ; 
    assert( myrow == 0 ) ; 
    assert( mycol == iam_ ) ; 
  cout << "Amesos_Scalapack.cpp::212" << endl ;
    DESCINIT_F77(DescA_, 
		 &NumGlobalElements_, 
		 &NumGlobalElements_, 
		 &m_per_p_,            //  Irrelevant as nprow_ = 1, but may affect blocking
		 &m_per_p_,
		 &zero,
		 &zero,
		 &ictxt_,
		 &NumGlobalElements_,
		 &info) ;
  cout << "Amesos_Scalapack.cpp::223" << endl ;
    assert( info == 0 ) ; 
    }
  

  cout << "Amesos_Scalapack.cpp::210" << endl ;

  return 0;
}


int Amesos_Scalapack::ConvertToScalapack(){
  
  //
  //  Convert matrix and vector to the form that Scalapack expects
  //  ScaLAPACK accepts the matrix to be in any 2D block-cyclic form
  //  We have chosen the simplest 2D block-cyclic form, a 1D blocked (not-cyclic)
  //  data distribution, for the matrix A.
  //  We use the same distribution for the multivectors X and B.  However, 
  //  except for very large numbers of right hand sides, this places all of X and B
  //  on process 0, making it effectively a serial matrix.  
  //  
  //  For now, we simply treat X and B as serial matrices (as viewed from epetra)
  //  though ScaLAPACK treats them as distributed matrices. 
  //
  //  int MyActualFirstElement = ScaLAPACK1DMatrix_->RowMatrixRowMap().MinMyGID() ; 
  int NumMyElements = ScaLAPACK1DMatrix_->NumMyRows() ; 
  //  int nnz_loc = ScaLAPACK1DMatrix_->NumMyNonzeros() ;
  assert( NumGlobalElements_ ==ScaLAPACK1DMatrix_->NumGlobalRows());
  assert( NumGlobalElements_ ==ScaLAPACK1DMatrix_->NumGlobalCols());
  DenseA_.resize( NumGlobalElements_ * NumMyElements ) ;
  for ( int i = 0 ; i < DenseA_.size() ; i++ ) DenseA_[i] = 0 ; 
  
  int NzThisRow ;
  int MyRow;
  int num_my_cols = ScaLAPACK1DMatrix_->NumMyCols() ; 
  if ( ( iam_ < nprow_ * npcol_ - 1 ) && ( m_per_p_ >= 1 ) ) 
    assert( num_my_cols == m_per_p_ ) ; 
  double *RowValues;
  int *ColIndices;
  int MaxNumEntries = ScaLAPACK1DMatrix_->MaxNumEntries();

  cout << "Amesos_Scalapack.cpp::247" << endl ;
  //
  //  
  //
  //  vector <int>Global_Columns( num_my_cols ) ; 

  for ( int i = 0 ; i < num_my_cols ; i ++ ) { 
    //    Global_Columns[i] = ScaLAPACK1DMatrix_->RowMatrixColMap().GID( i ) ; 
    //    assert( Global_Columns[i] == i + iam_ * m_per_p_ ); 
    assert( ScaLAPACK1DMatrix_->RowMatrixColMap().GID( i ) == i + iam_ * m_per_p_ ); 
  }
  assert( DescA_[7] == NumGlobalElements_ ) ; //  Double check Lda

  vector<int>ColIndicesV(MaxNumEntries);
  vector<double>RowValuesV(MaxNumEntries);

  Epetra_CrsMatrix *SuperluCrs = dynamic_cast<Epetra_CrsMatrix *>(ScaLAPACK1DMatrix_);
  for ( MyRow = 0; MyRow < NumMyElements ; MyRow++ ) {
    if ( SuperluCrs != 0 ) {
      EPETRA_CHK_ERR( SuperluCrs->
		      ExtractMyRowView( MyRow, NzThisRow, RowValues, 
					ColIndices ) != 0 ) ;
    }
    else {
      //      ColIndicesV.resize(MaxNumEntries);
      //      RowValuesV.resize(MaxNumEntries);
      EPETRA_CHK_ERR( ScaLAPACK1DMatrix_->
		      ExtractMyRowCopy( MyRow, MaxNumEntries, 
					NzThisRow, &RowValuesV[0], 
					&ColIndicesV[0] ) != 0 );
      RowValues =  &RowValuesV[0];
      ColIndices = &ColIndicesV[0];
    }
    for ( int j = 0; j < NzThisRow; j++ ) { 
      DenseA_[ MyRow + ( ColIndices[j] ) * NumGlobalElements_ ] = RowValues[j] ; 
      //      Ai_[Ai_index] = Global_Columns[ColIndices[j]] ; 
      //      Aval_[Ai_index] = RowValues[j] ; 
      //      Ai_index++;
    }
  }
  //  assert( NumMyElements == MyRow );
  //  Ap_[ NumMyElements ] = Ai_index ; 

  cout << "Amesos_Scalapack.cpp::290" << endl ;

  //
  //  Create a map to allow us to redistribute the vectors X and B 
  //
  Epetra_RowMatrix *RowMatrixA = dynamic_cast<Epetra_RowMatrix *>(Problem_->GetOperator());
  const Epetra_Map &OriginalMap = RowMatrixA->RowMatrixRowMap() ; 
  assert( NumGlobalElements_ == OriginalMap.NumGlobalElements() ) ;
  int NumMyElements_ = 0 ;
  if (iam_==0) NumMyElements_ = NumGlobalElements_;

#if 0
  IsLocal_ = ( OriginalMap.NumMyElements() == 
	       OriginalMap.NumGlobalElements() )?1:0;
  Comm().Broadcast( &IsLocal_, 1, 0 ) ; 
#endif

  if (SerialMap_) { delete SerialMap_ ; SerialMap_ = 0 ; } 
  SerialMap_ = new Epetra_Map( NumGlobalElements_, NumMyElements_, 0, Comm() );

  cout << "Amesos_Scalapack.cpp::310" << endl ;

  return 0;
}   


int Amesos_Scalapack::ReadParameterList() {

  //
  //  We have to set these to their defaults here because user codes 
  //  are not guaranteed to have a "Scalapack" parameter list.
  //
  MaxProcesses_ = - 1; 

  if (ParameterList_->isParameterSublist("Scalapack") ) {
    Teuchos::ParameterList ScalapackParams = ParameterList_->sublist("Scalapack") ;
    MaxProcesses_ = ScalapackParams.getParameter("MaxProcesses",MaxProcesses_);

  }  
  return 0;
}

//
//  Note:  The return value is only valid on the processes
//  which participate in the factorization.
//

int Amesos_Scalapack::PerformNumericFactorization( ) {

  Epetra_RowMatrix *RowMatrixA = dynamic_cast<Epetra_RowMatrix *>(Problem_->GetOperator());
  const Epetra_Map &OriginalMap = RowMatrixA->RowMatrixRowMap() ; 
  NumGlobalElements_ = OriginalMap.NumGlobalElements();
  Ipiv_.resize(NumGlobalElements_) ;

  cout << "Amesos_Scalapack.cpp::344" << endl ;

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

  cout << "Amesos_Scalapack.cpp::360" << endl ;

  //  All processes should return the same error code
  if ( nprow_ * npcol_ < Comm().NumProc() ) 
    Comm().Broadcast( Ierr, 1, 0 ) ; 

  return Ierr[0];
}




bool Amesos_Scalapack::MatrixShapeOK() const { 
  bool OK ;

  if ( GetProblem()->GetOperator()->OperatorRangeMap().NumGlobalPoints() != 
       GetProblem()->GetOperator()->OperatorDomainMap().NumGlobalPoints() ) OK = false;
  return OK; 
}


int Amesos_Scalapack::SymbolicFactorization() {

  return 0;
}

int Amesos_Scalapack::NumericFactorization() {
  
  iam_ = Comm().MyPID();
  
  RedistributeA();
  ConvertToScalapack();


  return PerformNumericFactorization( );
}


int Amesos_Scalapack::Solve() { 

  Epetra_MultiVector   *vecX = Problem_->GetLHS() ; 
  Epetra_MultiVector   *vecB = Problem_->GetRHS() ; 

  cout << "Amesos_Scalapack.cpp::401" << endl ;

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

  Epetra_MultiVector *SerialB =0;
  Epetra_MultiVector *SerialX =0;
  //
  //  Extract Serial versions of X and B 
  //
  double *SerialXvalues ;

  Epetra_RowMatrix *RowMatrixA = dynamic_cast<Epetra_RowMatrix *>(Problem_->GetOperator());
  Epetra_CrsMatrix *CastCrsMatrixA = dynamic_cast<Epetra_CrsMatrix*>(RowMatrixA) ; 
  //  Epetra_MultiVector *SerialXextract = 0;
  //  Epetra_MultiVector *SerialBextract = 0;
    
  //
  //  Copy B to the serial version of B
  //
  const Epetra_Map &OriginalMap = CastCrsMatrixA->RowMap();
  Epetra_MultiVector *SerialXextract = new Epetra_MultiVector( *SerialMap_, nrhs ) ; 
  Epetra_MultiVector *SerialBextract = new Epetra_MultiVector( *SerialMap_, nrhs ) ; 
  
  Epetra_Import ImportToSerial( *SerialMap_, OriginalMap );
  SerialBextract->Import( *vecB, ImportToSerial, Insert ) ;
  SerialB = SerialBextract ; 
  SerialX = SerialXextract ; 

  cout << "Amesos_Scalapack.cpp::440" << endl ;

  //
  //  Call SCALAPACK to perform the solve
  //

  int DescX[10];  
  
  SerialX->Scale(1.0, *SerialB) ;  

  int SerialXlda ; 

  //
  //  Setup DescX 
  //

  assert( nrhs <= m_per_p_ ) ;   // We may be able to remove this by change DescX( nb ) to nrhs

  int Ierr[1] ; 
  Ierr[0] = 0 ; 
  const int zero = 0 ; 
  const int one = 1 ; 
  if ( iam_ < nprow_ * npcol_ ) {
    assert( SerialX->ExtractView( &SerialXvalues, &SerialXlda ) == 0 ) ; 

    assert( SerialXlda == NumGlobalElements_ ) ; 
    
    DESCINIT_F77(DescX, 
		 &NumGlobalElements_, 
		 &nrhs, 
		 &m_per_p_,            //  Irrelevant as nprow_ = 1
		 &m_per_p_,            //  nrhs would be a better choice if legal
		 &zero,
		 &zero,
		 &ictxt_,
		 &NumGlobalElements_,
		 Ierr ) ;
    assert( Ierr[0] == 0 ) ; 
		  
    PDGETRS_F77(&NumGlobalElements_,  
		&nrhs, 
		&DenseA_[0],
		&one,
		&one, 
		&Ipiv_[0],
		DescA_,
		SerialXvalues,
		&one,
		&one, 
		DescX,
		Ierr ) ;
  }

  cout << "Amesos_Scalapack.cpp::494" << endl ;

   
  //
  //  Copy X back to the original vector
  // 

  Epetra_Import ImportFromSerial( OriginalMap, *SerialMap_ );
  vecX->Import( *SerialX, ImportFromSerial, Insert ) ;
  delete SerialBextract ;
  delete SerialXextract ;

  //  All processes should return the same error code
  if ( nprow_ * npcol_ < Comm().NumProc() ) 
    Comm().Broadcast( Ierr, 1, 0 ) ; 

  cout << "Amesos_Scalapack.cpp::510" << endl ;

  return Ierr[0];

}
