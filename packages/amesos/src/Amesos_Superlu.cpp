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

#include "Amesos_Superlu.h"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_Util.h"
//  #define USE_DGSTRF

namespace SLU
{
extern "C" {
#include "dsp_defs.h"
}
}
struct SLUData
{
  SLU::SuperMatrix A, B, X, L, U;
#ifdef USE_DGSTRF
  SLU::SuperMatrix AC;
#endif
  SLU::superlu_options_t SLU_options;
  SLU::mem_usage_t mem_usage;
  SLU::fact_t refactor_option ;         //  SamePattern or SamePattern_SameRowPerm 
};



  //=============================================================================
  Amesos_Superlu::Amesos_Superlu(const Epetra_LinearProblem &prob ):
    SerialCrsMatrixA_(0), 
    SerialMap_(0), 
    SerialMatrix_(0), 
    FactorizationDone_(false),
    FactorizationOK_(false),
    UseTranspose_(false),
    iam_(-1) {

  using namespace SLU;

  data_ = new SLUData();

  Problem_ = &prob ; 
  DestroyBandX_ = true;

  dCreate_Dense_Matrix( &(data_->X), 
			0, 
			0, 
			&DummyArray[0],
			0, 
			SLU_DN, SLU_D, SLU_GE);
    
  dCreate_Dense_Matrix( &(data_->B), 
			0, 
			0, 
			&DummyArray[0],
			0, 
			SLU_DN, SLU_D, SLU_GE);
    

  Teuchos::ParameterList ParamList ;
  SetParameters( ParamList ) ; 

}

//=============================================================================
Amesos_Superlu::~Amesos_Superlu(void) {

  using namespace SLU;

  if ( SerialMap_ ) delete SerialMap_ ; 
  if ( SerialCrsMatrixA_ ) delete SerialCrsMatrixA_ ; 

  Destroy_SuperMatrix_Store(&data_->B);
  Destroy_SuperMatrix_Store(&data_->X);
  if ( iam_ == 0 ) { 
    if ( FactorizationDone_ ) { 
      Destroy_SuperMatrix_Store(&data_->A);
      Destroy_SuperNode_Matrix(&data_->L);
      Destroy_CompCol_Matrix(&data_->U);
    }
  }
      

  delete data_; 

}
//  See  pre and post conditions in Amesos_Superlu.h

int Amesos_Superlu::ConvertToSerial() { 
  
  Epetra_RowMatrix *RowMatrixA = dynamic_cast<Epetra_RowMatrix *>(Problem_->GetOperator());
  EPETRA_CHK_ERR( RowMatrixA == 0 ) ; 

  Epetra_CrsMatrix *CastCrsMatrixA = dynamic_cast<Epetra_CrsMatrix*>(RowMatrixA) ; 
  EPETRA_CHK_ERR( CastCrsMatrixA == 0 ) ; 

  iam_ = Comm().MyPID() ;

  const Epetra_Map &OriginalMap = CastCrsMatrixA->RowMap() ; 

  NumGlobalElements_ = CastCrsMatrixA->NumGlobalRows();
  numentries_ = CastCrsMatrixA->NumGlobalNonzeros();
  assert( NumGlobalElements_ == CastCrsMatrixA->NumGlobalCols() );

  //
  //  Create a serial matrix 
  //
  assert( NumGlobalElements_ == OriginalMap.NumGlobalElements() ) ;
  int NumMyElements_ = 0 ;
  if (iam_==0) NumMyElements_ = NumGlobalElements_;


  IsLocal_ = ( OriginalMap.NumMyElements() == 
	       OriginalMap.NumGlobalElements() )?1:0;
  Comm().Broadcast( &IsLocal_, 1, 0 ) ; 

  //
  //  Convert Original Matrix to Serial (if it is not already) 
  //
  if (SerialMap_) { delete SerialMap_ ; SerialMap_ = 0 ; } 
  if ( SerialCrsMatrixA_ ) { delete SerialCrsMatrixA_ ; SerialCrsMatrixA_ = 0 ; } 
  if ( false && IsLocal_==1 ) {  //FIXME - remove the false here
     SerialMatrix_ = CastCrsMatrixA ;
  } else {
    assert( SerialMap_ == 0 ) ; 
    SerialMap_ = new Epetra_Map( NumGlobalElements_, NumMyElements_, 0, Comm() );

    Epetra_Export export_to_serial( OriginalMap, *SerialMap_);

    if ( SerialCrsMatrixA_ ) delete SerialCrsMatrixA_ ; 
    SerialCrsMatrixA_ = new Epetra_CrsMatrix(Copy, *SerialMap_, 0);
    SerialCrsMatrixA_->Export( *RowMatrixA, export_to_serial, Add ); 
    
    double zero = 0.0;
    SerialCrsMatrixA_->SetTracebackMode(0);
    for ( int i = 0 ; i < NumGlobalElements_; i++ ) {
      if ( SerialCrsMatrixA_->LRID(i) >= 0 ) 
	SerialCrsMatrixA_->InsertGlobalValues( i, 1, &zero, &i ) ;
      SerialCrsMatrixA_->SetTracebackMode(1);
    }

    SerialCrsMatrixA_->FillComplete();
    SerialMatrix_ = SerialCrsMatrixA_ ;
  }
  
  return 0;
}

int Amesos_Superlu::ReFactor(){
  
  using namespace SLU;
  //
  //  Convert matrix to the form that Superlu expects (Ap, Ai, Aval) 
  //

  if ( iam_==0 ) {
    assert( NumGlobalElements_ == SerialMatrix_->NumGlobalRows());
    assert( NumGlobalElements_ == SerialMatrix_->NumGlobalCols());
    assert( NumGlobalElements_ == SerialMatrix_->NumMyRows());
    assert( NumGlobalElements_ == SerialMatrix_->NumMyCols());

    //    assert( numentries_ == SerialMatrix_->NumGlobalNonzeros());
    numentries_ = SerialMatrix_->NumGlobalNonzeros();
    Ap_.resize( NumGlobalElements_+1 );
    Ai_.resize( EPETRA_MAX( NumGlobalElements_, numentries_) ) ; 
    Aval_.resize( EPETRA_MAX( NumGlobalElements_, numentries_) ) ; 


    int NzThisRow ;
    int Ai_index = 0 ; 
    int MyRow;
    double *RowValues;
    int *ColIndices;
    int MaxNumEntries_ = SerialMatrix_->MaxNumEntries();

    Epetra_CrsMatrix *SuperluCrs = dynamic_cast<Epetra_CrsMatrix *>(SerialMatrix_);
    if ( SuperluCrs == 0 ) {
      ColIndicesV_.resize(MaxNumEntries_);
      RowValuesV_.resize(MaxNumEntries_);
    }
    for ( MyRow = 0; MyRow < NumGlobalElements_ ; MyRow++ ) {
      if ( SuperluCrs != 0 ) {
	EPETRA_CHK_ERR( SuperluCrs->
			ExtractMyRowView( MyRow, NzThisRow, RowValues, 
					  ColIndices ) != 0 ) ;
      }
      else {
	EPETRA_CHK_ERR( SerialMatrix_->
			ExtractMyRowCopy( MyRow, MaxNumEntries_, 
					  NzThisRow, &RowValuesV_[0], 
					  &ColIndicesV_[0] ) != 0 );
	RowValues =  &RowValuesV_[0];
	ColIndices = &ColIndicesV_[0];
      }
      if ( Ap_[MyRow] != Ai_index ) EPETRA_CHK_ERR(-4);
      for ( int j = 0; j < NzThisRow; j++ ) { 
	assert( Ai_[Ai_index] == ColIndices[j] ) ;   // FIXME this may not work.   
	Aval_[Ai_index] = RowValues[j] ; 
	Ai_index++;
      }
    }
    assert( NumGlobalElements_ == MyRow );
    Ap_[ NumGlobalElements_ ] = Ai_index ; 
    
    
    assert ( FactorizationDone_ ) ; 
    Destroy_SuperMatrix_Store(&data_->A);
    Destroy_SuperNode_Matrix(&data_->L);
    Destroy_CompCol_Matrix(&data_->U);

    /* Create matrix A in the format expected by SuperLU. */
    dCreate_CompCol_Matrix( &(data_->A), NumGlobalElements_, NumGlobalElements_,
			    numentries_, &Aval_[0],
			    &Ai_[0], &Ap_[0], SLU_NR, SLU_D, SLU_GE );
  }
  return 0;
}
//
//  See also pre and post conditions in Amesos_Superlu.h
//  Preconditions:  
//    SerialMatrix_ points to the matrix to be factored and solved
//    NumGlobalElements_ has been set to the dimension of the matrix
//    numentries_ has been set to the number of non-zeros in the matrix
//      (i.e. ConvertToSerial() has been called) 
//    FactorizationDone_ and FactorizationOK_  must be accurate
//
//  Postconditions:
//    data->A 
//    FactorizationDone_ = true
//
//  Issues:
//
//
int Amesos_Superlu::Factor(){
  
  using namespace SLU;
  //
  //  Convert matrix to the form that Superlu expects (Ap, Ai, Aval) 
  //

  if ( iam_==0 ) {
    assert( NumGlobalElements_ == SerialMatrix_->NumGlobalRows());
    assert( NumGlobalElements_ == SerialMatrix_->NumGlobalCols());
    assert( NumGlobalElements_ == SerialMatrix_->NumMyRows());
    assert( NumGlobalElements_ == SerialMatrix_->NumMyCols());

    //    assert( numentries_ == SerialMatrix_->NumGlobalNonzeros());
    numentries_ = SerialMatrix_->NumGlobalNonzeros();
#ifdef NOVEC
    Ap_= new int[ NumGlobalElements_+1 ];
    int EMAX = EPETRA_MAX( NumGlobalElements_, numentries_) ; 
    Ai_= new int[ EMAX ];
    Aval_= new double[ EMAX ] ; 
#else
    Ap_.resize( NumGlobalElements_+1 );
    Ai_.resize( EPETRA_MAX( NumGlobalElements_, numentries_) ) ; 
    Aval_.resize( EPETRA_MAX( NumGlobalElements_, numentries_) ) ; 
#endif

    int NzThisRow ;
    int Ai_index = 0 ; 
    int MyRow;
    double *RowValues;
    int *ColIndices;
    int MaxNumEntries_ = SerialMatrix_->MaxNumEntries();

    Epetra_CrsMatrix *SuperluCrs = dynamic_cast<Epetra_CrsMatrix *>(SerialMatrix_);
    if ( SuperluCrs == 0 ) {
#ifdef NOVEC
      ColIndicesV_= new int[MaxNumEntries_];
      RowValuesV_= new double[MaxNumEntries_];
#else
      ColIndicesV_.resize(MaxNumEntries_);
      RowValuesV_.resize(MaxNumEntries_);
#endif
    }
    for ( MyRow = 0; MyRow < NumGlobalElements_ ; MyRow++ ) {
      if ( SuperluCrs != 0 ) {
	EPETRA_CHK_ERR( SuperluCrs->
			ExtractMyRowView( MyRow, NzThisRow, RowValues, 
					  ColIndices ) != 0 ) ;
      }
      else {
	EPETRA_CHK_ERR( SerialMatrix_->
			ExtractMyRowCopy( MyRow, MaxNumEntries_, 
					  NzThisRow, &RowValuesV_[0], 
					  &ColIndicesV_[0] ) != 0 );
	RowValues =  &RowValuesV_[0];
	ColIndices = &ColIndicesV_[0];
      }
      Ap_[MyRow] = Ai_index ; 
      for ( int j = 0; j < NzThisRow; j++ ) { 
	Ai_[Ai_index] = ColIndices[j] ; 
	Aval_[Ai_index] = RowValues[j] ; 
	Ai_index++;
      }
    }
    assert( NumGlobalElements_ == MyRow );
    Ap_[ NumGlobalElements_ ] = Ai_index ; 

    if ( FactorizationDone_ ) { 
      Destroy_SuperMatrix_Store(&data_->A);
      Destroy_SuperNode_Matrix(&data_->L);
      Destroy_CompCol_Matrix(&data_->U);
    }
      
    /* Create matrix A in the format expected by SuperLU. */
    dCreate_CompCol_Matrix( &(data_->A), NumGlobalElements_, NumGlobalElements_,
			    numentries_, &Aval_[0],
			    &Ai_[0], &Ap_[0], SLU_NR, SLU_D, SLU_GE );
  }
  return 0;
}   

int Amesos_Superlu::SetParameters( Teuchos::ParameterList &ParameterList ) {

  //  Some compilers reject the following cast:
  //  if( &ParameterList == 0 ) return 0;

  if (ParameterList.isSublist("Superlu") ) {
    Teuchos::ParameterList SuperluParams = ParameterList.sublist("Superlu") ;
  }  
  return 0;
}


bool Amesos_Superlu::MatrixShapeOK() const { 
  bool OK ;

  if ( GetProblem()->GetOperator()->OperatorRangeMap().NumGlobalPoints() != 
       GetProblem()->GetOperator()->OperatorDomainMap().NumGlobalPoints() ) OK = false;
  return OK; 
}


int Amesos_Superlu::SymbolicFactorization() {

  FactorizationOK_ = false ; 
  return 0;
}

int Amesos_Superlu::NumericFactorization() {
  
  using namespace SLU;

  Epetra_RowMatrix *RowMatrixA = dynamic_cast<Epetra_RowMatrix *>(Problem_->GetOperator());
  assert( RowMatrixA != 0 ) ;
  NumGlobalElements_ = RowMatrixA->NumGlobalRows();

#ifdef NOVEC
  perm_r_= new int[ NumGlobalElements_ ];
  perm_c_= new int[ NumGlobalElements_ ];
  etree_= new int[ NumGlobalElements_ ];
  R_= new double[ NumGlobalElements_ ];
  C_= new double[ NumGlobalElements_ ];
#else    
  perm_r_.resize( NumGlobalElements_ );
  perm_c_.resize( NumGlobalElements_ );
  etree_.resize( NumGlobalElements_ );
  R_.resize( NumGlobalElements_ );
  C_.resize( NumGlobalElements_ );
#endif    
  ConvertToSerial() ; 

  SLU::superlu_options_t& SLUopt =  data_->SLU_options ; 

  //  set_default_options( &(data_->SLU_options) ) ; 
  set_default_options( &SLUopt ) ; 
  if ( FactorizationOK_ ) {
    ReFactor() ; 
    SLUopt.Fact = data_->refactor_option ;
  }  else { 
    Factor() ; 
    FactorizationOK_ =true ;
    SLUopt.Fact = DOFACT ;
  }

  int Ierr[1];
  if ( iam_ == 0 ) { 
#if 0
    vector<double> DummyArray( NumGlobalElements_ );
    for ( int i=0; i <NumGlobalElements_; i++ ) DummyArray[i] = 1.0; // FIXME
#endif
    double rpg, rcond;
    equed_ = 'N';

#if 0
    if ( ! UseTranspose() )        // FIXME - I doubt we need this here.
      assert( SLUopt.Trans == NOTRANS ) ; 
    else
      SLUopt.Trans = TRANS ; 
    

    //    SLUopt.ColPerm  = COLAMD ;

    cout << " SLUopt.ColPerm  = " << SLUopt.ColPerm  << endl ; 
    cout << " SLUopt.Equil  = " << SLUopt.Equil  << endl ; 
    cout << " SLUopt.Fact  = " << SLUopt.Fact  << endl ; 
    cout << " SLUopt.IterRefine  = " << SLUopt.IterRefine  << endl ; 
    cout << " data_->A.Stype  = " << data_->A.Stype  
	 << " SLU_NC = " << SLU_NC 
	 << " SLU_NR = " << SLU_NR 
	 << endl ; 
    cout << " SLUopt.ColPerm  = " << SLUopt.ColPerm  << endl ; 
#endif

#ifdef NOVEC
    ferr_ = new double[1];
    berr_ = new double[1];
#else
    if ( 1 > ferr_.size() ) {
      ferr_.resize( 1 ) ; 
      berr_.resize( 1 ) ; 
    }
#endif


    data_->B.nrow = NumGlobalElements_; 
    data_->B.ncol = 0;
    SLU::DNformat* Bstore = (SLU::DNformat *) (data_->B.Store) ; 
    Bstore->lda = NumGlobalElements_; 
    Bstore->nzval = &DummyArray[0];
    data_->X.nrow = NumGlobalElements_; 
    data_->X.ncol = 0;
    SLU::DNformat* Xstore = (SLU::DNformat *) (data_->X.Store) ; 
    Xstore->lda = NumGlobalElements_; 
    Xstore->nzval = &DummyArray[0];

    SuperLUStat_t SLU_stat ;
    StatInit( &SLU_stat ) ; 
    assert( SLUopt.Fact == DOFACT); 
    dgssvx( &(SLUopt), &(data_->A), 
	    &perm_c_[0], &perm_r_[0], &etree_[0], &equed_, &R_[0], 
	    &C_[0], &(data_->L), &(data_->U), NULL, 0, 
	    &(data_->B), &(data_->X), &rpg, &rcond, &ferr_[0], 
	    &berr_[0], &(data_->mem_usage), &SLU_stat,
	    &Ierr[0] );
    StatFree( &SLU_stat ) ; 
#if 0
    if (m<200) dPrint_CompCol_Matrix("A", &(data_->A));
    if (m<200) dPrint_CompCol_Matrix("U", &(data_->U));
    if (m<200) dPrint_SuperNode_Matrix("L", &(data_->L));
#endif


    FactorizationDone_ = true ; 
  }
  

  return 0;
}


int Amesos_Superlu::Solve() { 

  using namespace SLU;

  Epetra_MultiVector   *vecX = Problem_->GetLHS() ; 
  Epetra_MultiVector   *vecB = Problem_->GetRHS() ; 
  int Ierr[1];

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

  Epetra_MultiVector *SerialB =0;
  Epetra_MultiVector *SerialX =0;
  //
  //  Extract Serial versions of X and B 
  //
  double *SerialXvalues ;
  double *SerialBvalues ;

  Epetra_RowMatrix *RowMatrixA = dynamic_cast<Epetra_RowMatrix *>(Problem_->GetOperator());
  Epetra_CrsMatrix *CastCrsMatrixA = dynamic_cast<Epetra_CrsMatrix*>(RowMatrixA) ; 
  Epetra_MultiVector *SerialXextract = 0;
  Epetra_MultiVector *SerialBextract = 0;
    
  Epetra_MultiVector CopyB = *vecB ; 

  //
  //  Copy B to the serial version of B
  //
  if ( IsLocal_ ==1 ) { 
    SerialB = &CopyB ; 
    SerialX = vecX ; 
  } else { 
    assert( IsLocal_ == 0 ) ;
    const Epetra_Map &OriginalMap = CastCrsMatrixA->RowMap();
    Epetra_MultiVector *SerialXextract = new Epetra_MultiVector( *SerialMap_, nrhs ) ; 
    Epetra_MultiVector *SerialBextract = new Epetra_MultiVector( *SerialMap_, nrhs ) ; 

    Epetra_Import ImportToSerial( *SerialMap_, OriginalMap );
    SerialBextract->Import( *vecB, ImportToSerial, Insert ) ;
    SerialB = SerialBextract ; 
    SerialX = SerialXextract ; 
  } 

  //
  //  Call SUPERLU's dgssvx to perform the solve
  //

  
  int SerialXlda ; 
  int SerialBlda ; 

  if ( iam_ == 0 ) {
    assert( SerialX->ExtractView( &SerialXvalues, &SerialXlda ) == 0 ) ; 
    assert( SerialXlda == NumGlobalElements_ ) ; 

    assert( SerialB->ExtractView( &SerialBvalues, &SerialBlda ) == 0 ) ; 
    assert( SerialBlda == NumGlobalElements_ ) ; 
    
    SLU::SuperMatrix& dataX = (data_->X) ;
    dataX.nrow =   NumGlobalElements_; 
    dataX.ncol =   nrhs ; 
    data_->X.nrow = NumGlobalElements_; 
    SLU::DNformat* Xstore = (SLU::DNformat *) (data_->X.Store) ; 
    Xstore->lda = SerialXlda; 
    Xstore->nzval = &SerialXvalues[0]; 

    SLU::SuperMatrix& dataB = (data_->B) ;
    dataB.nrow =   NumGlobalElements_; 
    dataB.ncol =   nrhs ; 
    data_->B.nrow = NumGlobalElements_; 
    SLU::DNformat* Bstore = (SLU::DNformat *) (data_->B.Store) ; 
    Bstore->lda = SerialBlda; 
    Bstore->nzval = &SerialBvalues[0]; 



#if 0
    dCreate_Dense_Matrix( &(data_->X), 
			  NumGlobalElements_, 
			  nrhs, 
			  &SerialXvalues[0], 
			  SerialXlda, 
			  SLU_DN, SLU_D, SLU_GE);

    dCreate_Dense_Matrix( &(data_->B), 
			  NumGlobalElements_, 
			  nrhs, 
			  &SerialBvalues[0], 
			  SerialBlda, 
			  SLU_DN, SLU_D, SLU_GE);
#endif

    double rpg, rcond;
#if 0
    char fact, trans, refact;
    fact = 'F';
    refact = 'N';
    trans = 'N' ;    //  This will allow us to try trans and not trans - see if either works.
    //    equed = 'N';
#endif
#ifdef NOVEC
    ferr_ = new double[ nrhs ] ; 
    berr_ = new double[ nrhs ] ; 
#else
    assert( berr_.size() == ferr_.size());
    if ( true || nrhs > ferr_.size() ) {
      ferr_.resize( nrhs ) ; 
      berr_.resize( nrhs ) ; 
    }
#endif
#if 0
    dgssvx( &fact, &trans, &refact, &(data_->A), &(data_->iparam), &perm_c_[0],
	    &perm_r_[0], &etree_[0], &equed_, &R_[0], &C_[0], &(data_->L), &(data_->U),
	    NULL, 0, &(data_->B), &(data_->X), &rpg, &rcond,
	    &ferr_[0], &berr_[0], &(data_->mem_usage), &Ierr[0] );
#endif

    SuperLUStat_t SLU_stat ;
    StatInit( &SLU_stat ) ;//    Copy the scheme used in dgssvx1.c 
    data_->SLU_options.Fact = FACTORED ; 
    SLU::superlu_options_t& SLUopt =  data_->SLU_options ; 
    if ( UseTranspose() ) 
      SLUopt.Trans = TRANS; 
    else
      assert( SLUopt.Trans == NOTRANS ); 


    //#ifdef USE_DGSTRF
    //    assert( equed_ == 'N' ) ; 
    dgssvx( &(SLUopt), &(data_->A), 
	    &perm_c_[0], &perm_r_[0], &etree_[0], &equed_, &R_[0], 
	    &C_[0], &(data_->L), &(data_->U), NULL, 0, 
	    &(data_->B), &(data_->X), &rpg, &rcond, &ferr_[0], 
	    &berr_[0], &(data_->mem_usage), &SLU_stat,
	    &Ierr[0] );
    //    assert( equed_ == 'N' ) ; 
    StatFree( &SLU_stat ) ; 

#if 0
    CopyB.Update( -1.0, *SerialB, 1.0 ) ; 
    double difference;
    CopyB.Norm2( &difference ); 
    assert( difference == 0.0 ) ; 
#endif


  }
  //
  //  Copy X back to the original vector
  // 

  if ( IsLocal_ == 0 ) { 
    const Epetra_Map &OriginalMap = CastCrsMatrixA->RowMap() ; 
    Epetra_Import ImportFromSerial( OriginalMap, *SerialMap_ );
    vecX->Import( *SerialX, ImportFromSerial, Insert ) ;
    delete SerialBextract ;
    delete SerialXextract ;
  } else {
    assert( SerialBextract == 0 ) ;
    assert( SerialXextract == 0 ) ;
  }  
  //  All processes should return the same error code
  if ( 1 < Comm().NumProc() ) {
    Comm().Broadcast( Ierr, 1, 0 ) ; 
  }

  return Ierr[0];
}
