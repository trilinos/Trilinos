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

#include "Amesos_Klu.h"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_Util.h"
#include "Amesos_Support.h"

extern "C" {
  // #include "amd.h"
#include "klu.h"
#if 0
#include "klu_dump.h"
#endif
}

namespace {

using Teuchos::RefCountPtr;

template<class T, class DeleteFunctor>
class DeallocFunctorDeleteWithCommon
{
public:
  DeallocFunctorDeleteWithCommon(
				 const RefCountPtr<klu_common>  &common
				 ,DeleteFunctor                        deleteFunctor
				 )
    : common_(common), deleteFunctor_(deleteFunctor)
    {}
  typedef T ptr_t;
  void free( T* ptr ) {
    if(ptr) deleteFunctor_(&ptr,&*common_);
  }
private:
  Teuchos::RefCountPtr<klu_common> common_;
  DeleteFunctor deleteFunctor_;
  DeallocFunctorDeleteWithCommon(); // Not defined and not to be called!
};

template<class T, class DeleteFunctor>
DeallocFunctorDeleteWithCommon<T,DeleteFunctor>
deallocFunctorDeleteWithCommon(
			       const RefCountPtr<klu_common>  &common
			       ,DeleteFunctor                        deleteFunctor
			       )
{
  return DeallocFunctorDeleteWithCommon<T,DeleteFunctor>(common,deleteFunctor);
}


} // namespace 

class Amesos_Klu_Pimpl 
{
public:
  Teuchos::RefCountPtr<klu_symbolic> Symbolic_;
  Teuchos::RefCountPtr<klu_numeric> Numeric_;
  Teuchos::RefCountPtr<klu_common> common_;
};

//=============================================================================
Amesos_Klu::Amesos_Klu(const Epetra_LinearProblem &prob ) :
  PrivateKluData_( rcp( new Amesos_Klu_Pimpl() ) ),
  CrsMatrixA_(0),
  TrustMe_(false),
  UseTranspose_(false),
  Problem_(&prob)
{
  // MS // move declaration of Problem_ above because I need it
  // MS // set up before calling Comm()
  Teuchos::ParameterList ParamList ;
  SetParameters( ParamList ) ;

}

//=============================================================================
Amesos_Klu::~Amesos_Klu(void) {

  // print out some information if required by the user
  if( (verbose_ && PrintTiming_) || verbose_ == 2 ) PrintTiming();
  if( (verbose_ && PrintStatus_) || verbose_ == 2 ) PrintStatus();
}

//=============================================================================
int Amesos_Klu::ExportToSerial() 
{

  if ( AddZeroToDiag_==0 && numentries_ != RowMatrixA_->NumGlobalNonzeros()) { 
    cerr << " The number of non zero entries in the matrix has changed since the last call to SymbolicFactorization().  " ;
    AMESOS_CHK_ERR( -2 );
  }
  if (UseDataInPlace_ != 1) {
    assert ( RowMatrixA_ != 0 ) ; 
    if (  AddZeroToDiag_==0 && numentries_ != RowMatrixA_->NumGlobalNonzeros()) { 
      cerr << " The number of non zero entries in the matrix has changed since the last call to SymbolicFactorization().  " ;
      AMESOS_CHK_ERR( -3 );
    }
    assert ( ImportToSerial_.get() != 0 ) ; 
    AMESOS_CHK_ERR(SerialCrsMatrixA_->Import(*StdIndexMatrix_, 
					     *ImportToSerial_, Insert ));

    if (AddZeroToDiag_ ) { 
      int OriginalTracebackMode = SerialCrsMatrixA_->GetTracebackMode() ; 
      SerialCrsMatrixA_->SetTracebackMode( EPETRA_MIN( OriginalTracebackMode, 0) ) ; // ExportToSerial is called both by PerformSymbolicFactorization() and PerformNumericFactorization().  When called by the latter, the call to insertglobalvalues is both unnecessary and illegal.  Fortunately, Epetra allows us to ignore the error message by setting the traceback mode to 0.

      //
      //  Add 0.0 to each diagonal entry to avoid empty diagonal entries;
      //
      double zero = 0.0;
      for ( int i = 0 ; i < SerialMap_->NumGlobalElements(); i++ ) 
	if ( SerialCrsMatrixA_->LRID(i) >= 0 ) 
	  SerialCrsMatrixA_->InsertGlobalValues( i, 1, &zero, &i ) ;
      SerialCrsMatrixA_->SetTracebackMode( OriginalTracebackMode ) ; 
    }
    AMESOS_CHK_ERR(SerialCrsMatrixA_->FillComplete());
    AMESOS_CHK_ERR(SerialCrsMatrixA_->OptimizeStorage());

    if( ! AddZeroToDiag_ &&  numentries_ != SerialMatrix_->NumGlobalNonzeros()) {
      cerr << " Amesos_Klu cannot handle this matrix.  " ;
      if ( Reindex_ ) {
	cerr << "Unknown error" << endl ; 
	AMESOS_CHK_ERR( -4 );
      } else {
	cerr << " Try setting the Reindex parameter to true. " << endl ; 
#ifndef HAVE_AMESOS_EPETRAEXT
	cerr << " You will need to rebuild the Amesos library with the EpetraExt library to use the reindex feature" << endl ; 
	cerr << " To rebuild Amesos with EpetraExt, add --enable-epetraext to your configure invocation" << endl ;
#endif
	AMESOS_CHK_ERR( -5 );
      }
    }
  }
  
  return 0;
}
//=============================================================================
//
// CreateLocalMatrixAndExporters() is called only by SymbolicFactorization() 
// for CrsMatrix objects.  All objects should be created here.  No assumptions 
// are made about the input operator.  I.e. it can be completely different from 
// the matrix at the time of the previous call to SymbolicFactorization(). 
//

int Amesos_Klu::CreateLocalMatrixAndExporters() 
{
  ResetTime(0);

  RowMatrixA_ = Problem_->GetMatrix(); // MS, 18-Apr-06 //
  if (RowMatrixA_ == 0) AMESOS_CHK_ERR(-1);

  const Epetra_Map &OriginalMatrixMap = RowMatrixA_->RowMatrixRowMap() ;
  const Epetra_Map &OriginalDomainMap = 
    UseTranspose()?GetProblem()->GetOperator()->OperatorRangeMap():
    GetProblem()->GetOperator()->OperatorDomainMap();
  const Epetra_Map &OriginalRangeMap = 
    UseTranspose()?GetProblem()->GetOperator()->OperatorDomainMap():
    GetProblem()->GetOperator()->OperatorRangeMap();

  NumGlobalElements_ = RowMatrixA_->NumGlobalRows();
  numentries_ = RowMatrixA_->NumGlobalNonzeros();
  assert( NumGlobalElements_ == RowMatrixA_->NumGlobalCols() );

  //
  //  Create a serial matrix
  //
  assert( NumGlobalElements_ == OriginalMatrixMap.NumGlobalElements() ) ;
  int NumMyElements_ = 0 ;
  if (MyPID_==0) NumMyElements_ = NumGlobalElements_;
  //
  //  UseDataInPlace_ is set to 1 (true) only if everything is perfectly
  //  normal.  Anything out of the ordinary reverts to the more expensive
  //  path. 
  //
  UseDataInPlace_ = ( OriginalMatrixMap.NumMyElements() ==
	       OriginalMatrixMap.NumGlobalElements() )?1:0;
  if ( ! OriginalRangeMap.SameAs( OriginalMatrixMap ) ) UseDataInPlace_ = 0 ; 
  if ( ! OriginalDomainMap.SameAs( OriginalMatrixMap ) ) UseDataInPlace_ = 0 ; 
  if ( AddZeroToDiag_ ) UseDataInPlace_ = 0 ; 
  Comm().Broadcast( &UseDataInPlace_, 1, 0 ) ;

  //
  //  Reindex matrix if necessary (and possible - i.e. CrsMatrix)
  //
  //  For now, since I don't know how to determine if we need to reindex the matrix, 
  //  we allow the user to control re-indexing through the "Reindex" parameter.
  //
  CrsMatrixA_ = dynamic_cast<Epetra_CrsMatrix *>(Problem_->GetOperator());

  if  ( Reindex_ ) {
    if ( CrsMatrixA_ == 0 ) AMESOS_CHK_ERR(-7);
#ifdef HAVE_AMESOS_EPETRAEXT
    const Epetra_Map& OriginalMap = CrsMatrixA_->RowMap();
    StdIndex_ = rcp( new Amesos_StandardIndex( OriginalMap  ) );
    //const Epetra_Map& OriginalColMap = CrsMatrixA_->RowMap();
    StdIndexDomain_ = rcp( new Amesos_StandardIndex( OriginalDomainMap  ) );
    StdIndexRange_ = rcp( new Amesos_StandardIndex( OriginalRangeMap  ) );

    StdIndexMatrix_ = StdIndex_->StandardizeIndex( CrsMatrixA_ );
#else
    cerr << "Amesos_Klu requires EpetraExt to reindex matrices." << endl 
	 <<  " Please rebuild with the EpetraExt library by adding --enable-epetraext to your configure invocation" << endl ;
    AMESOS_CHK_ERR(-8);
#endif
  } else { 
    StdIndexMatrix_ = RowMatrixA_ ;
  }

  //
  //  At this point, StdIndexMatrix_ points to a matrix with 
  //  standard indexing.  
  //

  //
  //  Convert Original Matrix to Serial (if it is not already)
  //
  if (UseDataInPlace_ == 1) {
    SerialMatrix_ = StdIndexMatrix_;
  } else {
    SerialMap_ = rcp(new Epetra_Map(NumGlobalElements_, NumMyElements_, 0, Comm()));
    
    if ( Problem_->GetRHS() ) 
      NumVectors_ = Problem_->GetRHS()->NumVectors() ; 
    else
      NumVectors_ = 1 ; 
    SerialXextract_ = rcp( new Epetra_MultiVector(*SerialMap_,NumVectors_));
    SerialBextract_ = rcp (new Epetra_MultiVector(*SerialMap_,NumVectors_));

    ImportToSerial_ = rcp(new Epetra_Import ( *SerialMap_, StdIndexMatrix_->RowMatrixRowMap() ) );
    if (ImportToSerial_.get() == 0) AMESOS_CHK_ERR(-9);

    SerialCrsMatrixA_ = rcp( new Epetra_CrsMatrix(Copy, *SerialMap_, 0) ) ;
    SerialMatrix_ = &*SerialCrsMatrixA_ ;
  }

  AddTime("matrix redistribution", 0);

  return(0);
}

//=============================================================================
//
//  See also pre and post conditions in Amesos_Klu.h
//  Preconditions:
//    firsttime specifies that this is the first time that 
//    ConertToKluCrs has been called, i.e. in symbolic factorization.  
//    No data allocation should happen unless firsttime=true.
//    SerialMatrix_ points to the matrix to be factored and solved
//    NumGlobalElements_ has been set to the dimension of the matrix
//    numentries_ has been set to the number of non-zeros in the matrix
//      (i.e. CreateLocalMatrixAndExporters() has been callded)
//
//  Postconditions:
//    Ap, VecAi, VecAval contain the matrix as Klu needs it
//
//
int Amesos_Klu::ConvertToKluCRS(bool firsttime)
{
  ResetTime(0);

  //
  //  Convert matrix to the form that Klu expects (Ap, VecAi, VecAval)
  //

  if (MyPID_==0) {
    assert( NumGlobalElements_ == SerialMatrix_->NumGlobalRows());
    assert( NumGlobalElements_ == SerialMatrix_->NumGlobalCols());
    if ( ! AddZeroToDiag_ ) {
      assert( numentries_ == SerialMatrix_->NumGlobalNonzeros()) ;
    } else {
      numentries_ = SerialMatrix_->NumGlobalNonzeros() ;
    }

    Epetra_CrsMatrix *CrsMatrix = dynamic_cast<Epetra_CrsMatrix *>(SerialMatrix_);
    bool StorageOptimized = ( CrsMatrix != 0 && CrsMatrix->StorageOptimized() );

    if ( AddToDiag_ != 0.0 ) StorageOptimized = false ;

    if ( firsttime ) { 
      Ap.resize( NumGlobalElements_+1 );
      if ( ! StorageOptimized ) { 
	VecAi.resize( EPETRA_MAX( NumGlobalElements_, numentries_) ) ;
	VecAval.resize( EPETRA_MAX( NumGlobalElements_, numentries_) ) ;
	Ai = &VecAi[0];
	Aval = &VecAval[0];
      }
    }

    double *RowValues;
    int *ColIndices;
    int NumEntriesThisRow;

    if( StorageOptimized ) {
      if ( firsttime ) {
	Ap[0] = 0;
	for ( int MyRow = 0; MyRow <NumGlobalElements_; MyRow++ ) {
	  if( CrsMatrix->
	      ExtractMyRowView( MyRow, NumEntriesThisRow, RowValues,
				ColIndices ) != 0 ) 
	    AMESOS_CHK_ERR( -10 );
	  if ( MyRow == 0 ) {
	    Ai = ColIndices ;
	    Aval = RowValues ;
	  }
	  Ap[MyRow+1] = Ap[MyRow] + NumEntriesThisRow ;
	}
      }
    } else { 
      
      int Ai_index = 0 ;
      int MyRow;
      
      int MaxNumEntries_ = SerialMatrix_->MaxNumEntries();
      if ( firsttime && CrsMatrix == 0 ) {
	ColIndicesV_.resize(MaxNumEntries_);
	RowValuesV_.resize(MaxNumEntries_);
      }
      
      for ( MyRow = 0; MyRow <NumGlobalElements_; MyRow++ ) {
	if ( CrsMatrix != 0 ) {
	  if( CrsMatrix->
	      ExtractMyRowView( MyRow, NumEntriesThisRow, RowValues,
				ColIndices ) != 0 ) 
	    AMESOS_CHK_ERR( -11 );

	} else {
	  if( SerialMatrix_->
			  ExtractMyRowCopy( MyRow, MaxNumEntries_,
					    NumEntriesThisRow, &RowValuesV_[0],
					    &ColIndicesV_[0] ) != 0 ) 
	    AMESOS_CHK_ERR( -12 );

	  RowValues =  &RowValuesV_[0];
	  ColIndices = &ColIndicesV_[0];
	}
	
	if ( firsttime ) {
	  Ap[MyRow] = Ai_index ;
	  for ( int j = 0; j < NumEntriesThisRow; j++ ) {
	    VecAi[Ai_index] = ColIndices[j] ;
	    //	  assert( VecAi[Ai_index] == Ai[Ai_index] ) ; 
	    VecAval[Ai_index] = RowValues[j] ;      //  We have to do this because of the hacks to get aroun bug #1502 
	    if (ColIndices[j] == MyRow) {
	      VecAval[Ai_index] += AddToDiag_;    
	    }
	    Ai_index++;
	  }
	} else { 
	  for ( int j = 0; j < NumEntriesThisRow; j++ ) {
	    VecAval[Ai_index] = RowValues[j] ;     
	    if (ColIndices[j] == MyRow) {
	      VecAval[Ai_index] += AddToDiag_;   
	    }
	    Ai_index++;
	  }
	}
      }
      Ap[MyRow] = Ai_index ;
    }

  }

  AddTime("matrix conversion", 0);

  return 0;
}

//=============================================================================
int Amesos_Klu::SetParameters( Teuchos::ParameterList &ParameterList ) {

  // ========================================= //
  // retrive KLU's parameters from list.       //
  // default values defined in the constructor //
  // ========================================= //

  // retrive general parameters

  SetStatusParameters( ParameterList );
  SetControlParameters( ParameterList );

  if (ParameterList.isParameter("TrustMe") ) 
    TrustMe_ = ParameterList.get<bool>( "TrustMe" );

#if 0

  unused for now

  if (ParameterList.isSublist("Klu") ) {
    Teuchos::ParameterList KluParams = ParameterList.sublist("Klu") ;
  }
#endif

  return 0;
}


//=============================================================================
int Amesos_Klu::PerformSymbolicFactorization() 
{
  ResetTime(0);

  if (MyPID_ == 0) {

    PrivateKluData_->common_ = rcp(new klu_common());
    klu_defaults(&*PrivateKluData_->common_) ;
    PrivateKluData_->Symbolic_ =
      rcp(
	  klu_analyze (NumGlobalElements_, &Ap[0], Ai,  &*PrivateKluData_->common_ ) 
	  ,deallocFunctorDeleteWithCommon<klu_symbolic>(PrivateKluData_->common_,klu_free_symbolic)
	  ,true
	  );

	const  bool symbolic_ok =  (PrivateKluData_->Symbolic_.get() != NULL) 
	  && PrivateKluData_->common_->status == KLU_OK ;
	if ( ! symbolic_ok ) {
	  cout << __FILE__ << "::" << __LINE__ <<
	    " NumericallySingularMatrixError MyPID_ == 0 " ;
	  AMESOS_CHK_ERR( NumericallySingularMatrixError ) ;
	}

  }

  AddTime("symbolic", 0);

  return 0;
}

//=============================================================================
int Amesos_Klu::PerformNumericFactorization( ) 
{

  // Changed this; it was "if (!TrustMe)...
  // The behavior is not intuitive. Maybe we should introduce a new
  // parameter, FastSolvers or something like that, that does not perform
  // any AddTime, ResetTime, GetTime.
  
  ResetTime(0);

  if (MyPID_ == 0) {

    bool factor_with_pivoting = true ;

    // set the default parameters
    PrivateKluData_->common_->scale = ScaleMethod_ ;

    const bool NumericNonZero =  PrivateKluData_->Numeric_.get() != 0 ; 
    
    // see if we can "refactorize"
    if ( refactorize_ && NumericNonZero ) { 
      // refactorize using the existing Symbolic and Numeric objects, and
      // using the identical pivot ordering as the prior klu_factor.
      // No partial pivoting is done.
      int result = klu_refactor (&Ap[0], Ai, Aval,
				 &*PrivateKluData_->Symbolic_, 
				 &*PrivateKluData_->Numeric_, &*PrivateKluData_->common_) ;
      // Did it work?
      const  bool refactor_ok = result == 1 && PrivateKluData_->common_->status == KLU_OK ;
      if ( refactor_ok ) { 
	
	double rcond ;
	
	klu_rcond (&*PrivateKluData_->Symbolic_,
		   &*PrivateKluData_->Numeric_,
		   &rcond, &*PrivateKluData_->common_) ;
	
	if ( rcond > rcond_threshold_ ) {
	  // factorizing without pivot worked fine.  We are done.
	  factor_with_pivoting = false ;
	}
      }
    }
    
    if ( factor_with_pivoting ) {
      
      // factor with partial pivoting:
      // either this is the first time we are factoring the matrix, or the
      // refactorize parameter is false, or we tried to refactorize and
      // found it to be too inaccurate.
      
      // factor the matrix using partial pivoting
      PrivateKluData_->Numeric_ =
	rcp( klu_factor(&Ap[0], Ai, Aval,
			&*PrivateKluData_->Symbolic_, &*PrivateKluData_->common_),
	     deallocFunctorDeleteWithCommon<klu_numeric>(PrivateKluData_->common_,klu_free_numeric)
	     ,true
	     );
      
      const  bool numeric_ok =  PrivateKluData_->Numeric_.get()!=NULL 
	&& PrivateKluData_->common_->status == KLU_OK ;
      if ( ! numeric_ok ) AMESOS_CHK_ERR( NumericallySingularMatrixError );  
      
    }
    
  }

  AddTime("numeric", 0);

  return 0;
}

//=============================================================================
bool Amesos_Klu::MatrixShapeOK() const {
  bool OK = true;

  // Comment by Tim:  The following code seems suspect.  The variable "OK"
  // is not set if the condition is true.
  // Does the variable "OK" default to true?
  if ( GetProblem()->GetOperator()->OperatorRangeMap().NumGlobalPoints() !=
       GetProblem()->GetOperator()->OperatorDomainMap().NumGlobalPoints() ) {
    OK = false;
  }
  return OK;
}

//=============================================================================
int Amesos_Klu::SymbolicFactorization() 
{
  MyPID_    = Comm().MyPID();
  NumProcs_ = Comm().NumProc();

  IsSymbolicFactorizationOK_ = false;
  IsNumericFactorizationOK_ = false;
  
  InitTime(Comm(), 2);

  ResetTime(1);

  // "overhead" time for the following method is considered here
  AMESOS_CHK_ERR( CreateLocalMatrixAndExporters() ) ;
  assert( NumGlobalElements_ == RowMatrixA_->NumGlobalCols() );

  //
  //  Perform checks in SymbolicFactorization(), but none in 
  //  NumericFactorization() or Solve()
  //
  if ( TrustMe_ ) { 
     if ( CrsMatrixA_ == 0 ) AMESOS_CHK_ERR( -15 );
     if( UseDataInPlace_ != 1 ) AMESOS_CHK_ERR( -15 ) ;
     if( Reindex_ )  AMESOS_CHK_ERR( -15 ) ;
     if( ! Problem_->GetLHS() )  AMESOS_CHK_ERR( -15 ) ;
     if( ! Problem_->GetRHS() )  AMESOS_CHK_ERR( -15 ) ;
     if( ! Problem_->GetLHS()->NumVectors() ) AMESOS_CHK_ERR( -15 ) ;
     if( ! Problem_->GetRHS()->NumVectors() ) AMESOS_CHK_ERR( -15 ) ; 
     SerialB_ = Problem_->GetRHS() ;
     SerialX_ = Problem_->GetLHS() ;
     NumVectors_ = SerialX_->NumVectors();
     if (MyPID_ == 0) {
       AMESOS_CHK_ERR(SerialX_->ExtractView(&SerialXBvalues_,&SerialXlda_ ));
       AMESOS_CHK_ERR(SerialB_->ExtractView(&SerialBvalues_,&SerialXlda_ ));
     }
  }

  // "overhead" time for the following two methods is considered here
  AMESOS_CHK_ERR( ExportToSerial() );

  AMESOS_CHK_ERR( ConvertToKluCRS(true) );

  AddTime("overhead", 1);

  // All this time if KLU time
  AMESOS_CHK_ERR( PerformSymbolicFactorization() );

  NumSymbolicFact_++;

  IsSymbolicFactorizationOK_ = true;
  
  return 0;
}

//=============================================================================
int Amesos_Klu::NumericFactorization() 
{
  if ( !TrustMe_ ) 
  { 
    IsNumericFactorizationOK_ = false;
    if (IsSymbolicFactorizationOK_ == false)
      AMESOS_CHK_ERR(SymbolicFactorization());
    
    ResetTime(1); // "overhead" time

    Epetra_CrsMatrix *CrsMatrixA = dynamic_cast<Epetra_CrsMatrix *>(RowMatrixA_);
    if ( CrsMatrixA == 0 )   // hack to get around Bug #1502 
      AMESOS_CHK_ERR( CreateLocalMatrixAndExporters() ) ;
    assert( NumGlobalElements_ == RowMatrixA_->NumGlobalCols() );
    if ( AddZeroToDiag_ == 0 ) 
      assert( numentries_ == RowMatrixA_->NumGlobalNonzeros() );

    AMESOS_CHK_ERR( ExportToSerial() );
    
    if ( CrsMatrixA == 0 ) {  // continuation of hack to avoid bug #1502
      AMESOS_CHK_ERR( ConvertToKluCRS(true) );
    }  else {
      AMESOS_CHK_ERR( ConvertToKluCRS(false) );
    }

    AddTime("overhead", 1);
  }

  // this time is all for KLU
  AMESOS_CHK_ERR( PerformNumericFactorization() );

  NumNumericFact_++;

  IsNumericFactorizationOK_ = true;
  
  return 0;
}

//=============================================================================
int Amesos_Klu::Solve() 
{
  Epetra_MultiVector* vecX = 0 ;
  Epetra_MultiVector* vecB = 0 ;


#ifdef Bug_8212
  //  This demonstrates Bug #2812 - Valgrind does not catch this
  //  memory leak
  lose_this_ = (int *) malloc( 300 ) ;

  
#ifdef Bug_8212_B
  //  This demonstrates Bug #2812 - Valgrind does catch this
  //  use of unitialized data - but only in TestOptions.exe 
  //  not in Test_Basic/amesos_test.exe 	
  //  		
    if ( lose_this_[0] == 12834 ) { 
	     cout << __FILE__ << "::"  << __LINE__ << " very unlikely to happen " << endl ; 
    }
#endif
#endif
  if ( !TrustMe_  ) { 

    SerialB_ = Problem_->GetRHS() ;
    SerialX_ = Problem_->GetLHS() ;

    Epetra_MultiVector* OrigVecX ;
    Epetra_MultiVector* OrigVecB ;

    if (IsNumericFactorizationOK_ == false)
      AMESOS_CHK_ERR(NumericFactorization());
    
    ResetTime(1);

    //
    //  Reindex the LHS and RHS 
    //
    OrigVecX = Problem_->GetLHS() ;
    OrigVecB = Problem_->GetRHS() ;
    
    if ( Reindex_ ) { 
#ifdef HAVE_AMESOS_EPETRAEXT
      vecX = StdIndexDomain_->StandardizeIndex( OrigVecX ) ;
      vecB = StdIndexRange_->StandardizeIndex( OrigVecB ) ;
#else
      AMESOS_CHK_ERR( -13 ) ; // Amesos_Klu can't handle non-standard indexing without EpetraExt 
#endif
    } else {
      vecX = OrigVecX ;
      vecB = OrigVecB ;
    } 
    
    if ((vecX == 0) || (vecB == 0))
      AMESOS_CHK_ERR(-1); // something wrong in input
    
    //  Extract Serial versions of X and B
    
    ResetTime(0);

    //  Copy B to the serial version of B
    //
    if (UseDataInPlace_ == 1) {
      SerialB_ = vecB;
      SerialX_ = vecX;
      NumVectors_ = Problem_->GetRHS()->NumVectors() ; 
    } else {
      assert (UseDataInPlace_ == 0);
      
      if( NumVectors_ != Problem_->GetRHS()->NumVectors() ) {
	NumVectors_ = Problem_->GetRHS()->NumVectors() ; 
	SerialXextract_ = rcp( new Epetra_MultiVector(*SerialMap_,NumVectors_));
	SerialBextract_ = rcp (new Epetra_MultiVector(*SerialMap_,NumVectors_));
      }
      if (NumVectors_ != vecB->NumVectors())
	AMESOS_CHK_ERR(-1); // internal error 
      
      ImportRangeToSerial_ = rcp(new Epetra_Import ( *SerialMap_, vecB->Map() ) );
      if ( SerialBextract_->Import(*vecB,*ImportRangeToSerial_,Insert) )
	AMESOS_CHK_ERR( -1 ) ; // internal error
      
      SerialB_ = &*SerialBextract_ ;
      SerialX_ = &*SerialXextract_ ;
    }
    AddTime("vector redistribution", 0);

    //  Call KLU to perform the solve
    
    ResetTime(0);
    if (MyPID_ == 0) {
      AMESOS_CHK_ERR(SerialB_->ExtractView(&SerialBvalues_,&SerialXlda_ ));
      AMESOS_CHK_ERR(SerialX_->ExtractView(&SerialXBvalues_,&SerialXlda_ ));
      if (SerialXlda_ != NumGlobalElements_)
	AMESOS_CHK_ERR(-1);
    }

    AddTime("overhead", 1);
  }
  if ( MyPID_ == 0) {
    if ( NumVectors_ == 1 ) {
      for ( int i = 0 ; i < NumGlobalElements_ ; i++ ) 
	SerialXBvalues_[i] = SerialBvalues_[i] ;
    } else {
      SerialX_->Scale(1.0, *SerialB_ ) ;    // X = B (Klu overwrites B with X)
    }
    if (UseTranspose()) {
      klu_solve( &*PrivateKluData_->Symbolic_, &*PrivateKluData_->Numeric_,
		     SerialXlda_, NumVectors_, &SerialXBvalues_[0], &*PrivateKluData_->common_ );
    } else {
      klu_tsolve( &*PrivateKluData_->Symbolic_, &*PrivateKluData_->Numeric_,
		      SerialXlda_, NumVectors_, &SerialXBvalues_[0], &*PrivateKluData_->common_ );
    }
  }

  AddTime("solve", 0);

  //  Copy X back to the original vector

  ResetTime(0);
  ResetTime(1);

  if (UseDataInPlace_ == 0) {
    ImportDomainToSerial_ = rcp(new Epetra_Import ( *SerialMap_, vecX->Map() ) );
    vecX->Export( *SerialX_, *ImportDomainToSerial_, Insert ) ;

  } // otherwise we are already in place.

  AddTime("vector redistribution", 0);

#if 0
  //
  //  ComputeTrueResidual causes TestOptions to fail on my linux box 
  //  Bug #1417
  if (ComputeTrueResidual_)
    ComputeTrueResidual(*SerialMatrix_, *vecX, *vecB, UseTranspose(), "Amesos_Klu");
#endif

  if (ComputeVectorNorms_)
    ComputeVectorNorms(*vecX, *vecB, "Amesos_Klu");

  AddTime("overhead", 1);

  ++NumSolve_;


  return(0) ;
}

// ================================================ ====== ==== ==== == =

void Amesos_Klu::PrintStatus() const
{

  if (MyPID_) return;

  PrintLine();

  cout << "Amesos_Klu : Matrix has " << NumGlobalElements_ << " rows"
       << " and " << numentries_ << " nonzeros" << endl;
  cout << "Amesos_Klu : Nonzero elements per row = "
       << 1.0*numentries_/NumGlobalElements_ << endl;
  cout << "Amesos_Klu : Percentage of nonzero elements = "
       << 100.0*numentries_/(pow(double(NumGlobalElements_),double(2.0))) << endl;
  cout << "Amesos_Klu : Use transpose = " << UseTranspose_ << endl;

  PrintLine();

  return;

}

// ================================================ ====== ==== ==== == =

void Amesos_Klu::PrintTiming() const
{
  if (MyPID_) return;

  double ConTime = GetTime("matrix conversion");
  double MatTime = GetTime("matrix redistribution");
  double VecTime = GetTime("vector redistribution");
  double SymTime = GetTime("symbolic");
  double NumTime = GetTime("numeric");
  double SolTime = GetTime("solve");
  double OveTime = GetTime("overhead");

  if (NumSymbolicFact_)
    SymTime /= NumSymbolicFact_;

  if (NumNumericFact_)
    NumTime /= NumNumericFact_;

  if (NumSolve_)
    SolTime /= NumSolve_;

  string p = "Amesos_Klu : ";
  PrintLine();

  cout << p << "Time to convert matrix to Klu format = "
       << ConTime << " (s)" << endl;
  cout << p << "Time to redistribute matrix = "
       << MatTime << " (s)" << endl;
  cout << p << "Time to redistribute vectors = "
       << VecTime << " (s)" << endl;
  cout << p << "Number of symbolic factorizations = "
       << NumSymbolicFact_ << endl;
  cout << p << "Time for sym fact = "
       << SymTime * NumSymbolicFact_ << " (s), avg = " << SymTime << " (s)" << endl;
  cout << p << "Number of numeric factorizations = "
       << NumNumericFact_ << endl;
  cout << p << "Time for num fact = "
       << NumTime * NumNumericFact_ << " (s), avg = " << NumTime << " (s)" << endl;
  cout << p << "Number of solve phases = "
       << NumSolve_ << endl;
  cout << p << "Time for solve = "
       << SolTime * NumSolve_ << " (s), avg = " << SolTime << " (s)" << endl;

  double tt = SymTime * NumSymbolicFact_ + NumTime * NumNumericFact_ + SolTime * NumSolve_;
  if (tt != 0)
  {
    cout << p << "Total time spent in Amesos = " << tt << " (s) " << endl;
    cout << p << "Total time spent in the Amesos interface = " << OveTime << " (s)" << endl;
    cout << p << "(the above time does not include KLU time)" << endl;
    cout << p << "Amesos interface time / total time = " << OveTime / tt << endl;
  }

  PrintLine();

  return;
}
