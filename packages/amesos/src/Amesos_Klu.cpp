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
#include "klu_btf.h"
#if 0
#include "klu_dump.h"
#endif
}

class Amesos_Klu_Pimpl {
public:
   klu_symbolic *Symbolic_ ;
   klu_numeric *Numeric_ ;

  Amesos_Klu_Pimpl():
    Symbolic_(0),
    Numeric_(0)
  {}

  ~Amesos_Klu_Pimpl(void){

    if ( Symbolic_ ) klu_btf_free_symbolic (&Symbolic_) ;
    if ( Numeric_ ) klu_btf_free_numeric (&Numeric_) ;
  }

} ;




//=============================================================================
Amesos_Klu::Amesos_Klu(const Epetra_LinearProblem &prob ) :
  PrivateKluData_( new Amesos_Klu_Pimpl() ),
  CrsMatrixA_(0),
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

  delete PrivateKluData_;

  // print out some information if required by the user
  if( (verbose_ && PrintTiming_) || verbose_ == 2 ) PrintTiming();
  if( (verbose_ && PrintStatus_) || verbose_ == 2 ) PrintStatus();
}

//=============================================================================
int Amesos_Klu::ExportToSerial() 
{
  if ( debug_ ) cout << __FILE__ << "::" << __LINE__ 
       << " UseDataInPlace_ = " << UseDataInPlace_ 
       << " iam = " << iam 
       << endl ; 

  if ( numentries_ != RowMatrixA_->NumGlobalNonzeros()) { 
    cerr << " The number of non zero entries in the matrix has changed since the last call to SymbolicFactorization().  " ;
    AMESOS_CHK_ERR( -2 );
  }
  if (UseDataInPlace_ != 1) {
    if ( debug_ ) cout << __FILE__ << "::" << __LINE__ << " iam = " << iam << endl ; 
    assert ( RowMatrixA_ != 0 ) ; 
    assert ( ImportToSerial_.get() != 0 ) ; 
    if ( debug_ ) cout << __FILE__ << "::" << __LINE__ << " iam = " << iam << endl ; 
    AMESOS_CHK_ERR(SerialCrsMatrixA_->Import(*StdIndexMatrix_, 
					     *ImportToSerial_, Insert ));

    Comm().Barrier();
    if ( debug_ ) cout << __FILE__ << "::" << __LINE__ << " iam = " << iam << endl ; 
    
    AMESOS_CHK_ERR(SerialCrsMatrixA_->FillComplete());
    if ( debug_ ) cout << __FILE__ << "::" << __LINE__ << " iam = " << iam << endl ; 

    if( numentries_ != SerialMatrix_->NumGlobalNonzeros()) {
      cerr << " Amesos_Klu cannot handle this matrix.  " ;
      if ( Reindex_ ) {
	cerr << "Unknown error" << endl ; 
	AMESOS_CHK_ERR( -5 );
      } else {
	cerr << " Try setting the Reindex parameter to true. " << endl ; 
	AMESOS_CHK_ERR( -3 );
      }
    }

  }
  
  return 0;
}
//=============================================================================
int Amesos_Klu::CreateLocalMatrixAndExporters() 
{

  ResetTime();

  RowMatrixA_ = dynamic_cast<Epetra_RowMatrix *>(Problem_->GetOperator());
  if (RowMatrixA_ == 0) AMESOS_CHK_ERR(-1);
  iam = Comm().MyPID() ;

  const Epetra_Map &OriginalMatrixMap = RowMatrixA_->RowMatrixRowMap() ;
  const Epetra_Map &OriginalDomainMap = 
    (UseTranspose())?GetProblem()->GetOperator()->OperatorRangeMap():GetProblem()->GetOperator()->OperatorDomainMap();
  const Epetra_Map &OriginalRangeMap = 
    UseTranspose()?GetProblem()->GetOperator()->OperatorDomainMap():
    GetProblem()->GetOperator()->OperatorRangeMap();

  NumGlobalElements_ = RowMatrixA_->NumGlobalRows();
  numentries_ = RowMatrixA_->NumGlobalNonzeros();
  assert( NumGlobalElements_ == RowMatrixA_->NumGlobalCols() );

  //
  //  Create a serial matrix
  //
  if ( debug_ ) cout << __FILE__ << "::" << __LINE__ << " iam = " << iam << endl ; 
  assert( NumGlobalElements_ == OriginalMatrixMap.NumGlobalElements() ) ;
  int NumMyElements_ = 0 ;
  if (iam==0) NumMyElements_ = NumGlobalElements_;
  //
  //  UseDataInPlace_ is set to 1 (true) only if everything is perfectly
  //  normal.  Anything out of the ordinary reverts to the more expensive
  //  path. 
  //
  if ( debug_ ) cout << __FILE__ << "::" << __LINE__ << " iam = " << iam << endl ; 
  UseDataInPlace_ = ( OriginalMatrixMap.NumMyElements() ==
	       OriginalMatrixMap.NumGlobalElements() )?1:0;
  if ( ! OriginalRangeMap.SameAs( OriginalMatrixMap ) ) UseDataInPlace_ = 0 ; 
  if ( ! OriginalDomainMap.SameAs( OriginalMatrixMap ) ) UseDataInPlace_ = 0 ; 
  Comm().Broadcast( &UseDataInPlace_, 1, 0 ) ;

  if ( debug_ ) cout << __FILE__ << "::" << __LINE__ << " iam = " << iam << endl ; 
  //
  //  Reindex matrix if necessary (and possible - i.e. CrsMatrix for now)
  //
  //  For now, since I don't know how to determine if we need to reindex the matrix, 
  //  I will reindex them all - and then either allow the user to choose or figure it out. 
  //
  CrsMatrixA_ = dynamic_cast<Epetra_CrsMatrix *>(Problem_->GetOperator());

  if (Reindex_) 
  {
    if (CrsMatrixA_ == 0)
      AMESOS_CHK_ERR(-1);
  }
  if ( debug_ ) cout << __FILE__ << "::" << __LINE__ << " iam = " << iam << endl ; 
  if  ( Reindex_ ) {
#ifdef HAVE_AMESOS_EPETRAEXT
    const Epetra_Map& OriginalMap = CrsMatrixA_->RowMap();
    StdIndex_ = rcp( new Amesos_StandardIndex( OriginalMap  ) );
    const Epetra_Map& OriginalColMap = CrsMatrixA_->RowMap();
    StdIndexDomain_ = rcp( new Amesos_StandardIndex( OriginalDomainMap  ) );
    StdIndexRange_ = rcp( new Amesos_StandardIndex( OriginalRangeMap  ) );

    if ( debug_ ) cout << __FILE__ << "::" << __LINE__ << " iam = " << iam << endl ; 
    //    StdIndexMatrix_ = &((*MatTrans_)( *CrsMatrixA_ ));
    StdIndexMatrix_ = StdIndex_->StandardizeIndex( CrsMatrixA_ );
#else
    cout << "Amesos_Klu requires EpetraExt to reindex matrices." << endl 
	 <<  " Please rebuild with the EpetraExt library by adding --enable-epetraext to your configure invocation" << endl ;
    return 13 ; 
#endif
  } else { 
    StdIndexMatrix_ = RowMatrixA_ ;
  }

  //
  //  Convert Original Matrix to Serial (if it is not already)
  //
  if ( debug_ ) cout << __FILE__ << "::" << __LINE__ << " iam = " << iam << endl ; 
  if (UseDataInPlace_ == 1) {
    SerialMatrix_ = StdIndexMatrix_;
  } else {
    SerialMap_ = rcp(new Epetra_Map(NumGlobalElements_, NumMyElements_, 0, Comm()));
    
  if ( debug_ ) cout << __FILE__ << "::" << __LINE__ << " iam = " << iam << endl ; 
    ImportToSerial_ = rcp(new Epetra_Import ( *SerialMap_, StdIndexMatrix_->RowMatrixRowMap() ) );
    
  if ( debug_ ) cout << __FILE__ << "::" << __LINE__ << " iam = " << iam << endl ; 
    if (ImportToSerial_.get() == 0) AMESOS_CHK_ERR(-1);

#if 0    
    assert ( OriginalRangeMap.SameAs( OriginalMatrixMap ) ) ;
    assert ( OriginalDomainMap.SameAs( OriginalMatrixMap ) );
#endif    
   
  if ( debug_ ) cout << __FILE__ << "::" << __LINE__ << " iam = " << iam << endl ; 
    SerialCrsMatrixA_ = rcp( new Epetra_CrsMatrix(Copy, *SerialMap_, 0) )
      ;
    SerialMatrix_ = &*SerialCrsMatrixA_ ;
  }
AddTime("matrix redistribution");

  if ( debug_ ) cout << __FILE__ << "::" << __LINE__ << " iam = " << iam << endl ; 
return 0;
}

//=============================================================================
int Amesos_Klu::CreateSerialMap()
{

assert( false ) ; 
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
//    Ap, Ai, Aval contain the matrix as Klu needs it
//
//
int Amesos_Klu::ConvertToKluCRS(bool firsttime)
{

  ResetTime();

  //
  //  Convert matrix to the form that Klu expects (Ap, Ai, Aval)
  //

  if (iam==0) {
    assert( NumGlobalElements_ == SerialMatrix_->NumGlobalRows());
    assert( NumGlobalElements_ == SerialMatrix_->NumGlobalCols());

#if 0    
    cout << __FILE__ << "::" << __LINE__ 
	 << " numentries_ = " <<  numentries_ 
	 << " SerialMatrix_->NumGlobalNonzeros() = " <<  SerialMatrix_->NumGlobalNonzeros() 
	 << endl ; 
#endif

    assert( numentries_ == SerialMatrix_->NumGlobalNonzeros()) ;
    if ( firsttime ) { 
      Ap.resize( NumGlobalElements_+1 );
      Ai.resize( EPETRA_MAX( NumGlobalElements_, numentries_) ) ;
      Aval.resize( EPETRA_MAX( NumGlobalElements_, numentries_) ) ;
    }

    int NumEntriesThisRow;
    int Ai_index = 0 ;
    int MyRow;
    Epetra_CrsMatrix *CrsMatrix = dynamic_cast<Epetra_CrsMatrix *>(SerialMatrix_);
    //    Epetra_CrsMatrix *CrsMatrix = 0 ;  //  Uncomment this and comment the above line to test how we do with Row Matrices

    int MaxNumEntries_ = SerialMatrix_->MaxNumEntries();
    if ( firsttime && CrsMatrix == 0 ) {
      ColIndicesV_.resize(MaxNumEntries_);
      RowValuesV_.resize(MaxNumEntries_);
    }
    double *RowValues;
    int *ColIndices;

    for ( MyRow = 0; MyRow <NumGlobalElements_; MyRow++ ) {
      if ( CrsMatrix != 0 ) {
	EPETRA_CHK_ERR( CrsMatrix->
			ExtractMyRowView( MyRow, NumEntriesThisRow, RowValues,
					  ColIndices ) != 0 ) ;
      } else {
	EPETRA_CHK_ERR( SerialMatrix_->
			ExtractMyRowCopy( MyRow, MaxNumEntries_,
					  NumEntriesThisRow, &RowValuesV_[0],
					  &ColIndicesV_[0] ) != 0 ) ;
	RowValues =  &RowValuesV_[0];
	ColIndices = &ColIndicesV_[0];
      }

      if ( firsttime ) {
	Ap[MyRow] = Ai_index ;
	for ( int j = 0; j < NumEntriesThisRow; j++ ) {
	  Ai[Ai_index] = ColIndices[j] ;
	  Ai_index++;
	}
      } else { 
	for ( int j = 0; j < NumEntriesThisRow; j++ ) {
	  Aval[Ai_index] = RowValues[j] ;     
          if (ColIndices[j] == MyRow) {
            Aval[Ai_index] += AddToDiag_;     // Bug #1405   - this fails if the matrix is missing diagonal entries 
	  }
	  Ai_index++;
	}
      }
    }
    Ap[MyRow] = Ai_index ;
  }

  AddTime("matrix conversion");

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

  // solve problem with transpose
  if( ParameterList.isParameter("UseTranspose") )
    SetUseTranspose(ParameterList.get("UseTranspose",UseTranspose()));

  // MS // now comment it out, if we have parameters for KLU sublist
  // MS // uncomment it
  /*
  if (ParameterList.isSublist("Klu") ) {
    Teuchos::ParameterList KluParams = ParameterList.sublist("Klu") ;
  }
  */

  return 0;
}


//=============================================================================
int Amesos_Klu::PerformSymbolicFactorization() 
{
  ResetTime();

  if (iam == 0) {
    if (PrivateKluData_->Symbolic_) {
	klu_btf_free_symbolic (&(PrivateKluData_->Symbolic_)) ;
    }

    PrivateKluData_->Symbolic_ =
	klu_btf_analyze (NumGlobalElements_, &Ap[0], &Ai[0], (klu_control *) 0);
    if ( PrivateKluData_->Symbolic_ == 0 ) EPETRA_CHK_ERR( 1 ) ;
  }

  AddTime("symbolic");

  return 0;
}

//=============================================================================
int Amesos_Klu::PerformNumericFactorization( ) 
{
  ResetTime();

  if (iam == 0) {

    bool factor_with_pivoting = true ;

    // set the default parameters
    klu_control control ;
    klu_btf_defaults (&control) ;
    control.scale = ScaleMethod_ ;

    // see if we can "refactorize"
    if ( refactorize_ && PrivateKluData_->Numeric_ ) {

	// refactorize using the existing Symbolic and Numeric objects, and
	// using the identical pivot ordering as the prior klu_btf_factor.
	// No partial pivoting is done.
	int result = klu_btf_refactor (&Ap[0], &Ai[0], &Aval[0],
		    PrivateKluData_->Symbolic_, &control,
		    PrivateKluData_->Numeric_) ;

	// Did it work?
	if ( result == KLU_OK) {

	    // Get the largest and smallest entry on the diagonal of U
	    double umin = (PrivateKluData_->Numeric_)->umin ;
	    double umax = (PrivateKluData_->Numeric_)->umax ;

	    // compute a crude estimate of the reciprocal of
	    // the condition number
	    double rcond = umin / umax ;

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

	// destroy the existing Numeric object, if it exists
	if ( PrivateKluData_->Numeric_ ) {
	    klu_btf_free_numeric (&(PrivateKluData_->Numeric_)) ;
	}

	// factor the matrix using partial pivoting
	PrivateKluData_->Numeric_ =
	    klu_btf_factor (&Ap[0], &Ai[0], &Aval[0],
		    PrivateKluData_->Symbolic_, &control) ;
	if ( PrivateKluData_->Numeric_ == 0 ) EPETRA_CHK_ERR( 2 ) ;
    }

  }

  AddTime("numeric");

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

  iam = Comm().MyPID();
  if ( debug_ ) cout << __FILE__ << "::" << __LINE__ << " iam = " << iam  << " Entering SymbolicFactorization()" << endl ; 

  IsSymbolicFactorizationOK_ = false;
  IsNumericFactorizationOK_ = false;
  
  InitTime(Comm());

  NumSymbolicFact_++;

  if ( debug_ ) cout << __FILE__ << "::" << __LINE__ << " iam = " << iam << endl ; 
  AMESOS_CHK_ERR( CreateLocalMatrixAndExporters() ) ;
  assert( NumGlobalElements_ == RowMatrixA_->NumGlobalCols() );
  if ( debug_ ) cout << __FILE__ << "::" << __LINE__ << " iam = " << iam << endl ; 

  AMESOS_CHK_ERR( ExportToSerial() );
  if ( debug_ ) cout << __FILE__ << "::" << __LINE__ << " iam = " << iam << endl ; 

  AMESOS_CHK_ERR( ConvertToKluCRS(true) );
  if ( debug_ ) cout << __FILE__ << "::" << __LINE__ << " iam = " << iam << endl ; 

  AMESOS_CHK_ERR( PerformSymbolicFactorization() );
  if ( debug_ ) cout << __FILE__ << "::" << __LINE__ << " iam = " << iam << endl ; 

  IsSymbolicFactorizationOK_ = true;
  
  if ( debug_ ) cout << __FILE__ << "::" << __LINE__ << " iam = " << iam  << " Leaving SymbolicFactorization()" << endl ; 
  return 0;
}

//=============================================================================
int Amesos_Klu::NumericFactorization() 
{
 
  if ( debug_ ) cout << __FILE__ << "::" << __LINE__ << " iam = " << iam  << " Entering NumericFactorization()" << endl ; 
  IsNumericFactorizationOK_ = false;
  if (IsSymbolicFactorizationOK_ == false)
    AMESOS_CHK_ERR(SymbolicFactorization());

  NumNumericFact_++;

  Epetra_CrsMatrix *CrsMatrixA = dynamic_cast<Epetra_CrsMatrix *>(RowMatrixA_);
  if ( CrsMatrixA == 0 )   // hack to get around Bug #1502 
    AMESOS_CHK_ERR( CreateLocalMatrixAndExporters() ) ;
  assert( NumGlobalElements_ == RowMatrixA_->NumGlobalCols() );
  assert( numentries_ == RowMatrixA_->NumGlobalNonzeros() );
  AMESOS_CHK_ERR( ExportToSerial() );

  AMESOS_CHK_ERR( ConvertToKluCRS(false) );

  AMESOS_CHK_ERR( PerformNumericFactorization() );

  IsNumericFactorizationOK_ = true;
  
  if ( debug_ ) cout << __FILE__ << "::" << __LINE__ << " iam = " << iam  << " Leaving NumericFactorization()" << endl ; 
  return 0;
}

//=============================================================================
int Amesos_Klu::Solve() 
{

 if ( debug_ )     cout << __FILE__ << "::" << __LINE__ << " iam = " << iam  << " Entering Solve()" << endl ; 
  if (IsNumericFactorizationOK_ == false)
    AMESOS_CHK_ERR(NumericFactorization());
  
  ++NumSolve_;

  //
  //  Reindex the LHS and RHS 
  //
  Epetra_MultiVector* OrigVecX = Problem_->GetLHS() ;
  Epetra_MultiVector* OrigVecB = Problem_->GetRHS() ;
  Epetra_MultiVector* vecX ;
  Epetra_MultiVector* vecB ;

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

  if ( debug_ ) cout << __FILE__ << "::" << __LINE__ << endl ; 
  if ((vecX == 0) || (vecB == 0))
    AMESOS_CHK_ERR(-1); // something wrong in input
  
  if ( debug_ ) cout << __FILE__ << "::" << __LINE__ << endl ; 
  int NumVectors = vecX->NumVectors();
  if (NumVectors != vecB->NumVectors())
    AMESOS_CHK_ERR(-1); // something wrong in input

  if ( debug_ ) cout << __FILE__ << "::" << __LINE__ << endl ; 
  // vectors with SerialMap_
  Epetra_MultiVector* SerialB = 0;
  Epetra_MultiVector* SerialX = 0;

  //  Extract Serial versions of X and B
  double *SerialXvalues ;

  Epetra_MultiVector* SerialXextract = 0;
  Epetra_MultiVector* SerialBextract = 0;

  ResetTime();

  //  Copy B to the serial version of B
  //
  if (UseDataInPlace_ == 1) {
    SerialB = vecB;
    SerialX = vecX;
  } else {
    assert (UseDataInPlace_ == 0);
    assert ( SerialBextract == 0 ) ; 
    assert ( SerialXextract == 0 ) ; 

    if ( debug_ ) cout << __FILE__ << "::" << __LINE__ << endl ; 
    SerialXextract = new Epetra_MultiVector(*SerialMap_,NumVectors);
    SerialBextract = new Epetra_MultiVector(*SerialMap_,NumVectors);
    
    if ( debug_ ) cout << __FILE__ << "::" << __LINE__ << endl ; 
    ImportRangeToSerial_ = rcp(new Epetra_Import ( *SerialMap_, vecB->Map() ) );
    SerialBextract->Import(*vecB,*ImportRangeToSerial_,Insert);
    if ( debug_ ) cout << __FILE__ << "::" << __LINE__ << endl ; 

    SerialB = SerialBextract ;
    SerialX = SerialXextract ;
  }

  if ( debug_ ) cout << __FILE__ << "::" << __LINE__ << endl ; 
  AddTime("vector redistribution");

  //  Call KLU to perform the solve

  if ( debug_ ) cout << __FILE__ << "::" << __LINE__ << endl ; 
  SerialX->Scale(1.0, *SerialB) ;    // X = B (Klu overwrites B with X)
  ResetTime();

  if ( debug_ ) cout << __FILE__ << "::" << __LINE__ << endl ; 
  int SerialXlda ;
  if (iam == 0) {
    AMESOS_CHK_ERR(SerialX->ExtractView(&SerialXvalues,&SerialXlda ));

    if (SerialXlda != NumGlobalElements_)
      AMESOS_CHK_ERR(-1);

    if (UseTranspose()) {
      klu_btf_solve( PrivateKluData_->Symbolic_, PrivateKluData_->Numeric_,
		     SerialXlda, NumVectors, &SerialXvalues[0] );
    } else {
      klu_btf_tsolve( PrivateKluData_->Symbolic_, PrivateKluData_->Numeric_,
		      SerialXlda, NumVectors, &SerialXvalues[0] );
    }
  }
  if ( debug_ ) cout << __FILE__ << "::" << __LINE__ << endl ; 

  AddTime("solve");

  //  Copy X back to the original vector

  ResetTime();

  if ( debug_ ) cout << __FILE__ << "::" << __LINE__ << endl ; 
  if (UseDataInPlace_ == 0) {
    ImportDomainToSerial_ = rcp(new Epetra_Import ( *SerialMap_, vecX->Map() ) );
    vecX->Export( *SerialX, *ImportDomainToSerial_, Insert ) ;
    delete SerialBextract ;
    delete SerialXextract ;

  } // otherwise we are already in place.

  if ( debug_ ) cout << __FILE__ << "::" << __LINE__ << endl ; 
  AddTime("vector redistribution");

#if 0
  //
  //  ComputeTrueResidual causes TestOptions to fail on my linux box 
  //  Bug #1147
  if (ComputeTrueResidual_)
    ComputeTrueResidual(*SerialMatrix_, *vecX, *vecB, UseTranspose(), "Amesos_Klu");
#endif

  if (ComputeVectorNorms_)
    ComputeVectorNorms(*vecX, *vecB, "Amesos_Klu");

  if ( debug_ ) cout << __FILE__ << "::" << __LINE__ << " iam = " << iam  << " Leaving Solve()" << endl ; 
  return(0) ;
}

// ================================================ ====== ==== ==== == =

void Amesos_Klu::PrintStatus() const
{

  if (iam) return;

  PrintLine();

  cout << "Amesos_Klu : Matrix has " << NumGlobalElements_ << " rows"
       << " and " << numentries_ << " nonzeros" << endl;
  cout << "Amesos_Klu : Nonzero elements per row = "
       << 1.0*numentries_/NumGlobalElements_ << endl;
  cout << "Amesos_Klu : Percentage of nonzero elements = "
       << 100.0*numentries_/(pow(NumGlobalElements_,2.0)) << endl;
  cout << "Amesos_Klu : Use transpose = " << UseTranspose_ << endl;

  PrintLine();

  return;

}

// ================================================ ====== ==== ==== == =

void Amesos_Klu::PrintTiming() const
{
  if (iam) return;

  double ConTime = GetTime("matrix conversion");
  double MatTime = GetTime("matrix redistribution");
  double VecTime = GetTime("vector redistribution");
  double SymTime = GetTime("symbolic");
  double NumTime = GetTime("numeric");
  double SolTime = GetTime("solve");

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
       << SymTime << " (s), avg = " << SymTime << " (s)" << endl;
  cout << p << "Number of numeric factorizations = "
       << NumNumericFact_ << endl;
  cout << p << "Time for num fact = "
       << NumTime << " (s), avg = " << NumTime << " (s)" << endl;
  cout << p << "Number of solve phases = "
       << NumSolve_ << endl;
  cout << p << "Time for solve = "
       << SolTime << " (s), avg = " << SolTime << " (s)" << endl;

  PrintLine();

  return;
}
