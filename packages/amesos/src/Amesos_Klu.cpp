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
#include "Epetra_Export.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_Util.h"
#include "CrsMatrixTranspose.h"
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
  SerialMap_(0),
  SerialCrsMatrixA_(0),
  SerialMatrix_(0),
  Matrix_(0),
  UseTranspose_(false),
  Problem_(&prob),
  PrintTiming_(false),
  PrintStatus_(false),
  ComputeVectorNorms_(false),
  ComputeTrueResidual_(false),
  verbose_(1),
  IsSymbolicFactorizationOK_(false),
  IsNumericFactorizationOK_(false),
  refactorize_(false),
  rcond_threshold_(1e-12),
  ScaleMethod_(1),
  ConTime_(0.0),
  SymTime_(0.0),
  NumTime_(0.0),
  SolTime_(0.0),
  VecTime_(0.0),
  MatTime_(0.0),
  NumSymbolicFact_(0),
  NumNumericFact_(0),
  NumSolve_(0),
  Time_(0),
  ImportToSerial_(0)
{
  // MS // move declaration of Problem_ above because I need it
  // MS // set up before calling Comm()
  Teuchos::ParameterList ParamList ;
  SetParameters( ParamList ) ;

}

//=============================================================================
Amesos_Klu::~Amesos_Klu(void) {

  if (SerialMap_) 
    delete SerialMap_;
  if (SerialCrsMatrixA_) 
    delete SerialCrsMatrixA_;

  delete PrivateKluData_;

  if (Time_) 
    delete Time_;

  if (ImportToSerial_) 
    delete ImportToSerial_;

  // print out some information if required by the user
  if( (verbose_ && PrintTiming_) || verbose_ == 2 ) PrintTiming();
  if( (verbose_ && PrintStatus_) || verbose_ == 2 ) PrintStatus();

}

//=============================================================================
int Amesos_Klu::ConvertToSerial() {

  Time_->ResetStartTime();

  Epetra_RowMatrix *RowMatrixA = dynamic_cast<Epetra_RowMatrix *>(Problem_->GetOperator());
  EPETRA_CHK_ERR( RowMatrixA == 0 ) ;

  iam = Comm().MyPID() ;

  const Epetra_Map &OriginalMap = RowMatrixA->RowMatrixRowMap() ;

  NumGlobalElements_ = RowMatrixA->NumGlobalRows();
  numentries_ = RowMatrixA->NumGlobalNonzeros();
  assert( NumGlobalElements_ == RowMatrixA->NumGlobalCols() );

  //
  //  Create a serial matrix
  //
  assert( NumGlobalElements_ == OriginalMap.NumGlobalElements() ) ;
  int NumMyElements_ = 0 ;
  if (iam==0) NumMyElements_ = NumGlobalElements_;


  IsLocal_ = ( OriginalMap.NumMyElements() ==
	       OriginalMap.NumGlobalElements() )?1:0;
  Comm().Broadcast( &IsLocal_, 1, 0 ) ;

  //
  //  KEN:  Consider giving Epetra_RowMatrix_Transposer a shot
  //  I am not confident that  Epetra_RowMatrix_Transposer works,
  //  but it is worth a try.
  //
  //  Convert Original Matrix to Serial (if it is not already)
  //
  if (SerialMap_) { delete SerialMap_ ; SerialMap_ = 0 ; }
  if ( SerialCrsMatrixA_ ) { delete SerialCrsMatrixA_ ; SerialCrsMatrixA_ = 0 ; }
  if ( IsLocal_==1 ) {
     SerialMatrix_ = RowMatrixA ;
  } else {
    if ( SerialMap_ ) delete SerialMap_ ;
    SerialMap_ = new Epetra_Map( NumGlobalElements_, NumMyElements_, 0, Comm() );

    // FIXME: can I be ImportToSerial??
    Epetra_Export export_to_serial( OriginalMap, *SerialMap_);

    if ( SerialCrsMatrixA_ ) delete SerialCrsMatrixA_ ;
    SerialCrsMatrixA_ = new Epetra_CrsMatrix(Copy, *SerialMap_, 0);
    SerialCrsMatrixA_->Export( *RowMatrixA, export_to_serial, Add );

    SerialCrsMatrixA_->TransformToLocal() ;
    SerialMatrix_ = SerialCrsMatrixA_ ;
  }

  MatTime_ += Time_->ElapsedTime();

  return 0;
}

//=============================================================================
int Amesos_Klu::CreateSerialMap()
{
  Epetra_RowMatrix *RowMatrixA;
  RowMatrixA = dynamic_cast<Epetra_RowMatrix *>(Problem_->GetOperator());
  if (RowMatrixA == 0)
    AMESOS_CHK_ERR(-1);

  iam = Comm().MyPID();
  NumGlobalElements_ = RowMatrixA->NumGlobalRows();
  numentries_ = RowMatrixA->NumGlobalNonzeros();
  int NumMyElements = 0;
  if (iam == 0) 
    NumMyElements = NumGlobalElements_;

  SerialMap_ = new Epetra_Map(NumGlobalElements_,NumMyElements,0,Comm());
  if (SerialMap_ == 0)
    AMESOS_CHK_ERR(-1);
  
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
//      (i.e. ConvertToSerial() has been callded)
//
//  Postconditions:
//    Ap, Ai, Aval contain the matrix as Klu needs it
//
//
int Amesos_Klu::ConvertToKluCRS(bool firsttime){

  Time_->ResetStartTime();

  Matrix_ = SerialMatrix_ ;
  //
  //  Convert matrix to the form that Klu expects (Ap, Ai, Aval)
  //

  if (iam==0) {
    assert( NumGlobalElements_ == Matrix_->NumGlobalRows());
    assert( NumGlobalElements_ == Matrix_->NumGlobalCols());
    assert( numentries_ == Matrix_->NumGlobalNonzeros());
    if ( firsttime ) { 
      Ap.resize( NumGlobalElements_+1 );
      Ai.resize( EPETRA_MAX( NumGlobalElements_, numentries_) ) ;
      Aval.resize( EPETRA_MAX( NumGlobalElements_, numentries_) ) ;
    }

    int NumEntriesThisRow;
    int Ai_index = 0 ;
    int MyRow;
    Epetra_CrsMatrix *CrsMatrix = dynamic_cast<Epetra_CrsMatrix *>(Matrix_);
    //    Epetra_CrsMatrix *CrsMatrix = 0 ;  //  Uncomment this and comment the above line to test how we do with Row Matrices

    int MaxNumEntries_ = Matrix_->MaxNumEntries();
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
	EPETRA_CHK_ERR( Matrix_->
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
	  Ai_index++;
	}
      }
    }
    Ap[MyRow] = Ai_index ;
  }

  ConTime_ += Time_->ElapsedTime();

  return 0;
}

//=============================================================================
int Amesos_Klu::SetParameters( Teuchos::ParameterList &ParameterList ) {

  //  Some compilers reject the following cast:
  //  if( &ParameterList == 0 ) return 0;

  // ========================================= //
  // retrive KLU's parameters from list.       //
  // default values defined in the constructor //
  // ========================================= //

  // retrive general parameters

  // solve problem with transpose
  if( ParameterList.isParameter("UseTranspose") )
    SetUseTranspose(ParameterList.get("UseTranspose",false));

  // print some timing information (on process 0)
  if( ParameterList.isParameter("PrintTiming") )
    PrintTiming_ = ParameterList.get("PrintTiming", false);

  // print some statistics (on process 0). Do not include timing
  if( ParameterList.isParameter("PrintStatus") )
    PrintStatus_ = ParameterList.get("PrintStatus", false);

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
    verbose_ = ParameterList.get("OutputLevel",1);

  // refactorize
  if( ParameterList.isParameter("Refactorize") )
    refactorize_ = ParameterList.get("Refactorize", false);

  // threshold for determining if refactorize worked OK
  // UNUSED at present - KSS June 2004
  if( ParameterList.isParameter("RcondThreshold") )
    rcond_threshold_ = ParameterList.get("RcondThreshold", 1e-12);

  // scaling method: 0: none, 1: use method's default, 2: use
  // the method's 1st alternative, 3: etc.
  if( ParameterList.isParameter("ScaleMethod") )
    ScaleMethod_ = ParameterList.get("ScaleMethod", 1);

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
int Amesos_Klu::PerformSymbolicFactorization() {

  Time_->ResetStartTime();

  if (iam == 0) {
    if (PrivateKluData_->Symbolic_) {
	klu_btf_free_symbolic (&(PrivateKluData_->Symbolic_)) ;
    }

    PrivateKluData_->Symbolic_ =
	klu_btf_analyze (NumGlobalElements_, &Ap[0], &Ai[0], (klu_control *) 0);
    if ( PrivateKluData_->Symbolic_ == 0 ) EPETRA_CHK_ERR( 1 ) ;
  }

  SymTime_ += Time_->ElapsedTime();

  return 0;
}

//=============================================================================
int Amesos_Klu::PerformNumericFactorization( ) {

  Time_->ResetStartTime();

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

  NumTime_ += Time_->ElapsedTime();

  return 0;
}

//=============================================================================
bool Amesos_Klu::MatrixShapeOK() const {
  bool OK ;

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

  IsSymbolicFactorizationOK_ = false;
  IsNumericFactorizationOK_ = false;
  
  if( Time_ == 0 ) Time_ = new Epetra_Time( Comm() );

  NumSymbolicFact_++;

  ConvertToSerial() ;

  ConvertToKluCRS(true);

  PerformSymbolicFactorization();

  IsSymbolicFactorizationOK_ = true;
  
  return 0;
}

//=============================================================================
int Amesos_Klu::NumericFactorization() 
{
 
  IsNumericFactorizationOK_ = false;
  if (IsSymbolicFactorizationOK_ == false)
    AMESOS_CHK_ERR(SymbolicFactorization());

  if( Time_ == 0 ) Time_ = new Epetra_Time( Comm() );

  NumNumericFact_++;

  ConvertToSerial() ;
  ConvertToKluCRS(false);

  PerformNumericFactorization();

  IsNumericFactorizationOK_ = true;
  
  return 0;
}

//=============================================================================
int Amesos_Klu::Solve() 
{

  if (IsNumericFactorizationOK_ == false)
    AMESOS_CHK_ERR(NumericFactorization());
  
  if( Time_ == 0 ) Time_ = new Epetra_Time( Comm() );

  ++NumSolve_;

  Epetra_MultiVector* vecX = Problem_->GetLHS() ;
  Epetra_MultiVector* vecB = Problem_->GetRHS() ;

  if ((vecX == 0) || (vecB == 0))
    AMESOS_CHK_ERR(-1); // something wrong in input
  
  int NumVectors = vecX->NumVectors();
  if (NumVectors != vecB->NumVectors())
    AMESOS_CHK_ERR(-1); // something wrong in input

  // vectors with SerialMap_
  Epetra_MultiVector* SerialB = 0;
  Epetra_MultiVector* SerialX = 0;

  //  Extract Serial versions of X and B
  double *SerialXvalues ;

  Epetra_RowMatrix *RowMatrixA = dynamic_cast<Epetra_RowMatrix *>(Problem_->GetOperator());
  Epetra_MultiVector* SerialXextract = 0;
  Epetra_MultiVector* SerialBextract = 0;

  Time_->ResetStartTime(); // track time to broadcast vectors

  //  Copy B to the serial version of B
  //
  if (IsLocal_ == 1) {
    SerialB = vecB;
    SerialX = vecX;
  } else {
    assert (IsLocal_ == 0);
    const Epetra_Map &OriginalMap = RowMatrixA->RowMatrixRowMap();

    // check whether the stored ImportToSerial_ (if allocated) 
    // is still valid or not.
    if (ImportToSerial_ != 0) {
      if (!(OriginalMap.SameAs(ImportToSerial_->TargetMap()))) {
	delete SerialMap_;
	AMESOS_CHK_ERR(CreateSerialMap());
	delete ImportToSerial_;
	ImportToSerial_ = 0;
      }
    }

    if (ImportToSerial_ == 0) {
      ImportToSerial_ = new Epetra_Import(*SerialMap_,OriginalMap);
      assert (ImportToSerial_ != 0);
    }

    Epetra_MultiVector *SerialXextract = 
      new Epetra_MultiVector(*SerialMap_,NumVectors);
    Epetra_MultiVector *SerialBextract = 
      new Epetra_MultiVector(*SerialMap_,NumVectors);
    
    SerialBextract->Import(*vecB,*ImportToSerial_,Insert);
    SerialB = SerialBextract ;
    SerialX = SerialXextract ;
  }

  VecTime_ += Time_->ElapsedTime();

  //  Call KLU to perform the solve

  SerialX->Scale(1.0, *SerialB) ;

  Time_->ResetStartTime(); // tract time to solve

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

  SolTime_ += Time_->ElapsedTime();

  //  Copy X back to the original vector

  Time_->ResetStartTime();  // track time to broadcast vectors

  if (IsLocal_ == 0) {
    vecX->Export( *SerialX, *ImportToSerial_, Insert ) ;
    delete SerialBextract ;
    delete SerialXextract ;
  } // otherwise we are already in place.

  VecTime_ += Time_->ElapsedTime();

  // MS // compute vector norms
  if( ComputeVectorNorms_ == true || verbose_ == 2 ) {
    double NormLHS, NormRHS;
    for( int i=0 ; i<NumVectors ; ++i ) {
      assert((*vecX)(i)->Norm2(&NormLHS)==0);
      assert((*vecB)(i)->Norm2(&NormRHS)==0);
      if( verbose_ && Comm().MyPID() == 0 ) {
	cout << "Amesos_Klu : vector " << i << ", ||x|| = " << NormLHS
	     << ", ||b|| = " << NormRHS << endl;
      }
    }
  }

  // MS // compute true residual
  if( ComputeTrueResidual_ == true || verbose_ == 2  ) {
    double Norm;
    Epetra_MultiVector Ax(vecB->Map(),NumVectors);
    for( int i=0 ; i<NumVectors ; ++i ) {
      (Problem_->GetMatrix()->Multiply(UseTranspose(), *((*vecX)(i)), Ax));
      (Ax.Update(1.0, *((*vecB)(i)), -1.0));
      (Ax.Norm2(&Norm));

      if( verbose_ && Comm().MyPID() == 0 ) {
	cout << "Amesos_Klu : vector " << i << ", ||Ax - b|| = " << Norm << endl;
      }
    }
  }

  return(0) ;
}

// ================================================ ====== ==== ==== == =

void Amesos_Klu::PrintStatus()
{

  if( iam != 0  ) return;

  cout << "----------------------------------------------------------------------------" << endl;
  cout << "Amesos_Klu : Matrix has " << NumGlobalElements_ << " rows"
       << " and " << numentries_ << " nonzeros" << endl;
  cout << "Amesos_Klu : Nonzero elements per row = "
       << 1.0*numentries_/NumGlobalElements_ << endl;
  cout << "Amesos_Klu : Percentage of nonzero elements = "
       << 100.0*numentries_/(pow(NumGlobalElements_,2.0)) << endl;
  cout << "Amesos_Klu : Use transpose = " << UseTranspose_ << endl;
  cout << "----------------------------------------------------------------------------" << endl;

  return;

}

// ================================================ ====== ==== ==== == =

void Amesos_Klu::PrintTiming()
{
  if( iam ) return;

  cout << "----------------------------------------------------------------------------" << endl;
  cout << "Amesos_Klu : Time to convert matrix to KLU format = "
       << ConTime_ << " (s)" << endl;
  cout << "Amesos_Klu : Time to redistribute matrix = "
       << MatTime_ << " (s)" << endl;
  cout << "Amesos_Klu : Time to redistribute vectors = "
       << VecTime_ << " (s)" << endl;
  cout << "Amesos_Klu : Number of symbolic factorizations = "
       << NumSymbolicFact_ << endl;
  cout << "Amesos_Klu : Time for sym fact = "
       << SymTime_ << " (s), avg = " << SymTime_/NumSymbolicFact_
       << " (s)" << endl;
  cout << "Amesos_Klu : Number of numeric factorizations = "
       << NumNumericFact_ << endl;
  cout << "Amesos_Klu : Time for num fact = "
       << NumTime_ << " (s), avg = " << NumTime_/NumNumericFact_
       << " (s)" << endl;
  cout << "Amesos_Klu : Number of solve phases = "
       << NumSolve_ << endl;
  cout << "Amesos_Klu : Time for solve = "
       << SolTime_ << " (s), avg = " << SolTime_/NumSolve_
       << " (s)" << endl;
  cout << "----------------------------------------------------------------------------" << endl;

  return;
}
