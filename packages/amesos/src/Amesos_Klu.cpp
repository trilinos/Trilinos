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
  TransposeMatrix_(0),
  Matrix_(0),
  UseTranspose_(false),
  Problem_(&prob), 
  PrintTiming_(false),
  PrintStatus_(false),
  ComputeVectorNorms_(false),
  ComputeTrueResidual_(false),
  verbose_(1),
  debug_(0),
  ConTime_(0.0),
  SymTime_(0.0),
  NumTime_(0.0),
  SolTime_(0.0),
  VecTime_(0.0),
  MatTime_(0.0),
  NumSymbolicFact_(0),
  NumNumericFact_(0),
  NumSolve_(0),
  Time_(0)
{
  // MS // move declaration of Problem_ above because I need it
  // MS // set up before calling Comm()
  Teuchos::ParameterList ParamList ;
  SetParameters( ParamList ) ;

}

//=============================================================================
Amesos_Klu::~Amesos_Klu(void) {

  if ( SerialMap_ ) delete SerialMap_ ; 
  if ( SerialCrsMatrixA_ ) delete SerialCrsMatrixA_ ; 
  if ( TransposeMatrix_ ) delete TransposeMatrix_ ; 
  delete PrivateKluData_; 

  // print out some information if required by the user
  if( (verbose_ && PrintTiming_) || verbose_ == 2 ) PrintTiming();
  if( (verbose_ && PrintStatus_) || verbose_ == 2 ) PrintStatus();

  if( Time_ ) { delete Time_; Time_ = 0; }
  
}
//  See  pre and post conditions in Amesos_Klu.h

int Amesos_Klu::ConvertToSerial() { 

  if( debug_ == 1 ) cout << "Entering `ConvertToSerial()'" << endl;

  Time_->ResetStartTime();
  
  Epetra_RowMatrix *RowMatrixA = dynamic_cast<Epetra_RowMatrix *>(Problem_->GetOperator());
  EPETRA_CHK_ERR( RowMatrixA == 0 ) ; 

  Epetra_CrsMatrix *CastCrsMatrixA = dynamic_cast<Epetra_CrsMatrix*>(RowMatrixA) ; 
  EPETRA_CHK_ERR( CastCrsMatrixA == 0 ) ; 

  iam = Comm().MyPID() ;

  const Epetra_Map &OriginalMap = CastCrsMatrixA->RowMap() ; 

  NumGlobalElements_ = CastCrsMatrixA->NumGlobalRows();
  numentries_ = CastCrsMatrixA->NumGlobalNonzeros();
  assert( NumGlobalElements_ == CastCrsMatrixA->NumGlobalCols() );

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
     SerialMatrix_ = CastCrsMatrixA ;
  } else {
    if ( SerialMap_ ) delete SerialMap_ ; 
    SerialMap_ = new Epetra_Map( NumGlobalElements_, NumMyElements_, 0, Comm() );

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

//
//  See also pre and post conditions in Amesos_Klu.h
//  Preconditions:  
//    SerialMatrix_ points to the matrix to be factored and solved
//    NumGlobalElements_ has been set to the dimension of the matrix
//    numentries_ has been set to the number of non-zeros in the matrix
//      (i.e. ConvertToSerial() has been callded) 
//
//  Postconditions:
//    Ap, Ai, Aval contain the matrix as Klu needs it (transposed 
//     iff UseTranspose() is NOT true)
//
//  Issues:
//    Is there a way to avoid deleting and recreating TransposeMatrix_?
//    The straight forward attempt at a solution, i.e. re-using 
//    TransposeMatrix_, fails in the call to CrsMatrixTranspose.
//    However, since we own CrsMatrixTranspose, we can fix it.
//
//    Mike has mentioned a desire to have a Matrix Transpose operation
//    in epetra (or epetraext).
//
int Amesos_Klu::ConvertToKluCRS(bool firsttime){

  if( debug_ == 1 ) cout << "Entering `ConvertToKluCRSt()'" << endl;

  Time_->ResetStartTime();
  
  if ( UseTranspose() ) { 
    Matrix_ = SerialMatrix_ ; 
  } else { 
#if 1
    if ( TransposeMatrix_ ) delete TransposeMatrix_ ; 
    TransposeMatrix_ =  new Epetra_CrsMatrix(Copy, (Epetra_Map &)SerialMatrix_->Map(), 0);
#else
    //
    //  This fails in CrsMatrixTranspose
    //
    if ( firsttime && TransposeMatrix_ ) delete TransposeMatrix_ ; 
    if ( firsttime ) TransposeMatrix_ =  new Epetra_CrsMatrix(Copy, (Epetra_Map &)SerialMatrix_->Map(), 0);
#endif
    EPETRA_CHK_ERR( CrsMatrixTranspose( SerialMatrix_, TransposeMatrix_ ) ) ; 
    Matrix_ = TransposeMatrix_ ; 
  }
  //
  //  Convert matrix to the form that Klu expects (Ap, Ai, Aval) 
  //

  if ( iam==0 ) {
    assert( NumGlobalElements_ == Matrix_->NumGlobalRows());
    assert( NumGlobalElements_ == Matrix_->NumGlobalCols());
    assert( numentries_ == Matrix_->NumGlobalNonzeros());
    Ap.resize( NumGlobalElements_+1 );
    Ai.resize( EPETRA_MAX( NumGlobalElements_, numentries_) ) ; 
    Aval.resize( EPETRA_MAX( NumGlobalElements_, numentries_) ) ; 

    int NumEntriesThisRow;
    int Ai_index = 0 ; 
    int MyRow;
#ifdef USE_VIEW
    double *RowValues;
    int *ColIndices;
#else
    int MaxNumEntries_ = Matrix_->MaxNumEntries();
    vector<int>ColIndices(MaxNumEntries_);
    vector<double>RowValues(MaxNumEntries_);
#endif
    for ( MyRow = 0; MyRow <NumGlobalElements_; MyRow++ ) {
#ifdef USE_VIEW
      EPETRA_CHK_ERR( Matrix_->
		      ExtractMyRowView( MyRow, NumEntriesThisRow, RowValues, 
					ColIndices ) != 0 ) ;
#else
      EPETRA_CHK_ERR( Matrix_->
		      ExtractMyRowCopy( MyRow, MaxNumEntries_, 
					NumEntriesThisRow, &RowValues[0], 
					&ColIndices[0] ) != 0 ) ;
#endif
      Ap[MyRow] = Ai_index ; 
      for ( int j = 0; j < NumEntriesThisRow; j++ ) { 
	Ai[Ai_index] = ColIndices[j] ; 
	Aval[Ai_index] = RowValues[j] ; 
	Ai_index++;
      }
    }
    Ap[MyRow] = Ai_index ; 
  }

  ConTime_ += Time_->ElapsedTime(); 

  return 0;
}   


int Amesos_Klu::SetParameters( Teuchos::ParameterList &ParameterList ) {

  //  Some compilers reject the following cast:
  //  if( &ParameterList == 0 ) return 0;

  if( debug_ == 1 ) cout << "Entering `SetParameters()' ..." << endl;
  
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

  // possible debug statements
  // 0 - no debug
  // 1 - debug
  if( ParameterList.isParameter("DebugLevel") )
    debug_ = ParameterList.get("DebugLevel",0);
  
  // MS // now comment it out, if we have parameters for KLU sublist
  // MS // uncomment it
  /*
  if (ParameterList.isSublist("Klu") ) {
    Teuchos::ParameterList KluParams = ParameterList.sublist("Klu") ;
  }
  */
  
  return 0;
}


int Amesos_Klu::PerformSymbolicFactorization() {

  if( debug_ == 1 ) cout << "Entering `PerformSymbolicFactorization()'" << endl;

  Time_->ResetStartTime();  

  if ( iam == 0 ) { 
    if ( PrivateKluData_->Symbolic_ ) klu_btf_free_symbolic (&(PrivateKluData_->Symbolic_)) ;

    PrivateKluData_->Symbolic_ = klu_btf_analyze (NumGlobalElements_, &Ap[0], &Ai[0] ) ;
    if ( PrivateKluData_->Symbolic_ == 0 ) EPETRA_CHK_ERR( 1 ) ; 
  }

  SymTime_ += Time_->ElapsedTime();

  return 0;
}

int Amesos_Klu::PerformNumericFactorization( ) {
  
  if( debug_ == 1 ) cout << "Entering `PerformNumericFactorization()'" << endl;

  Time_->ResetStartTime();
  
  if ( iam == 0 ) {

    const double tol = 0.001; //  At some point we need to expose this to the user 

    if ( PrivateKluData_->Numeric_ ) klu_btf_free_numeric (&(PrivateKluData_->Numeric_)) ;
    PrivateKluData_->Numeric_ = klu_btf_factor (&Ap[0], &Ai[0], &Aval[0], tol, PrivateKluData_->Symbolic_) ;
    if ( PrivateKluData_->Numeric_ == 0 ) EPETRA_CHK_ERR( 2 ) ; 
    
  }

  NumTime_ += Time_->ElapsedTime();

  return 0;
}




bool Amesos_Klu::MatrixShapeOK() const { 
  bool OK ;

  if ( GetProblem()->GetOperator()->OperatorRangeMap().NumGlobalPoints() != 
       GetProblem()->GetOperator()->OperatorDomainMap().NumGlobalPoints() ) OK = false;
  return OK; 
}


int Amesos_Klu::SymbolicFactorization() {
  
  if( debug_ == 1 ) cout << "Entering `SymbolicFactorization()'" << endl;
  if( Time_ == 0 ) Time_ = new Epetra_Time( Comm() );
  
  NumSymbolicFact_++;  
  
  ConvertToSerial() ; 
  
  ConvertToKluCRS(true);
  
  PerformSymbolicFactorization();

  return 0;
}

int Amesos_Klu::NumericFactorization() {
  
  if( debug_ == 1 ) cout << "Entering `NumericFactorization()'" << endl;
  if( Time_ == 0 ) Time_ = new Epetra_Time( Comm() );
  
  NumNumericFact_++;

  ConvertToSerial() ; 
  ConvertToKluCRS(false);

  PerformNumericFactorization( );
  return 0;
}


int Amesos_Klu::Solve() { 

  if( debug_ == 1 ) cout << "Entering `Solve()'" << endl;
  if( Time_ == 0 ) Time_ = new Epetra_Time( Comm() );
  
  NumSolve_++;
  
  Epetra_MultiVector   *vecX = Problem_->GetLHS() ; 
  Epetra_MultiVector   *vecB = Problem_->GetRHS() ; 

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

  Epetra_RowMatrix *RowMatrixA = dynamic_cast<Epetra_RowMatrix *>(Problem_->GetOperator());
  Epetra_CrsMatrix *CastCrsMatrixA = dynamic_cast<Epetra_CrsMatrix*>(RowMatrixA) ; 
  Epetra_MultiVector *SerialXextract = 0;
  Epetra_MultiVector *SerialBextract = 0;

  Time_->ResetStartTime(); // track time to broadcast vectors
  //
  //  Copy B to the serial version of B
  //
  if ( IsLocal_ ==1 ) { 
    SerialB = vecB ; 
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

  VecTime_ += Time_->ElapsedTime();
    
  //
  //  Call KLU to perform the solve
  //

  SerialX->Scale(1.0, *SerialB) ;

  Time_->ResetStartTime(); // tract time to solve
  
  int SerialXlda ; 
  if ( iam == 0 ) {
    assert( SerialX->ExtractView( &SerialXvalues, &SerialXlda ) == 0 ) ; 

    assert( SerialXlda == NumGlobalElements_ ) ; 
    
    for ( int j =0 ; j < nrhs; j++ ) { 

      vector<double> workspace( NumGlobalElements_ ) ; 
      klu_btf_solve( PrivateKluData_->Symbolic_, PrivateKluData_->Numeric_,
		     &SerialXvalues[j*SerialXlda], &workspace[0] );
    }
  }
    
  SolTime_ += Time_->ElapsedTime();
    
   
  //
  //  Copy X back to the original vector
  // 

  Time_->ResetStartTime();  // track time to broadcast vectors

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

  VecTime_ += Time_->ElapsedTime();

  // MS // compute vector norms
  if( ComputeVectorNorms_ == true || verbose_ == 2 ) {
    double NormLHS, NormRHS;
    for( int i=0 ; i<nrhs ; ++i ) {
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
    Epetra_MultiVector Ax(vecB->Map(),nrhs);
    for( int i=0 ; i<nrhs ; ++i ) {
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
