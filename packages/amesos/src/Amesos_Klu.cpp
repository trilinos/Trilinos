  /* Copyright (2003) Sandia Corportation. Under the terms of Contract 
   * DE-AC04-94AL85000, there is a non-exclusive license for use of this 
   * work by or on behalf of the U.S. Government.  Export of this program
   * may require a license from the United States Government. */


  /* NOTICE:  The United States Government is granted for itself and others
   * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
   * license in ths data to reproduce, prepare derivative works, and
   * perform publicly and display publicly.  Beginning five (5) years from
   * July 25, 2001, the United States Government is granted for itself and
   * others acting on its behalf a paid-up, nonexclusive, irrevocable
   * worldwide license in this data to reproduce, prepare derivative works,
   * distribute copies to the public, perform publicly and display
   * publicly, and to permit others to do so.
   * 
   * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
   * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
   * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
   * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
   * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
   * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */

#include "Amesos_Klu.h"
#if 0
extern "C" {
#include "klu_dump.h"
}
#endif
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_Util.h"
#include "CrsMatrixTranspose.h"

  //=============================================================================
  Amesos_Klu::Amesos_Klu(const Epetra_LinearProblem &prob, 
				 const AMESOS::Parameter::List &ParameterList ) :  
    SerialCrsMatrixA_(0), 
    SerialMap_(0), 
    SerialMatrix_(0), 
    TransposeMatrix_(0),
    Symbolic_(0),
    Numeric_(0),
    Matrix_(0) {


  Problem_ = &prob ; 
  ParameterList_ = &ParameterList ; 
}

//=============================================================================
Amesos_Klu::~Amesos_Klu(void) {

  if ( SerialMap_ ) delete SerialMap_ ; 
  if ( SerialCrsMatrixA_ ) delete SerialCrsMatrixA_ ; 
  if ( TransposeMatrix_ ) delete TransposeMatrix_ ; 
  if ( Symbolic_ ) klu_btf_free_symbolic (&Symbolic_) ;
  if ( Numeric_ ) klu_btf_free_numeric (&Numeric_) ;

}

int Amesos_Klu::ConvertToSerial() { 
  
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
    SerialCrsMatrixA_->Export( *CastCrsMatrixA, export_to_serial, Add ); 
    
    SerialCrsMatrixA_->TransformToLocal() ; 
    SerialMatrix_ = SerialCrsMatrixA_ ;
  }

  return 0;
} 

int Amesos_Klu::ConvertToKluCRS(){
  
  if ( UseTranspose() ) { 
    Matrix_ = SerialMatrix_ ; 
  } else { 
    if ( TransposeMatrix_ ) delete TransposeMatrix_ ; 
    TransposeMatrix_ =  new Epetra_CrsMatrix(Copy, (Epetra_Map &)SerialMatrix_->Map(), 0);
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
    double *RowValues;
    int *ColIndices;
    int Ai_index = 0 ; 
    int MyRow;
    for ( MyRow = 0; MyRow <NumGlobalElements_; MyRow++ ) {
      assert( Matrix_->ExtractMyRowView( MyRow, NumEntriesThisRow, RowValues, ColIndices ) == 0 ) ;
      Ap[MyRow] = Ai_index ; 
      for ( int j = 0; j < NumEntriesThisRow; j++ ) { 
	Ai[Ai_index] = ColIndices[j] ; 
	Aval[Ai_index] = RowValues[j] ; 
	Ai_index++;
      }
    }
    Ap[MyRow] = Ai_index ; 
  }

  
  return 0;
}   


int Amesos_Klu::PerformSymbolicFactorization() {

  if ( iam == 0 ) { 
    Symbolic_ = klu_btf_analyze (NumGlobalElements_, &Ap[0], &Ai[0] ) ;
    if ( Symbolic_ == 0 ) EPETRA_CHK_ERR( 1 ) ; 
  }
  return 0;
}

int Amesos_Klu::PerformNumericFactorization( ) {

  if ( iam == 0 ) {

    const double tol = 0.001; //  At some point we need to expose this to the user 

    Numeric_ = klu_btf_factor (&Ap[0], &Ai[0], &Aval[0], tol, Symbolic_) ;
    if ( Numeric_ == 0 ) EPETRA_CHK_ERR( 2 ) ; 
    
  }
  
  return 0;
}




bool Amesos_Klu::MatrixShapeOK() const { 
  bool OK ;

  if ( GetProblem()->GetOperator()->OperatorRangeMap().NumGlobalPoints() != 
       GetProblem()->GetOperator()->OperatorDomainMap().NumGlobalPoints() ) OK = false;
  return OK; 
}


int Amesos_Klu::SymbolicFactorization() {

  //  cout << " SymbolicFactorization() B = " << *(Problem_->GetRHS()) << endl ; 

  ConvertToSerial() ; 
  
  ConvertToKluCRS();
  
  PerformSymbolicFactorization();

  return 0;
}

int Amesos_Klu::NumericFactorization() {
  
  //  cout << " NumericFactorization() B = " << *(Problem_->GetRHS()) << endl ; 

  PerformNumericFactorization( );
  return 0;
}


int Amesos_Klu::Solve() { 

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

  //
  //  Call KLU to perform the solve
  //

  
  SerialX->Scale(1.0, *SerialB) ;

  int SerialXlda ; 
  if ( iam == 0 ) {
    assert( SerialX->ExtractView( &SerialXvalues, &SerialXlda ) == 0 ) ; 

    assert( SerialXlda == NumGlobalElements_ ) ; 
    
    for ( int j =0 ; j < nrhs; j++ ) { 
      double *Control = (double *) NULL, *Info = (double *) NULL ;

      int status = 0 ; 
	//  kludge debugxx knen fix this:  
      vector<double> workspace( 100* (10+NumGlobalElements_) ) ; 
      klu_btf_solve( Symbolic_, Numeric_,  
		     &SerialXvalues[j*SerialXlda], &workspace[0] );
    }
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
  return(0) ; 
}
