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
    Matrix_(0) {


  Problem_ = &prob ; 
  ParameterList_ = &ParameterList ; 
}

//=============================================================================
Amesos_Klu::~Amesos_Klu(void) {

  if ( SerialMap_ ) delete SerialMap_ ; 
  if ( SerialCrsMatrixA_ ) delete SerialCrsMatrixA_ ; 
  if ( TransposeMatrix_ ) delete TransposeMatrix_ ; 
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

  assert( NumGlobalElements_ == Matrix_->NumGlobalRows());
  assert( NumGlobalElements_ == Matrix_->NumGlobalCols());
  assert( numentries_ == Matrix_->NumGlobalNonzeros());
  Ap.resize( NumGlobalElements_+1 );
  Ai.resize( EPETRA_MAX( NumGlobalElements_, numentries_) ) ; 
  Aval.resize( EPETRA_MAX( NumGlobalElements_, numentries_) ) ; 

  if ( iam==0 ) {
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

#ifdef USE_AMD
  double *Control = (double *) NULL, *Info = (double *) NULL ;
  
  if ( iam== 0 ) {
    amd_defaults( Control ) ; 
    permuation.resize(NumGlobalElements_);
    int status = amd_order (NumGlobalElements_, &Ap[0], 
			    &Ai[0], &Permutation_[0],  
				Control, Info) ;
    EPETRA_CHK_ERR( status ) ; 
  }
#endif 
  return 0;
}

int Amesos_Klu::PerformNumericFactorization( ) {

  if ( iam == 0 ) {
    double Klu_Control[KLU_CONTROL];
    
    klu_defaults( Klu_Control ) ; 
    int status = klu( NumGlobalElements_, &Ap[0], &Ai[0], 
		      &Aval[0], Klu_Control, 
		      &Lp, &Li, &Lx, &Up, &Ui, &Ux, &P) ;

#if 0
    if (Numeric) umfpack_di_free_numeric (&Numeric) ;
    int status = umfpack_di_numeric (&Ap[0], &Ai[0], 
				 &Aval[0], Symbolic, 
				 &Numeric, Control, Info) ;
#endif
    EPETRA_CHK_ERR( status ) ; 
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

  double *SerialBvalues ;
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


  int SerialBlda, SerialXlda ; 
  if ( iam == 0 ) {
    assert( SerialB->ExtractView( &SerialBvalues, &SerialBlda ) == 0 ) ; 
    assert( SerialX->ExtractView( &SerialXvalues, &SerialXlda ) == 0 ) ; 
    assert( SerialBlda == NumGlobalElements_ ) ; 
    assert( SerialXlda == NumGlobalElements_ ) ; 
    
    for ( int j =0 ; j < nrhs; j++ ) { 
      double *Control = (double *) NULL, *Info = (double *) NULL ;

      int status = 0 ; 
#if 0
      cout << " B = " << 
	SerialBvalues[0] << " " <<
	SerialBvalues[1] << " " <<
	SerialBvalues[2] << " " <<
	SerialBvalues[3] << " " <<
	SerialBvalues[4] << " at the top" << endl ; 
      cout << " P = " << P[0] << " " << 
	P[1] << " " << 
	P[2] << " " << 
	P[3] << " " << 
	P[4] << endl ; 
#endif
#if 0
      defined in klu_dump
      assert( klu_valid( NumGlobalElements_, Lp, Li, Lx ) ) ; 
      assert( klu_valid( NumGlobalElements_, Up, Ui, Ux ) ) ; 
#endif

      klu_permute( NumGlobalElements_, P, 
		   &SerialBvalues[j*SerialBlda],
		   &SerialXvalues[j*SerialXlda] ) ;
#if 0
      cout << " X = " << 
	SerialXvalues[0] << " " <<
	SerialXvalues[1] << " " <<
	SerialXvalues[2] << " " <<
	SerialXvalues[3] << " " <<
	SerialXvalues[4] << " before the left solve" << endl ; 
#endif
      klu_lsolve(  NumGlobalElements_, 
		   Lp, Li, Lx, &SerialXvalues[j*SerialXlda] );
#if 0
      cout << " X = " << 
	SerialXvalues[0] << " " <<
	SerialXvalues[1] << " " <<
	SerialXvalues[2] << " " <<
	SerialXvalues[3] << " " <<
	SerialXvalues[4] << " before the right solve" << endl ; 
#endif
      klu_usolve(  NumGlobalElements_, 
		   Up, Ui, Ux, &SerialXvalues[j*SerialXlda] );
	
#if 0			    
      cout << " X = " << 
	SerialXvalues[0] << " " <<
	SerialXvalues[1] << " " <<
	SerialXvalues[2] << " " <<
	SerialXvalues[3] << " " <<
	SerialXvalues[4] << " after the right solve" << endl ; 
#endif
#if 0
      int status = umfpack_di_solve (KluRequest, &Ap[0], 
				     &Ai[0], &Aval[0], 
				     &SerialXvalues[j*SerialXlda], 
				     &SerialBvalues[j*SerialBlda], 
				     Numeric, Control, Info) ;
#endif
#if 0
      for ( int k=0; k < nrhs ; k++ ) {
	for ( int i =0; i < NumGlobalElements_ ; i++ ) {
	  cout << "h" << j <<"X( " << i+1 << "," << k+1 << ") = " << SerialXvalues[i+SerialXlda*k] << ";" << endl ; 
	}
      }
      for ( int k=0; k < nrhs ; k++ ) {
	for ( int i =0; i < NumGlobalElements_ ; i++ ) {
	  cout << "h" << j <<"B( " << i+1 << "," << k+1 << ") = " << SerialBvalues[i+SerialBlda*k] << ";" << endl ; 
	}
      }
#endif
      EPETRA_CHK_ERR( status ) ; 
    }
  }
    

#if 0
  Comm().Barrier();
  if  (iam == 0 ) { 
    cout << " SerialXlda = " << SerialXlda << endl ; 
    cout << " SerialBlda = " << SerialBlda << endl ; 
  //  Print for matlab 
  //
  for (int i = 0; i < NumGlobalElements_ ; i++ ) { 
    for ( int j = Ap[i]; j < Ap[i+1] ; j++ ) { 
      cout << "A(" << i +1  << "," << Ai[j]+1 << " ) = " << Aval[j] << "; % iam = " << iam <<endl ; 
    } 
  }
  for ( int j=0; j < nrhs ; j++ ) {
    for ( int i =0; i < NumGlobalElements_ ; i++ ) {
      cout << "X( " << i+1 << "," << j+1 << ") = " << SerialXvalues[i+SerialXlda*j] << ";" << endl ; 
    }
  }
  for ( int j=0; j < nrhs ; j++ ) {
    for ( int i =0; i < NumGlobalElements_ ; i++ ) {
      cout << "B( " << i+1 << "," << j+1 << ") = " << SerialBvalues[i+SerialBlda*j] << ";" << endl ; 
    }
  }
  }
#endif
   
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
#if 0
  cout << " Here is SerialB " << endl ; 
  SerialB->Print(cout ) ; 
  cout << " There was SerialB " << endl ; 
  
  cout << " Here is SerialX " << endl ; 
  SerialX->Print(cout ) ; 
  cout << " There was SerialX " << endl ; 
  
  cout << " Here is VecX " << endl ; 
  vecX->Print(cout ) ; 
  cout << " There was VecX " << endl ; 
#endif  

  return(0) ; 
}
