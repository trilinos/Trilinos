//  As of July 1st, USE_STL_SORT and USE_LOCAL both work together or separatelys
//  (But you have to set at least one)
//  #define USE_STL_SORT
#define USE_LOCAL

#ifndef USE_LOCAL
#ifndef USE_STL_SORT
At present, either USE_LOCAL or USE_STL_SORT is required
#endif
#endif

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

#include "Amesos_Umfpack.h"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_Util.h"
#ifdef USE_STL_SORT
#include <algorithm>
#endif

  //=============================================================================
  Amesos_Umfpack::Amesos_Umfpack(const Epetra_LinearProblem &prob, 
				 const AMESOS::Parameter::List &ParameterList ) :  
    Symbolic(0),
    Numeric(0),
    SerialCrsMatrixA_(0), 
    SerialMap_(0), 
    SerialMatrix_(0), 
    SymbolicFactorizationOK_(false), 
    NumericFactorizationOK_(false)  {


  Problem_ = &prob ; 
  ParameterList_ = &ParameterList ; 
}

//=============================================================================
Amesos_Umfpack::~Amesos_Umfpack(void) {

  if ( SerialMap_ ) delete SerialMap_ ; 
  if ( SerialCrsMatrixA_ ) delete SerialCrsMatrixA_ ; 
  if ( Symbolic ) umfpack_di_free_symbolic (&Symbolic) ;
  if ( Numeric ) umfpack_di_free_numeric (&Numeric) ;
}

int Amesos_Umfpack::ConvertToSerial() { 
  
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
    SerialMap_ = new Epetra_Map( NumGlobalElements_, NumMyElements_, 0, Comm() );

    Epetra_Export export_to_serial( OriginalMap, *SerialMap_);

    SerialCrsMatrixA_ = new Epetra_CrsMatrix(Copy, *SerialMap_, 0);
    SerialCrsMatrixA_->Export( *CastCrsMatrixA, export_to_serial, Add ); 
    
    SerialCrsMatrixA_->TransformToLocal() ; 
    SerialMatrix_ = SerialCrsMatrixA_ ;
  }

  return 0;
} 

int Amesos_Umfpack::ConvertToUmfpackCRS(){
  
  //
  //  Convert matrix to the form that Umfpack expects (Ap, Ai, Aval) 
  //

  assert( NumGlobalElements_ == SerialMatrix_->NumGlobalRows());
  assert( NumGlobalElements_ == SerialMatrix_->NumGlobalCols());
  assert( numentries_ == SerialMatrix_->NumGlobalNonzeros());
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
      assert( SerialMatrix_->ExtractMyRowView( MyRow, NumEntriesThisRow, RowValues, ColIndices ) == 0 ) ;
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


int Amesos_Umfpack::PerformSymbolicFactorization() {

  double *Control = (double *) NULL, *Info = (double *) NULL ;
  
  //    if ( Symbolic ) umfpack_di_free_symbolic (&Symbolic) ;
  if ( iam== 0 ) {
    (void) umfpack_di_symbolic (NumGlobalElements_, NumGlobalElements_, &Ap[0], 
				&Ai[0], &Aval[0], 
				&Symbolic, Control, Info) ;
  }
  SymbolicFactorizationOK_ = true ; 
  return 0;
}

int Amesos_Umfpack::PerformNumericFactorization( ) {

  if ( iam == 0 ) {
    double *Control = (double *) NULL, *Info = (double *) NULL ;
    if (Numeric) umfpack_di_free_numeric (&Numeric) ;
    int status = umfpack_di_numeric (&Ap[0], &Ai[0], 
				 &Aval[0], Symbolic, 
				 &Numeric, Control, Info) ;
    assert( status == 0 ) ; 
  }
  
  NumericFactorizationOK_ = true ; 
  return 0;
}




bool Amesos_Umfpack::MatrixShapeOK() const { 
  bool OK ;

  if ( GetProblem()->GetOperator()->OperatorRangeMap().NumGlobalPoints() != 
       GetProblem()->GetOperator()->OperatorDomainMap().NumGlobalPoints() ) OK = false;
  return OK; 
}


int Amesos_Umfpack::SymbolicFactorization() {

  ConvertToSerial() ; 
  
  ConvertToUmfpackCRS();
  
  PerformSymbolicFactorization();

  NumericFactorizationOK_ = false; 
  return 0;
}

int Amesos_Umfpack::NumericFactorization() {
  
  ConvertToSerial() ; 
  
  ConvertToUmfpackCRS();
  
  if ( ! SymbolicFactorizationOK_ ) {
    PerformSymbolicFactorization();
  }

  PerformNumericFactorization( );
  return 0;
}


int Amesos_Umfpack::Solve() { 

  if ( ! ( SymbolicFactorizationOK_ &&  NumericFactorizationOK_ ) ) {
    ConvertToSerial() ; 
  
    ConvertToUmfpackCRS();
  }

  if ( ! SymbolicFactorizationOK_ ) {
    PerformSymbolicFactorization();
    assert( ! NumericFactorizationOK_ );  // Can't redo Symbolic Phase and not the Numeric
  }

  if ( ! NumericFactorizationOK_ ) PerformNumericFactorization( );

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

  Epetra_MultiVector *SerialB, *SerialX; 
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
  //  Call UMFPACK to perform the solve
  //  Note:  UMFPACK uses a Compressed Column Storage instead of compressed row storage, 
  //  Hence to compute A X = B, we ask UMFPACK to perform A^T X = B and vice versa
  //


  int SerialBlda, SerialXlda ; 
  int UmfpackRequest = UseTranspose()?UMFPACK_A:UMFPACK_At ;
  if ( iam == 0 ) {
    assert( SerialB->ExtractView( &SerialBvalues, &SerialBlda ) == 0 ) ; 
    assert( SerialX->ExtractView( &SerialXvalues, &SerialXlda ) == 0 ) ; 
    assert( SerialBlda == NumGlobalElements_ ) ; 
    assert( SerialXlda == NumGlobalElements_ ) ; 
    
#if 0
    for ( int j=0; j < 1 ; j++ ) {
      for ( int i =0; i < NumGlobalElements_ ; i++ ) {
	SerialBvalues[i+SerialBlda*j] = 10*j + i ; 
	SerialXvalues[i+SerialXlda*j] = 100*j + 10*i ; 
      }
    }
#endif
    for ( int j =0 ; j < nrhs; j++ ) { 
      double *Control = (double *) NULL, *Info = (double *) NULL ;
      
      int status = umfpack_di_solve (UmfpackRequest, &Ap[0], 
				     &Ai[0], &Aval[0], 
				     &SerialXvalues[j*SerialXlda], 
				     &SerialBvalues[j*SerialBlda], 
				     Numeric, Control, Info) ;

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
