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

  /*
    Information about MUMPS:

    MUMPS will not perform and column permutation on a matrix provided
    in distributed form.  Amesos_Mumps will match this, allowing
    column permutation only if the matrix is provided in serial form.
    This is unfortunate because it is an exception to the general rule
    that the capability (given adequate memory) of any class
    implementing the Amesos_BaseSolver base class does not depend on
    the distribution of the input matrix.  However, neither of the
    other options are attractive.  Coalescing the matrix to a single
    process independent of whether column permutation is requested
    unnecessarily limits the size problem that can be solved.
    Coalescing the matrix to a single process only when column
    permutation is requested would cause some problems to run out of memory
    when column permutation is requested.

    QUESTION:  Should ICNTL(18) be set to 2 or to 3?  

    I think that ICNTL(18) == 2 is the right choice.  That ends up requiring the
    most data redistribution, and it requires that the matrix structure be stored
    on a single process.  But, it allows the symbolic factorization to be reused.

    Plan of attack:
    1)  First make it work with ICNTL(18) == 0 
    2)  Then move on to ICNTL(18) == 2 

   */


#include "Amesos_Mumps.h"
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
  Amesos_Mumps::Amesos_Mumps(const Epetra_LinearProblem &prob, 
				 const AMESOS::Parameter::List &ParameterList ) :  
    SerialCrsMatrixA_(0), 
    SerialMap_(0), 
    SerialMatrix_(0), 
    SymbolicFactorizationOK_(false), 
    NumericFactorizationOK_(false)  {


  Problem_ = &prob ; 
  ParameterList_ = &ParameterList ; 
}

//=============================================================================
Amesos_Mumps::~Amesos_Mumps(void) {

  if ( SerialMap_ ) delete SerialMap_ ; 
  if ( SerialCrsMatrixA_ ) delete SerialCrsMatrixA_ ; 
}

int Amesos_Mumps::ConvertToSerial() { 
  
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

int Amesos_Mumps::ConvertToTriplet(){
  
  //
  //  Convert matrix to the form that MUMPS expects:  Row, Col, Val
  //

  assert( NumGlobalElements_ == SerialMatrix_->NumGlobalRows());
  assert( NumGlobalElements_ == SerialMatrix_->NumGlobalCols());
  assert( numentries_ == SerialMatrix_->NumGlobalNonzeros());
  Row.resize( numentries_ ) ;  
  Col.resize( numentries_ ) ;  
  Val.resize( numentries_ ) ;  

  if ( iam==0 ) {
    int NumEntriesThisRow;
    double *RowValues;
    int *ColIndices;
    int Ai_index = 0 ; 
    for ( int MyRow = 0; MyRow <NumGlobalElements_; MyRow++ ) {
      assert( SerialMatrix_->ExtractMyRowView( MyRow, NumEntriesThisRow, RowValues, ColIndices ) == 0 ) ;
      for ( int j = 0; j < NumEntriesThisRow; j++ ) { 
	Row[Ai_index] = MyRow   + 1  ;   // "+1" converts to Fortran indices
	Col[Ai_index] = ColIndices[j]   +1 ;   
	Val[Ai_index] = RowValues[j] ; 
	Ai_index++;
      }
    }
  }

  
  return 0;
}   


int Amesos_Mumps::PerformSymbolicFactorization() {

  //
  Epetra_RowMatrix *RowMatrixA = dynamic_cast<Epetra_RowMatrix *>(Problem_->GetOperator());
  EPETRA_CHK_ERR( RowMatrixA == 0 ) ; 
  Epetra_CrsMatrix *CastCrsMatrixA = dynamic_cast<Epetra_CrsMatrix*>(RowMatrixA) ; 
  EPETRA_CHK_ERR( CastCrsMatrixA == 0 ) ; 
  const Epetra_Comm &Comm = RowMatrixA->Comm();
  const Epetra_MpiComm & comm1 = dynamic_cast<const Epetra_MpiComm &> (Comm);
//  MPI_Comm MPIC = comm1.Comm() ;
//  MDS.comm_fortran = (F_INT) MPI_Comm_c2f( MPIC ) ;   // Compiled on cygwin but not on Atlantis

//  MDS.comm_fortran = (F_INT) MPIR_FromPointer( MPIC ) ;  // Other recommendation from the MUMPS manual - did not compile on Atlantis either
  MDS.job = -1  ;     //  Initialization
  MDS.par = 1 ;       //  Host IS involved in computations
  MDS.sym = 0 ;       //  Matrix is not symmetric
  dmumps_c( &MDS ) ;   //  Initialize MUMPS 





  MDS.n = NumGlobalElements_ ; 
  MDS.nz = numentries_ ; 
  MDS.irn = &Row[0]; 
  MDS.jcn = &Col[0]; 
  MDS.a = &Val[0]; 
#define ICNTL(I) icntl[(I)-1]
  MDS.ICNTL(1)= -1;  // Turn off error messages
  MDS.ICNTL(2)= -1;  // Turn off diagnostic printing
  MDS.ICNTL(3)= -1;  // Turn off global information messages
  MDS.ICNTL(4)= 0;   // No messages output
  assert(MDS.ICNTL(5)== 0);  // Matrix is given in elemental (i.e. triplet) from 
  assert(MDS.ICNTL(6)== 7);   // Choose column permutation automatically
  assert(MDS.ICNTL(7)== 7);   // Choose ordering method automatically
  assert(MDS.ICNTL(8)== 7);   // Choose scaling automatically
  assert(MDS.ICNTL(9)== 1);   // Compute A^T x = b - changed in Solve()
  assert(MDS.ICNTL(10)==0);   // Maximum steps of iterative refinement
  assert(MDS.ICNTL(11)== 0);  // Do not collect statistics
  assert(MDS.ICNTL(12)== 0);  // Use Node level parallelism
  assert(MDS.ICNTL(13)== 0);  // Use ScaLAPACK for root node 
  assert(MDS.ICNTL(14)== 20); // Increase memory allocation 20% at a time 
  assert(MDS.ICNTL(15)== 0);  // Minimize memory use (not flop count)
  assert(MDS.ICNTL(16)== 0);  // Do not perform null space detection
  assert(MDS.ICNTL(17)== 0);  // Unused (null space dimension)
  assert(MDS.ICNTL(18)== 0);  // Pass matrix as a serial matrix
  assert(MDS.ICNTL(19)== 0);  // Do not compute the Schur complement
  MDS.job = 1  ;     // Request symbolic factorization
  dmumps_c( &MDS ) ;  // Perform symbolic factorization

  SymbolicFactorizationOK_ = true ; 
  return 0;
}

int Amesos_Mumps::PerformNumericFactorization( ) {

  MDS.job = 2  ;     // Request numeric factorization
  dmumps_c( &MDS ) ;  // Perform numeric factorization

  NumericFactorizationOK_ = true ; 
  return 0;
}




bool Amesos_Mumps::MatrixShapeOK() const { 
  bool OK ;

  if ( GetProblem()->GetOperator()->OperatorRangeMap().NumGlobalPoints() != 
       GetProblem()->GetOperator()->OperatorDomainMap().NumGlobalPoints() ) OK = false;
  return OK; 
}


int Amesos_Mumps::SymbolicFactorization() {

  ConvertToSerial() ; 
  
  ConvertToTriplet();
  
  PerformSymbolicFactorization();

  NumericFactorizationOK_ = false; 
  return 0;
}

int Amesos_Mumps::NumericFactorization() {
  
  ConvertToSerial() ; 
  
  ConvertToTriplet();
  
  if ( ! SymbolicFactorizationOK_ ) {
    PerformSymbolicFactorization();
  }

  PerformNumericFactorization( );
  return 0;
}


int Amesos_Mumps::Solve() { 

  if ( ! ( SymbolicFactorizationOK_ &&  NumericFactorizationOK_ ) ) {
    ConvertToSerial() ; 
  
    ConvertToTriplet();
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
  //  Call MUMPS to perform the solve
  //


  int SerialBlda, SerialXlda ; 
  if ( UseTranspose() ) 
    MDS.ICNTL(9) = 0 ; 
  else
    MDS.ICNTL(9) = 1 ; 

  assert( SerialB->ExtractView( &SerialBvalues, &SerialBlda ) == 0 ) ; 
  assert( SerialX->ExtractView( &SerialXvalues, &SerialXlda ) == 0 ) ; 
  assert( iam != 0 ||  SerialBlda == NumGlobalElements_ ) ; 
  assert(  iam != 0 || SerialXlda == NumGlobalElements_ ) ; 
  MDS.job = 3  ;     // Request solve
  
  for ( int j =0 ; j < nrhs; j++ ) { 
	if ( iam == 0 ) { 
      for ( int i=0; i< NumGlobalElements_ ; i++ ) 
        SerialXvalues[j*SerialXlda+i] = SerialBvalues[j*SerialXlda+i];
      MDS.rhs =  &SerialXvalues[j*SerialXlda];
}
    dmumps_c( &MDS ) ;  // Perform solve
    
//  #define PRINTALL
#ifdef PRINTALL
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
  }

#ifdef PRINTALL
  Comm().Barrier();
  if  (iam == 0 ) { 
    cout << " SerialXlda = " << SerialXlda << endl ; 
    cout << " SerialBlda = " << SerialBlda << endl ; 
    //  Print for matlab 
    //
    for (int i = 0; i < numentries_ ; i++ ) { 
	cout << "A(" << Row[i] << "," << Col[i] << " ) = " << Val[i] << "; % iam = " << iam <<endl ; 
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
  
#ifdef PRINTALL
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
