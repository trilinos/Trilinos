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



    Ken, I changed a little bit this class. Now almost all the MUMPS functionality are
    supported. It is not supported:
    - matrices in elemental format;
    - ICNTL values of 1 and 2. In fact:
      * ICNTL(18) = 3 is the default choice. I think it is the best for paralle
        applications;
      * ICTNL(0) can be used as follows: the matrix is given in distributed format
        by the user, but then he/she sets `SetKeepMatrixDistributed(false)'. The
	matrix will be converted to COO format, then those vectors are shipped to
	processor 0, and MUMPS takes are of using the processes.
    - timing can be obtained using functions GetSymFactTime, GetNumFactTime,
      and GetSolFactTime()

	
    NOTES:
    - This class should be used with MUMPS 4.3 or 4.3.1 (never tested with
      older versions of MUMPS, and developed with 4.3.1)
    - Never tested the reuse of symbolic factorization and change of
      numerical values. In particular the code assumes that the user is
      doing the right thing (i.e., it is not changing the structure of the
      matrix)
    - Still to do: Schur complement should be returned as a VBR matrix
      if the input matrix is VBR
    - return MFLOPS? required memory?

    Marzio Sala, SNL 9214, 24-Jan-03

   */


#include "Amesos_Mumps.h"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_Util.h"
#include "Epetra_Time.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_IntVector.h"
#include "Epetra_Vector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Amesos_EpetraInterface.h"
#include "Amesos_EpetraRedistributor.h"

#define ICNTL(I) icntl[(I)-1]
#define CNTL(I)  cntl[(I)-1]
#define INFOG(I) infog[(I)-1]
#define RINFOG(I) rinfog[(I)-1]

#define DEF_VALUE_INT -123456789
#define DEF_VALUE_DOUBLE -123456.789
  
//=============================================================================
Amesos_Mumps::Amesos_Mumps(const Epetra_LinearProblem &prob,
			   const AMESOS::Parameter::List &ParameterList ) :  
  Amesos_EpetraRedistributor(&prob),
  SymbolicFactorizationOK_(false), 
  NumericFactorizationOK_(false),
  KeepMatrixDistributed_(true) ,
  Map_(0),
  NumMUMPSNonzeros_(0),
  IsConvertToTripletOK_(false),
  Col(0),
  Row(0),
  Val(0),
  ErrorMsgLevel_(-1),
  NumSchurComplementRows_(-1),
  SchurComplementRows_(0),
  SchurComplement_(0),
  IsComputeSchurComplementOK_(false),
  RowSca_(0),
  ColSca_(0),
  PermIn_(0),
  Maxis_(DEF_VALUE_INT),
  Maxs_(DEF_VALUE_INT)
{
  ParameterList_ = &ParameterList ; 
  MyPID = Comm().MyPID() ;

  if( Comm().NumProc() == 1 ) {
    KeepMatrixDistributed_ = false;
  }
  
  // set to -1 icntl_ and cntl_. The use can override default values by using
  // SetICNTL(pos,value) and SetCNTL(pos,value).
  for( int i=0 ; i<40 ; ++i ) icntl_[i] = DEF_VALUE_INT;
  for( int i=0 ; i<5 ; ++i ) cntl_[i] = DEF_VALUE_DOUBLE;
}

//=============================================================================
Amesos_Mumps::~Amesos_Mumps(void)
{

  // destroy instance of the package
  MDS.job = -2;
  dmumps_c(&MDS);
  
  if( Row ) delete Row;
  if( Col ) delete Col;
  if( Val ) delete Val;

  if( IsComputeSchurComplementOK_ && MyPID == 0 ) {
    delete [] MDS.schur;
  }
  
}

//=============================================================================
int Amesos_Mumps::SetICNTL(int pos, int value)
{
  // NOTE: suppose first position is 1 (as in FORTRAN)
  if( pos>0 && pos<41 ) {
    icntl_[pos-1] = value;
    return 0;
  } else {
    return -1;
  }
}

//=============================================================================
int Amesos_Mumps::SetCNTL(int pos, double value)
{
  // NOTE: suppose first position is 1 (as in FORTRAN)
  if( pos>0 && pos<6 ) {
    cntl_[pos-1] = value;
    return 0;
  } else {
    return -1;
  }
}

int Amesos_Mumps::ConvertToTriplet()
{

  if( IsConvertToTripletOK_ ) return 0;
  
  NumMUMPSNonzeros_ = 0;
  
  // MS // convert to MUMPS format, keeping in distributed form.
  // MS // This doesn't require the matrix to be shipped to proc 0,
  // MS // then decomposed and reshipped by MUMPS

  // be sure that those vectors are empty. Otherwise when the user
  // call UpdateValues we will reallocate all of them
  if( Row ) delete Row;
  if( Col ) delete Col;
  if( Val ) delete Val;

  Row = new Epetra_IntSerialDenseVector( NumMyNonzeros() ) ;  
  Col = new Epetra_IntSerialDenseVector( NumMyNonzeros() ) ;  
  Val = new Epetra_SerialDenseVector( NumMyNonzeros() ) ;  

  // MS // retrive already allocated pointers for rows, cols, and vals
  // MS // Those arrays are allocated and freed by the Amesos_EpetraInterface
  // MS // class. Note that, if a View method is implemented, those pointers
  // MS // will change to reflect the action of the View. Otherwise,
  // MS // the allocated vectors will be used to copy data  

  int NumIndices;
  
  int * RowIndices = GetRowIndices();
  int * ColIndices = GetColIndices();
  double * MatrixValues = GetValues();
  
  // MS // do all in terms of BlockRows. This way we can support VBR matrices
  
  for( int LocalBlockRow=0; LocalBlockRow<NumMyBlockRows() ; ++LocalBlockRow ) {

    GetRow(LocalBlockRow,NumIndices,RowIndices,ColIndices,MatrixValues);

    for( int i=0 ; i<NumIndices ; ++i ) {
      //      if( MatrixValues[i] != 0.0 ) {
	(*Row)[NumMUMPSNonzeros_] = RowIndices[i]+1;
	(*Col)[NumMUMPSNonzeros_] = ColIndices[i]+1;
	(*Val)[NumMUMPSNonzeros_] = MatrixValues[i];
	NumMUMPSNonzeros_++;
	//      }
    }
  }

  assert(NumMUMPSNonzeros_<=NumMyNonzeros());

  // MS // bring matrix to proc zero if required
  
  if( KeepMatrixDistributed_ == false ) {

    Epetra_IntSerialDenseVector * OldRow = Row;
    Epetra_IntSerialDenseVector * OldCol = Col;
    Epetra_SerialDenseVector * OldVal = Val;

    int OldNumMUMPSNonzeros = NumMUMPSNonzeros_;
    Comm().SumAll(&OldNumMUMPSNonzeros,&NumMUMPSNonzeros_,1);

    if( MyPID == 0 ) {
      Row = new Epetra_IntSerialDenseVector(NumMUMPSNonzeros_);
      Col = new Epetra_IntSerialDenseVector(NumMUMPSNonzeros_);
      Val = new Epetra_SerialDenseVector(NumMUMPSNonzeros_);
    } else {
      NumMUMPSNonzeros_ = 0;
    } 

    Epetra_Map OldNnz(-1,OldNumMUMPSNonzeros,0,Comm());
    
    Epetra_Map NewNnz(-1,NumMUMPSNonzeros_,0,Comm());
    Epetra_Import OldToNew(NewNnz,OldNnz);

    Epetra_IntVector GOldRow(View,OldNnz,OldRow->Values());
    Epetra_IntVector GOldCol(View,OldNnz,OldCol->Values());
    Epetra_Vector GOldVal(View,OldNnz,OldVal->Values());

    Epetra_IntVector GRow(View,NewNnz,Row->Values());
    Epetra_IntVector GCol(View,NewNnz,Col->Values());
    Epetra_Vector GVal(View,NewNnz,Val->Values());
    
    GRow.Import(GOldRow,OldToNew,Insert);
    GCol.Import(GOldCol,OldToNew,Insert);
    GVal.Import(GOldVal,OldToNew,Insert);

    delete OldRow;
    delete OldCol;
    delete OldVal;

    if( MyPID != 0 ) {
      delete Row; Row = 0;
      delete Col; Col = 0;
      delete Val; Val = 0;
    }
  }

  IsConvertToTripletOK_ = true;
  
  return 0;

}   


int Amesos_Mumps::UpdateMatrixValues()
{

  CreateSerialMap();
  CreateImportAndExport();

  return( ConvertToTriplet() );

  NumericFactorizationOK_ = false;
  
}

int Amesos_Mumps::ReadParameterList() {
  if (ParameterList_->isParameterSublist("Mumps") ) {
    AMESOS::Parameter::List MumpsParams = ParameterList_->sublist("Mumps") ;
  }  
  return 0;
}

int Amesos_Mumps::PerformSymbolicFactorization()
{

  if( IsLocal() ) {
#ifdef EPETRA_MPI
    MPI_Comm MPIC = MPI_COMM_SELF ;
    MDS.comm_fortran = (F_INT) MPI_Comm_c2f( MPIC ) ;   // Compiled on cygwin but not on Atlantis
#else
    MDS.comm_fortran = -987654 ;  // Needed by MUMPS 4.3 
#endif
  }
  else
    MDS.comm_fortran = -987654 ;  // Needed by MUMPS 4.3 
  
  //  MDS.comm_fortran = (F_INT) MPIR_FromPointer( MPIC ) ;  // Other recommendation from the MUMPS manual - did not compile on Atlantis either

  MDS.job = -1  ;     //  Initialization
  MDS.par = 1 ;       //  Host IS involved in computations
  MDS.sym = MatrixProperty();
  dmumps_c( &MDS ) ;   //  Initialize MUMPS 

  // check the error flag. MUMPS is siccessful only if
  // infog(1) is zereo
  //  if( MDS.INFOG(1) ) return( MDS.INFOG(1) );
  
  MDS.n = NumGlobalRows() ;
  if( KeepMatrixDistributed_ ) {
    MDS.nz_loc = NumMUMPSNonzeros_;
    MDS.irn_loc = Row->Values(); 
    MDS.jcn_loc = Col->Values(); 
    MDS.a_loc = Val->Values();
  } else {
    if( MyPID == 0 ) {
      MDS.nz = NumMUMPSNonzeros_;
      MDS.irn = Row->Values(); 
      MDS.jcn = Col->Values(); 
      MDS.a = Val->Values();
    }
  }

  // scaling if provided by the user
  if( RowSca_ != 0 ) {
    MDS.rowsca = RowSca_;
    MDS.colsca = ColSca_;
  }

  // given ordering if provided by the user
  if( PermIn_ != 0 ) {
    MDS.perm_in = PermIn_;
  }

  //  if( Maxis_ != DEF_VALUE_INT ) MDS.maxis = Maxis_;
  //  if( Maxs_ != DEF_VALUE_INT ) MDS.maxs = Maxs_;
  
  MDS.job = 1  ;     // Request symbolic factorization

  SetICNTLandCNTL(); // initialize icntl and cntl. NOTE: I initialize those vectors
                     // here, and I don't change them anymore. This is not that much
                     // limitation, though

  dmumps_c( &MDS ) ;  // Perform symbolic factorization
  // check the error flag. MUMPS is siccessful only if
  // infog(1) is zereo
  //  if( MDS.INFOG(1) ) return( MDS.INFOG(1) );
  
  SymbolicFactorizationOK_ = true ; 
  return 0;

}

void Amesos_Mumps::SetICNTLandCNTL()
{
  
  MDS.ICNTL(1)  = -1;  // Turn off error messages
  MDS.ICNTL(2)  = -1;  // Turn off diagnostic printing
  MDS.ICNTL(3)  = -1;  // Turn off global information messages
  MDS.ICNTL(4)  = ErrorMsgLevel_;
  
  MDS.ICNTL(6)  = 7;   // Choose column permutation automatically
  MDS.ICNTL(7)  = 7;   // Choose ordering method automatically
  MDS.ICNTL(8)  = 7;   // Choose scaling automatically
  MDS.ICNTL(9)  = 1;   // Compute A^T x = b - changed in Solve()
  MDS.ICNTL(10) = 0;   // Maximum steps of iterative refinement
  MDS.ICNTL(11) = 0;   // Do not collect statistics
  MDS.ICNTL(12) = 0;   // Use Node level parallelism
  MDS.ICNTL(13) = 0;   // Use ScaLAPACK for root node 
  MDS.ICNTL(14) = 20;  // Increase memory allocation 20% at a time 
  MDS.ICNTL(15) = 0;   // Minimize memory use (not flop count)
  MDS.ICNTL(16) = 0;   // Do not perform null space detection
  MDS.ICNTL(17) = 0;   // Unused (null space dimension)
  
  if( KeepMatrixDistributed_ ) MDS.ICNTL(18)= 3;
  else                         MDS.ICNTL(18)= 0;

  if ( UseTranspose() )  MDS.ICNTL(9) = 0 ; 
  else                   MDS.ICNTL(9) = 1 ; 

  if( IsComputeSchurComplementOK_ ) MDS.ICNTL(19) = 1;
  else                              MDS.ICNTL(19) = 0;

  // something to do if the Schur complement is required.
  if( IsComputeSchurComplementOK_ && MyPID == 0 ) {
    MDS.size_schur = NumSchurComplementRows_;
    MDS.listvar_schur = SchurComplementRows_;
    MDS.schur = new double[NumSchurComplementRows_*NumSchurComplementRows_];
  }
  
  // retrive user's specified options
  for( int i=0 ; i<40 ; ++i ) {
    if( icntl_[i] != DEF_VALUE_INT )
      MDS.ICNTL(i+1) = icntl_[i];
  }
  for( int i=0 ; i<5 ; ++i ) {
    if( icntl_[i] != DEF_VALUE_DOUBLE )
      MDS.CNTL(i+1) = cntl_[i];
  }
  
  // take care that required options are not overwritten
  assert(MDS.ICNTL(5)== 0);  // Matrix is given in elemental (i.e. triplet) from
  
  return;
  
}

int Amesos_Mumps::PerformNumericFactorization( )
{
  
  MDS.job = 2  ;     // Request numeric factorization
  dmumps_c( &MDS ) ;  // Perform numeric factorization
  // check the error flag. MUMPS is siccessful only if
  // infog(1) is zereo
  //  if( MDS.INFOG(1) ) return( MDS.INFOG(1) );
  
  NumericFactorizationOK_ = true ;

  return 0;
}


bool Amesos_Mumps::MatrixShapeOK() const { 
  bool OK ;

  if ( GetProblem()->GetOperator()->OperatorRangeMap().NumGlobalPoints() != 
       GetProblem()->GetOperator()->OperatorDomainMap().NumGlobalPoints() ) OK = false;
  return OK; 
}


int Amesos_Mumps::SymbolicFactorization()
{

  CreateSerialMap();
  CreateImportAndExport();
  ConvertToTriplet();

  Epetra_Time Time(Comm());
  
  PerformSymbolicFactorization();

  NumericFactorizationOK_ = false;

  AddToSymFactTime(Time.ElapsedTime());
  
  return 0;
}

int Amesos_Mumps::NumericFactorization() {

  CreateSerialMap();
  CreateImportAndExport();
  ConvertToTriplet();

  if( SymbolicFactorizationOK_ == false )  PerformSymbolicFactorization();

  Epetra_Time Time(Comm());
  PerformNumericFactorization( );
  AddToNumFactTime(Time.ElapsedTime());
  
  return 0;
}

int Amesos_Mumps::Solve()
{ 

  CreateSerialMap();
  CreateImportAndExport();
  ConvertToTriplet();

  if( SymbolicFactorizationOK_ == false )  PerformSymbolicFactorization();
  if( NumericFactorizationOK_  == false )  PerformSymbolicFactorization();
  
  Epetra_Time Time(Comm());
  
  Epetra_MultiVector   *vecX = GetProblem()->GetLHS() ; 
  Epetra_MultiVector   *vecB = GetProblem()->GetRHS() ;
      
  //
  //  Compute the number of right hands sides (and check that X and B have the same shape) 
  //

  if( vecB == NULL || vecX == NULL ) {
    if( MyPID == 0 ) cout << "Amesos ERROR : LHS or RHS not set" << endl;
    return -1;
  }

  int nrhs = vecX->NumVectors() ;
  if( nrhs != vecB->NumVectors() ) {
    if( MyPID == 0 ) cout << "Amesos ERROR : multivectors LHS or RHS have different sizes" << endl;
    return -1;
  }

  // will be used only if IsLocal == false
  CreateSerialRHS( nrhs );
  
  // MS // Extract Serial versions of X and B. 
  // MS // NOTE: both serial and distributed version of MUMPS require
  // MS // rhs to be entirely stored on processor 0. This is because we
  // MS // still need SerialMap_ also for distributed input

  if( IsLocal() ) {
    if( MyPID == 0 ) {
      for ( int j =0 ; j < nrhs; j++ ) {
	for( int i=0 ; i<NumMyRows() ; ++i ) (*vecX)[j][i] = (*vecB)[j][i];
	MDS.rhs = (*vecX)[j];
	MDS.job = 3  ;     // Request solve
	dmumps_c( &MDS ) ;  // Perform solve
	// check the error flag. MUMPS is siccessful only if
	// infog(1) is zereo
	//	if( MDS.INFOG(1) ) return( MDS.INFOG(1) );

      }
    }
  } else {

    SerialRHS()->Import( *vecB, *(ImportToProcZero()), Insert ) ;

    for ( int j =0 ; j < nrhs; j++ ) { 
      if ( MyPID == 0 ) {
	MDS.rhs = (*SerialRHS())[j];
      }
      MDS.job = 3  ;     // Request solve
      dmumps_c( &MDS ) ;  // Perform solve
      // check the error flag. MUMPS is siccessful only if
      // infog(1) is zereo
      //      if( MDS.INFOG(1) ) return( MDS.INFOG(1) );
    }

    vecX->Import( *(SerialRHS()), *(ImportFromProcZero()), Insert ) ;
    
  }
  
  // Retrive the SC and put it in an Epetra_CrsMatrix on host only
  if( IsComputeSchurComplementOK_ ) {
    
    if( MyPID != 0 ) NumSchurComplementRows_ = 0;
    
    Epetra_Map SCMap(-1,NumSchurComplementRows_, 0, Comm());

    SchurComplement_ = new Epetra_CrsMatrix(Copy,SCMap,NumSchurComplementRows_);

    if( MyPID == 0 )
      for( int i=0 ; i<NumSchurComplementRows_ ; ++i ) {
	for( int j=0 ; j<NumSchurComplementRows_ ; ++j ) {
	  int pos = i+ j *NumSchurComplementRows_;
	  SchurComplement_->InsertGlobalValues(i,1,&(MDS.schur[pos]),&j);
	}
      }
    
    SchurComplement_->FillComplete();
  }
  
  AddToSolFactTime(Time.ElapsedTime());
  
  return(0) ; 
}


int Amesos_Mumps::ComputeSchurComplement(bool flag, int NumSchurComplementRows,
					 int * SchurComplementRows)
{
  NumSchurComplementRows_ = NumSchurComplementRows;
  
  SchurComplementRows_ = SchurComplementRows;

  // modify because MUMPS is fortran-driven
  if( MyPID == 0 )
    for( int i=0 ; i<NumSchurComplementRows ; ++i ) SchurComplementRows_[i]++;
  
  IsComputeSchurComplementOK_ = flag;

  return 0;
}

int Amesos_Mumps::PrintInformation() 
{
  if( Comm().MyPID() ) return 0;

  cout << "------------ : ---------------------------------" << endl;
  cout << "Amesos_Mumps : global information for all phases" << endl;
  cout << "------------ : ---------------------------------" << endl;
  cout << "Amesos_Mumps : Matrix has " << NumGlobalRows() << " rows"
       << " and " << NumGlobalNonzeros() << " nonzeros" << endl;

  cout << "Amesos_Mumps : estimated FLOPS for elimination = "
       << MDS.RINFOG(1) << endl;
  cout << "Amesos_Mumps : total FLOPS for assembly = "
       << MDS.RINFOG(2) << endl;
  cout << "Amesos_Mumps : total FLOPS for elimination = "
       << MDS.RINFOG(3) << endl;
  
  cout << "Amesos_Mumps : total real space to store the LU factors = "
       << MDS.INFOG(9) << endl;
  cout << "Amesos_Mumps : total integer space to store the LU factors = "
       << MDS.INFOG(10) << endl;
  cout << "Amesos_Mumps : total number of iterative steps refinement = "
       << MDS.INFOG(15) << endl;
  cout << "Amesos_Mumps : estimated size of MUMPS internal data\n"
       << "Amesos_Mumps : for running factorization = "
       << MDS.INFOG(16) << " Mbytes" << endl;
  cout << "Amesos_Mumps : estimated size of MUMPS internal data\n"
       << "Amesos_Mumps : for running factorization = "
       << MDS.INFOG(17) << " Mbytes" << endl;
  cout << "Amesos_Mumps : estimated size of MUMPS internal data\n"
       << "Amesos_Mumps : allocated during factorization = "
       << MDS.INFOG(19) << " Mbytes" << endl;
  cout << "------------ : ---------------------------------" << endl;

  return 0;
  
}
