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
#include "Amesos_EpetraBaseSolver.h"
#include "EpetraExt_Redistor.h"

#define ICNTL(I) icntl[(I)-1]
#define CNTL(I)  cntl[(I)-1]
#define INFOG(I) infog[(I)-1]
#define RINFOG(I) rinfog[(I)-1]

#define DEF_VALUE_INT -123456789
#define DEF_VALUE_DOUBLE -123456.789
  
//=============================================================================

Amesos_Mumps::Amesos_Mumps(const Epetra_LinearProblem &prob,
			   const Teuchos::ParameterList &ParameterList ) :
  Amesos_EpetraBaseSolver(prob,ParameterList),
  SymbolicFactorizationOK_(false), 
  NumericFactorizationOK_(false),
  KeepMatrixDistributed_(true) ,
  UseMpiCommSelf_(false),
  Map_(0),
  NumMUMPSNonzeros_(0),
  NumMyMUMPSNonzeros_(0),
  Col(0),
  Row(0),
  Val(0),
  ErrorMsgLevel_(-1),
  NumSchurComplementRows_(-1),
  SchurComplementRows_(0),
  CrsSchurComplement_(0),
  DenseSchurComplement_(0),
  IsComputeSchurComplementOK_(false),
  RowSca_(0),
  ColSca_(0),
  PermIn_(0),
  Maxis_(DEF_VALUE_INT),
  Maxs_(DEF_VALUE_INT),
  OldMatrix_(0),
  Redistor_(0),
  TargetVector_(0),
  verbose_(0),
  AddZeroToDiag_(false),
  PrintTiming_(true),
  PrintStatistics_(true),
  Threshold_(0.0)
{
  // -777 is for me. It means : never called MUMPS, so
  // SymbolicFactorization will not call Destroy();
  MDS.job = -777;
  
  // set to -1 icntl_ and cntl_. The use can override default values by using
  // SetICNTL(pos,value) and SetCNTL(pos,value).
  for( int i=0 ; i<40 ; ++i ) icntl_[i] = DEF_VALUE_INT;
  for( int i=0 ; i<5 ; ++i ) cntl_[i] = DEF_VALUE_DOUBLE;
}

//=============================================================================

int Amesos_Mumps::Destroy()
{
  if( verbose_ == 3 ) cout << "Calling MUMPS with job=-2..." << endl;
  
  // destroy instance of the package
  MDS.job = -2;

  dmumps_c(&MDS);

  if( Redistor_ ) delete Redistor_;
  if( TargetVector_ ) delete TargetVector_;
  
  return 0;
}

//=============================================================================

Amesos_Mumps::~Amesos_Mumps(void)
{

  Destroy();
  
  if( Row ) delete Row;
  if( Col ) delete Col;
  if( Val ) delete Val;

  if( IsComputeSchurComplementOK_ && Comm().MyPID() == 0 ) {
    delete [] MDS.schur;
  }
  
}

//=============================================================================

int Amesos_Mumps::ConvertToTriplet()
{

  if( verbose_ == 3 ) cout << "Entering `ConvertToTriplet'" << endl;
  
  Epetra_Time Time(Comm());

  // MS // convert to MUMPS format, keeping in distributed form.
  // MS // This doesn't require the matrix to be shipped to proc 0,
  // MS // then decomposed and reshipped by MUMPS

  // MS // be sure that those vectors are empty. Allocate them
  
  if( Row ) delete Row;
  if( Col ) delete Col;
  if( Val ) delete Val;

  Row = new Epetra_IntSerialDenseVector( NumMyNonzeros() ) ;  
  Col = new Epetra_IntSerialDenseVector( NumMyNonzeros() ) ;  
  Val = new Epetra_SerialDenseVector( NumMyNonzeros() ) ;  

  NumMyMUMPSNonzeros_ = 0;
  
  // MS // retrive already allocated pointers for rows, cols, and vals
  // MS // Those arrays are allocated and freed by the Amesos_EpetraBaseSolver
  // MS // class. Note that, if a View method is implemented, those pointers
  // MS // will change to reflect the action of the View. Otherwise,
  // MS // the allocated vectors will be used to copy data  

  int NumIndices;
  
  int * RowIndices = GetRowIndices();
  int * ColIndices = GetColIndices();
  double * MatrixValues = GetValues();
  
  // MS // for CrsA and VbrA, handle the case of IndexBase different from 0
  // MS // diff is the difference between what we have and what we want 
  // MS // (1, FORTRAN style)
  int diff = 1-IndexBase();
  
  // MS // do all in terms of BlockRows. This way we can support VBR matrices

  for( int LocalBlockRow=0; LocalBlockRow<NumMyBlockRows() ; ++LocalBlockRow ) {

    bool FoundDiagonal = false;
  
    GetRow(LocalBlockRow,NumIndices,RowIndices,ColIndices,MatrixValues);

    for( int i=0 ; i<NumIndices ; ++i ) {
      if( abs(MatrixValues[i]) >= Threshold_ ) {
	if( RowIndices[i] == ColIndices[i] ) FoundDiagonal = true;
	(*Row)[NumMyMUMPSNonzeros_] = RowIndices[i]+diff;
	(*Col)[NumMyMUMPSNonzeros_] = ColIndices[i]+diff;
	(*Val)[NumMyMUMPSNonzeros_] = MatrixValues[i];
	NumMyMUMPSNonzeros_++;
      }
      // to add: zero on diagonal
    }
  }

  assert(NumMyMUMPSNonzeros_<=NumMyNonzeros());

  // MS // bring matrix to proc zero if required. Note that I first
  // MS // convert to COO format (locally), then I redistribute COO vectors.
  
  if( KeepMatrixDistributed_ == false ) {

    Epetra_IntSerialDenseVector * OldRow = Row;
    Epetra_IntSerialDenseVector * OldCol = Col;
    Epetra_SerialDenseVector * OldVal = Val;

    Comm().SumAll(&NumMyMUMPSNonzeros_,&NumMUMPSNonzeros_,1);

    if( Comm().MyPID() != 0 ) NumMUMPSNonzeros_ = 0;
      
    Row = new Epetra_IntSerialDenseVector(NumMUMPSNonzeros_);
    Col = new Epetra_IntSerialDenseVector(NumMUMPSNonzeros_);
    Val = new Epetra_SerialDenseVector(NumMUMPSNonzeros_);

    Epetra_Map OldNnz(-1,NumMyMUMPSNonzeros_,0,Comm());
    
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

    delete OldRow; OldRow = 0;
    delete OldCol; OldCol = 0;
    delete OldVal; OldVal = 0;

    if( Comm().MyPID() != 0 ) {
      delete Row; Row = 0;
      delete Col; Col = 0;
      delete Val; Val = 0;
    }

  } else {

    NumMUMPSNonzeros_ = NumMyMUMPSNonzeros_;

  }

  if( PrintTiming_ && Comm().MyPID() == 0 )
    cout << "Amesos_Mumps : Time to convert matrix to MUMPS format = "
	 << Time.ElapsedTime() << " (s)" << endl;
  
  AddToConvTime(Time.ElapsedTime());  
  
  return 0;

}

//=============================================================================

int Amesos_Mumps::ConvertToTripletValues()
{

  if( verbose_ == 3 ) cout << "Entering `ConvertToTripletValues'" << endl;

  Epetra_Time Time(Comm());

  // MS // convert to MUMPS format, keeping in distributed form.
  // MS // This doesn't require the matrix to be shipped to proc 0,
  // MS // then decomposed and reshipped by MUMPS

  // MS // be sure that those vectors are empty. Allocate them
  
  if( Val ) delete Val;

  Val = new Epetra_SerialDenseVector( NumMyNonzeros() ) ;  

  int NumMUMPSNonzerosValues = 0;
  
  // MS // retrive already allocated pointers for rows, cols, and vals
  // MS // Those arrays are allocated and freed by the Amesos_EpetraBaseSolver
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
      // following if : please test me................................
      //      if( MatrixValues[i] != 0.0 ) {
	(*Val)[NumMUMPSNonzerosValues] = MatrixValues[i];
	NumMUMPSNonzerosValues++;
	//      }
    }
  }

  if( NumMUMPSNonzerosValues != NumMyMUMPSNonzeros_ ) {
     cerr << "Amesos_Mumps ERROR : it appears that matrix has changed\n"
          << "Amesos_Mumps ERROR : its structure since SymbolicFactorization() has\n"
	  << "Amesos_Mumps ERROR : called. # nonzeros (sym) = " << NumMyMUMPSNonzeros_ << endl
	  << "Amesos_Mumps ERROR : called. # nonzeros (num) = " << NumMUMPSNonzerosValues << endl;
     return -2;
  }  
    
  // MS // bring matrix to proc zero if required
  
  if( KeepMatrixDistributed_ == false ) {

    Epetra_SerialDenseVector * OldVal = Val;

    Comm().SumAll(&NumMyMUMPSNonzeros_,&NumMUMPSNonzeros_,1);

    if( Comm().MyPID() != 0 ) NumMUMPSNonzeros_ = 0;
      
    Val = new Epetra_SerialDenseVector(NumMUMPSNonzeros_);

    Epetra_Map OldNnz(-1,NumMyMUMPSNonzeros_,0,Comm());
    
    Epetra_Map NewNnz(-1,NumMUMPSNonzeros_,0,Comm());
    Epetra_Import OldToNew(NewNnz,OldNnz);

    Epetra_Vector GOldVal(View,OldNnz,OldVal->Values());

    Epetra_Vector GVal(View,NewNnz,Val->Values());
    
    GVal.Import(GOldVal,OldToNew,Insert);

    delete OldVal; OldVal = 0;

    if( Comm().MyPID() != 0 ) {
      delete Val; Val = 0;
    }
  } else {

    NumMUMPSNonzeros_ = NumMyMUMPSNonzeros_;

  }

  if( PrintTiming_ && Comm().MyPID() == 0 )
    cout << "Amesos_Mumps : Time to convert matrix values to MUMPS format = "
	 << Time.ElapsedTime() << " (s)" << endl;
  
  AddToConvTime(Time.ElapsedTime());  
 
  return 0;

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

//=============================================================================

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
  if( IsComputeSchurComplementOK_ && Comm().MyPID() == 0 ) {
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

//=============================================================================

int Amesos_Mumps::ReadParameterList()
{

  if( ParameterList_ == NULL ) return 0;

  // retrive general parameters

  if( ParameterList_->getParameter("UseTranspose",false) )
    SetUseTranspose(true);

  Threshold_ = ParameterList_->getParameter("Threshold", 0.0);

  AddZeroToDiag_ = ParameterList_->getParameter("AddZeroToDiag", false);
  if( AddZeroToDiag_ && Comm().MyPID() == 0 ) {
    cerr << "Amesos_Mumps Warning: AddZeroToDiag not yet implemented." << endl;
  }
  
  PrintTiming_ = ParameterList_->getParameter("PrintTiming", true);

  PrintStatistics_ = ParameterList_->getParameter("PrintStatistics", true);

  // retrive MUMPS' parameters
  
  if (ParameterList_->isParameterSublist("mumps") ) {
    Teuchos::ParameterList MumpsParams = ParameterList_->sublist("mumps") ;
    // integer array of parameters
    if( MumpsParams.isParameter("ICNTL") ) {
      int * ICNTL = NULL;
      ICNTL = MumpsParams.getParameter("ICNTL", ICNTL);
      SetICNTL(ICNTL);
    }
    // double array of parameters
    if( MumpsParams.isParameter("CNTL") ) {
      double * CNTL = NULL;
      CNTL = MumpsParams.getParameter("CNTL", CNTL);
      SetCNTL(CNTL);
    }
    // ordering
     if( MumpsParams.isParameter("PermIn") ) {
      int * PermIn = NULL;
      PermIn = MumpsParams.getParameter("PermIn", PermIn);
      SetOrdering(PermIn);
    }
     // Maxis
     if( MumpsParams.isParameter("Maxis") ) {
       int Maxis = 0;
       Maxis = MumpsParams.getParameter("Maxis", Maxis);
       SetMaxis(Maxis);
     }
     // Maxs
     if( MumpsParams.isParameter("Maxs") ) {
       int Maxs = 0;
       Maxs = MumpsParams.getParameter("Maxs", Maxs);
       SetMaxs(Maxs);
     }
     // Col scaling
     if( MumpsParams.isParameter("ColScaling") ) {
       double * ColSca = NULL;
       ColSca = MumpsParams.getParameter("ColScaling", ColSca);
       SetColScaling(ColSca);
     }
     // Row scaling
     if( MumpsParams.isParameter("RowScaling") ) {
       double * RowSca = NULL;
       RowSca = MumpsParams.getParameter("RowScaling", RowSca);
       SetRowScaling(RowSca);
     }
     // that's all folks
  }  
  return 0;
}

//=============================================================================

int Amesos_Mumps::PerformSymbolicFactorization()
{

  if( verbose_ == 3 ) cout << "Entering `PerformSymbolicFactorization'" << endl;

  // erase data if present
  if( SymbolicFactorizationOK_ && MDS.job != -777 ) {
   Destroy();
  }

  if( Comm().NumProc() == 1 ) {
      KeepMatrixDistributed_ = false;
  }
  
  Epetra_Time Time(Comm());

  if( IsLocal() || UseMpiCommSelf_ ) {
#ifdef EPETRA_MPI
#ifndef TFLOP
    MPI_Comm MPIC = MPI_COMM_SELF ;
    MDS.comm_fortran = (F_INT) MPI_Comm_c2f( MPIC ) ;   // Compiled on cygwin but not on Atlantis
#else
    MDS.comm_fortran = -987654 ;  // Needed by MUMPS 4.3 
#endif
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
  //  infog(1) is zereo
  MDS.n = NumGlobalRows() ;
  
  // fix pointers for nonzero pattern of A. Numerical values
  // will be entered in PerformNumericalFactorization()
  if( KeepMatrixDistributed_ ) {
    MDS.nz_loc = NumMUMPSNonzeros_;
    MDS.irn_loc = Row->Values(); 
    MDS.jcn_loc = Col->Values(); 
  } else {
    if( Comm().MyPID() == 0 ) {
      MDS.nz = NumMUMPSNonzeros_;
      MDS.irn = Row->Values(); 
      MDS.jcn = Col->Values(); 
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
                     // a limitation, though
  
  dmumps_c( &MDS ) ;  // Perform symbolic factorization
  // check the error flag. MUMPS is successful only if
  // infog(1) is zereo
  if( MDS.INFOG(1)!=0 && verbose_ ) 
  EPETRA_CHK_ERR( MDS.INFOG(1) );
  if( MDS.INFOG(1) ) return( MDS.INFOG(1) );

  if( PrintTiming_ && Comm().MyPID() == 0 )
    cout << "Amesos_Mumps : Time for symbolic factorization = "
	 << Time.ElapsedTime() << " (s)" << endl;
  
  AddToSymFactTime(Time.ElapsedTime());  

  SymbolicFactorizationOK_ = true ;

  return 0;

}

//=============================================================================

int Amesos_Mumps::PerformNumericFactorization( )
{
  
  Epetra_Time Time(Comm());
  
  if( verbose_ == 3 ) cout << "Entering `PerformNumericFactorization'" << endl;

  // set vector for matrix entries. User may have changed it
  // since PerformSymbolicFactorization() has been called
  if( KeepMatrixDistributed_ ) MDS.a_loc = Val->Values();
  else {
    if( Comm().MyPID() == 0 ) {
      MDS.a = Val->Values();
    }
  }
  
  MDS.job = 2  ;     // Request numeric factorization
  dmumps_c( &MDS ) ;  // Perform numeric factorization
  // check the error flag. MUMPS is siccessful only if
  // infog(1) is zereo
  EPETRA_CHK_ERR( MDS.INFOG(1) );
  if( MDS.INFOG(1) ) return( MDS.INFOG(1) );

  if( PrintTiming_ && Comm().MyPID() == 0 )
    cout << "Amesos_Mumps : Time for numeric factorization = "
	 << Time.ElapsedTime() << " (s)" << endl;
  
  AddToNumFactTime(Time.ElapsedTime());
  
  NumericFactorizationOK_ = true;

  return 0;
}

//=============================================================================

int Amesos_Mumps::SymbolicFactorization()
{

  if( verbose_ == 3 ) cout << "Entering `SymbolicFactorization'" << endl;
  
  ReadParameterList();  
  
  // first, gather matrix property using base class
  // Amesos_EpetraBaseSolver

  Epetra_RowMatrix * Matrix = GetMatrix();
  assert( Matrix != NULL );
  if( Matrix != OldMatrix_ ) {
    // initialize redistor, target map is on proc 0.
    // Redistor will be used to bring rhs to proc 0, and solution to
    // all procs
    EPETRA_CHK_ERR( SetInterface(Matrix) );
    OldMatrix_ = Matrix;
    SymbolicFactorizationOK_ = false;
    NumericFactorizationOK_ = false;
  }

  // now, create a map, in which all nodes are stored on processor 0.
  // Also creates import objects. This is done using the base class
  // Amesos_EpetraRedistributor. NOTE: the following classes return
  // if they have already been called.  

  EPETRA_CHK_ERR(ConvertToTriplet());
  EPETRA_CHK_ERR(PerformSymbolicFactorization());

  if( verbose_ == 3 ) cout << "Exiting `SymbolicFactorization'" << endl;

  return 0;
}

//=============================================================================

int Amesos_Mumps::NumericFactorization()
{

  if( verbose_ == 3 ) cout << "Entering `NumericFactorization'" << endl;

  ReadParameterList();

  // MS // As in SymbolicFactorization. All those functions
  // MS // returns if they have already been called.
  
  Epetra_RowMatrix * Matrix = GetMatrix();
  assert( Matrix != NULL );  
  if( Matrix != OldMatrix_ ) {
    EPETRA_CHK_ERR( SetInterface(Matrix) );
    OldMatrix_ = Matrix;
    SymbolicFactorizationOK_ = false;
    NumericFactorizationOK_ = false;
  }

  // MS // reship the matrix. User can update numerical values
  // MS // Should the values remain the same, it is better to use
  // MS // Solve() and NOT SymbolicFactorization() + NumericFactorization()
  // MS // + Solve()
  
  if( SymbolicFactorizationOK_ == false ) {
    EPETRA_CHK_ERR( ConvertToTriplet() );
    EPETRA_CHK_ERR( PerformSymbolicFactorization() );
  } else {
    EPETRA_CHK_ERR( ConvertToTripletValues() );
  }  
  
  EPETRA_CHK_ERR( PerformNumericFactorization() );
  
  if( verbose_ == 3 ) cout << "Exiting `NumericFactorization'" << endl;
  
  return 0;
}
 
//=============================================================================

int Amesos_Mumps::Solve()
{ 

  if( verbose_ == 3 ) cout << "Entering `Solve'" << endl;

  ReadParameterList();

  double RedistorTime = 0.0;
  Epetra_Time TimeForRedistor(Comm());
  
  Epetra_RowMatrix * Matrix = GetMatrix();
  assert( Matrix != NULL );
  if( Matrix != OldMatrix_ ) {
    EPETRA_CHK_ERR( SetInterface(Matrix) );
    OldMatrix_ = Matrix;
    SymbolicFactorizationOK_ = false;
    NumericFactorizationOK_ = false;
  }
  
  if( SymbolicFactorizationOK_ == false && NumericFactorizationOK_ == false ) {

    // MS // it is the first time the user call something. The
    // MS // matrix is neither sym fact not num fact
    // MS // We convert it to COO format, then sym+fact
    EPETRA_CHK_ERR( ConvertToTriplet() );
    EPETRA_CHK_ERR( PerformSymbolicFactorization() );
    EPETRA_CHK_ERR( PerformNumericFactorization() );

  }  else if( NumericFactorizationOK_  == false )  {

    // MS // User has already called SymFact, but maybe he changed
    // MS // the matrix entries. Ww reship the matrix
    EPETRA_CHK_ERR( ConvertToTripletValues() );
    EPETRA_CHK_ERR( PerformNumericFactorization() );

  } else {

    // MS // just for comment it out: This case, the user has already
    // MS // numerically factorized the matrix. We can solve without
    // MS // shipping the matrix
    
  }
  
  Epetra_Time Time(Comm());

  Epetra_MultiVector * vecX = GetLHS() ; 
  Epetra_MultiVector * vecB = GetRHS() ;

  // create redistor the first time enters here. I suppose that LHS and RHS
  // (all of them) have the same map.
  if( Redistor_ == 0 ) {
    Redistor_ = new EpetraExt_Redistor(vecX->Map(),1);
    assert( Redistor_ != 0 );
  }
  
  //  Compute the number of right hands sides (and check that X and B have the same shape) 

  if( vecB == NULL || vecX == NULL ) {
    if( Comm().MyPID() == 0 ) cout << "Amesos_Mumps ERROR : LHS or RHS not set" << endl;
    return -1;
  }

  int nrhs = vecX->NumVectors() ;
  if( nrhs != vecB->NumVectors() ) {
    if( Comm().MyPID() == 0 ) cout << "Amesos_Mumps ERROR : multivectors LHS or RHS have different sizes" << endl;
    return -2;
  }

  // will be used only if IsLocal == false
  if( IsLocal() == false && UseMpiCommSelf_ == false && TargetVector_ == 0 ) {
    TargetVector_ = new Epetra_MultiVector(*(Redistor_->TargetMap()),nrhs);
    assert( TargetVector_ != 0 );
  }
  
  // MS // Extract Target versions of X and B. 
  // MS // NOTE: both Target and distributed version of MUMPS require
  // MS // rhs to be entirely stored on processor 0. This is because we
  // MS // still need TargetMap_ also for distributed input

  if( IsLocal() || UseMpiCommSelf_ ) {
    if( Comm().MyPID() == 0 ) {
      for ( int j =0 ; j < nrhs; j++ ) {
	for( int i=0 ; i<NumMyRows() ; ++i ) (*vecX)[j][i] = (*vecB)[j][i];
	MDS.rhs = (*vecX)[j];
	MDS.job = 3  ;     // Request solve
	dmumps_c( &MDS ) ;  // Perform solve
	// check the error flag. MUMPS is siccessful only if
	// infog(1) is zereo
        EPETRA_CHK_ERR( MDS.INFOG(1) );
  	if( MDS.INFOG(1) ) return( MDS.INFOG(1) );

      }
    }
  } else {

    TimeForRedistor.ResetStartTime();
    Redistor_->TargetImport( *vecB, *TargetVector_ ) ;
    RedistorTime += TimeForRedistor.ElapsedTime();
    
    for ( int j =0 ; j < nrhs; j++ ) { 
      if ( Comm().MyPID() == 0 ) {
	MDS.rhs = (*TargetVector_)[j];
      }
      MDS.job = 3  ;     // Request solve
      dmumps_c( &MDS ) ;  // Perform solve
      // check the error flag. MUMPS is siccessful only if
      // infog(1) is zereo
      EPETRA_CHK_ERR( MDS.INFOG(1) );
      if( MDS.INFOG(1) ) return( MDS.INFOG(1) );
    }

    TimeForRedistor.ResetStartTime();
    Redistor_->SourceImport( *vecX, *TargetVector_ ) ;
    RedistorTime += TimeForRedistor.ElapsedTime();

  }
  
  EPETRA_CHK_ERR( UpdateLHS() );

  if( PrintTiming_ && Comm().MyPID() == 0 ) {
     cout << "Amesos_Mumps : Time for gather B and scatter X = "
	 << RedistorTime << " (s)" << endl;
    cout << "Amesos_Mumps : Time for solve = "
	 << Time.ElapsedTime() << " (s)" << endl;
  }
  
  AddToSolTime(Time.ElapsedTime());

  // compute vector norms
  if( ParameterList_->getParameter("ComputeVectorNorms",false) ) {
    double NormLHS, NormRHS;
    for( int i=0 ; i<nrhs ; ++i ) {
      assert((*vecX)(i)->Norm2(&NormLHS)==0);
      assert((*vecB)(i)->Norm2(&NormRHS)==0);
      if( Comm().MyPID() == 0 ) {
	cout << "Amesos_Mumps : vector " << i << ", ||x|| = " << NormLHS
	     << ", ||b|| = " << NormRHS << endl;
      }
    }
  }

  // compute true residual
  if( ParameterList_->getParameter("ComputeTrueResidual",false) ) {
    double Norm;
    Epetra_Vector Ax(vecB->Map());
    for( int i=0 ; i<nrhs ; ++i ) {
      assert(GetMatrix()->Multiply(UseTranspose(), *((*vecX)(i)), Ax)==0);
      assert(Ax.Update(1.0, *((*vecB)(i)), -1.0)==0);
      assert(Ax.Norm2(&Norm)==0);
      if( Comm().MyPID() == 0 ) {
	cout << "Amesos_Mumps : vector " << i << ", ||Ax - b|| = " << Norm << endl;
      }
    }
  }

  //  ParameterList_->setParameter("solution time",GetSolTime());
			      
  if( ParameterList_->getParameter("PrintStatistics",true) ) PrintInformation();
  
  return(0) ; 
}

// ================================================ ====== ==== ==== == =

Epetra_CrsMatrix * Amesos_Mumps::GetCrsSchurComplement() 
{

  if( IsComputeSchurComplementOK_ ) {
    
    if( Comm().MyPID() != 0 ) NumSchurComplementRows_ = 0;
    
    Epetra_Map SCMap(-1,NumSchurComplementRows_, 0, Comm());

    CrsSchurComplement_ = new Epetra_CrsMatrix(Copy,SCMap,NumSchurComplementRows_);

    if( Comm().MyPID() == 0 )
      for( int i=0 ; i<NumSchurComplementRows_ ; ++i ) {
	for( int j=0 ; j<NumSchurComplementRows_ ; ++j ) {
	  int pos = i+ j *NumSchurComplementRows_;
	  CrsSchurComplement_->InsertGlobalValues(i,1,&(MDS.schur[pos]),&j);
	}
      }
    
    CrsSchurComplement_->FillComplete();

    return CrsSchurComplement_;
  }

  return 0;
  
}

// ================================================ ====== ==== ==== == =

Epetra_SerialDenseMatrix * Amesos_Mumps::GetDenseSchurComplement() 
{
  
  if( IsComputeSchurComplementOK_ ) {
    
    if( Comm().MyPID() != 0 ) return 0;
    
    DenseSchurComplement_ = new Epetra_SerialDenseMatrix(NumSchurComplementRows_,
							NumSchurComplementRows_);
    
    for( int i=0 ; i<NumSchurComplementRows_ ; ++i ) {
      for( int j=0 ; j<NumSchurComplementRows_ ; ++j ) {
	int pos = i+ j *NumSchurComplementRows_;
	(*DenseSchurComplement_)(i,j) = MDS.schur[pos];
      }
    }
    
    return DenseSchurComplement_;
    
  }
  
  return 0;
  
}

// ================================================ ====== ==== ==== == =

int Amesos_Mumps::ComputeSchurComplement(bool flag, int NumSchurComplementRows,
					 int * SchurComplementRows)
{
  NumSchurComplementRows_ = NumSchurComplementRows;
  
  SchurComplementRows_ = SchurComplementRows;

  // modify because MUMPS is fortran-driven
  if( Comm().MyPID() == 0 )
    for( int i=0 ; i<NumSchurComplementRows ; ++i ) SchurComplementRows_[i]++;
  
  IsComputeSchurComplementOK_ = flag;

  return 0;
}

// ================================================ ====== ==== ==== == =

int Amesos_Mumps::PrintInformation() 
{
  if( Comm().MyPID() ) return 0;

  cout << "------------ : ---------------------------------------------------------" << endl;
  cout << "Amesos_Mumps : global information for all phases" << endl;
  cout << "------------ : ---------------------------------------------------------" << endl;
  cout << "Amesos_Mumps : Matrix has " << NumGlobalRows() << " rows"
       << " and " << NumGlobalNonzeros() << " nonzeros" << endl;
  cout << "Amesos_Mumps : Time to convert matrix into MUMPS format = "
       << GetConvTime() << " (s) " << endl;
  cout << "Amesos_Mumps : Symbolic factorization time = "
       << GetSymFactTime() << " (s) " << endl;
  cout << "Amesos_Mumps : Numerical factorization time = "
       << GetSymFactTime() << " (s) " << endl;
  cout << "Amesos_Mumps : Solution time = "
       << GetSymFactTime() << " (s) " << endl;
  cout << "Amesos_Mumps : Estimated FLOPS for elimination = "
       << MDS.RINFOG(1) << endl;
  cout << "Amesos_Mumps : Total FLOPS for assembly = "
       << MDS.RINFOG(2) << endl;
  cout << "Amesos_Mumps : Total FLOPS for elimination = "
       << MDS.RINFOG(3) << endl;
  
  cout << "Amesos_Mumps : Total real space to store the LU factors = "
       << MDS.INFOG(9) << endl;
  cout << "Amesos_Mumps : Total integer space to store the LU factors = "
       << MDS.INFOG(10) << endl;
  cout << "Amesos_Mumps : Total number of iterative steps refinement = "
       << MDS.INFOG(15) << endl;
  cout << "Amesos_Mumps : Estimated size of MUMPS internal data\n"
       << "Amesos_Mumps : for running factorization = "
       << MDS.INFOG(16) << " Mbytes" << endl;
  cout << "Amesos_Mumps : Estimated size of MUMPS internal data\n"
       << "Amesos_Mumps : for running factorization = "
       << MDS.INFOG(17) << " Mbytes" << endl;
  cout << "Amesos_Mumps : Estimated size of MUMPS internal data\n"
       << "Amesos_Mumps : Allocated during factorization = "
       << MDS.INFOG(19) << " Mbytes" << endl;
 cout << "------------ : ---------------------------------------------------------" << endl;
 
  return 0;
  
}

int Amesos_Mumps::SetICNTL(int * icntl)
{
  for( int i=0 ; i<40 ; ++i ) icntl_[i] = icntl[i];
  return 0;
}

int Amesos_Mumps::SetCNTL(double * cntl)  
{
  for( int i=0 ; i<5 ; ++i ) cntl_[i] = cntl[i];
  return 0;
}
 
