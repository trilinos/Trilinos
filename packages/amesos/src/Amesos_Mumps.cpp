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
#include "Epetra_Util.h"

#define ICNTL(I) icntl[(I)-1]
#define CNTL(I)  cntl[(I)-1]
#define INFOG(I) infog[(I)-1]
#define INFO(I) info[(I)-1]
#define RINFOG(I) rinfog[(I)-1]

#ifndef HAVE_AMESOS_SMUMPS
#define MUMPS_INTERFACE dmumps_c
#else
#define MUMPS_INTERFACE smumps_c
#endif

const int DEF_VALUE_INT = -123456789;
const double DEF_VALUE_DOUBLE = -123456.789;
  
//=============================================================================

Amesos_Mumps::Amesos_Mumps(const Epetra_LinearProblem &prob ) :
  Problem_(&prob),
  SymbolicFactorizationOK_(false), 
  NumericFactorizationOK_(false),
  MaxProcs_(-1),
  RedistrMap_(0),
  RedistrMatrix_(0),
  RedistrImporter_(0),
  SerialMap_(0),
  SerialImporter_(0),
  NumSchurComplementRows_(-1),
  SchurComplementRows_(0),
  CrsSchurComplement_(0),
  DenseSchurComplement_(0),
  IsComputeSchurComplementOK_(false),
  ComputeVectorNorms_(false),
  ComputeTrueResidual_(false),
  MatrixProperty_(0),
  RowSca_(0),
  ColSca_(0),
  PermIn_(0),
  Maxis_(DEF_VALUE_INT),
  Maxs_(DEF_VALUE_INT),
  verbose_(1),
  AddToDiag_(0.0),
  AddDiagElement_(false),
  PrintTiming_(false),
  PrintStatus_(false),
  Threshold_(0.0),
  MUMPSComm_(0),
  UseTranspose_(false),
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
  MatrixType_(0)

{
  // -777 is for me. It means : never called MUMPS, so
  // SymbolicFactorization will not call Destroy();
  MDS.job = -777;
  
  // set to -1 icntl_ and cntl_. The use can override default values by using
  // SetICNTL(pos,value) and SetCNTL(pos,value).
  for( int i=0 ; i<40 ; ++i ) icntl_[i] = DEF_VALUE_INT;
  for( int i=0 ; i<5 ; ++i ) cntl_[i] = DEF_VALUE_DOUBLE;

  Teuchos::ParameterList ParamList;
  SetParameters( ParamList );
  
}

//=============================================================================

void Amesos_Mumps::Destroy()
{
  
  // destroy instance of the package
  MDS.job = -2;
  
  if (Comm().MyPID() < MaxProcs_) 
    MUMPS_INTERFACE(&MDS);
    
  if (RedistrMap_) {
    delete RedistrMap_;
    RedistrMap_ = 0;
  }

  if (RedistrMatrix_) {
    delete RedistrMatrix_;
    RedistrMatrix_ = 0;
  }

  if (RedistrImporter_) {
    delete RedistrImporter_;
    RedistrImporter_ = 0;
  }

  if (SerialMap_) {
    delete SerialMap_;
    SerialMap_ = 0;
  }

  if (SerialImporter_) {
    delete SerialImporter_;
    SerialImporter_ = 0;
  }

  if (IsComputeSchurComplementOK_ && (Comm().MyPID() == 0)
      && MDS.schur) {
    delete [] MDS.schur;
    MDS.schur = 0;
  }

  if (MUMPSComm_) {
    MPI_Comm_free( &MUMPSComm_ );
    MUMPSComm_ = 0;
  }

  if( (verbose_ && PrintTiming_) || verbose_ == 2) 
    PrintTiming();
  if( (verbose_ && PrintStatus_) || verbose_ == 2) 
    PrintStatus();

  if (Time_) { 
    delete Time_; 
    Time_ = 0; 
  }
     
  return;
}

//=============================================================================

Amesos_Mumps::~Amesos_Mumps(void)
{

  Destroy();
 
}

//=============================================================================

int Amesos_Mumps::ConvertToTriplet(const bool OnlyValues)
{

  Epetra_RowMatrix* ptr;
  if (Comm().NumProc() == 1)
    ptr = &Matrix();
  else {
    ptr = &RedistrMatrix(true);
  }

  Time_->ResetStartTime();


  Row.resize(ptr->NumMyNonzeros());
  Col.resize(ptr->NumMyNonzeros());
  Val.resize(ptr->NumMyNonzeros());

  int MaxNumEntries = ptr->MaxNumEntries();
  vector<int> Indices;
  vector<double> Values;
  Indices.resize(MaxNumEntries);
  Values.resize(MaxNumEntries);

  int count = 0;

  for (int i = 0; i < ptr->NumMyRows() ; ++i) {

    int GlobalRow = ptr->RowMatrixRowMap().GID(i);

    int NumEntries = 0;
    int ierr;
    ierr = ptr->ExtractMyRowCopy(i, MaxNumEntries,
				   NumEntries, &Values[0],
				   &Indices[0]);
    AMESOS_CHK_ERR(ierr);

    for (int j = 0 ; j < NumEntries ; ++j) {
      if (OnlyValues == false) {
	Row[count] = GlobalRow + 1;
	Col[count] = ptr->RowMatrixColMap().GID(Indices[j]) + 1;
      }
      Val[count] = Values[j];
      count++;
    }
  }

  ConTime_ == Time_->ElapsedTime(); 
  
  assert (count <= ptr->NumMyNonzeros());

  return(0);

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
  MDS.ICNTL(4)  = -1;
  MDS.ICNTL(5)  = 0;   // Matrix is given in elemental (i.e. triplet) from
  MDS.ICNTL(6)  = 7;   // Choose column permutation automatically
  MDS.ICNTL(7)  = 7;   // Choose ordering method automatically
  MDS.ICNTL(8)  = 7;   // Choose scaling automatically
  MDS.ICNTL(9)  = 1;   // Compute A x = b 
  MDS.ICNTL(10) = 0;   // Maximum steps of iterative refinement
  MDS.ICNTL(11) = 0;   // Do not collect statistics
  MDS.ICNTL(12) = 0;   // Use Node level parallelism
  MDS.ICNTL(13) = 0;   // Use ScaLAPACK for root node 
  MDS.ICNTL(14) = 20;  // Increase memory allocation 20% at a time 
  MDS.ICNTL(15) = 0;   // Minimize memory use (not flop count)
  MDS.ICNTL(16) = 0;   // Do not perform null space detection
  MDS.ICNTL(17) = 0;   // Unused (null space dimension)
  
  if (Comm().NumProc() != 1) 
    MDS.ICNTL(18)= 3;
  else
    MDS.ICNTL(18)= 0;
  
  if ( UseTranspose() )  MDS.ICNTL(9) = 0 ; 
  else                   MDS.ICNTL(9) = 1 ;
  
  if( IsComputeSchurComplementOK_ ) MDS.ICNTL(19) = 1;
  else                              MDS.ICNTL(19) = 0;

  // something to do if the Schur complement is required.
  if( IsComputeSchurComplementOK_ && Comm().MyPID() == 0 ) {
    MDS.size_schur = NumSchurComplementRows_;
    MDS.listvar_schur = SchurComplementRows_;
    MDS.schur = new AMESOS_TYPE[NumSchurComplementRows_*NumSchurComplementRows_];
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
  assert(MDS.ICNTL(5)== 0); 
  
  return;
  
}

//=============================================================================

int Amesos_Mumps::SetParameters( Teuchos::ParameterList & ParameterList)
{

  // ========================================= //
  // retrive MUMPS' parameters from list.      //
  // default values defined in the constructor //
  // ========================================= //
  
  // retrive general parameters

  // solve problem with transpose
  if( ParameterList.isParameter("UseTranspose") )
    SetUseTranspose(ParameterList.get("UseTranspose",false));
  
  // ignore all elements below given tolerance
  if( ParameterList.isParameter("Threshold") )
    Threshold_ = ParameterList.get("Threshold", 0.0);

  // add zero to diagonal if diagonal element is not present
  if( ParameterList.isParameter("AddZeroToDiag") )
    AddDiagElement_ = ParameterList.get("AddZeroToDiag", false);

  // add this value to diagonal
  if( ParameterList.isParameter("AddToDiag") )
    AddToDiag_ = ParameterList.get("AddToDiag", 0.0);

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

  // define on how many processes matrix should be converted into MUMPS
  // format. (this value must be less than available procs)
  if( ParameterList.isParameter("MaxProcs") )
    MaxProcs_ = ParameterList.get("MaxProcs",-1);

  // Matrix property, defined internally in Amesos_Mumps as an integer,
  // whose value can be:
  // - 0 : general unsymmetric matrix;
  // - 1 : SPD;
  // - 2 : general symmetric matrix.
  if( ParameterList.isParameter("MatrixType") ) {
    string MatrixType;
    MatrixType = ParameterList.get("MatrixType",MatrixType);
    if( MatrixType == "SPD" )
      MatrixType_ = 1;
    else if( MatrixType == "symmetric" ) 
      MatrixType_ = 2;
    else if( MatrixType == "general" )   
      MatrixType_ = 0;
    else {
      if( Comm().MyPID() == 0 ) {
	cerr << "Amesos_Mumps : ERROR" << endl 
	     << "Amesos_Mumps : MatrixType value not recognized ("
	     << MatrixType << ")" << endl;
      }
    }
  }

  // some verbose output:
  // 0 - no output at all
  // 1 - output as specified by other parameters
  // 2 - all possible output
  if( ParameterList.isParameter("OutputLevel") )
    verbose_ = ParameterList.get("OutputLevel",1);

  // retrive MUMPS' specific parameters
  
  if (ParameterList.isSublist("mumps") ) {
    Teuchos::ParameterList MumpsParams = ParameterList.sublist("mumps") ;
    // integer array of parameters
    if( MumpsParams.isParameter("ICNTL") ) {
      int * ICNTL = 0;
      ICNTL = MumpsParams.get("ICNTL", ICNTL);
      if( ICNTL ) SetICNTL(ICNTL);
    }
    // double array of parameters
    if( MumpsParams.isParameter("CNTL") ) {
      double * CNTL = 0;
      CNTL = MumpsParams.get("CNTL", CNTL);
      if( CNTL ) SetCNTL(CNTL);
    }
    // ordering
     if( MumpsParams.isParameter("PermIn") ) {
      int * PermIn = 0;
      PermIn = MumpsParams.get("PermIn", PermIn);
      if( PermIn ) SetOrdering(PermIn);
    }
     // Maxis
     if( MumpsParams.isParameter("Maxis") ) {
       int Maxis = 0;
       Maxis = MumpsParams.get("Maxis", Maxis);
       SetMaxis(Maxis);
     }
     // Maxs
     if( MumpsParams.isParameter("Maxs") ) {
       int Maxs = 0;
       Maxs = MumpsParams.get("Maxs", Maxs);
       SetMaxs(Maxs);
     }
     // Col scaling
     if( MumpsParams.isParameter("ColScaling") ) {
       AMESOS_TYPE * ColSca = 0;
       ColSca = MumpsParams.get("ColScaling", ColSca);
       if( ColSca ) SetColScaling(ColSca);
     }
     // Row scaling
     if( MumpsParams.isParameter("RowScaling") ) {
       AMESOS_TYPE * RowSca = 0;
       RowSca = MumpsParams.get("RowScaling", RowSca);
       if( RowSca ) SetRowScaling(RowSca);
     }
     // that's all folks
  }  

  return 0;
}

//=============================================================================

void Amesos_Mumps::CheckParameters() 
{

#ifndef HAVE_AMESOS_MPI_C2F
  MaxProcs_ = -3;
#endif

  // check parameters and fix values of MaxProcs_

  int NumGlobalNonzeros, NumRows;
  
  NumGlobalNonzeros = Matrix().NumGlobalNonzeros(); 
  NumRows = Matrix().NumGlobalRows(); 

  // optimal value for MaxProcs == -1
  
  int OptNumProcs1 = 1 + EPETRA_MAX(NumRows/10000, NumGlobalNonzeros/1000000);
  OptNumProcs1 = EPETRA_MIN(Comm().NumProc(),OptNumProcs1);

  // optimal value for MaxProcs == -2

  int OptNumProcs2 = (int)sqrt(1.0 *  Comm().NumProc());
  if( OptNumProcs2 < 1 ) OptNumProcs2 = 1;

  // fix the value of MaxProcs

  switch (MaxProcs_) {
  case -1:
    MaxProcs_ = OptNumProcs1;
    break;
  case -2:
    MaxProcs_ = OptNumProcs2;
    break;
  case -3:
    MaxProcs_ = Comm().NumProc();
    break;
  }

  // few checks
  if (MaxProcs_ > Comm().NumProc()) 
    MaxProcs_ = Comm().NumProc();

  return;
  
}

//=============================================================================

int Amesos_Mumps::SymbolicFactorization()
{

  // erase data if present. 
  if (SymbolicFactorizationOK_ && MDS.job != -777) {
   Destroy();
  }

  SymbolicFactorizationOK_ = false;
  NumericFactorizationOK_ = false;

  if (Time_ == 0) 
    Time_ = new Epetra_Time(Comm());

  CheckParameters();
  AMESOS_CHK_ERR(ConvertToTriplet(false));

#if defined(HAVE_MPI) && defined(HAVE_AMESOS_MPI_C2F)
  if (MaxProcs_ != Comm().NumProc()) {

    if(MUMPSComm_) 
      MPI_Comm_free(&MUMPSComm_);

    int * ProcsInGroup = new int[MaxProcs_];
    for (int i = 0 ; i < MaxProcs_ ; ++i) 
      ProcsInGroup[i] = i;

    MPI_Group OrigGroup, MumpsGroup;
    MPI_Comm_group(MPI_COMM_WORLD, &OrigGroup);
    MPI_Group_incl(OrigGroup, MaxProcs_, ProcsInGroup, &MumpsGroup);
    MPI_Comm_create(MPI_COMM_WORLD, MumpsGroup, &MUMPSComm_);
    MDS.comm_fortran = (F_INT) MPI_Comm_c2f( MUMPSComm_);

    delete [] ProcsInGroup;

  } else {
    const Epetra_MpiComm* MpiComm = dynamic_cast<const Epetra_MpiComm*>(&Comm());
    assert (MpiComm != 0);
    MDS.comm_fortran = (F_INT) MPI_Comm_c2f(MpiComm->GetMpiComm());
  }
#else
  // only thing I can do, use MPI_COMM_WORLD.
  MDS.comm_fortran = -987654;
#endif
  
  MDS.job = -1  ;     //  Initialization
  MDS.par = 1 ;       //  Host IS involved in computations
  MDS.sym = MatrixProperty_;

  RedistrMatrix(true);

  if (Comm().MyPID() < MaxProcs_) {
    MUMPS_INTERFACE(&MDS);   //  Initialize MUMPS

    CheckError();
  }

  MDS.n = Matrix().NumGlobalRows() ;

  // fix pointers for nonzero pattern of A. Numerical values
  // will be entered in PerformNumericalFactorization()
  if (Comm().NumProc() != 1) {

    MDS.nz_loc = RedistrMatrix().NumMyNonzeros();

    if (Comm().MyPID() < MaxProcs_) {
      MDS.irn_loc = &Row[0]; 
      MDS.jcn_loc = &Col[0];

    }
  } else {
    if (Comm().MyPID() == 0) {
      MDS.nz = Matrix().NumMyNonzeros();
      MDS.irn = &Row[0]; 
      MDS.jcn = &Col[0]; 
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

  SetICNTLandCNTL();

  // Perform symbolic factorization

  Time_->ResetStartTime();
  if (Comm().MyPID() < MaxProcs_) 
    MUMPS_INTERFACE(&MDS);

  SymTime_ += Time_->ElapsedTime();

  CheckError();

  SymbolicFactorizationOK_ = true ;
  NumSymbolicFact_++;  
  return 0;

}

//=============================================================================

int Amesos_Mumps::NumericFactorization()
{

  NumericFactorizationOK_ = false;
  
  if (SymbolicFactorizationOK_ == false)
    AMESOS_CHK_ERR(SymbolicFactorization());

  RedistrMatrix(true);
  AMESOS_CHK_ERR(ConvertToTriplet(true));

#ifndef HAVE_AMESOS_SMUMPS
  if (Comm().NumProc() != 1) {
    if (Comm().MyPID() < MaxProcs_) 
      MDS.a_loc = &Val[0];
  } else {
    MDS.a = &Val[0];
  }
#else
  SVal.resize(Val.size());

  if (Comm().MyPID() < MaxProcs_) 
    for (int i = 0 ; i < Val.size() ; ++i) 
      SVal[i] = (float)(Val[i]);

  if (Comm().NumProc() != 1) {
    if (Comm().MyPID() < MaxProcs_) 
      MDS.a_loc = &SVal[0];
  } else {
    MDS.a = &SVal[0];
  }
#endif

  // Request numeric factorization 
  MDS.job = 2;
  // Perform numeric factorization
  Time_->ResetStartTime();
  if (Comm().MyPID() < MaxProcs_) {
    MUMPS_INTERFACE(&MDS);
  }

  NumTime_ += Time_->ElapsedTime();
  
  CheckError();

  NumericFactorizationOK_ = true;
  NumNumericFact_++;  
  return 0;
}

//=============================================================================

int Amesos_Mumps::Solve()
{ 

  if (NumericFactorizationOK_ == false)
    AMESOS_CHK_ERR(NumericFactorization());

  NumSolve_++;  
  
  if( Time_ == 0 ) Time_ = new Epetra_Time( Comm() );

  Epetra_MultiVector* vecX = Problem_->GetLHS() ; 
  Epetra_MultiVector* vecB = Problem_->GetRHS() ;
  int NumVectors = vecX->NumVectors();

  if ((vecX == 0) || (vecB == 0))
    AMESOS_CHK_ERR(-1);

  if (NumVectors != vecB->NumVectors())
    AMESOS_CHK_ERR(-1);

#ifdef HAVE_AMESOS_SMUMPS
  // need to allocate space for single-precision solution/rhs
  if (Comm().MyPID() == 0) {
    SVector.resize(Matrix().NumGlobalRows());
  }
#endif
  
  if (Comm().NumProc() == 1) {

    for (int j = 0 ; j < NumVectors; j++) {

      Time_->ResetStartTime();
      MDS.job = 3;     // Request solve

#ifndef HAVE_AMESOS_SMUMPS      
      for (int i = 0 ; i < Matrix().NumMyRows() ; ++i) 
	(*vecX)[j][i] = (*vecB)[j][i];
      MDS.rhs = (*vecX)[j];
      MUMPS_INTERFACE(&MDS) ;  // Perform solve
#else
      for (int i = 0 ; i < Matrix().NumMyRows() ; ++i) 
	SVector[i] = (float)(*vecB)[j][i];
      MDS.rhs = &SVector[0];
      MUMPS_INTERFACE(&MDS) ;  // Perform solve
      for (int i = 0 ; i < Matrix().NumMyRows() ; ++i) 
	(*vecX)[j][i] = (double)SVector[i];
#endif
      SolTime_ += Time_->ElapsedTime();
      
      CheckError();
    }

  } else {

    Epetra_MultiVector SerialVector(SerialMap(),NumVectors);

    Time_->ResetStartTime();
    AMESOS_CHK_ERR(SerialVector.Import(*vecB,SerialImporter(),Insert));
    VecTime_ += Time_->ElapsedTime();
    
    for (int j = 0 ; j < NumVectors; j++) { 

      if (Comm().MyPID() == 0) {
#ifndef HAVE_AMESOS_SMUMPS
	MDS.rhs = SerialVector[j];
#else
	// copy the double-precision into the single-precision SVector
	for (int i = 0 ; i < Matrix().NumGlobalRows() ; ++i) 
	  SVector[i] = (float)(SerialVector[j][i]);
	MDS.rhs = &SVector[0];
#endif
      }
      // solve the linear system and take time
      MDS.job = 3;     
      Time_->ResetStartTime();
      if (Comm().MyPID() < MaxProcs_) 
	MUMPS_INTERFACE(&MDS) ;  // Perform solve
#ifdef HAVE_AMESOS_SMUMPS  
      // copy back MUMPS' solution into TargetVector, that will be later
      // redistributed using Epetra utilities
      if (Comm().MyPID() == 0)
	for (int i = 0 ; i < Matrix().NumGlobalRows() ; ++i) 
	  SerialVector[j][i] = (double)SVector[i];
#endif
      SolTime_ += Time_->ElapsedTime();
      
      CheckError();

    }

    // ship solution back and take timing
    Time_->ResetStartTime();
    AMESOS_CHK_ERR(vecX->Export(SerialVector,SerialImporter(),Insert));
    VecTime_ += Time_->ElapsedTime();

  }

  // compute vector norms
  if( ComputeVectorNorms_ == true || verbose_ == 2 ) {
    double NormLHS, NormRHS;
    for( int i=0 ; i<NumVectors ; ++i ) {
      assert((*vecX)(i)->Norm2(&NormLHS)==0);
      assert((*vecB)(i)->Norm2(&NormRHS)==0);
      if( verbose_ && Comm().MyPID() == 0 ) {
	cout << "Amesos_Mumps : vector " << i << ", ||x|| = " << NormLHS
	     << ", ||b|| = " << NormRHS << endl;
      }
    }
  }
  
  // compute true residual
  if( ComputeTrueResidual_ == true || verbose_ == 2  ) {
    double Norm;
    Epetra_MultiVector Ax(vecB->Map(),NumVectors);
    for( int i=0 ; i<NumVectors ; ++i ) {
      (Matrix().Multiply(UseTranspose(), *vecX, Ax));
      (Ax.Update(1.0, *vecB, -1.0));
      (Ax.Norm2(&Norm));
      
      if( verbose_ && Comm().MyPID() == 0 ) {
	cout << "Amesos_Mumps : vector " << i << ", ||Ax - b|| = " << Norm << endl;
      }
    }
  }

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
#ifndef HAVE_AMESOS_SMUMPS
	  CrsSchurComplement_->InsertGlobalValues(i,1,&(MDS.schur[pos]),&j);
#else
	  double val = (double)MDS.schur[pos];
	  CrsSchurComplement_->InsertGlobalValues(i,1,&val,&j);
#endif
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

void Amesos_Mumps::PrintStatus() 
{

  if( Comm().MyPID() != 0  ) return;

  cout << "----------------------------------------------------------------------------" << endl;
#ifdef HAVE_AMESOS_SMUMPS
  cout << "Amesos_Mumps : Using single precision" << endl;
#endif
  cout << "Amesos_Mumps : Matrix has " << Matrix().NumGlobalRows() << " rows"
       << " and " << Matrix().NumGlobalNonzeros() << " nonzeros" << endl;
  cout << "Amesos_Mumps : Nonzero elements per row = "
       << 1.0*Matrix().NumGlobalNonzeros()/Matrix().NumGlobalRows() << endl;
  cout << "Amesos_Mumps : Percentage of nonzero elements = "
       << 100.0*Matrix().NumGlobalNonzeros()/(pow(Matrix().NumGlobalRows(),2.0)) << endl;
  cout << "Amesos_Mumps : Use transpose = " << UseTranspose_ << endl;
  if( MatrixProperty_ == 0 ) cout << "Amesos_Mumps : Matrix is general unsymmetric" << endl;
  if( MatrixProperty_ == 2 ) cout << "Amesos_Mumps : Matrix is general symmetric" << endl;
  if( MatrixProperty_ == 1 ) cout << "Amesos_Mumps : Matrix is SPD" << endl;
  cout << "Amesos_Mumps : Available process(es) = " << Comm().NumProc() << endl;
  cout << "Amesos_Mumps : Using " << MaxProcs_ << " process(es)" << endl;
  
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
  cout << "Amesos_Mumps : for running factorization = "
       << MDS.INFOG(17) << " Mbytes" << endl;
  cout << "Amesos_Mumps : Allocated during factorization = "
       << MDS.INFOG(19) << " Mbytes" << endl;
  cout << "----------------------------------------------------------------------------" << endl;
 
  return;
  
}

// ================================================ ====== ==== ==== == =

int Amesos_Mumps::SetICNTL(int * icntl)
{
  for( int i=0 ; i<40 ; ++i ) icntl_[i] = icntl[i];
  return 0;

}

// ================================================ ====== ==== ==== == =

int Amesos_Mumps::SetCNTL(double * cntl)  
{
  for( int i=0 ; i<5 ; ++i ) cntl_[i] = cntl[i];
  return 0;
}
 
// ================================================ ====== ==== ==== == =

void Amesos_Mumps::CheckError() 
{
  
  bool Wrong = ((MDS.INFOG(1) != 0) || (MDS.INFO(1) != 0))
               && (Comm().MyPID() < MaxProcs_);
  
  // an error occurred in MUMPS. Print out information and quit.

  if (Comm().MyPID() == 0 && Wrong) {
    cerr << "Amesos_Mumps : ERROR" << endl;
    cerr << "Amesos_Mumps : INFOG(1) = " << MDS.INFOG(1) << endl;
    cerr << "Amesos_Mumps : INFOG(2) = " << MDS.INFOG(2) << endl;
  }
  
  if (MDS.INFO(1) != 0 && Wrong) {
    cerr << "Amesos_Mumps : On process " << Comm().MyPID()
	 << ", INFO(1) = " << MDS.INFO(1) << endl;
    cerr << "Amesos_Mumps : On process " << Comm().MyPID()
	 << ", INFO(2) = " << MDS.INFO(2) << endl;
  }

  if (Wrong) 
    exit(EXIT_FAILURE);
  
}

// ================================================ ====== ==== ==== == =
void Amesos_Mumps::PrintTiming()
{

  if( Comm().MyPID() ) return;

  cout << "----------------------------------------------------------------------------" << endl;
  cout << "Amesos_Mumps : Time to convert matrix to MUMPS format = "
    << ConTime_ << " (s)" << endl;
  if( MaxProcs_ != Comm().NumProc() ) 
    cout << "Amesos_Mumps : Time to redistribute matrix = "
      << MatTime_ << " (s)" << endl;
  cout << "Amesos_Mumps : Time to redistribute vectors = "
    << VecTime_ << " (s)" << endl;
  cout << "Amesos_Mumps : Number of symbolic factorizations = "
    << NumSymbolicFact_ << endl;
  cout << "Amesos_Mumps : Time for sym fact = "
    << SymTime_ << " (s), avg = " << SymTime_/NumSymbolicFact_
    << " (s)" << endl;
  cout << "Amesos_Mumps : Number of numeric factorizations = "
    << NumNumericFact_ << endl;
  cout << "Amesos_Mumps : Time for num fact = "
    << NumTime_ << " (s), avg = " << NumTime_/NumNumericFact_
    << " (s)" << endl;
  cout << "Amesos_Mumps : Number of solve phases = "
    << NumSolve_ << endl;
  cout << "Amesos_Mumps : Time for solve = "
    << SolTime_ << " (s), avg = " << SolTime_/NumSolve_
    << " (s)" << endl;
  cout << "----------------------------------------------------------------------------" << endl;

  return;
}

// ================================================ ====== ==== ==== == =
Epetra_RowMatrix& Amesos_Mumps::Matrix() 
{
  Epetra_RowMatrix* Matrix = dynamic_cast<Epetra_RowMatrix*>(Problem_->GetOperator());
  assert (Matrix != 0);
  return(*Matrix);
}

// ================================================ ====== ==== ==== == =
Epetra_Map& Amesos_Mumps::RedistrMap() 
{
  assert (Comm().NumProc() != 1);
  if (RedistrMap_ == 0) {
    int i = Matrix().NumGlobalRows() / MaxProcs_;
    if (Comm().MyPID() == 0)
      i += Matrix().NumGlobalRows() % MaxProcs_;
    else if (Comm().MyPID() >= MaxProcs_)
      i = 0;

    RedistrMap_ = new Epetra_Map(Matrix().NumGlobalRows(),i,0,Comm());
    assert (RedistrMap_ != 0);
  }
  return(*RedistrMap_);
}

// ================================================ ====== ==== ==== == =
Epetra_Import& Amesos_Mumps::RedistrImporter()
{
  assert (Comm().NumProc() != 1);

  if (RedistrImporter_ == 0) {
    RedistrImporter_ = new Epetra_Import(RedistrMap(),Matrix().RowMatrixRowMap());
    assert (RedistrImporter_ != 0);
  }
  return(*RedistrImporter_);
}

// ================================================ ====== ==== ==== == =
Epetra_RowMatrix& Amesos_Mumps::RedistrMatrix(const bool ImportMatrix)
{
  if (Comm().NumProc() == 1)
    return(Matrix());

  if (ImportMatrix) {
    delete RedistrMatrix_;
    RedistrMatrix_ = 0;
  }

  if (RedistrMatrix_ == 0) {
    RedistrMatrix_ = new Epetra_CrsMatrix(Copy,RedistrMap(),0);
    if (ImportMatrix) {
      int ierr =RedistrMatrix_->Import(Matrix(),RedistrImporter(),Insert);
      assert (ierr == 0);
      assert(RedistrMatrix_->FillComplete() == 0);
    }
  }

  return(*RedistrMatrix_);
}

// ================================================ ====== ==== ==== == =
Epetra_Map& Amesos_Mumps::SerialMap()
{
  if (SerialMap_ == 0) {
    int i = Matrix().NumGlobalRows();
    if (Comm().MyPID()) i = 0;
    SerialMap_ = new Epetra_Map(-1,i,0,Comm());
    assert (SerialMap_ != 0);
  }
  return(*SerialMap_);
}

// ================================================ ====== ==== ==== == =
Epetra_Import& Amesos_Mumps::SerialImporter()
{ 
  if (SerialImporter_ == 0) {
    SerialImporter_ = new Epetra_Import(SerialMap(),Matrix().OperatorDomainMap());
    assert (SerialImporter_ != 0);
  }
  return(*SerialImporter_);
}
