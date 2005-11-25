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
extern "C" {
#ifndef HAVE_AMESOS_SMUMPS
#include "dmumps_c.h"
#else
#include "smumps_c.h"
#endif
}
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_Util.h"
#include "Epetra_Time.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_IntVector.h"
#include "Epetra_Vector.h"
#include "Epetra_SerialDenseMatrix.h"
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
  
class Amesos_Mumps_Pimpl {
public:
#ifndef HAVE_AMESOS_SMUMPS  
  //! Mumps data structure for double-precision
  DMUMPS_STRUC_C MDS;
#else
  //! Mumps data structure for single-precision
  SMUMPS_STRUC_C MDS;
#endif
} ;
//=============================================================================

Amesos_Mumps::Amesos_Mumps(const Epetra_LinearProblem &prob ) :
  PrivateMumpsData_( rcp( new Amesos_Mumps_Pimpl() ) ),
  Problem_(&prob),
  NoDestroy_(false),
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
  RowSca_(0),
  ColSca_(0),
  PermIn_(0),
  Maxis_(DEF_VALUE_INT),
  Maxs_(DEF_VALUE_INT),
  //  Threshold_(0.0),
  MUMPSComm_(0),
  UseTranspose_(false)

{
  // -777 is for me. It means : never called MUMPS, so
  // SymbolicFactorization will not call Destroy();
  PrivateMumpsData_->MDS.job = -777;
  
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
  if (!NoDestroy_) 
  { 
    // destroy instance of the package
    PrivateMumpsData_->MDS.job = -2;
    
    if (Comm().MyPID() < MaxProcs_) 
      MUMPS_INTERFACE(&(PrivateMumpsData_->MDS));
    
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
	&& PrivateMumpsData_->MDS.schur) {
      delete [] PrivateMumpsData_->MDS.schur;
      PrivateMumpsData_->MDS.schur = 0;
    }
    
    if (MUMPSComm_) {
      MPI_Comm_free( &MUMPSComm_ );
      MUMPSComm_ = 0;
    }
    
    if( (verbose_ && PrintTiming_) || verbose_ == 2) 
      PrintTiming();
    if( (verbose_ && PrintStatus_) || verbose_ == 2) 
      PrintStatus();
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

  ResetTime();


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
      
      // MS // Added on 15-Mar-05.
      if (AddToDiag_ && Indices[j] == i)
        Values[j] += AddToDiag_;

      Val[count] = Values[j];
      count++;
    }
  }

  AddTime("matrix conversion");
  
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
  
  PrivateMumpsData_->MDS.ICNTL(1)  = -1;  // Turn off error messages
  PrivateMumpsData_->MDS.ICNTL(2)  = -1;  // Turn off diagnostic printing
  PrivateMumpsData_->MDS.ICNTL(3)  = -1;  // Turn off global information messages
  PrivateMumpsData_->MDS.ICNTL(4)  = -1;
  PrivateMumpsData_->MDS.ICNTL(5)  = 0;   // Matrix is given in elemental (i.e. triplet) from
  PrivateMumpsData_->MDS.ICNTL(6)  = 7;   // Choose column permutation automatically
  PrivateMumpsData_->MDS.ICNTL(7)  = 7;   // Choose ordering method automatically
  PrivateMumpsData_->MDS.ICNTL(8)  = 7;   // Choose scaling automatically
  PrivateMumpsData_->MDS.ICNTL(9)  = 1;   // Compute A x = b 
  PrivateMumpsData_->MDS.ICNTL(10) = 0;   // Maximum steps of iterative refinement
  PrivateMumpsData_->MDS.ICNTL(11) = 0;   // Do not collect statistics
  PrivateMumpsData_->MDS.ICNTL(12) = 0;   // Use Node level parallelism
  PrivateMumpsData_->MDS.ICNTL(13) = 0;   // Use ScaLAPACK for root node 
  PrivateMumpsData_->MDS.ICNTL(14) = 20;  // Increase memory allocation 20% at a time 
  PrivateMumpsData_->MDS.ICNTL(15) = 0;   // Minimize memory use (not flop count)
  PrivateMumpsData_->MDS.ICNTL(16) = 0;   // Do not perform null space detection
  PrivateMumpsData_->MDS.ICNTL(17) = 0;   // Unused (null space dimension)
  
  if (Comm().NumProc() != 1) 
    PrivateMumpsData_->MDS.ICNTL(18)= 3;
  else
    PrivateMumpsData_->MDS.ICNTL(18)= 0;
  
  if ( UseTranspose() )  PrivateMumpsData_->MDS.ICNTL(9) = 0 ; 
  else                   PrivateMumpsData_->MDS.ICNTL(9) = 1 ;
  
  if( IsComputeSchurComplementOK_ ) PrivateMumpsData_->MDS.ICNTL(19) = 1;
  else                              PrivateMumpsData_->MDS.ICNTL(19) = 0;

  // something to do if the Schur complement is required.
  if( IsComputeSchurComplementOK_ && Comm().MyPID() == 0 ) {
    PrivateMumpsData_->MDS.size_schur = NumSchurComplementRows_;
    PrivateMumpsData_->MDS.listvar_schur = SchurComplementRows_;
    PrivateMumpsData_->MDS.schur = new AMESOS_TYPE[NumSchurComplementRows_*NumSchurComplementRows_];
  }
  
  // retrive user's specified options
  for( int i=0 ; i<40 ; ++i ) {
    if( icntl_[i] != DEF_VALUE_INT )
      PrivateMumpsData_->MDS.ICNTL(i+1) = icntl_[i];
  }
  for( int i=0 ; i<5 ; ++i ) {
    if( icntl_[i] != DEF_VALUE_DOUBLE )
      PrivateMumpsData_->MDS.CNTL(i+1) = cntl_[i];
  }
  
  // take care that required options are not overwritten
  assert(PrivateMumpsData_->MDS.ICNTL(5)== 0); 
  
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

  SetStatusParameters( ParameterList );

  SetControlParameters( ParameterList );

  //  Tthreshold_ is unused at the moment
  //  // ignore all elements below given tolerance
  //  if( ParameterList.isParameter("Threshold") )
  //    Threshold_ = ParameterList.get("Threshold", 0.0);

  if( ParameterList.isParameter("NoDestroy") )
    NoDestroy_ = ParameterList.get("NoDestroy", false);
  
  if ( debug_ ) 
    PrivateMumpsData_->MDS.ICNTL(2)=6; // Turn on Mumps verbose output 


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
  if (IsSymbolicFactorizationOK_ && PrivateMumpsData_->MDS.job != -777) {
   Destroy();
  }

  IsSymbolicFactorizationOK_ = false;
  IsNumericFactorizationOK_ = false;

  InitTime(Comm());
  
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
    PrivateMumpsData_->MDS.comm_fortran = (F_INT) MPI_Comm_c2f( MUMPSComm_);

    delete [] ProcsInGroup;

  } else {
    const Epetra_MpiComm* MpiComm = dynamic_cast<const Epetra_MpiComm*>(&Comm());
    assert (MpiComm != 0);
    PrivateMumpsData_->MDS.comm_fortran = (F_INT) MPI_Comm_c2f(MpiComm->GetMpiComm());
  }
#else
  // only thing I can do, use MPI_COMM_WORLD.
  PrivateMumpsData_->MDS.comm_fortran = -987654;
#endif
  
  PrivateMumpsData_->MDS.job = -1  ;     //  Initialization
  PrivateMumpsData_->MDS.par = 1 ;       //  Host IS involved in computations
  PrivateMumpsData_->MDS.sym = MatrixProperty_;

  RedistrMatrix(true);

  if (Comm().MyPID() < MaxProcs_) {
    MUMPS_INTERFACE(&(PrivateMumpsData_->MDS));   //  Initialize MUMPS

    CheckError();
  }

  PrivateMumpsData_->MDS.n = Matrix().NumGlobalRows() ;

  // fix pointers for nonzero pattern of A. Numerical values
  // will be entered in PerformNumericalFactorization()
  if (Comm().NumProc() != 1) {

    PrivateMumpsData_->MDS.nz_loc = RedistrMatrix().NumMyNonzeros();

    if (Comm().MyPID() < MaxProcs_) {
      PrivateMumpsData_->MDS.irn_loc = &Row[0]; 
      PrivateMumpsData_->MDS.jcn_loc = &Col[0];

    }
  } else {
    if (Comm().MyPID() == 0) {
      PrivateMumpsData_->MDS.nz = Matrix().NumMyNonzeros();
      PrivateMumpsData_->MDS.irn = &Row[0]; 
      PrivateMumpsData_->MDS.jcn = &Col[0]; 
    }
  }

  // scaling if provided by the user
  if( RowSca_ != 0 ) {
    PrivateMumpsData_->MDS.rowsca = RowSca_;
    PrivateMumpsData_->MDS.colsca = ColSca_;
  }

  // given ordering if provided by the user
  if( PermIn_ != 0 ) {
    PrivateMumpsData_->MDS.perm_in = PermIn_;
  }

  //  if( Maxis_ != DEF_VALUE_INT ) PrivateMumpsData_->MDS.maxis = Maxis_;
  //  if( Maxs_ != DEF_VALUE_INT ) PrivateMumpsData_->MDS.maxs = Maxs_;

  PrivateMumpsData_->MDS.job = 1  ;     // Request symbolic factorization

  SetICNTLandCNTL();

  // Perform symbolic factorization

  ResetTime();

  if (Comm().MyPID() < MaxProcs_) 
    MUMPS_INTERFACE(&(PrivateMumpsData_->MDS));

  AddTime("symbolic");

  CheckError();

  IsSymbolicFactorizationOK_ = true ;
  NumSymbolicFact_++;  
  return 0;

}

//=============================================================================

int Amesos_Mumps::NumericFactorization()
{

  IsNumericFactorizationOK_ = false;
  
  if (IsSymbolicFactorizationOK_ == false)
    AMESOS_CHK_ERR(SymbolicFactorization());

  RedistrMatrix(true);
  AMESOS_CHK_ERR(ConvertToTriplet(true));

#ifndef HAVE_AMESOS_SMUMPS
  if (Comm().NumProc() != 1) {
    if (Comm().MyPID() < MaxProcs_) 
      PrivateMumpsData_->MDS.a_loc = &Val[0];
  } else {
    PrivateMumpsData_->MDS.a = &Val[0];
  }
#else
  SVal.resize(Val.size());


  //
  //  Pass the numeric values
  //
  if (Comm().MyPID() < MaxProcs_) 
    for (int i = 0 ; i < Val.size() ; ++i) 
      SVal[i] = (float)(Val[i]);

  if (Comm().NumProc() != 1) {
    if (Comm().MyPID() < MaxProcs_) 
      PrivateMumpsData_->MDS.a_loc = &SVal[0];
  } else {
    PrivateMumpsData_->MDS.a = &SVal[0];
  }
#endif

  // Request numeric factorization 
  PrivateMumpsData_->MDS.job = 2;
  // Perform numeric factorization
  ResetTime();

  if (Comm().MyPID() < MaxProcs_) {
    MUMPS_INTERFACE(&(PrivateMumpsData_->MDS));
  }

  AddTime("numeric");
  
  CheckError();

  IsNumericFactorizationOK_ = true;
  NumNumericFact_++;  
  return 0;
}

//=============================================================================

int Amesos_Mumps::Solve()
{ 

  if (IsNumericFactorizationOK_ == false)
    AMESOS_CHK_ERR(NumericFactorization());

  NumSolve_++;  
  
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

      ResetTime();

      PrivateMumpsData_->MDS.job = 3;     // Request solve

#ifndef HAVE_AMESOS_SMUMPS      
      for (int i = 0 ; i < Matrix().NumMyRows() ; ++i) 
	(*vecX)[j][i] = (*vecB)[j][i];
      PrivateMumpsData_->MDS.rhs = (*vecX)[j];
      MUMPS_INTERFACE(&(PrivateMumpsData_->MDS)) ;  // Perform solve
#else
      for (int i = 0 ; i < Matrix().NumMyRows() ; ++i) 
	SVector[i] = (float)(*vecB)[j][i];
      PrivateMumpsData_->MDS.rhs = &SVector[0];
      MUMPS_INTERFACE(&(PrivateMumpsData_->MDS)) ;  // Perform solve
      for (int i = 0 ; i < Matrix().NumMyRows() ; ++i) 
	(*vecX)[j][i] = (double)SVector[i];
#endif
      AddTime("solve");
      
      CheckError();
    }

  } else {

    Epetra_MultiVector SerialVector(SerialMap(),NumVectors);

    ResetTime();
    AMESOS_CHK_ERR(SerialVector.Import(*vecB,SerialImporter(),Insert));
    AddTime("vector redistribution");
    
    for (int j = 0 ; j < NumVectors; j++) { 

      if (Comm().MyPID() == 0) {
#ifndef HAVE_AMESOS_SMUMPS
	PrivateMumpsData_->MDS.rhs = SerialVector[j];
#else
	// copy the double-precision into the single-precision SVector
	for (int i = 0 ; i < Matrix().NumGlobalRows() ; ++i) 
	  SVector[i] = (float)(SerialVector[j][i]);
	PrivateMumpsData_->MDS.rhs = &SVector[0];
#endif
      }
      // solve the linear system and take time
      PrivateMumpsData_->MDS.job = 3;     
      ResetTime();
      if (Comm().MyPID() < MaxProcs_) 
	MUMPS_INTERFACE(&(PrivateMumpsData_->MDS)) ;  // Perform solve
#ifdef HAVE_AMESOS_SMUMPS  
      // copy back MUMPS' solution into TargetVector, that will be later
      // redistributed using Epetra utilities
      if (Comm().MyPID() == 0)
	for (int i = 0 ; i < Matrix().NumGlobalRows() ; ++i) 
	  SerialVector[j][i] = (double)SVector[i];
#endif
      AddTime("solve");
      
      CheckError();

    }

    // ship solution back and take timing
    ResetTime();
    AMESOS_CHK_ERR(vecX->Export(SerialVector,SerialImporter(),Insert));
    AddTime("vector redistribution");
  }

  if (ComputeTrueResidual_)
    ComputeTrueResidual(Matrix(), *vecX, *vecB, UseTranspose(), "Amesos_Mumps");

  if (ComputeVectorNorms_)
    ComputeVectorNorms(*vecX, *vecB, "Amesos_Mumps");

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
	  CrsSchurComplement_->InsertGlobalValues(i,1,&(PrivateMumpsData_->MDS.schur[pos]),&j);
#else
	  double val = (double)PrivateMumpsData_->MDS.schur[pos];
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
	(*DenseSchurComplement_)(i,j) = PrivateMumpsData_->MDS.schur[pos];
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

void Amesos_Mumps::PrintStatus() const
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
       << PrivateMumpsData_->MDS.RINFOG(1) << endl;
  cout << "Amesos_Mumps : Total FLOPS for assembly = "
       << PrivateMumpsData_->MDS.RINFOG(2) << endl;
  cout << "Amesos_Mumps : Total FLOPS for elimination = "
       << PrivateMumpsData_->MDS.RINFOG(3) << endl;
  
  cout << "Amesos_Mumps : Total real space to store the LU factors = "
       << PrivateMumpsData_->MDS.INFOG(9) << endl;
  cout << "Amesos_Mumps : Total integer space to store the LU factors = "
       << PrivateMumpsData_->MDS.INFOG(10) << endl;
  cout << "Amesos_Mumps : Total number of iterative steps refinement = "
       << PrivateMumpsData_->MDS.INFOG(15) << endl;
  cout << "Amesos_Mumps : Estimated size of MUMPS internal data\n"
       << "Amesos_Mumps : for running factorization = "
       << PrivateMumpsData_->MDS.INFOG(16) << " Mbytes" << endl;
  cout << "Amesos_Mumps : for running factorization = "
       << PrivateMumpsData_->MDS.INFOG(17) << " Mbytes" << endl;
  cout << "Amesos_Mumps : Allocated during factorization = "
       << PrivateMumpsData_->MDS.INFOG(19) << " Mbytes" << endl;
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
  
  bool Wrong = ((PrivateMumpsData_->MDS.INFOG(1) != 0) || (PrivateMumpsData_->MDS.INFO(1) != 0))
               && (Comm().MyPID() < MaxProcs_);
  
  // an error occurred in MUMPS. Print out information and quit.

  if (Comm().MyPID() == 0 && Wrong) {
    cerr << "Amesos_Mumps : ERROR" << endl;
    cerr << "Amesos_Mumps : INFOG(1) = " << PrivateMumpsData_->MDS.INFOG(1) << endl;
    cerr << "Amesos_Mumps : INFOG(2) = " << PrivateMumpsData_->MDS.INFOG(2) << endl;
  }
  
  if (PrivateMumpsData_->MDS.INFO(1) != 0 && Wrong) {
    cerr << "Amesos_Mumps : On process " << Comm().MyPID()
	 << ", INFO(1) = " << PrivateMumpsData_->MDS.INFO(1) << endl;
    cerr << "Amesos_Mumps : On process " << Comm().MyPID()
	 << ", INFO(2) = " << PrivateMumpsData_->MDS.INFO(2) << endl;
  }

  if (Wrong) 
    exit(EXIT_FAILURE);
  
}

// ================================================ ====== ==== ==== == =
void Amesos_Mumps::PrintTiming() const
{
  if (Problem_->GetOperator() == 0 || Comm().MyPID() != 0)
    return;

  double ConTime = GetTime("conversion");
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

  string p = "Amesos_Mumps : ";
  PrintLine();

  cout << p << "Time to convert matrix to MUMPS format = "
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
      ierr = RedistrMatrix_->FillComplete();
      assert (ierr == 0);
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
AMESOS_TYPE * Amesos_Mumps::GetRINFO() 
{
  return ( PrivateMumpsData_->MDS.rinfo);
}

int * Amesos_Mumps::GetINFO() 
{
  return (PrivateMumpsData_->MDS.info);
}

AMESOS_TYPE * Amesos_Mumps::GetRINFOG()
{
  return (PrivateMumpsData_->MDS.rinfog);
}

int * Amesos_Mumps::GetINFOG()
{
  return (PrivateMumpsData_->MDS.infog);
}

