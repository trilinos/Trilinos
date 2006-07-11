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

//=============================================================================
Amesos_Mumps::Amesos_Mumps(const Epetra_LinearProblem &prob ) :
  Problem_(&prob),
  NoDestroy_(false),
  MaxProcs_(-1),
  NumSchurComplementRows_(-1),
  IsComputeSchurComplementOK_(false),
  RowSca_(0),
  ColSca_(0),
  PermIn_(0),
#ifdef HAVE_MPI  
  MUMPSComm_(0),
#endif
  UseTranspose_(false)
{
  // -777 is for me. It means: I never called MUMPS, so
  // SymbolicFactorization should not call Destroy() and ask MUMPS to
  // free its space.
  MDS.job = -777;
  
  // load up my default parameters
  ICNTL[1]  = -1;  // Turn off error messages  6=on, -1 =off
  ICNTL[2]  = -1;  // Turn off diagnostic printing  6=on, -1=off
  ICNTL[3]  = -1;  // Turn off global information messages   6=on, -1=off
  ICNTL[4]  = -1;  // 3 = most msgs; -1= none  
  ICNTL[5]  = 0;   // Matrix is given in elemental (i.e. triplet) from
  ICNTL[6]  = 7;   // Choose column permutation automatically
  ICNTL[7]  = 7;   // Choose ordering method automatically
  ICNTL[8]  = 7;   // Choose scaling automatically
  ICNTL[9]  = 1;   // Compute A x = b 
  ICNTL[10] = 0;   // Maximum steps of iterative refinement
  ICNTL[11] = 0;   // Do not collect statistics
  ICNTL[12] = 0;   // Use Node level parallelism
  ICNTL[13] = 0;   // Use ScaLAPACK for root node 
  ICNTL[14] = 20;  // Increase memory allocation 20% at a time 
  ICNTL[15] = 0;   // Minimize memory use (not flop count)
  ICNTL[16] = 0;   // Do not perform null space detection
  ICNTL[17] = 0;   // Unused (null space dimension)

  Teuchos::ParameterList ParamList;
  SetParameters(ParamList);
}

//=============================================================================
void Amesos_Mumps::Destroy()
{
  if (!NoDestroy_) 
  { 
    // destroy instance of the package
    MDS.job = -2;

    if (Comm().MyPID() < MaxProcs_) dmumps_c(&(MDS));
    
#if 0
    if (IsComputeSchurComplementOK_ && (Comm().MyPID() == 0)
	&& MDS.schur) {
      delete [] MDS.schur;
      MDS.schur = 0;
    }
#endif
    
#ifdef HAVE_MPI
    if (MUMPSComm_) 
    {
      MPI_Comm_free( &MUMPSComm_ );
      MUMPSComm_ = 0;
    }
#endif
    
    if( (verbose_ && PrintTiming_) || verbose_ == 2) PrintTiming();
    if( (verbose_ && PrintStatus_) || verbose_ == 2) PrintStatus();
  }
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
  
#ifdef EXTRA_DEBUG_INFO
  Epetra_CrsMatrix* Eptr = dynamic_cast<Epetra_CrsMatrix*>( ptr );
  if ( ptr->NumGlobalNonzeros() < 300 ) SetICNTL(4,3 );  // Enable more debug info for small matrices
  if ( ptr->NumGlobalNonzeros() < 42 && Eptr ) { 
      cout << " Matrix = " << endl ; 
      Eptr->Print( cout ) ; 
  } else {
      assert( Eptr );
  }
#endif

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
// NOTE: suppose first position is 1 (as in FORTRAN)
void Amesos_Mumps::SetICNTL(int pos, int value)
{
  ICNTL[pos] = value;
}

//=============================================================================
// NOTE: suppose first position is 1 (as in FORTRAN)
void Amesos_Mumps::SetCNTL(int pos, double value)
{
  CNTL[pos] = value;
}

//=============================================================================
void Amesos_Mumps::SetICNTLandCNTL()
{
  map<int,int>::iterator i_iter;
  for (i_iter = ICNTL.begin() ; i_iter != ICNTL.end() ; ++i_iter)
  {
    int pos = i_iter->first;
    int val = i_iter->second;
    if (pos < 1 || pos > 40) continue;
    MDS.ICNTL(pos) = val;
  }

  map<int,double>::iterator d_iter;
  for (d_iter = CNTL.begin() ; d_iter != CNTL.end() ; ++d_iter)
  {
    int pos = d_iter->first;
    double val = d_iter->second;
    if (pos < 1 || pos > 5) continue;
    MDS.CNTL(pos) = val;
  }

  // fix some options
  
  if (Comm().NumProc() != 1) MDS.ICNTL(18)= 3;
  else                       MDS.ICNTL(18)= 0;
  
  if (UseTranspose())  MDS.ICNTL(9) = 0; 
  else                 MDS.ICNTL(9) = 1;
  
  MDS.ICNTL(5) = 0;

#if 0
  if (IsComputeSchurComplementOK_) MDS.ICNTL(19) = 1;
  else                             MDS.ICNTL(19) = 0;

  // something to do if the Schur complement is required.
  if (IsComputeSchurComplementOK_ && Comm().MyPID() == 0) 
  {
    MDS.size_schur = NumSchurComplementRows_;
    MDS.listvar_schur = SchurComplementRows_;
    MDS.schur = new double[NumSchurComplementRows_*NumSchurComplementRows_];
  }
#endif
}

//=============================================================================
int Amesos_Mumps::SetParameters( const Teuchos::ParameterList & ParameterList)
{
  // ========================================= //
  // retrive MUMPS' parameters from list.      //
  // default values defined in the constructor //
  // ========================================= //
  
  // retrive general parameters

  SetStatusParameters(ParameterList);

  SetControlParameters(ParameterList);

  if (ParameterList.isParameter("NoDestroy"))
    NoDestroy_ = ParameterList.get<bool>("NoDestroy");
  
  // can be use using mumps sublist, via "ICNTL(2)"
  //  if (debug_) 
  //    MDS.ICNTL(2) = 6; // Turn on Mumps verbose output 

  // retrive MUMPS' specific parameters
  
  if (ParameterList.isSublist("mumps")) 
  {
    Teuchos::ParameterList MumpsParams = ParameterList.sublist("mumps") ;

    // ICNTL
    for (int i = 1 ; i <= 40 ; ++i)
    {
      char what[80];
      sprintf(what, "ICNTL(%d)", i);
      if (MumpsParams.isParameter(what)) 
        SetICNTL(i, MumpsParams.get<int>(what));
    }

    // CNTL
    for (int i = 1 ; i <= 5 ; ++i)
    {
      char what[80];
      sprintf(what, "CNTL(%d)", i);
      if (MumpsParams.isParameter(what)) 
        SetCNTL(i, MumpsParams.get<double>(what));
    }

    // ordering
     if (MumpsParams.isParameter("PermIn")) 
     {
      int* PermIn = 0;
      PermIn = MumpsParams.get<int*>("PermIn");
      if (PermIn) SetOrdering(PermIn);
     }
     // Col scaling
     if (MumpsParams.isParameter("ColScaling")) 
     {
       double * ColSca = 0;
       ColSca = MumpsParams.get<double *>("ColScaling");
       if (ColSca) SetColScaling(ColSca);
     }
     // Row scaling
     if (MumpsParams.isParameter("RowScaling")) 
     {
       double * RowSca = 0;
       RowSca = MumpsParams.get<double *>("RowScaling");
       if (RowSca) SetRowScaling(RowSca);
     }
     // that's all folks
  }  

  return(0);
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
  
  int OptNumProcs1 = 1 + EPETRA_MAX(NumRows/10000, NumGlobalNonzeros/100000);
  OptNumProcs1 = EPETRA_MIN(Comm().NumProc(),OptNumProcs1);

  // optimal value for MaxProcs == -2

  int OptNumProcs2 = (int)sqrt(1.0 *  Comm().NumProc());
  if (OptNumProcs2 < 1) OptNumProcs2 = 1;

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
  if (MaxProcs_ > Comm().NumProc()) MaxProcs_ = Comm().NumProc();
//  if ( MaxProcs_ > 1 ) MaxProcs_ = Comm().NumProc();     // Bug - bogus kludge here  - didn't work anyway
}

//=============================================================================
int Amesos_Mumps::SymbolicFactorization()
{

  // erase data if present. 
  if (IsSymbolicFactorizationOK_ && MDS.job != -777)
   Destroy();

  IsSymbolicFactorizationOK_ = false;
  IsNumericFactorizationOK_ = false;

  InitTime(Comm());
  
  CheckParameters();
  AMESOS_CHK_ERR(ConvertToTriplet(false));

#if defined(HAVE_MPI) && defined(HAVE_AMESOS_MPI_C2F)
  if (MaxProcs_ != Comm().NumProc()) 
  {
    if(MUMPSComm_) 
      MPI_Comm_free(&MUMPSComm_);

    vector<int> ProcsInGroup(MaxProcs_);
    for (int i = 0 ; i < MaxProcs_ ; ++i) 
      ProcsInGroup[i] = i;

    MPI_Group OrigGroup, MumpsGroup;
    MPI_Comm_group(MPI_COMM_WORLD, &OrigGroup);
    MPI_Group_incl(OrigGroup, MaxProcs_, &ProcsInGroup[0], &MumpsGroup);
    MPI_Comm_create(MPI_COMM_WORLD, MumpsGroup, &MUMPSComm_);
    MDS.comm_fortran = (F_INT) MPI_Comm_c2f( MUMPSComm_);
  } 
  else 
  {
    const Epetra_MpiComm* MpiComm = dynamic_cast<const Epetra_MpiComm*>(&Comm());
    assert (MpiComm != 0);
    MDS.comm_fortran = (F_INT) MPI_Comm_c2f(MpiComm->GetMpiComm());
  }
#else
  // only thing I can do, use MPI_COMM_WORLD. This will work in serial as well
  MDS.comm_fortran = -987654;
#endif
  
  MDS.job = -1  ;     //  Initialization
  MDS.par = 1 ;       //  Host IS involved in computations
//  MDS.sym = MatrixProperty_;
  MDS.sym =  0;       //  MatrixProperty_ is unititalized.  Furthermore MUMPS 
                      //  expects only half of the matrix to be provided for
                      //  symmetric matrices.  Hence setting MDS.sym to be non-zero
                      //  indicating that the matrix is symmetric will only work
                      //  if we change ConvertToTriplet to pass only half of the 
                      //  matrix.  Bug #2331 and Bug #2332 


  RedistrMatrix(true);

  if (Comm().MyPID() < MaxProcs_) 
  {
    dmumps_c(&(MDS));   //  Initialize MUMPS
    static_cast<void>( CheckError( ) );  
  }

  MDS.n = Matrix().NumGlobalRows();

  // fix pointers for nonzero pattern of A. Numerical values
  // will be entered in PerformNumericalFactorization()
  if (Comm().NumProc() != 1) 
  {
    MDS.nz_loc = RedistrMatrix().NumMyNonzeros();

    if (Comm().MyPID() < MaxProcs_) 
    {
      MDS.irn_loc = &Row[0]; 
      MDS.jcn_loc = &Col[0];
    }
  } 
  else 
  {
    if (Comm().MyPID() == 0) 
    {
      MDS.nz = Matrix().NumMyNonzeros();
      MDS.irn = &Row[0]; 
      MDS.jcn = &Col[0]; 
    }
  }

  // scaling if provided by the user
  if (RowSca_ != 0) 
  {
    MDS.rowsca = RowSca_;
    MDS.colsca = ColSca_;
  }

  // given ordering if provided by the user
  if (PermIn_ != 0) 
  {
    MDS.perm_in = PermIn_;
  }

  MDS.job = 1;     // Request symbolic factorization

  SetICNTLandCNTL();

  // Perform symbolic factorization

  ResetTime();

  if (Comm().MyPID() < MaxProcs_) 
    dmumps_c(&(MDS));

  AddTime("symbolic");

  int IntWrong = CheckError()?1:0 ; 
  int AnyWrong;
  Comm().SumAll( &IntWrong, &AnyWrong, 1 ) ; 
  bool Wrong = AnyWrong > 0 ; 


  if ( Wrong ) {
      AMESOS_CHK_ERR( StructurallySingularMatrixError ) ; 
  }

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

  if (Comm().NumProc() != 1) 
  {
    if (Comm().MyPID() < MaxProcs_) 
      MDS.a_loc = &Val[0];
  } 
  else 
    MDS.a = &Val[0];

  // Request numeric factorization 
  MDS.job = 2;
  // Perform numeric factorization
  ResetTime();

  if (Comm().MyPID() < MaxProcs_) {
    dmumps_c(&(MDS));
  }

  AddTime("numeric");
  
  int IntWrong = CheckError()?1:0 ; 
  int AnyWrong;
  Comm().SumAll( &IntWrong, &AnyWrong, 1 ) ; 
  bool Wrong = AnyWrong > 0 ; 


  if ( Wrong ) {
      AMESOS_CHK_ERR( NumericallySingularMatrixError ) ; 
  }

  IsNumericFactorizationOK_ = true;
  NumNumericFact_++;  
  return(0);
}

//=============================================================================

int Amesos_Mumps::Solve()
{ 
  if (IsNumericFactorizationOK_ == false)
    AMESOS_CHK_ERR(NumericFactorization());

  Epetra_MultiVector* vecX = Problem_->GetLHS(); 
  Epetra_MultiVector* vecB = Problem_->GetRHS();
  int NumVectors = vecX->NumVectors();

  if ((vecX == 0) || (vecB == 0))
    AMESOS_CHK_ERR(-1);

  if (NumVectors != vecB->NumVectors())
    AMESOS_CHK_ERR(-1);

  if (Comm().NumProc() == 1) 
  {
    // do not import any data
    for (int j = 0 ; j < NumVectors; j++) 
    {
      ResetTime();

      MDS.job = 3;     // Request solve

      for (int i = 0 ; i < Matrix().NumMyRows() ; ++i) 
	(*vecX)[j][i] = (*vecB)[j][i];
      MDS.rhs = (*vecX)[j];

      dmumps_c(&(MDS)) ;  // Perform solve
      static_cast<void>( CheckError( ) );   // Can hang 
      AddTime("solve");
    }
  } 
  else 
  {
    Epetra_MultiVector SerialVector(SerialMap(),NumVectors);

    ResetTime();
    AMESOS_CHK_ERR(SerialVector.Import(*vecB,SerialImporter(),Insert));
    AddTime("vector redistribution");
    
    for (int j = 0 ; j < NumVectors; j++) 
    {
      if (Comm().MyPID() == 0)
	MDS.rhs = SerialVector[j];

      // solve the linear system and take time
      MDS.job = 3;     
      ResetTime();
      if (Comm().MyPID() < MaxProcs_) 
	dmumps_c(&(MDS)) ;  // Perform solve
      static_cast<void>( CheckError( ) );   // Can hang 

      AddTime("solve");
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

  NumSolve_++;  
  
  return(0) ; 
}

#if 0
//=============================================================================
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

//=============================================================================
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
  
  return(0);
}

//=============================================================================
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
#endif

//=============================================================================
void Amesos_Mumps::PrintStatus() const
{

  if (Comm().MyPID() != 0 ) return;

  //  The following lines are commented out to deal with bug #1887 - kss
#ifndef IRIX64
  PrintLine();
  cout << "Amesos_Mumps : Matrix has " << Matrix().NumGlobalRows() << " rows"
       << " and " << Matrix().NumGlobalNonzeros() << " nonzeros" << endl;
  cout << "Amesos_Mumps : Nonzero elements per row = "
       << 1.0*Matrix().NumGlobalNonzeros()/Matrix().NumGlobalRows() << endl;
  cout << "Amesos_Mumps : Percentage of nonzero elements = "
       << 100.0*Matrix().NumGlobalNonzeros()/(pow(Matrix().NumGlobalRows(),2.0)) << endl;
  cout << "Amesos_Mumps : Use transpose = " << UseTranspose_ << endl;
//  MatrixProperty_ is unused - see bug #2331 and bug #2332 in this file and bugzilla
  if (MatrixProperty_ == 0) cout << "Amesos_Mumps : Matrix is general unsymmetric" << endl;
  if (MatrixProperty_ == 2) cout << "Amesos_Mumps : Matrix is general symmetric" << endl;
  if (MatrixProperty_ == 1) cout << "Amesos_Mumps : Matrix is SPD" << endl;
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
  PrintLine();
#endif
}

//=============================================================================
int Amesos_Mumps::CheckError() 
{
  bool Wrong = ((MDS.INFOG(1) != 0) || (MDS.INFO(1) != 0))
               && (Comm().MyPID() < MaxProcs_);
  
  // an error occurred in MUMPS. Print out information and quit.

  if (Comm().MyPID() == 0 && Wrong) 
  {
    cerr << "Amesos_Mumps : ERROR" << endl;
    cerr << "Amesos_Mumps : INFOG(1) = " << MDS.INFOG(1) << endl;
    cerr << "Amesos_Mumps : INFOG(2) = " << MDS.INFOG(2) << endl;
  }
  
  if (MDS.INFO(1) != 0 && Wrong) 
  {
    cerr << "Amesos_Mumps : On process " << Comm().MyPID()
	 << ", INFO(1) = " << MDS.INFO(1) << endl;
    cerr << "Amesos_Mumps : On process " << Comm().MyPID()
	 << ", INFO(2) = " << MDS.INFO(2) << endl;
  }

  if (Wrong) 
      return 1 ;
  else
      return 0 ;
}

// ======================================================================
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

  if (NumSymbolicFact_) SymTime /= NumSymbolicFact_;
  if (NumNumericFact_)  NumTime /= NumNumericFact_;
  if (NumSolve_)        SolTime /= NumSolve_;

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
}

// ===========================================================================
Epetra_RowMatrix& Amesos_Mumps::Matrix() 
{
  Epetra_RowMatrix* Matrix = dynamic_cast<Epetra_RowMatrix*>(Problem_->GetOperator());
  assert (Matrix != 0);
  return(*Matrix);
}

// ===========================================================================
const Epetra_RowMatrix& Amesos_Mumps::Matrix() const
{
  Epetra_RowMatrix* Matrix = dynamic_cast<Epetra_RowMatrix*>(Problem_->GetOperator());
  assert (Matrix != 0);
  return(*Matrix);
}

// ===========================================================================
Epetra_Map& Amesos_Mumps::RedistrMap() 
{
  assert (Comm().NumProc() != 1);

  if (RedistrMap_ == Teuchos::null) {
    int i = Matrix().NumGlobalRows() / MaxProcs_;
    if (Comm().MyPID() == 0)
      i += Matrix().NumGlobalRows() % MaxProcs_;
    else if (Comm().MyPID() >= MaxProcs_)
      i = 0;

    RedistrMap_ = rcp(new Epetra_Map(Matrix().NumGlobalRows(),i,0,Comm()));
    assert (RedistrMap_.get() != 0);
  }
  return(*RedistrMap_);
}

// ===========================================================================
Epetra_Import& Amesos_Mumps::RedistrImporter()
{
  assert (Comm().NumProc() != 1);

  if (RedistrImporter_ == null) {
    RedistrImporter_ = rcp(new Epetra_Import(RedistrMap(),Matrix().RowMatrixRowMap()));
    assert (RedistrImporter_.get() != 0);
  }
  return(*RedistrImporter_);
}

// ===========================================================================
Epetra_RowMatrix& Amesos_Mumps::RedistrMatrix(const bool ImportMatrix)
{
  if (Comm().NumProc() == 1) return(Matrix());

  if (ImportMatrix) RedistrMatrix_ = null;

  if (RedistrMatrix_ == null) 
  {
    RedistrMatrix_ = rcp(new Epetra_CrsMatrix(Copy,RedistrMap(),0));
    if (ImportMatrix) {
      int ierr = RedistrMatrix_->Import(Matrix(),RedistrImporter(),Insert);
      assert (ierr == 0);
      ierr = RedistrMatrix_->FillComplete();
      assert (ierr == 0);
    }
  }

  return(*RedistrMatrix_);
}

// ===========================================================================
Epetra_Map& Amesos_Mumps::SerialMap()
{
  if (SerialMap_ == null) 
  {
    int i = Matrix().NumGlobalRows();
    if (Comm().MyPID()) i = 0;
    SerialMap_ = rcp(new Epetra_Map(-1,i,0,Comm()));
    assert (SerialMap_ != null);
  }
  return(*SerialMap_);
}

// ===========================================================================
Epetra_Import& Amesos_Mumps::SerialImporter()
{ 
  if (SerialImporter_ == null) {
    SerialImporter_ = rcp(new Epetra_Import(SerialMap(),Matrix().OperatorDomainMap()));
    assert (SerialImporter_ != null);
  }
  return(*SerialImporter_);
}

// ===========================================================================
double* Amesos_Mumps::GetRINFO() 
{
  return ( MDS.rinfo);
}

// ===========================================================================
int* Amesos_Mumps::GetINFO() 
{
  return (MDS.info);
}

// ===========================================================================
double* Amesos_Mumps::GetRINFOG()
{
  return (MDS.rinfog);
}

// ===========================================================================
int* Amesos_Mumps::GetINFOG()
{
  return (MDS.infog);
}

