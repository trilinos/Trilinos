// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Ifpack2_ConfigDefs.hpp"
#include "Ifpack2_Preconditioner.hpp"
#include "Ifpack2_IC.hpp"
#include "Ifpack2_IC_Utils.hpp"
#include "Ifpack2_Condest.hpp"
#include "Tpetra_Comm.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_RowMatrix.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Util.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"
using Teuchos::RCP;
using Teuchos::rcp;

//==============================================================================
Ifpack2_IC::Ifpack2_IC(Tpetra_RowMatrix* A) :
  A_(rcp(A,false)),
  Comm_(A->Comm()),
  UseTranspose_(false),
  Condest_(-1.0),
  Athresh_(0.0),
  Rthresh_(1.0),
  Droptol_(0.0),
  Lfil_(0),
  Aict_(0),
  Lict_(0),
  Ldiag_(0),
  IsInitialized_(false),
  IsComputed_(false),
  NumInitialize_(0),
  NumCompute_(0),
  NumApplyInverse_(0),
  InitializeTime_(0.0),
  ComputeTime_(0),
  ApplyInverseTime_(0),
  ComputeFlops_(0.0),
  ApplyInverseFlops_(0.0)
{
  Teuchos::ParameterList List;
  SetParameters(List);

}
//==============================================================================
Ifpack2_IC::~Ifpack2_IC()
{
  if (Lict_ != 0) {
    Ifpack2_AIJMatrix * Lict = (Ifpack2_AIJMatrix *) Lict_;
    delete [] Lict->ptr;
    delete [] Lict->col;
    delete [] Lict->val;
    delete Lict;
  }
  if (Aict_ != 0) {
    Ifpack2_AIJMatrix * Aict = (Ifpack2_AIJMatrix *) Aict_;
    delete Aict;
  }
  if (Ldiag_ != 0) delete [] Ldiag_;

  IsInitialized_ = false;
  IsComputed_ = false;
}

//==========================================================================
int Ifpack2_IC::SetParameters(Teuchos::ParameterList& List)
{

  Lfil_ = List.get("fact: level-of-fill",Lfil_);
  Athresh_ = List.get("fact: absolute threshold", Athresh_);
  Rthresh_ = List.get("fact: relative threshold", Rthresh_);
  Droptol_ = List.get("fact: drop tolerance", Droptol_);

  // set label
  sprintf(Label_, "TIFPACK IC (fill=%d, drop=%f)",
	  Lfil_, Droptol_);
  return(0);
}

//==========================================================================
int Ifpack2_IC::Initialize()
{

  IsInitialized_ = false;
  // FIXME: construction of ILUK graph must be here
  
  IsInitialized_ = true;
  return(0);
}

//==========================================================================
int Ifpack2_IC::ComputeSetup()
{
  // (re)allocate memory for ICT factors
  U_ = rcp(new Tpetra_CrsMatrix(Copy, Matrix().RowMatrixRowMap(), 
                                Matrix().RowMatrixRowMap(), 0));
  D_ = rcp(new Tpetra_Vector(Matrix().RowMatrixRowMap()));
  
  if (U_.get() == 0 || D_.get() == 0)
    IFPACK2_CHK_ERR(-5); // memory allocation error

  int ierr = 0;
  int i, j;
  int NumIn, NumL, NumU;
  bool DiagFound;
  int NumNonzeroDiags = 0;

  // Get Maximun Row length
  int MaxNumEntries = Matrix().MaxNumEntries();

  vector<int> InI(MaxNumEntries); // Allocate temp space
  vector<int> UI(MaxNumEntries);
  vector<double> InV(MaxNumEntries);
  vector<double> UV(MaxNumEntries);

  double *DV;
  ierr = D_->ExtractView(&DV); // Get view of diagonal

  // First we copy the user's matrix into diagonal vector and U, regardless of fill level

  int NumRows = Matrix().NumMyRows();

  for (i=0; i< NumRows; i++) {

    Matrix().ExtractMyRowCopy(i, MaxNumEntries, NumIn, &InV[0], &InI[0]); // Get Values and Indices
    
    // Split into L and U (we don't assume that indices are ordered).
    NumL = 0; 
    NumU = 0; 
    DiagFound = false;
    
    for (j=0; j< NumIn; j++) {
      int k = InI[j];

      if (k==i) {
	DiagFound = true;
	// Store perturbed diagonal in Tpetra_Vector D_
	DV[i] += Rthresh_ * InV[j] + EPETRA_SGN(InV[j]) * Athresh_; 
      }

      else if (k < 0) return(-1); // Out of range
      else if (i<k && k<NumRows) {
	UI[NumU] = k;
	UV[NumU] = InV[j];
	NumU++;
      }
    }
    
    // Check in things for this row of L and U

    if (DiagFound) NumNonzeroDiags++;
    if (NumU) U_->InsertMyValues(i, NumU, &UV[0], &UI[0]);
    
  }

  U_->FillComplete(Matrix().OperatorDomainMap(), 
		   Matrix().OperatorRangeMap());
  
  int ierr1 = 0;
  if (NumNonzeroDiags<U_->NumMyRows()) ierr1 = 1;
  Matrix().Comm().MaxAll(&ierr1, &ierr, 1);
  IFPACK2_CHK_ERR(ierr);
  
  return(0);
}

//==========================================================================
int Ifpack2_IC::Compute() {

  if (!IsInitialized()) 
    IFPACK2_CHK_ERR(Initialize());

  IsComputed_ = false;
  
  // copy matrix into L and U.
  IFPACK2_CHK_ERR(ComputeSetup());
  
  int i;
  
  int m, n, nz, Nrhs, ldrhs, ldlhs;
  int * ptr=0, * ind;
  double * val, * rhs, * lhs;
  
  int ierr = Tpetra_Util_ExtractHbData(U_.get(), 0, 0, m, n, nz, ptr, ind,
				       val, Nrhs, rhs, ldrhs, lhs, ldlhs);
  if (ierr < 0) 
    IFPACK2_CHK_ERR(ierr);
  
  Ifpack2_AIJMatrix * Aict;
  if (Aict_==0) {
    Aict = new Ifpack2_AIJMatrix;
    Aict_ = (void *) Aict;
  }
  else Aict = (Ifpack2_AIJMatrix *) Aict_;
  Ifpack2_AIJMatrix * Lict;
  if (Lict_==0) {
    Lict = new Ifpack2_AIJMatrix;
    Lict_ = (void *) Lict;
  }
  else {
    Lict = (Ifpack2_AIJMatrix *) Lict_;
    Ifpack2_AIJMatrix_dealloc( Lict );  // Delete storage, crout_ict will allocate it again.
  }
  if (Ldiag_ != 0) delete [] Ldiag_; // Delete storage, crout_ict will allocate it again.
  Aict->val = val;
  Aict->col = ind;
  Aict->ptr = ptr;
  double *DV;
  EPETRA_CHK_ERR(D_->ExtractView(&DV)); // Get view of diagonal
    
  crout_ict(m, Aict, DV, Droptol_, Lfil_, Lict, &Ldiag_);

  // Get rid of unnecessary data
  delete [] ptr;
  
  // Create Epetra View of L from crout_ict
  U_ = rcp(new Tpetra_CrsMatrix(View, A_->RowMatrixRowMap(), A_->RowMatrixRowMap(),0));
  D_ = rcp(new Tpetra_Vector(View, A_->RowMatrixRowMap(), Ldiag_));
  
  ptr = Lict->ptr;
  ind = Lict->col;
  val = Lict->val;
  
  for (i=0; i< m; i++) {
    int NumEntries = ptr[i+1]-ptr[i];
    int * Indices = ind+ptr[i];
    double * Values = val+ptr[i];
    U_->InsertMyValues(i, NumEntries, Values, Indices);
  }
  
  U_->FillComplete(A_->OperatorDomainMap(), A_->OperatorRangeMap());
  D_->Reciprocal(*D_); // Put reciprocal of diagonal in this vector
  
  double current_flops = 2 * nz; // Just an estimate
  double total_flops = 0;
  
  A_->Comm().SumAll(&current_flops, &total_flops, 1); // Get total madds across all PEs
  
  ComputeFlops_ += total_flops; 
  // Now count the rest. NOTE: those flops are *global*
  ComputeFlops_ += (double) U_->NumGlobalNonzeros(); // Accounts for multiplier above
  ComputeFlops_ += (double) D_->GlobalLength(); // Accounts for reciprocal of diagonal
  
  IsComputed_ = true;
  
  return(0);
  
}

//=============================================================================
// This function finds Y such that LDU Y = X or U(trans) D L(trans) Y = X for multiple RHS
int Ifpack2_IC::ApplyInverse(const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, 
			    Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const
{
  
  if (!IsComputed())
    IFPACK2_CHK_ERR(-3); // compute preconditioner first
  
  if (X.NumVectors() != Y.NumVectors()) 
    IFPACK2_CHK_ERR(-2); // Return error: X and Y not the same size
  
  bool Upper = true;
  bool UnitDiagonal = true;
  
  // AztecOO gives X and Y pointing to the same memory location,
  // need to create an auxiliary vector, Xcopy
  RefCountPtr< const Tpetra_MultiVector > Xcopy;
  if (X.Pointers()[0] == Y.Pointers()[0])
    Xcopy = rcp( new Tpetra_MultiVector(X) );
  else
    Xcopy = rcp( &X, false );
  
  U_->Solve(Upper, true, UnitDiagonal, *Xcopy, Y);
  Y.Multiply(1.0, *D_, Y, 0.0); // y = D*y (D_ has inverse of diagonal)
  U_->Solve(Upper, false, UnitDiagonal, Y, Y); // Solve Uy = y
  
    ++NumApplyInverse_;
  ApplyInverseFlops_ += 4.0 * U_->NumGlobalNonzeros();
  ApplyInverseFlops_ += D_->GlobalLength();
  return(0);

}

//=============================================================================
// This function finds X such that LDU Y = X or U(trans) D L(trans) Y = X for multiple RHS
int Ifpack2_IC::Apply(const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, 
		      Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const 
{

  if (X.NumVectors() != Y.NumVectors()) 
    IFPACK2_CHK_ERR(-2); // Return error: X and Y not the same size

  Tpetra_MultiVector * X1 = (Tpetra_MultiVector *) &X;
  Tpetra_MultiVector * Y1 = (Tpetra_MultiVector *) &Y;
  
  U_->Multiply(false, *X1, *Y1);
  Y1->Update(1.0, *X1, 1.0); // Y1 = Y1 + X1 (account for implicit unit diagonal)
  Y1->ReciprocalMultiply(1.0, *D_, *Y1, 0.0); // y = D*y (D_ has inverse of diagonal)
  Tpetra_MultiVector Y1temp(*Y1); // Need a temp copy of Y1
  U_->Multiply(true, Y1temp, *Y1);
  Y1->Update(1.0, Y1temp, 1.0); // (account for implicit unit diagonal)
  return(0);
}

//=============================================================================
double Ifpack2_IC::Condest(const Ifpack2_CondestType CT, 
			  const int MaxIters, const double Tol,
			  Tpetra_RowMatrix* Matrix_in)
{
  if (!IsComputed()) // cannot compute right now
    return(-1.0);
  
  if (Condest_ == -1.0)
    Condest_ = Ifpack2_Condest(*this, CT, MaxIters, Tol, Matrix_in);
  
  return(Condest_);
}

//=============================================================================
std::ostream&
Ifpack2_IC::Print(std::ostream& os) const
{
  if (!Comm().MyPID()) {
    os << endl;
    os << "================================================================================" << endl;
    os << "Ifpack2_IC: " << Label() << endl << endl;
    os << "Level-of-fill      = " << LevelOfFill() << endl;
    os << "Absolute threshold = " << AbsoluteThreshold() << endl;
    os << "Relative threshold = " << RelativeThreshold() << endl;
    os << "Drop tolerance     = " << DropTolerance() << endl;
    os << "Condition number estimate = " << Condest() << endl;
    os << "Global number of rows            = " << A_->NumGlobalRows() << endl;
    if (IsComputed_) {
      os << "Number of nonzeros of H         = " << U_->NumGlobalNonzeros() << endl;
      os << "nonzeros / rows                 = " 
         << 1.0 * U_->NumGlobalNonzeros() / U_->NumGlobalRows() << endl;
    }
    os << endl;
    os << "Phase           # calls   Total Time (s)       Total MFlops     MFlops/s" << endl;
    os << "-----           -------   --------------       ------------     --------" << endl;
    os << "Initialize()    "   << std::setw(5) << NumInitialize() 
       << "  " << std::setw(15) << InitializeTime() 
       << "               0.0            0.0" << endl;
    os << "Compute()       "   << std::setw(5) << NumCompute() 
       << "  " << std::setw(15) << ComputeTime()
       << "  " << std::setw(15) << 1.0e-6 * ComputeFlops();
    if (ComputeTime() != 0.0)
      os << "  " << std::setw(15) << 1.0e-6 * ComputeFlops() / ComputeTime() << endl;
    else
      os << "  " << std::setw(15) << 0.0 << endl;
    os << "ApplyInverse()  "   << std::setw(5) << NumApplyInverse() 
       << "  " << std::setw(15) << ApplyInverseTime()
       << "  " << std::setw(15) << 1.0e-6 * ApplyInverseFlops();
    if (ApplyInverseTime() != 0.0) 
      os << "  " << std::setw(15) << 1.0e-6 * ApplyInverseFlops() / ApplyInverseTime() << endl;
    else
      os << "  " << std::setw(15) << 0.0 << endl;
    os << "================================================================================" << endl;
    os << endl;
  }

  
  return(os);
} 
