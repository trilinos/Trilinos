// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Ifpack_ConfigDefs.h"
#ifdef HAVE_IFPACK_ML
#include "Ifpack_Preconditioner.h"
#include "Ifpack_ML.h"
#include "Ifpack_Condest.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Map.h"
#include "Epetra_Comm.h"
#include "ml_include.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_Time.h"
#include "Teuchos_ParameterList.hpp"

static bool FirstTime = true;

//==============================================================================
Ifpack_ML::Ifpack_ML(Epetra_RowMatrix* Matrix) :
  Matrix_(Teuchos::rcp( Matrix, false )),
  Label_("ML"),
  IsInitialized_(false),
  IsComputed_(false),
  UseTranspose_(false),
  NumInitialize_(0),
  NumCompute_(0),
  NumApplyInverse_(0),
  InitializeTime_(0.0),
  ComputeTime_(0.0),
  ApplyInverseTime_(0.0),
  ComputeFlops_(0),
  ApplyInverseFlops_(0),
  Condest_(-1.0)
{
}

//==============================================================================
int Ifpack_ML::SetParameters(Teuchos::ParameterList& List)
{

  List_ = List;
  return(0);
}

//==============================================================================
int Ifpack_ML::Initialize()
{

  IsInitialized_ = false;
  IsComputed_ = false;

  if (Matrix_ == Teuchos::null)
    IFPACK_CHK_ERR(-1);

  // only square matrices
  if (Matrix_->NumGlobalRows() != Matrix_->NumGlobalCols())
    IFPACK_CHK_ERR(-1);

  // at least one nonzero
  if (Matrix_->NumMyNonzeros() == 0) 
    IFPACK_CHK_ERR(-1);

  if (Time_ == Teuchos::null)
    Time_ = Teuchos::rcp( new Epetra_Time(Comm()) );

  MLPrec_ = Teuchos::rcp( new ML_Epetra::MultiLevelPreconditioner(Matrix(),
                          List_, false) );
  
  // if empty, give up.
  if (MLPrec_ == Teuchos::null)
    IFPACK_CHK_ERR(-1);

  IsInitialized_ = true;
  ++NumInitialize_;
  InitializeTime_ += Time_->ElapsedTime();
  return(0);
}

//==============================================================================
int Ifpack_ML::Compute()
{

  if (!IsInitialized())
    IFPACK_CHK_ERR(Initialize());

  IsComputed_ = false;
  Time_->ResetStartTime();

  if (Matrix_ == Teuchos::null)
    IFPACK_CHK_ERR(-1);

  IFPACK_CHK_ERR(MLPrec_->ComputePreconditioner());

  IsComputed_ = true;
  ++NumCompute_;
  ComputeTime_ += Time_->ElapsedTime();
  return(0);
}

//==============================================================================
int Ifpack_ML::SetUseTranspose(bool UseTranspose)
{
  UseTranspose_ = UseTranspose;
  return(0);
}

//==============================================================================
int Ifpack_ML::
Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  // check for maps ???
  IFPACK_CHK_ERR(Matrix_->Apply(X,Y));
  return(0);
}

//==============================================================================
int Ifpack_ML::
ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{

  if (IsComputed() == false)
    IFPACK_CHK_ERR(-1);

  int numVectors = X.NumVectors();
  if (numVectors != Y.NumVectors())
    IFPACK_CHK_ERR(-1); // wrong input

  
  Time_->ResetStartTime();

  // AztecOO gives X and Y pointing to the same memory location,
  // need to create an auxiliary vector, Xcopy
  Teuchos::RefCountPtr<const Epetra_MultiVector> Xcopy;
  if (X.Pointers()[0] == Y.Pointers()[0])
    Xcopy = Teuchos::rcp( new Epetra_MultiVector(X) );
  else
    Xcopy = Teuchos::rcp( &X, false );

  for (int i=0; i<numVectors; i++) {
    
  IFPACK_CHK_ERR(MLPrec_->ApplyInverse(*Xcopy,Y));
  }

  ++NumApplyInverse_;
  ApplyInverseTime_ += Time_->ElapsedTime();

  return(0);
}

//==============================================================================
double Ifpack_ML::NormInf() const
{
  return(-1.0);
}

//==============================================================================
const char* Ifpack_ML::Label() const
{
  return((char*)Label_.c_str());
}

//==============================================================================
bool Ifpack_ML::UseTranspose() const
{
  return(UseTranspose_);
}

//==============================================================================
bool Ifpack_ML::HasNormInf() const
{
  return(false);
}

//==============================================================================
const Epetra_Comm & Ifpack_ML::Comm() const
{
  return(Matrix_->Comm());
}

//==============================================================================
const Epetra_Map & Ifpack_ML::OperatorDomainMap() const
{
  return(Matrix_->OperatorDomainMap());
}

//==============================================================================
const Epetra_Map & Ifpack_ML::OperatorRangeMap() const
{
  return(Matrix_->OperatorRangeMap());
}

//==============================================================================
double Ifpack_ML::Condest(const Ifpack_CondestType CT,
                              const int MaxIters, const double Tol,
			      Epetra_RowMatrix* Matrix)
{

  if (!IsComputed()) // cannot compute right now
    return(-1.0);

  if (Condest_ == -1.0)
    Condest_ = Ifpack_Condest(*this, CT, MaxIters, Tol, Matrix);

  return(Condest_);
}

//==============================================================================
std::ostream& Ifpack_ML::Print(std::ostream& os) const
{
  if (!Comm().MyPID()) {
    os << endl;
    os << "================================================================================" << endl;
    os << "Ifpack_ML: " << Label () << endl << endl;
    os << "Condition number estimate = " << Condest() << endl;
    os << "Global number of rows            = " << Matrix_->NumGlobalRows() << endl;
    os << endl;
    os << "Phase           # calls   Total Time (s)       Total MFlops     MFlops/s" << endl;
    os << "-----           -------   --------------       ------------     --------" << endl;
    os << "Initialize()    "   << std::setw(5) << NumInitialize_ 
       << "  " << std::setw(15) << InitializeTime_ 
       << "              0.0              0.0" << endl;
    os << "Compute()       "   << std::setw(5) << NumCompute_ 
       << "  " << std::setw(15) << ComputeTime_
       << "  " << std::setw(15) << 1.0e-6 * ComputeFlops_;
    if (ComputeTime_ != 0.0) 
      os << "  " << std::setw(15) << 1.0e-6 * ComputeFlops_ / ComputeTime_ << endl;
    else
      os << "  " << std::setw(15) << 0.0 << endl;
    os << "ApplyInverse()  "   << std::setw(5) << NumApplyInverse_ 
       << "  " << std::setw(15) << ApplyInverseTime_
       << "  " << std::setw(15) << 1.0e-6 * ApplyInverseFlops_;
    if (ApplyInverseTime_ != 0.0) 
      os << "  " << std::setw(15) << 1.0e-6 * ApplyInverseFlops_ / ApplyInverseTime_ << endl;
    else
      os << "  " << std::setw(15) << 0.0 << endl;
    os << "================================================================================" << endl;
    os << endl;
  }

  return(os);
}
#endif // HAVE_IFPACK_ML 
