/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
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
//@HEADER
*/

#include "Ifpack_ConfigDefs.h"
#include "Ifpack_Preconditioner.h"
#include "Ifpack_Amesos.h"
#include "Ifpack_Condest.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Map.h"
#include "Epetra_Comm.h"
#include "Amesos.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_Time.h"
#include "Teuchos_ParameterList.hpp"

static bool FirstTime = true;

//==============================================================================
Ifpack_Amesos::Ifpack_Amesos(Epetra_RowMatrix* Matrix_in) :
  Matrix_(Teuchos::rcp( Matrix_in, false )),
  Label_("Amesos_Klu"),
  IsEmpty_(false),
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
  Problem_ = Teuchos::rcp( new Epetra_LinearProblem );
}

//==============================================================================
Ifpack_Amesos::Ifpack_Amesos(const Ifpack_Amesos& rhs) :
  Matrix_(Teuchos::rcp( &rhs.Matrix(), false )),
  Label_(rhs.Label()),
  IsEmpty_(false),
  IsInitialized_(false),
  IsComputed_(false),
  NumInitialize_(rhs.NumInitialize()),
  NumCompute_(rhs.NumCompute()),
  NumApplyInverse_(rhs.NumApplyInverse()),
  InitializeTime_(rhs.InitializeTime()),
  ComputeTime_(rhs.ComputeTime()),
  ApplyInverseTime_(rhs.ApplyInverseTime()),
  ComputeFlops_(rhs.ComputeFlops()),
  ApplyInverseFlops_(rhs.ApplyInverseFlops()),
  Condest_(rhs.Condest())
{

  Problem_ = Teuchos::rcp( new Epetra_LinearProblem );

  // copy the RHS list in *this.List
  Teuchos::ParameterList RHSList(rhs.List());
  List_ = RHSList;

  // I do not have a copy constructor for Amesos,
  // so Initialize() and Compute() of this object 
  // are called if the rhs did so
  if (rhs.IsInitialized()) {
    IsInitialized_ = true;
    Initialize();
  }
  if (rhs.IsComputed()) {
    IsComputed_ = true;
    Compute();
  }

}
//==============================================================================
int Ifpack_Amesos::SetParameters(Teuchos::ParameterList& List_in)
{

  List_ = List_in;
  Label_ = List_in.get("amesos: solver type", Label_);
  return(0);
}

//==============================================================================
int Ifpack_Amesos::Initialize()
{

  IsEmpty_ = false;
  IsInitialized_ = false;
  IsComputed_ = false;

  if (Matrix_ == Teuchos::null)
    IFPACK_CHK_ERR(-1);

#if 0
  // better to avoid strange games with maps, this class should be
  // used for Ifpack_LocalFilter'd matrices only
  if (Comm().NumProc() != 1) {
    cout << "Class Ifpack_Amesos must be used for serial runs;" << endl;
    cout << "for parallel runs you should declare objects as:" << endl; 
    cout << "Ifpack_AdditiveSchwarz<Ifpack_Amesos> APrec(Matrix)" << endl;
    exit(EXIT_FAILURE);
  }
#endif

  // only square matrices
  if (Matrix_->NumGlobalRows() != Matrix_->NumGlobalCols())
    IFPACK_CHK_ERR(-1);

  // if the matrix has a dimension of 0, this is an empty preconditioning object.
  if (Matrix_->NumGlobalRows() == 0) {
    IsEmpty_ = true;
    IsInitialized_ = true;
    ++NumInitialize_;
    return(0);
  }

  Problem_->SetOperator(const_cast<Epetra_RowMatrix*>(Matrix_.get()));

  // create timer, which also starts it.
  if (Time_ == Teuchos::null)
    Time_ = Teuchos::rcp( new Epetra_Time(Comm()) );

  Amesos Factory;
  Solver_ = Teuchos::rcp( Factory.Create((char*)Label_.c_str(),*Problem_) );
  
  if (Solver_ == Teuchos::null) 
  {
    // try to create KLU, it is generally enabled
    Solver_ = Teuchos::rcp( Factory.Create("Amesos_Klu",*Problem_) );
  }
  if (Solver_ == Teuchos::null)
  {
    // finally try to create LAPACK, it is generally enabled
    // NOTE: the matrix is dense, so this should only be for
    // small problems...
    if (FirstTime)
    {
      cerr << "IFPACK WARNING: In class Ifpack_Amesos:" << endl;
      cerr << "IFPACK WARNING: Using LAPACK because other Amesos" << endl;
      cerr << "IFPACK WARNING: solvers are not available. LAPACK" << endl;
      cerr << "IFPACK WARNING: allocates memory to store the matrix as" << endl;
      cerr << "IFPACK WARNING: dense, I hope you have enough memory..." << endl;
      cerr << "IFPACK WARNING: (file " << __FILE__ << ", line " << __LINE__ 
           << ")" << endl;
      FirstTime = false;
    }
    Solver_ = Teuchos::rcp( Factory.Create("Amesos_Lapack",*Problem_) );
  }
  // if empty, give up.
  if (Solver_ == Teuchos::null)
    IFPACK_CHK_ERR(-1);

  IFPACK_CHK_ERR(Solver_->SetUseTranspose(UseTranspose_));
  Solver_->SetParameters(List_);
  IFPACK_CHK_ERR(Solver_->SymbolicFactorization());

  IsInitialized_ = true;
  ++NumInitialize_;
  InitializeTime_ += Time_->ElapsedTime();
  return(0);
}

//==============================================================================
int Ifpack_Amesos::Compute()
{

  if (!IsInitialized())
    IFPACK_CHK_ERR(Initialize());

  if (IsEmpty_) {
    IsComputed_ = true;
    ++NumCompute_;
    return(0);
  }

  IsComputed_ = false;
  Time_->ResetStartTime();

  if (Matrix_ == Teuchos::null)
    IFPACK_CHK_ERR(-1);

  IFPACK_CHK_ERR(Solver_->NumericFactorization());

  IsComputed_ = true;
  ++NumCompute_;
  ComputeTime_ += Time_->ElapsedTime();
  return(0);
}

//==============================================================================
int Ifpack_Amesos::SetUseTranspose(bool UseTranspose_in)
{
  // store the value in UseTranspose_. If we have the solver, we pass to it
  // right away, otherwise we wait till when it is created.
  UseTranspose_ = UseTranspose_in;
  if (Solver_ != Teuchos::null)
    IFPACK_CHK_ERR(Solver_->SetUseTranspose(UseTranspose_in));

  return(0);
}

//==============================================================================
int Ifpack_Amesos::
Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  // check for maps ???
  IFPACK_CHK_ERR(Matrix_->Apply(X,Y));
  return(0);
}

//==============================================================================
int Ifpack_Amesos::
ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  if (IsEmpty_) {
    ++NumApplyInverse_;
    return(0);
  }

  if (IsComputed() == false)
    IFPACK_CHK_ERR(-1);

  if (X.NumVectors() != Y.NumVectors())
    IFPACK_CHK_ERR(-1); // wrong input
  
  Time_->ResetStartTime();

  // AztecOO gives X and Y pointing to the same memory location,
  // need to create an auxiliary vector, Xcopy
  Teuchos::RefCountPtr<const Epetra_MultiVector> Xcopy;
  if (X.Pointers()[0] == Y.Pointers()[0])
    Xcopy = Teuchos::rcp( new Epetra_MultiVector(X) );
  else
    Xcopy = Teuchos::rcp( &X, false );
    
  Problem_->SetLHS(&Y);
  Problem_->SetRHS((Epetra_MultiVector*)Xcopy.get());
  IFPACK_CHK_ERR(Solver_->Solve());

  ++NumApplyInverse_;
  ApplyInverseTime_ += Time_->ElapsedTime();

  return(0);
}

//==============================================================================
double Ifpack_Amesos::NormInf() const
{
  return(-1.0);
}

//==============================================================================
const char* Ifpack_Amesos::Label() const
{
  return((char*)Label_.c_str());
}

//==============================================================================
bool Ifpack_Amesos::UseTranspose() const
{
  return(UseTranspose_);
}

//==============================================================================
bool Ifpack_Amesos::HasNormInf() const
{
  return(false);
}

//==============================================================================
const Epetra_Comm & Ifpack_Amesos::Comm() const
{
  return(Matrix_->Comm());
}

//==============================================================================
const Epetra_Map & Ifpack_Amesos::OperatorDomainMap() const
{
  return(Matrix_->OperatorDomainMap());
}

//==============================================================================
const Epetra_Map & Ifpack_Amesos::OperatorRangeMap() const
{
  return(Matrix_->OperatorRangeMap());
}

//==============================================================================
double Ifpack_Amesos::Condest(const Ifpack_CondestType CT,
                              const int MaxIters, const double Tol,
			      Epetra_RowMatrix* Matrix_in)
{

  if (!IsComputed()) // cannot compute right now
    return(-1.0);

  if (Condest_ == -1.0)
    Condest_ = Ifpack_Condest(*this, CT, MaxIters, Tol, Matrix_in);

  return(Condest_);
}

//==============================================================================
std::ostream& Ifpack_Amesos::Print(std::ostream& os) const
{
  if (!Comm().MyPID()) {
    os << endl;
    os << "================================================================================" << endl;
    os << "Ifpack_Amesos: " << Label () << endl << endl;
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
