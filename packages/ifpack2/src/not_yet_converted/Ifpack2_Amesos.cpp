/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
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

#include "Ifpack2_ConfigDefs.hpp"
#include "Ifpack2_Preconditioner.hpp"
#include "Ifpack2_Amesos.hpp"
#include "Ifpack2_Condest.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Comm.hpp"
#include "Amesos.hpp"
#include "Tpetra_LinearProblem.hpp"
#include "Tpetra_RowMatrix.hpp"
#include "Tpetra_Time.hpp"
#include "Teuchos_ParameterList.hpp"

static bool FirstTime = true;

//==============================================================================
Ifpack2_Amesos::Ifpack2_Amesos(Tpetra_RowMatrix* Matrix_in) :
  Matrix_(Teuchos::rcp( Matrix_in, false )),
  Label_("Amesos_Klu"),
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
  Problem_ = Teuchos::rcp( new Tpetra_LinearProblem );
}

//==============================================================================
Ifpack2_Amesos::Ifpack2_Amesos(const Ifpack2_Amesos& rhs) :
  Matrix_(Teuchos::rcp( &rhs.Matrix(), false )),
  Label_(rhs.Label()),
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

  Problem_ = Teuchos::rcp( new Tpetra_LinearProblem );

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
int Ifpack2_Amesos::SetParameters(Teuchos::ParameterList& List_in)
{

  List_ = List_in;
  Label_ = List_in.get("amesos: solver type", Label_);
  return(0);
}

//==============================================================================
int Ifpack2_Amesos::Initialize()
{

  IsInitialized_ = false;
  IsComputed_ = false;

  if (Matrix_ == Teuchos::null)
    IFPACK2_CHK_ERR(-1);

#if 0
  // better to avoid strange games with maps, this class should be
  // used for Ifpack2_LocalFilter'd matrices only
  if (Comm().NumProc() != 1) {
    cout << "Class Ifpack2_Amesos must be used for serial runs;" << endl;
    cout << "for parallel runs you should declare objects as:" << endl; 
    cout << "Ifpack2_AdditiveSchwarz<Ifpack2_Amesos> APrec(Matrix)" << endl;
    exit(EXIT_FAILURE);
  }
#endif

  // only square matrices
  if (Matrix_->NumGlobalRows() != Matrix_->NumGlobalCols())
    IFPACK2_CHK_ERR(-1);

  // at least one nonzero
  if (Matrix_->NumMyNonzeros() == 0) 
    IFPACK2_CHK_ERR(-1);

  Problem_->SetOperator(const_cast<Tpetra_RowMatrix*>(Matrix_.get()));

  if (Time_ == Teuchos::null)
    Time_ = Teuchos::rcp( new Tpetra_Time(Comm()) );

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
      cerr << "TIFPACK WARNING: In class Ifpack2_Amesos:" << endl;
      cerr << "TIFPACK WARNING: Using LAPACK because other Amesos" << endl;
      cerr << "TIFPACK WARNING: solvers are not available. LAPACK" << endl;
      cerr << "TIFPACK WARNING: allocates memory to store the matrix as" << endl;
      cerr << "TIFPACK WARNING: dense, I hope you have enough memory..." << endl;
      cerr << "TIFPACK WARNING: (file " << __FILE__ << ", line " << __LINE__ 
           << ")" << endl;
      FirstTime = false;
    }
    Solver_ = Teuchos::rcp( Factory.Create("Amesos_Lapack",*Problem_) );
  }
  // if empty, give up.
  if (Solver_ == Teuchos::null)
    IFPACK2_CHK_ERR(-1);

  IFPACK2_CHK_ERR(Solver_->SetUseTranspose(UseTranspose_));
  Solver_->SetParameters(List_);
  IFPACK2_CHK_ERR(Solver_->SymbolicFactorization());

  IsInitialized_ = true;
  ++NumInitialize_;
  InitializeTime_ += Time_->ElapsedTime();
  return(0);
}

//==============================================================================
int Ifpack2_Amesos::Compute()
{

  if (!IsInitialized())
    IFPACK2_CHK_ERR(Initialize());

  IsComputed_ = false;
  Time_->ResetStartTime();

  if (Matrix_ == Teuchos::null)
    IFPACK2_CHK_ERR(-1);

  IFPACK2_CHK_ERR(Solver_->NumericFactorization());

  IsComputed_ = true;
  ++NumCompute_;
  ComputeTime_ += Time_->ElapsedTime();
  return(0);
}

//==============================================================================
int Ifpack2_Amesos::SetUseTranspose(bool UseTranspose_in)
{
  // store the value in UseTranspose_. If we have the solver, we pass to it
  // right away, otherwise we wait till when it is created.
  UseTranspose_ = UseTranspose_in;
  if (Solver_ != Teuchos::null)
    IFPACK2_CHK_ERR(Solver_->SetUseTranspose(UseTranspose_in));

  return(0);
}

//==============================================================================
int Ifpack2_Amesos::
Apply(const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const
{
  // check for maps ???
  IFPACK2_CHK_ERR(Matrix_->Apply(X,Y));
  return(0);
}

//==============================================================================
int Ifpack2_Amesos::
ApplyInverse(const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const
{

  if (IsComputed() == false)
    IFPACK2_CHK_ERR(-1);

  if (X.NumVectors() != Y.NumVectors())
    IFPACK2_CHK_ERR(-1); // wrong input
  
  Time_->ResetStartTime();

  // AztecOO gives X and Y pointing to the same memory location,
  // need to create an auxiliary vector, Xcopy
  Teuchos::RCP<const Tpetra_MultiVector> Xcopy;
  if (X.Pointers()[0] == Y.Pointers()[0])
    Xcopy = Teuchos::rcp( new Tpetra_MultiVector(X) );
  else
    Xcopy = Teuchos::rcp( &X, false );
    
  Problem_->SetLHS(&Y);
  Problem_->SetRHS((Tpetra_MultiVector*)Xcopy.get());
  IFPACK2_CHK_ERR(Solver_->Solve());

  ++NumApplyInverse_;
  ApplyInverseTime_ += Time_->ElapsedTime();

  return(0);
}

//==============================================================================
double Ifpack2_Amesos::NormInf() const
{
  return(-1.0);
}

//==============================================================================
const char* Ifpack2_Amesos::Label() const
{
  return((char*)Label_.c_str());
}

//==============================================================================
bool Ifpack2_Amesos::UseTranspose() const
{
  return(UseTranspose_);
}

//==============================================================================
bool Ifpack2_Amesos::HasNormInf() const
{
  return(false);
}

//==============================================================================
const Tpetra_Comm & Ifpack2_Amesos::Comm() const
{
  return(Matrix_->Comm());
}

//==============================================================================
const Tpetra_Map & Ifpack2_Amesos::OperatorDomainMap() const
{
  return(Matrix_->OperatorDomainMap());
}

//==============================================================================
const Tpetra_Map & Ifpack2_Amesos::OperatorRangeMap() const
{
  return(Matrix_->OperatorRangeMap());
}

//==============================================================================
double Ifpack2_Amesos::Condest(const Ifpack2_CondestType CT,
                              const int MaxIters, const double Tol,
			      Tpetra_RowMatrix* Matrix_in)
{

  if (!IsComputed()) // cannot compute right now
    return(-1.0);

  if (Condest_ == -1.0)
    Condest_ = Ifpack2_Condest(*this, CT, MaxIters, Tol, Matrix_in);

  return(Condest_);
}

//==============================================================================
std::ostream& Ifpack2_Amesos::Print(std::ostream& os) const
{
  if (!Comm().MyPID()) {
    os << endl;
    os << "================================================================================" << endl;
    os << "Ifpack2_Amesos: " << Label () << endl << endl;
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
