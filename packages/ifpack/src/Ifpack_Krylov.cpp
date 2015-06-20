/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#include "Ifpack_ConfigDefs.h"
#include <iomanip>
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Time.h"
#include "Ifpack_Krylov.h"
#include "Ifpack_Utils.h"
#include "Ifpack_Condest.h"
#ifdef HAVE_IFPACK_AZTECOO
#include "AztecOO.h"
#endif

#ifdef HAVE_IFPACK_EPETRAEXT
#include "Epetra_CrsMatrix.h"
#include "EpetraExt_PointToBlockDiagPermute.h"
#endif

//==============================================================================
// NOTE: any change to the default values should be committed to the other
//       constructor as well.
Ifpack_Krylov::
Ifpack_Krylov(Epetra_Operator* Operator) :
  IsInitialized_(false),
  IsComputed_(false),
  NumInitialize_(0),
  NumCompute_(0),
  NumApplyInverse_(0),
  InitializeTime_(0.0),
  ComputeTime_(0.0),
  ApplyInverseTime_(0.0),
  ComputeFlops_(0.0),
  ApplyInverseFlops_(0.0),
  Iterations_(5),
  Tolerance_(1e-6),
  SolverType_(1),
  PreconditionerType_(1),
  NumSweeps_(0),
  BlockSize_(1),
  DampingParameter_(1.0),
  UseTranspose_(false),
  Condest_(-1.0),
  /* ComputeCondest_(false), (unused; commented out to avoid build warnings) */
  Label_(),
  NumMyRows_(0),
  NumMyNonzeros_(0),
  NumGlobalRows_(0),
  NumGlobalNonzeros_(0),
  Operator_(Teuchos::rcp(Operator,false)),
  IsRowMatrix_(false),
  ZeroStartingSolution_(true)
{
}

//==============================================================================
// NOTE: This constructor has been introduced because SWIG does not appear
//       to appreciate dynamic_cast. An instruction of type
//       Matrix_ = dynamic_cast<const Epetra_RowMatrix*> in the
//       other construction does not work in PyTrilinos -- of course
//       it does in any C++ code (for an Epetra_Operator that is also
//       an Epetra_RowMatrix).
//
// FIXME: move declarations into a separate method?
Ifpack_Krylov::
Ifpack_Krylov(Epetra_RowMatrix* Operator) :
  IsInitialized_(false),
  IsComputed_(false),
  NumInitialize_(0),
  NumCompute_(0),
  NumApplyInverse_(0),
  InitializeTime_(0.0),
  ComputeTime_(0.0),
  ApplyInverseTime_(0.0),
  ComputeFlops_(0.0),
  ApplyInverseFlops_(0.0),
  Iterations_(5),
  Tolerance_(1e-6),
  SolverType_(1),
  PreconditionerType_(1),
  NumSweeps_(0),
  BlockSize_(1),
  DampingParameter_(1.0),
  UseTranspose_(false),
  Condest_(-1.0),
  /* ComputeCondest_(false), (unused; commented out to avoid build warnings) */
  Label_(),
  NumMyRows_(0),
  NumMyNonzeros_(0),
  NumGlobalRows_(0),
  NumGlobalNonzeros_(0),
  Operator_(Teuchos::rcp(Operator,false)),
  Matrix_(Teuchos::rcp(Operator,false)),
  IsRowMatrix_(true),
  ZeroStartingSolution_(true)
{
}

//==============================================================================
int Ifpack_Krylov::SetParameters(Teuchos::ParameterList& List)
{
  Iterations_           = List.get("krylov: iterations",Iterations_);
  Tolerance_            = List.get("krylov: tolerance",Tolerance_);
  SolverType_           = List.get("krylov: solver",SolverType_);
  PreconditionerType_   = List.get("krylov: preconditioner",PreconditionerType_);
  NumSweeps_            = List.get("krylov: number of sweeps",NumSweeps_);
  BlockSize_            = List.get("krylov: block size",BlockSize_);
  DampingParameter_     = List.get("krylov: damping parameter",DampingParameter_);
  ZeroStartingSolution_ = List.get("krylov: zero starting solution",ZeroStartingSolution_);
  SetLabel();
  return(0);
}

//==============================================================================
const Epetra_Comm& Ifpack_Krylov::Comm() const
{
  return(Operator_->Comm());
}

//==============================================================================
const Epetra_Map& Ifpack_Krylov::OperatorDomainMap() const
{
  return(Operator_->OperatorDomainMap());
}

//==============================================================================
const Epetra_Map& Ifpack_Krylov::OperatorRangeMap() const
{
  return(Operator_->OperatorRangeMap());
}

//==============================================================================
int Ifpack_Krylov::
Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  if (IsComputed() == false)
    IFPACK_CHK_ERR(-3);

  if (X.NumVectors() != Y.NumVectors())
    IFPACK_CHK_ERR(-2);

  if (IsRowMatrix_)
  {
    IFPACK_CHK_ERR(Matrix_->Multiply(UseTranspose(),X,Y));
  }
  else
  {
    IFPACK_CHK_ERR(Operator_->Apply(X,Y));
  }

  return(0);
}

//==============================================================================
int Ifpack_Krylov::Initialize()
{
  IsInitialized_ = false;

  if (Operator_ == Teuchos::null)
    IFPACK_CHK_ERR(-2);

  if (Time_ == Teuchos::null)
    Time_ = Teuchos::rcp( new Epetra_Time(Comm()) );

  if (IsRowMatrix_)
  {
    if (Matrix().NumGlobalRows64() != Matrix().NumGlobalCols64())
      IFPACK_CHK_ERR(-2); // only square matrices

    NumMyRows_ = Matrix_->NumMyRows();
    NumMyNonzeros_ = Matrix_->NumMyNonzeros();
    NumGlobalRows_ = Matrix_->NumGlobalRows64();
    NumGlobalNonzeros_ = Matrix_->NumGlobalNonzeros64();
  }
  else
  {
    if (Operator_->OperatorDomainMap().NumGlobalElements64() !=
        Operator_->OperatorRangeMap().NumGlobalElements64())
      IFPACK_CHK_ERR(-2); // only square operators
  }

  ++NumInitialize_;
  InitializeTime_ += Time_->ElapsedTime();
  IsInitialized_ = true;
  return(0);
}

//==============================================================================
int Ifpack_Krylov::Compute()
{
  if (!IsInitialized())
    IFPACK_CHK_ERR(Initialize());

  Time_->ResetStartTime();

#ifdef HAVE_IFPACK_AZTECOO
  // setup Aztec solver
  AztecSolver_ = Teuchos::rcp( new AztecOO );
  if(IsRowMatrix_==true) {
    AztecSolver_ -> SetUserMatrix(&*Matrix_);
  }
  else {
    AztecSolver_ -> SetUserOperator(&*Operator_);
  }
  if(SolverType_==0) {
    AztecSolver_ -> SetAztecOption(AZ_solver, AZ_cg);
  }
  else {
    AztecSolver_ -> SetAztecOption(AZ_solver, AZ_gmres);
  }
  AztecSolver_ -> SetAztecOption(AZ_output, AZ_none);
  // setup preconditioner
  Teuchos::ParameterList List;
  List.set("relaxation: damping factor", DampingParameter_);
  List.set("relaxation: sweeps",NumSweeps_);
  if(PreconditionerType_==0)      { }
  else if(PreconditionerType_==1) { List.set("relaxation: type", "Jacobi"                ); }
  else if(PreconditionerType_==2) { List.set("relaxation: type", "Gauss-Seidel"          ); }
  else if(PreconditionerType_==3) { List.set("relaxation: type", "symmetric Gauss-Seidel"); }
  if(BlockSize_==1) {
    IfpackPrec_ = Teuchos::rcp( new Ifpack_PointRelaxation(&*Matrix_) );
  }
  else {
    IfpackPrec_ = Teuchos::rcp( new Ifpack_BlockRelaxation< Ifpack_DenseContainer > (&*Matrix_) );
    int NumRows;
    if(IsRowMatrix_==true) {
      NumRows = Matrix_->NumMyRows();
    }
    else {
      long long NumRows_LL = Operator_->OperatorDomainMap().NumGlobalElements64();
      if(NumRows_LL > std::numeric_limits<int>::max())
        throw "Ifpack_Krylov::Compute: NumGlobalElements don't fit an int";
      else
        NumRows = static_cast<int>(NumRows_LL);
    }
    List.set("partitioner: type", "linear");
    List.set("partitioner: local parts", NumRows/BlockSize_);
  }
  if(PreconditionerType_>0) {
    IfpackPrec_ -> SetParameters(List);
    IfpackPrec_ -> Initialize();
    IfpackPrec_ -> Compute();
    AztecSolver_ -> SetPrecOperator(&*IfpackPrec_);
  }
#else
  cout << "You need to configure IFPACK with support for AztecOO" << endl;
  cout << "to use this preconditioner. This may require --enable-aztecoo" << endl;
  cout << "in your configure script." << endl;
  IFPACK_CHK_ERR(-1);
#endif

  // reset values
  IsComputed_ = false;
  Condest_ = -1.0;
  ++NumCompute_;
  ComputeTime_ += Time_->ElapsedTime();
  IsComputed_ = true;

  return(0);
}

//==============================================================================
ostream& Ifpack_Krylov::Print(ostream & os) const
{

  if (!Comm().MyPID()) {
    os << endl;
    os << "================================================================================" << endl;
    os << "Ifpack_Krylov" << endl;
    os << "Number of iterations                 = " << Iterations_ << endl;
    os << "Residual Tolerance                   = " << Tolerance_ << endl;
    os << "Solver type (O for CG, 1 for GMRES)  = " << SolverType_ << endl;
    os << "Preconditioner type                  = " << PreconditionerType_ << endl;
    os << "(0 for none, 1 for Jacobi, 2 for GS, 3 for SGS )" << endl;
    os << "Condition number estimate            = " << Condest() << endl;
    os << "Global number of rows                = " << Operator_->OperatorRangeMap().NumGlobalElements64() << endl;
    if (ZeroStartingSolution_)
      os << "Using zero starting solution" << endl;
    else
      os << "Using input starting solution" << endl;
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

//==============================================================================
double Ifpack_Krylov::
Condest(const Ifpack_CondestType CT,
        const int MaxIters, const double Tol,
        Epetra_RowMatrix* Matrix_in)
{
  if (!IsComputed()) // cannot compute right now
    return(-1.0);

  // always computes it. Call Condest() with no parameters to get
  // the previous estimate.
  Condest_ = Ifpack_Condest(*this, CT, MaxIters, Tol, Matrix_in);

  return(Condest_);
}

//==============================================================================
void Ifpack_Krylov::SetLabel()
{
  Label_ = "IFPACK (Krylov smoother), iterations=" + Ifpack_toString(Iterations_);
}

//==============================================================================
int Ifpack_Krylov::
ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{

  if (!IsComputed())
    IFPACK_CHK_ERR(-3);

  if (Iterations_ == 0)
    return 0;

  int nVec = X.NumVectors();
  if (nVec != Y.NumVectors())
    IFPACK_CHK_ERR(-2);

  Time_->ResetStartTime();

  // AztecOO gives X and Y pointing to the same memory location,
  // need to create an auxiliary vector, Xcopy
  Teuchos::RCP<Epetra_MultiVector> Xcopy = Teuchos::rcp( new Epetra_MultiVector(X) );
  if(ZeroStartingSolution_==true) {
    Y.PutScalar(0.0);
  }

#ifdef HAVE_IFPACK_AZTECOO
  AztecSolver_ -> SetLHS(&Y);
  AztecSolver_ -> SetRHS(&*Xcopy);
  AztecSolver_ -> Iterate(Iterations_,Tolerance_);
#else
  cout << "You need to configure IFPACK with support for AztecOO" << endl;
  cout << "to use this preconditioner. This may require --enable-aztecoo" << endl;
  cout << "in your configure script." << endl;
  IFPACK_CHK_ERR(-1);
#endif

  // Flops are updated in each of the following.
  ++NumApplyInverse_;
  ApplyInverseTime_ += Time_->ElapsedTime();
  return(0);
}
