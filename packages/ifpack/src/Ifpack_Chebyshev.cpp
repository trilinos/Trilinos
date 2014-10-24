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
#include "Ifpack_Chebyshev.h"
#include "Ifpack_Utils.h"
#include "Ifpack_Condest.h"
#ifdef HAVE_IFPACK_AZTECOO
#include "Ifpack_DiagPreconditioner.h"
#include "AztecOO.h"
#endif

#ifdef HAVE_IFPACK_EPETRAEXT
#include "Epetra_CrsMatrix.h"
#include "EpetraExt_PointToBlockDiagPermute.h"
#endif


#define ABS(x) ((x)>0?(x):-(x))

// Helper function for normal equations
inline void Apply_Transpose(Teuchos::RCP<const Epetra_Operator> Operator_,const Epetra_MultiVector &X,Epetra_MultiVector &Y){
  Epetra_Operator * Operator=const_cast<Epetra_Operator*>(&*Operator_);
  Operator->SetUseTranspose(true);
  Operator->Apply(X,Y);
  Operator->SetUseTranspose(false);
}




//==============================================================================
// NOTE: any change to the default values should be committed to the other
//       constructor as well.
Ifpack_Chebyshev::
Ifpack_Chebyshev(const Epetra_Operator* Operator) :
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
  PolyDegree_(1),
  UseTranspose_(false),
  Condest_(-1.0),
  ComputeCondest_(false),
  EigRatio_(30.0),
  Label_(),
  LambdaMin_(0.0),
  LambdaMax_(-1.0),
  MinDiagonalValue_(0.0),
  NumMyRows_(0),
  NumMyNonzeros_(0),
  NumGlobalRows_(0),
  NumGlobalNonzeros_(0),
  Operator_(Teuchos::rcp(Operator,false)),
  UseBlockMode_(false),
  SolveNormalEquations_(false),
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
Ifpack_Chebyshev::
Ifpack_Chebyshev(const Epetra_RowMatrix* Operator) :
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
  PolyDegree_(1),
  UseTranspose_(false),
  Condest_(-1.0),
  ComputeCondest_(false),
  EigRatio_(30.0),
  EigMaxIters_(10),
  Label_(),
  LambdaMin_(0.0),
  LambdaMax_(-1.0),
  MinDiagonalValue_(0.0),
  NumMyRows_(0),
  NumMyNonzeros_(0),
  NumGlobalRows_(0),
  NumGlobalNonzeros_(0),
  Operator_(Teuchos::rcp(Operator,false)),
  Matrix_(Teuchos::rcp(Operator,false)),
  UseBlockMode_(false),
  SolveNormalEquations_(false),
  IsRowMatrix_(true),
  ZeroStartingSolution_(true)
{
}

//==============================================================================
int Ifpack_Chebyshev::SetParameters(Teuchos::ParameterList& List)
{

  EigRatio_             = List.get("chebyshev: ratio eigenvalue", EigRatio_);
  LambdaMin_            = List.get("chebyshev: min eigenvalue", LambdaMin_);
  LambdaMax_            = List.get("chebyshev: max eigenvalue", LambdaMax_);
  PolyDegree_           = List.get("chebyshev: degree",PolyDegree_);
  MinDiagonalValue_     = List.get("chebyshev: min diagonal value",
                                   MinDiagonalValue_);
  ZeroStartingSolution_ = List.get("chebyshev: zero starting solution",
                                   ZeroStartingSolution_);

  Epetra_Vector* ID     = List.get("chebyshev: operator inv diagonal",
                                   (Epetra_Vector*)0);
  EigMaxIters_          = List.get("chebyshev: eigenvalue max iterations",EigMaxIters_);

#ifdef HAVE_IFPACK_EPETRAEXT
  // This is *always* false if EpetraExt isn't enabled
  UseBlockMode_         = List.get("chebyshev: use block mode",UseBlockMode_);
  if(!List.isParameter("chebyshev: block list")) UseBlockMode_=false;
  else{
    BlockList_          = List.get("chebyshev: block list",BlockList_);

    // Since we know we're doing a matrix inverse, clobber the block list
    // w/"invert" if it's set to multiply
    Teuchos::ParameterList Blist;
    Blist=BlockList_.get("blockdiagmatrix: list",Blist);
    string dummy("invert");
    string ApplyMode=Blist.get("apply mode",dummy);
    if(ApplyMode==string("multiply")){
      Blist.set("apply mode","invert");
      BlockList_.set("blockdiagmatrix: list",Blist);
    }
  }
#endif

  SolveNormalEquations_ = List.get("chebyshev: solve normal equations",SolveNormalEquations_);

  if (ID != 0)
  {
    InvDiagonal_ = Teuchos::rcp( new Epetra_Vector(*ID) );
  }

  SetLabel();

  return(0);
}

//==============================================================================
const Epetra_Comm& Ifpack_Chebyshev::Comm() const
{
  return(Operator_->Comm());
}

//==============================================================================
const Epetra_Map& Ifpack_Chebyshev::OperatorDomainMap() const
{
  return(Operator_->OperatorDomainMap());
}

//==============================================================================
const Epetra_Map& Ifpack_Chebyshev::OperatorRangeMap() const
{
  return(Operator_->OperatorRangeMap());
}

//==============================================================================
int Ifpack_Chebyshev::
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
int Ifpack_Chebyshev::Initialize()
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
int Ifpack_Chebyshev::Compute()
{
  if (!IsInitialized())
    IFPACK_CHK_ERR(Initialize());

  Time_->ResetStartTime();

  // reset values
  IsComputed_ = false;
  Condest_ = -1.0;

  if (PolyDegree_ <= 0)
    IFPACK_CHK_ERR(-2); // at least one application

#ifdef HAVE_IFPACK_EPETRAEXT
  // Check to see if we can run in block mode
  if (IsRowMatrix_ && InvDiagonal_ == Teuchos::null && UseBlockMode_){
    const Epetra_CrsMatrix *CrsMatrix = dynamic_cast<const Epetra_CrsMatrix*>(&*Matrix_);

    // If we don't have CrsMatrix, we can't use the block preconditioner
    if (!CrsMatrix) {
      UseBlockMode_ = false;

    } else {
      int ierr;
      InvBlockDiagonal_ = Teuchos::rcp(new EpetraExt_PointToBlockDiagPermute(*CrsMatrix));
      if (InvBlockDiagonal_ == Teuchos::null)
        IFPACK_CHK_ERR(-6);

      ierr = InvBlockDiagonal_->SetParameters(BlockList_);
      if (ierr)
        IFPACK_CHK_ERR(-7);

      ierr = InvBlockDiagonal_->Compute();
      if (ierr)
        IFPACK_CHK_ERR(-8);
    }

    // Automatically Compute Eigenvalues
    double lambda_max = 0;
    PowerMethod(EigMaxIters_, lambda_max);
    LambdaMax_ = lambda_max;

    // Test for Exact Preconditioned case
    if (ABS(LambdaMax_-1) < 1e-6)
      LambdaMax_ = LambdaMin_ = 1.0;
    else
      LambdaMin_ = LambdaMax_/EigRatio_;
  }
#endif

  if (IsRowMatrix_ && InvDiagonal_ == Teuchos::null && !UseBlockMode_) {
    InvDiagonal_ = Teuchos::rcp(new Epetra_Vector(Matrix().Map()));

    if (InvDiagonal_ == Teuchos::null)
      IFPACK_CHK_ERR(-5);

    IFPACK_CHK_ERR(Matrix().ExtractDiagonalCopy(*InvDiagonal_));

    // Inverse diagonal elements
    // Replace zeros with 1.0
    for (int i = 0 ; i < NumMyRows_ ; ++i) {
      double diag = (*InvDiagonal_)[i];
      if (IFPACK_ABS(diag) < MinDiagonalValue_)
        (*InvDiagonal_)[i] = MinDiagonalValue_;
      else
        (*InvDiagonal_)[i] = 1.0 / diag;
    }
    // Automatically compute maximum eigenvalue estimate of D^{-1}A if user hasn't provided one
    double lambda_max=0;
    if (LambdaMax_ == -1) {
      PowerMethod(Matrix(), *InvDiagonal_, EigMaxIters_, lambda_max);
      LambdaMax_=lambda_max;
      // Test for Exact Preconditioned case
      if (ABS(LambdaMax_-1) < 1e-6) LambdaMax_=LambdaMin_=1.0;
      else                          LambdaMin_=LambdaMax_/EigRatio_;
    }
    // otherwise the inverse of the diagonal has been given by the user
  }
#ifdef IFPACK_FLOPCOUNTERS
  ComputeFlops_ += NumMyRows_;
#endif

  ++NumCompute_;
  ComputeTime_ += Time_->ElapsedTime();
  IsComputed_ = true;

  SetLabel();

  return(0);
}

//==============================================================================
ostream& Ifpack_Chebyshev::Print(ostream & os) const
{

  double MyMinVal, MyMaxVal;
  double MinVal, MaxVal;

  if (IsComputed_) {
    InvDiagonal_->MinValue(&MyMinVal);
    InvDiagonal_->MaxValue(&MyMaxVal);
    Comm().MinAll(&MyMinVal,&MinVal,1);
    Comm().MinAll(&MyMaxVal,&MaxVal,1);
  }

  if (!Comm().MyPID()) {
    os << endl;
    os << "================================================================================" << endl;
    os << "Ifpack_Chebyshev" << endl;
    os << "Degree of polynomial      = " << PolyDegree_ << endl;
    os << "Condition number estimate = " << Condest() << endl;
    os << "Global number of rows     = " << Operator_->OperatorRangeMap().NumGlobalElements64() << endl;
    if (IsComputed_) {
      os << "Minimum value on stored inverse diagonal = " << MinVal << endl;
      os << "Maximum value on stored inverse diagonal = " << MaxVal << endl;
    }
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
double Ifpack_Chebyshev::
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
void Ifpack_Chebyshev::SetLabel()
{
  std::ostringstream oss;
  oss << "\"Ifpack Chebyshev polynomial\": {"
      << "Initialized: " << (IsInitialized() ? "true" : "false")
      << ", Computed: " << (IsComputed() ? "true" : "false")
      << ", degree: " << PolyDegree_
      << ", lambdaMax: " << LambdaMax_
      << ", alpha: "  << EigRatio_
      << ", lambdaMin: " << LambdaMin_
      << "}";
  Label_ = oss.str();
}

//==============================================================================
int Ifpack_Chebyshev::
ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  if (!IsComputed())
    IFPACK_CHK_ERR(-3);

  if (PolyDegree_ == 0)
    return 0;

  int nVec = X.NumVectors();
  int len  = X.MyLength();
  if (nVec != Y.NumVectors())
    IFPACK_CHK_ERR(-2);

  Time_->ResetStartTime();

  // AztecOO gives X and Y pointing to the same memory location,
  // need to create an auxiliary vector, Xcopy
  Teuchos::RefCountPtr<const Epetra_MultiVector> Xcopy;
  if (X.Pointers()[0] == Y.Pointers()[0])
    Xcopy = Teuchos::rcp( new Epetra_MultiVector(X) );
  else
    Xcopy = Teuchos::rcp( &X, false );

  double **xPtr = 0, **yPtr = 0;
  Xcopy->ExtractView(&xPtr);
  Y.ExtractView(&yPtr);

#ifdef HAVE_IFPACK_EPETRAEXT
  EpetraExt_PointToBlockDiagPermute* IBD=0;
  if (UseBlockMode_)
    IBD = &*InvBlockDiagonal_;
#endif

  //--- Do a quick solve when the matrix is identity
  double *invDiag = 0;
  if (!UseBlockMode_)
    invDiag = InvDiagonal_->Values();

  if ((LambdaMin_ == 1.0) && (LambdaMax_ == LambdaMin_)) {
#ifdef HAVE_IFPACK_EPETRAEXT
    if (UseBlockMode_)
      IBD->ApplyInverse(*Xcopy, Y);
    else
#endif
    if (nVec == 1) {
      double *yPointer = yPtr[0], *xPointer = xPtr[0];
      for (int i = 0; i < len; ++i)
        yPointer[i] = xPointer[i] * invDiag[i];

    } else {
      for (int i = 0; i < len; ++i) {
        double coeff = invDiag[i];
        for (int k = 0; k < nVec; ++k)
          yPtr[k][i] = xPtr[k][i] * coeff;
      }
    } // if (nVec == 1)

    return 0;
  } // if ((LambdaMin_ == 1.0) && (LambdaMax_ == LambdaMin_))

  //--- Initialize coefficients
  // Note that delta stores the inverse of ML_Cheby::delta
  double alpha = LambdaMax_ / EigRatio_;
  double beta  = 1.1 * LambdaMax_;
  double delta = 2.0 / (beta - alpha);
  double theta = 0.5 * (beta + alpha);
  double s1    = theta * delta;

  // Temporary vectors
  // In ML_Cheby, V corresponds to pAux and W to dk
  // NOTE: we would like to move the construction to the Compute()
  // call, but that is not possible as we don't know how many
  // vectors are in the multivector
  bool               zeroOut = false;
  Epetra_MultiVector V(X.Map(), X.NumVectors(), zeroOut);
  Epetra_MultiVector W(X.Map(), X.NumVectors(), zeroOut);
#ifdef HAVE_IFPACK_EPETRAEXT
  Epetra_MultiVector Temp(X.Map(), X.NumVectors(), zeroOut);
#endif

  double *vPointer = V.Values(), *wPointer = W.Values();

  double oneOverTheta = 1.0/theta;

  //--- If solving normal equations, multiply RHS by A^T
  if (SolveNormalEquations_) {
    Apply_Transpose(Operator_, Y, V);
    Y = V;
  }

  // Do the smoothing when block scaling is turned OFF
  // --- Treat the initial guess
  if (ZeroStartingSolution_ == false) {
    Operator_->Apply(Y, V);

    // Compute W = invDiag * ( X - V )/ Theta
#ifdef HAVE_IFPACK_EPETRAEXT
    if (UseBlockMode_) {
      Temp.Update(oneOverTheta, X, -oneOverTheta, V, 0.0);
      IBD->ApplyInverse(Temp, W);

      // Perform additional matvecs for normal equations
      // CMS: Testing this only in block mode FOR NOW
      if (SolveNormalEquations_){
        IBD->ApplyInverse(W, Temp);
        Apply_Transpose(Operator_, Temp, W);
      }
    }
    else
#endif
    if (nVec == 1) {
      double *xPointer = xPtr[0];
      for (int i = 0; i < len; ++i)
        wPointer[i] = invDiag[i] * (xPointer[i] - vPointer[i]) * oneOverTheta;

    } else {
      for (int i = 0; i < len; ++i) {
        double coeff = invDiag[i]*oneOverTheta;
        double *wi = wPointer + i, *vi = vPointer + i;
        for (int k = 0; k < nVec; ++k) {
          *wi = (xPtr[k][i] - (*vi)) * coeff;
          wi = wi + len; vi = vi + len;
        }
      }
    } // if (nVec == 1)
    // Update the vector Y
    Y.Update(1.0, W, 1.0);

  } else { // if (ZeroStartingSolution_ == false)
    // Compute W = invDiag * X / Theta
#ifdef HAVE_IFPACK_EPETRAEXT
    if (UseBlockMode_) {
      IBD->ApplyInverse(X, W);

      // Perform additional matvecs for normal equations
      // CMS: Testing this only in block mode FOR NOW
      if (SolveNormalEquations_) {
        IBD->ApplyInverse(W, Temp);
        Apply_Transpose(Operator_, Temp, W);
      }

      W.Scale(oneOverTheta);
      Y.Update(1.0, W, 0.0);
    }
    else
#endif
    if (nVec == 1) {
      double *xPointer = xPtr[0];
      for (int i = 0; i < len; ++i)
        wPointer[i] = invDiag[i] * xPointer[i] * oneOverTheta;

      memcpy(yPtr[0], wPointer, len*sizeof(double));

    } else {
      for (int i = 0; i < len; ++i) {
        double coeff = invDiag[i]*oneOverTheta;
        double *wi = wPointer + i;
        for (int k = 0; k < nVec; ++k) {
          *wi = xPtr[k][i] * coeff;
          wi = wi + len;
        }
      }

      for (int k = 0; k < nVec; ++k)
        memcpy(yPtr[k], wPointer + k*len, len*sizeof(double));
    } // if (nVec == 1)
  } // if (ZeroStartingSolution_ == false)

  //--- Apply the polynomial
  double rhok = 1.0/s1, rhokp1;
  double dtemp1, dtemp2;
  int degreeMinusOne = PolyDegree_ - 1;
  if (nVec == 1) {
    double *xPointer = xPtr[0];
    for (int k = 0; k < degreeMinusOne; ++k) {
      Operator_->Apply(Y, V);
      rhokp1 = 1.0 / (2.0*s1 - rhok);
      dtemp1 = rhokp1 * rhok;
      dtemp2 = 2.0 * rhokp1 * delta;
      rhok   = rhokp1;

      // Compute W = dtemp1 * W
      W.Scale(dtemp1);

      // Compute W = W + dtemp2 * invDiag * ( X - V )
#ifdef HAVE_IFPACK_EPETRAEXT
      if (UseBlockMode_) {
        //NTS: We can clobber V since it will be reset in the Apply
        V.Update(dtemp2, X, -dtemp2);
        IBD->ApplyInverse(V, Temp);

        // Perform additional matvecs for normal equations
        // CMS: Testing this only in block mode FOR NOW
        if (SolveNormalEquations_) {
          IBD->ApplyInverse(V, Temp);
          Apply_Transpose(Operator_, Temp, V);
        }

        W.Update(1.0, Temp, 1.0);
      }
      else
#endif
      for (int i = 0; i < len; ++i)
        wPointer[i] += dtemp2* invDiag[i] * (xPointer[i] - vPointer[i]);

      // Update the vector Y
      Y.Update(1.0, W, 1.0);
    } // for (k = 0; k < degreeMinusOne; ++k)

  } else { // if (nVec == 1) {
    for (int k = 0; k < degreeMinusOne; ++k) {
      Operator_->Apply(Y, V);
      rhokp1 = 1.0 / (2.0*s1 - rhok);
      dtemp1 = rhokp1 * rhok;
      dtemp2 = 2.0 * rhokp1 * delta;
      rhok   = rhokp1;

      // Compute W = dtemp1 * W
      W.Scale(dtemp1);

      // Compute W = W + dtemp2 * invDiag * ( X - V )
#ifdef HAVE_IFPACK_EPETRAEXT
      if (UseBlockMode_) {
        // We can clobber V since it will be reset in the Apply
        V.Update(dtemp2, X, -dtemp2);
        IBD->ApplyInverse(V, Temp);

        // Perform additional matvecs for normal equations
        // CMS: Testing this only in block mode FOR NOW
        if (SolveNormalEquations_) {
          IBD->ApplyInverse(V,Temp);
          Apply_Transpose(Operator_,Temp,V);
        }

        W.Update(1.0, Temp, 1.0);
      }
      else
#endif
        for (int i = 0; i < len; ++i) {
          double coeff = invDiag[i]*dtemp2;
          double *wi = wPointer + i, *vi = vPointer + i;
          for (int j = 0; j < nVec; ++j) {
            *wi += (xPtr[j][i] - (*vi)) * coeff;
            wi = wi + len; vi = vi + len;
          }
        }

      // Update the vector Y
      Y.Update(1.0, W, 1.0);
    } // for (k = 0; k < degreeMinusOne; ++k)
  } // if (nVec == 1)


  // Flops are updated in each of the following.
  ++NumApplyInverse_;
  ApplyInverseTime_ += Time_->ElapsedTime();

  return(0);
}

//==============================================================================
int Ifpack_Chebyshev::
PowerMethod(const Epetra_Operator& Operator,
            const Epetra_Vector& InvPointDiagonal,
            const int MaximumIterations,
            double& lambda_max)
{
  // this is a simple power method
  lambda_max = 0.0;
  double RQ_top, RQ_bottom, norm;
  Epetra_Vector x(Operator.OperatorDomainMap());
  Epetra_Vector y(Operator.OperatorRangeMap());
  x.Random();
  x.Norm2(&norm);
  if (norm == 0.0) IFPACK_CHK_ERR(-1);

  x.Scale(1.0 / norm);

  for (int iter = 0; iter < MaximumIterations; ++iter)
  {
    Operator.Apply(x, y);
    IFPACK_CHK_ERR(y.Multiply(1.0, InvPointDiagonal, y, 0.0));
    IFPACK_CHK_ERR(y.Dot(x, &RQ_top));
    IFPACK_CHK_ERR(x.Dot(x, &RQ_bottom));
    lambda_max = RQ_top / RQ_bottom;
    IFPACK_CHK_ERR(y.Norm2(&norm));
    if (norm == 0.0) IFPACK_CHK_ERR(-1);
    IFPACK_CHK_ERR(x.Update(1.0 / norm, y, 0.0));
  }

  return(0);
}

//==============================================================================
int Ifpack_Chebyshev::
CG(const Epetra_Operator& Operator,
   const Epetra_Vector& InvPointDiagonal,
   const int MaximumIterations,
   double& lambda_min, double& lambda_max)
{
#ifdef HAVE_IFPACK_AZTECOO
  Epetra_Vector x(Operator.OperatorDomainMap());
  Epetra_Vector y(Operator.OperatorRangeMap());
  x.Random();
  y.PutScalar(0.0);

  Epetra_LinearProblem LP(const_cast<Epetra_Operator*>(&Operator), &x, &y);
  AztecOO solver(LP);
  solver.SetAztecOption(AZ_solver, AZ_cg_condnum);
  solver.SetAztecOption(AZ_output, AZ_none);

  Ifpack_DiagPreconditioner diag(Operator.OperatorDomainMap(),
                                 Operator.OperatorRangeMap(),
                                 InvPointDiagonal);
  solver.SetPrecOperator(&diag);
  solver.Iterate(MaximumIterations, 1e-10);

  const double* status = solver.GetAztecStatus();

  lambda_min = status[AZ_lambda_min];
  lambda_max = status[AZ_lambda_max];

  return(0);
#else
  cout << "You need to configure IFPACK with support for AztecOO" << endl;
  cout << "to use the CG estimator. This may require --enable-aztecoo" << endl;
  cout << "in your configure script." << endl;
  IFPACK_CHK_ERR(-1);
#endif
}

//==============================================================================
#ifdef HAVE_IFPACK_EPETRAEXT
int Ifpack_Chebyshev::
PowerMethod(const int MaximumIterations,  double& lambda_max)
{

  if(!UseBlockMode_) IFPACK_CHK_ERR(-1);
  // this is a simple power method
  lambda_max = 0.0;
  double RQ_top, RQ_bottom, norm;
  Epetra_Vector x(Operator_->OperatorDomainMap());
  Epetra_Vector y(Operator_->OperatorRangeMap());
  Epetra_Vector z(Operator_->OperatorRangeMap());
  x.Random();
  x.Norm2(&norm);
  if (norm == 0.0) IFPACK_CHK_ERR(-1);

  x.Scale(1.0 / norm);

  for (int iter = 0; iter < MaximumIterations; ++iter)
  {
    Operator_->Apply(x, z);
    InvBlockDiagonal_->ApplyInverse(z,y);
    if(SolveNormalEquations_){
      InvBlockDiagonal_->ApplyInverse(y,z);
      Apply_Transpose(Operator_,z, y);
    }

    IFPACK_CHK_ERR(y.Dot(x, &RQ_top));
    IFPACK_CHK_ERR(x.Dot(x, &RQ_bottom));
    lambda_max = RQ_top / RQ_bottom;
    IFPACK_CHK_ERR(y.Norm2(&norm));
    if (norm == 0.0) IFPACK_CHK_ERR(-1);
    IFPACK_CHK_ERR(x.Update(1.0 / norm, y, 0.0));
  }

  return(0);
}
#endif

//==============================================================================
#ifdef HAVE_IFPACK_EPETRAEXT
int Ifpack_Chebyshev::
CG(const int MaximumIterations,
   double& lambda_min, double& lambda_max)
{
  IFPACK_CHK_ERR(-1);// NTS: This always seems to yield errors in AztecOO, ergo,
                     // I turned it off.

  if(!UseBlockMode_) IFPACK_CHK_ERR(-1);

#ifdef HAVE_IFPACK_AZTECOO
  Epetra_Vector x(Operator_->OperatorDomainMap());
  Epetra_Vector y(Operator_->OperatorRangeMap());
  x.Random();
  y.PutScalar(0.0);
  Epetra_LinearProblem LP(const_cast<Epetra_RowMatrix*>(&*Matrix_), &x, &y);

  AztecOO solver(LP);
  solver.SetAztecOption(AZ_solver, AZ_cg_condnum);
  solver.SetAztecOption(AZ_output, AZ_none);

  solver.SetPrecOperator(&*InvBlockDiagonal_);
  solver.Iterate(MaximumIterations, 1e-10);

  const double* status = solver.GetAztecStatus();

  lambda_min = status[AZ_lambda_min];
  lambda_max = status[AZ_lambda_max];

  return(0);
#else
  cout << "You need to configure IFPACK with support for AztecOO" << endl;
  cout << "to use the CG estimator. This may require --enable-aztecoo" << endl;
  cout << "in your configure script." << endl;
  IFPACK_CHK_ERR(-1);
#endif
}
#endif
