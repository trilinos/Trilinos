/*@HEADER
// ***********************************************************************
//
//       Tifpack: Tempated Object-Oriented Algebraic Preconditioner Package
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

#include "Tifpack_ConfigDefs.hpp"
#include <iomanip>
#include "Tpetra_Operator.hpp"
#include "Tpetra_RowMatrix.hpp"
#include "Tpetra_Comm.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_Time.hpp"
#include "Tifpack_Chebyshev.hpp"
#include "Tifpack_Utils.hpp"
#include "Tifpack_Condest.hpp"
#ifdef HAVE_TIFPACK_AZTECOO
#include "Tifpack_DiagPreconditioner.hpp"
#include "AztecOO.hpp"
#endif

//==============================================================================
// NOTE: any change to the default values should be committed to the other
//       constructor as well.
Tifpack_Chebyshev::
Tifpack_Chebyshev(const Tpetra_Operator* Operator) :
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
  LambdaMax_(100.0),
  MinDiagonalValue_(0.0),
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
//       Matrix_ = dynamic_cast<const Tpetra_RowMatrix*> in the
//       other construction does not work in PyTrilinos -- of course
//       it does in any C++ code (for an Tpetra_Operator that is also
//       an Tpetra_RowMatrix).
//
// FIXME: move declarations into a separate method?
Tifpack_Chebyshev::
Tifpack_Chebyshev(const Tpetra_RowMatrix* Operator) :
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
  LambdaMax_(100.0),
  MinDiagonalValue_(0.0),
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
int Tifpack_Chebyshev::SetParameters(Teuchos::ParameterList& List)
{

  EigRatio_             = List.get("chebyshev: ratio eigenvalue", EigRatio_);
  LambdaMin_            = List.get("chebyshev: min eigenvalue", LambdaMin_);
  LambdaMax_            = List.get("chebyshev: max eigenvalue", LambdaMax_);
  PolyDegree_           = List.get("chebyshev: degree",PolyDegree_);
  MinDiagonalValue_     = List.get("chebyshev: min diagonal value", 
                                   MinDiagonalValue_);
  ZeroStartingSolution_ = List.get("chebyshev: zero starting solution", 
                                   ZeroStartingSolution_);

  Tpetra_Vector* ID     = List.get("chebyshev: operator inv diagonal", 
                                   (Tpetra_Vector*)0);

  if (ID != 0) 
  {
    InvDiagonal_ = Teuchos::rcp( new Tpetra_Vector(*ID) );
  }

  SetLabel();

  return(0);
}

//==============================================================================
const Tpetra_Comm& Tifpack_Chebyshev::Comm() const
{
  return(Operator_->Comm());
}

//==============================================================================
const Tpetra_Map& Tifpack_Chebyshev::OperatorDomainMap() const
{
  return(Operator_->OperatorDomainMap());
}

//==============================================================================
const Tpetra_Map& Tifpack_Chebyshev::OperatorRangeMap() const
{
  return(Operator_->OperatorRangeMap());
}

//==============================================================================
int Tifpack_Chebyshev::
Apply(const Tpetra_MultiVector& X, Tpetra_MultiVector& Y) const
{
  if (IsComputed() == false)
    TIFPACK_CHK_ERR(-3);

  if (X.NumVectors() != Y.NumVectors())
    TIFPACK_CHK_ERR(-2);

  if (IsRowMatrix_)
  {
    TIFPACK_CHK_ERR(Matrix_->Multiply(UseTranspose(),X,Y));
  }
  else
  {
    TIFPACK_CHK_ERR(Operator_->Apply(X,Y));
  }

  return(0);
}

//==============================================================================
int Tifpack_Chebyshev::Initialize()
{
  IsInitialized_ = false;

  if (Operator_ == Teuchos::null)
    TIFPACK_CHK_ERR(-2);

  if (Time_ == Teuchos::null)
    Time_ = Teuchos::rcp( new Tpetra_Time(Comm()) );

  if (IsRowMatrix_)
  {
    if (Matrix().NumGlobalRows() != Matrix().NumGlobalCols())
      TIFPACK_CHK_ERR(-2); // only square matrices

    NumMyRows_ = Matrix_->NumMyRows();
    NumMyNonzeros_ = Matrix_->NumMyNonzeros();
    NumGlobalRows_ = Matrix_->NumGlobalRows();
    NumGlobalNonzeros_ = Matrix_->NumGlobalNonzeros();
  }
  else
  {
    if (Operator_->OperatorDomainMap().NumGlobalElements() !=       
        Operator_->OperatorRangeMap().NumGlobalElements())
      TIFPACK_CHK_ERR(-2); // only square operators
  }

  ++NumInitialize_;
  InitializeTime_ += Time_->ElapsedTime();
  IsInitialized_ = true;
  return(0);
}

//==============================================================================
int Tifpack_Chebyshev::Compute()
{
  if (!IsInitialized())
    TIFPACK_CHK_ERR(Initialize());

  Time_->ResetStartTime();

  // reset values
  IsComputed_ = false;
  Condest_ = -1.0;

  if (PolyDegree_ <= 0)
    TIFPACK_CHK_ERR(-2); // at least one application
  
  if (IsRowMatrix_ && InvDiagonal_ == Teuchos::null)
  {
    InvDiagonal_ = Teuchos::rcp( new Tpetra_Vector(Matrix().Map()) );

    if (InvDiagonal_ == Teuchos::null)
      TIFPACK_CHK_ERR(-5);

    TIFPACK_CHK_ERR(Matrix().ExtractDiagonalCopy(*InvDiagonal_));

    // Inverse diagonal elements
    // Replace zeros with 1.0
    for (int i = 0 ; i < NumMyRows_ ; ++i) {
      double diag = (*InvDiagonal_)[i];
      if (TIFPACK_ABS(diag) < MinDiagonalValue_)
        (*InvDiagonal_)[i] = MinDiagonalValue_;
      else
        (*InvDiagonal_)[i] = 1.0 / diag;
    }
  }
  // otherwise the inverse of the diagonal has been given by the user

  ComputeFlops_ += NumMyRows_;

  ++NumCompute_;
  ComputeTime_ += Time_->ElapsedTime();
  IsComputed_ = true;

  return(0);
}

//==============================================================================
ostream& Tifpack_Chebyshev::Print(ostream & os) const
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
    os << "Tifpack_Chebyshev" << endl;
    os << "Degree of polynomial      = " << PolyDegree_ << endl;
    os << "Condition number estimate = " << Condest() << endl;
    os << "Global number of rows     = " << Operator_->OperatorRangeMap().NumGlobalElements() << endl;
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
double Tifpack_Chebyshev::
Condest(const Tifpack_CondestType CT, 
        const int MaxIters, const double Tol,
	Tpetra_RowMatrix* Matrix_in)
{
  if (!IsComputed()) // cannot compute right now
    return(-1.0);

  // always computes it. Call Condest() with no parameters to get
  // the previous estimate.
  Condest_ = Tifpack_Condest(*this, CT, MaxIters, Tol, Matrix_in);

  return(Condest_);
}

//==============================================================================
void Tifpack_Chebyshev::SetLabel()
{
  Label_ = "TIFPACK (Chebyshev polynomial), degree=" + Tifpack_toString(PolyDegree_);
}

//==============================================================================
int Tifpack_Chebyshev::
ApplyInverse(const Tpetra_MultiVector& X, Tpetra_MultiVector& Y) const
{
  if (!IsComputed())
    TIFPACK_CHK_ERR(-3);

  if (PolyDegree_ == 0)
    return 0;

  int nVec = X.NumVectors();
  int len = X.MyLength();
  if (nVec != Y.NumVectors())
    TIFPACK_CHK_ERR(-2);

  Time_->ResetStartTime();

  // AztecOO gives X and Y pointing to the same memory location,
  // need to create an auxiliary vector, Xcopy
  Teuchos::RefCountPtr<const Tpetra_MultiVector> Xcopy;
  if (X.Pointers()[0] == Y.Pointers()[0])
    Xcopy = Teuchos::rcp( new Tpetra_MultiVector(X) );
  else
    Xcopy = Teuchos::rcp( &X, false );

  double **xPtr = 0, **yPtr = 0;
  Xcopy->ExtractView(&xPtr);
  Y.ExtractView(&yPtr);

  //--- Do a quick solve when the matrix is identity
  double *invDiag = InvDiagonal_->Values();
  if ((LambdaMin_ == 1.0) && (LambdaMax_ == LambdaMin_)) {
    if (nVec == 1) {
      double *yPointer = yPtr[0], *xPointer = xPtr[0];
      for (int i = 0; i < len; ++i)
        yPointer[i] = xPointer[i]*invDiag[i];
    }
    else {
      int i, k;
      for (i = 0; i < len; ++i) {
        double coeff = invDiag[i];
        for (k = 0; k < nVec; ++k)
          yPtr[k][i] = xPtr[k][i] * coeff;
      }
    } // if (nVec == 1)
    return 0;
  } // if ((LambdaMin_ == 1.0) && (LambdaMax_ == LambdaMin_))

  //--- Initialize coefficients
  // Note that delta stores the inverse of ML_Cheby::delta
  double alpha = LambdaMax_ / EigRatio_;
  double beta = 1.1 * LambdaMax_;
  double delta = 2.0 / (beta - alpha);
  double theta = 0.5 * (beta + alpha);
  double s1 = theta * delta;

  //--- Define vectors
  // In ML_Cheby, V corresponds to pAux and W to dk
  Tpetra_MultiVector V(X);
  Tpetra_MultiVector W(X);

  double *vPointer = V.Values(), *wPointer = W.Values();

  double oneOverTheta = 1.0/theta;
  int i, j, k;

  // Do the smoothing when block scaling is turned OFF
  // --- Treat the initial guess
  if (ZeroStartingSolution_ == false) {
    Operator_->Apply(Y, V);
    // Compute W = invDiag * ( X - V )/ Theta
    if (nVec == 1) {
      double *xPointer = xPtr[0];
      for (i = 0; i < len; ++i)
        wPointer[i] = invDiag[i] * (xPointer[i] - vPointer[i]) * oneOverTheta;
    }
    else {
      for (i = 0; i < len; ++i) {
        double coeff = invDiag[i]*oneOverTheta;
        double *wi = wPointer + i, *vi = vPointer + i;
        for (k = 0; k < nVec; ++k) {
          *wi = (xPtr[k][i] - (*vi)) * coeff;
          wi = wi + len; vi = vi + len;
        }
      }
    } // if (nVec == 1)
    // Update the vector Y
    Y.Update(1.0, W, 1.0);
  }
  else {
    // Compute W = invDiag * X / Theta
    if (nVec == 1) {
      double *xPointer = xPtr[0];
      for (i = 0; i < len; ++i)
        wPointer[i] = invDiag[i] * xPointer[i] * oneOverTheta;
      memcpy(yPtr[0], wPointer, len*sizeof(double));
    }
    else {
      for (i = 0; i < len; ++i) {
        double coeff = invDiag[i]*oneOverTheta;
        double *wi = wPointer + i;
        for (k = 0; k < nVec; ++k) {
          *wi = xPtr[k][i] * coeff;
          wi = wi + len;
        }
      }
      for (k = 0; k < nVec; ++k)
        memcpy(yPtr[k], wPointer + k*len, len*sizeof(double));
    } // if (nVec == 1)
  } // if (ZeroStartingSolution_ == false)

  //--- Apply the polynomial
  double rhok = 1.0/s1, rhokp1;
  double dtemp1, dtemp2;
  int degreeMinusOne = PolyDegree_ - 1;
  if (nVec == 1) {
    double *xPointer = xPtr[0];
    for (k = 0; k < degreeMinusOne; ++k) {
      Operator_->Apply(Y, V);
      rhokp1 = 1.0 / (2.0*s1 - rhok);
      dtemp1 = rhokp1 * rhok;
      dtemp2 = 2.0 * rhokp1 * delta;
      rhok = rhokp1;
      // Compute W = dtemp1 * W
      W.Scale(dtemp1);
      // Compute W = W + dtemp2 * invDiag * ( X - V )
      for (i = 0; i < len; ++i)
        wPointer[i] += dtemp2* invDiag[i] * (xPointer[i] - vPointer[i]);
      // Update the vector Y
      Y.Update(1.0, W, 1.0);
    } // for (k = 0; k < degreeMinusOne; ++k)
  }
  else {
    for (k = 0; k < degreeMinusOne; ++k) {
      Operator_->Apply(Y, V);
      rhokp1 = 1.0 / (2.0*s1 - rhok);
      dtemp1 = rhokp1 * rhok;
      dtemp2 = 2.0 * rhokp1 * delta;
      rhok = rhokp1;
      // Compute W = dtemp1 * W
      W.Scale(dtemp1);
      // Compute W = W + dtemp2 * invDiag * ( X - V )
      for (i = 0; i < len; ++i) {
        double coeff = invDiag[i]*dtemp2;
        double *wi = wPointer + i, *vi = vPointer + i;
        for (j = 0; j < nVec; ++j) {
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
int Tifpack_Chebyshev::
PowerMethod(const Tpetra_Operator& Operator, 
            const Tpetra_Vector& InvPointDiagonal, 
            const int MaximumIterations, 
            double& lambda_max)
{
  // this is a simple power method
  lambda_max = 0.0;
  double RQ_top, RQ_bottom, norm;
  Tpetra_Vector x(Operator.OperatorDomainMap());
  Tpetra_Vector y(Operator.OperatorRangeMap());
  x.Random();
  x.Norm2(&norm);
  if (norm == 0.0) TIFPACK_CHK_ERR(-1);

  x.Scale(1.0 / norm);

  for (int iter = 0; iter < MaximumIterations; ++iter)
  {
    Operator.Apply(x, y);
    TIFPACK_CHK_ERR(y.Multiply(1.0, InvPointDiagonal, y, 0.0));
    TIFPACK_CHK_ERR(y.Dot(x, &RQ_top));
    TIFPACK_CHK_ERR(x.Dot(x, &RQ_bottom));
    lambda_max = RQ_top / RQ_bottom;
    TIFPACK_CHK_ERR(y.Norm2(&norm));
    if (norm == 0.0) TIFPACK_CHK_ERR(-1);
    TIFPACK_CHK_ERR(x.Update(1.0 / norm, y, 0.0));
  }

  return(0);
}

//==============================================================================
int Tifpack_Chebyshev::
CG(const Tpetra_Operator& Operator, 
   const Tpetra_Vector& InvPointDiagonal, 
   const int MaximumIterations, 
   double& lambda_min, double& lambda_max)
{
#ifdef HAVE_TIFPACK_AZTECOO
  Tpetra_Vector x(Operator.OperatorDomainMap());
  Tpetra_Vector y(Operator.OperatorRangeMap());
  x.Random();
  y.PutScalar(0.0);

  Tpetra_LinearProblem LP(const_cast<Tpetra_Operator*>(&Operator), &x, &y);
  AztecOO solver(LP);
  solver.SetAztecOption(AZ_solver, AZ_cg_condnum);
  solver.SetAztecOption(AZ_output, AZ_none);

  Tifpack_DiagPreconditioner diag(Operator.OperatorDomainMap(),
                                 Operator.OperatorRangeMap(),
                                 InvPointDiagonal);
  solver.SetPrecOperator(&diag);
  solver.Iterate(MaximumIterations, 1e-10);

  const double* status = solver.GetAztecStatus();

  lambda_min = status[AZ_lambda_min];
  lambda_max = status[AZ_lambda_max];

  return(0);
#else
  cout << "You need to configure TIFPACK with support for AztecOO" << endl;
  cout << "to use the CG estimator. This may require --enable-aztecoo" << endl;
  cout << "in your configure script." << endl;
  TIFPACK_CHK_ERR(-1);
#endif
}

