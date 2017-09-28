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
#include <cmath>
#include "Epetra_Operator.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Ifpack_Preconditioner.h"
#include "Ifpack_PointRelaxation.h"
#include "Ifpack_Utils.h"
#include "Ifpack_Condest.h"

static const int IFPACK_JACOBI = 0;
static const int IFPACK_GS = 1;
static const int IFPACK_SGS = 2;

//==============================================================================
Ifpack_PointRelaxation::
Ifpack_PointRelaxation(const Epetra_RowMatrix* Matrix_in) :
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
  NumSweeps_(1),
  DampingFactor_(1.0),
  UseTranspose_(false),
  Condest_(-1.0),
  /* ComputeCondest_(false), (unused; commented out to avoid build warnings) */
  PrecType_(IFPACK_JACOBI),
  MinDiagonalValue_(0.0),
  NumMyRows_(0),
  NumMyNonzeros_(0),
  NumGlobalRows_(0),
  NumGlobalNonzeros_(0),
  Matrix_(Teuchos::rcp(Matrix_in,false)),
  IsParallel_(false),
  ZeroStartingSolution_(true),
  DoBackwardGS_(false),
  DoL1Method_(false),
  L1Eta_(1.5),
  NumLocalSmoothingIndices_(Matrix_in->NumMyRows()),
  LocalSmoothingIndices_(0)
{
}

//==============================================================================
int Ifpack_PointRelaxation::SetParameters(Teuchos::ParameterList& List)
{
  using std::cout;
  using std::endl;

  std::string PT;
  if (PrecType_ == IFPACK_JACOBI)
    PT = "Jacobi";
  else if (PrecType_ == IFPACK_GS)
    PT = "Gauss-Seidel";
  else if (PrecType_ == IFPACK_SGS)
    PT = "symmetric Gauss-Seidel";

  PT = List.get("relaxation: type", PT);

  if (PT == "Jacobi")
    PrecType_ = IFPACK_JACOBI;
  else if (PT == "Gauss-Seidel")
    PrecType_ = IFPACK_GS;
  else if (PT == "symmetric Gauss-Seidel")
    PrecType_ = IFPACK_SGS;
  else {
    IFPACK_CHK_ERR(-2);
  }

  NumSweeps_            = List.get("relaxation: sweeps",NumSweeps_);
  DampingFactor_        = List.get("relaxation: damping factor",
                                   DampingFactor_);
  MinDiagonalValue_     = List.get("relaxation: min diagonal value",
                                   MinDiagonalValue_);
  ZeroStartingSolution_ = List.get("relaxation: zero starting solution",
                                   ZeroStartingSolution_);

  DoBackwardGS_         = List.get("relaxation: backward mode",DoBackwardGS_);

  DoL1Method_           = List.get("relaxation: use l1",DoL1Method_);

  L1Eta_                = List.get("relaxation: l1 eta",L1Eta_);


  // This (partial) reordering would require a very different implementation of Jacobi than we have no, so we're not
  // going to do it.
  NumLocalSmoothingIndices_= List.get("relaxation: number of local smoothing indices",NumLocalSmoothingIndices_);
  LocalSmoothingIndices_   = List.get("relaxation: local smoothing indices",LocalSmoothingIndices_);
  if(PrecType_ == IFPACK_JACOBI && LocalSmoothingIndices_) {
    NumLocalSmoothingIndices_=NumMyRows_;
    LocalSmoothingIndices_=0;
    if(!Comm().MyPID()) cout<<"Ifpack_PointRelaxation: WARNING: Reordered/Local Smoothing not implemented for Jacobi.  Defaulting to regular Jacobi"<<endl;
  }

  SetLabel();

  return(0);
}

//==============================================================================
const Epetra_Comm& Ifpack_PointRelaxation::Comm() const
{
  return(Matrix_->Comm());
}

//==============================================================================
const Epetra_Map& Ifpack_PointRelaxation::OperatorDomainMap() const
{
  return(Matrix_->OperatorDomainMap());
}

//==============================================================================
const Epetra_Map& Ifpack_PointRelaxation::OperatorRangeMap() const
{
  return(Matrix_->OperatorRangeMap());
}

//==============================================================================
int Ifpack_PointRelaxation::Initialize()
{
  IsInitialized_ = false;

  if (Matrix_ == Teuchos::null)
    IFPACK_CHK_ERR(-2);

  if (Time_ == Teuchos::null)
    Time_ = Teuchos::rcp( new Epetra_Time(Comm()) );

  if (Matrix().NumGlobalRows64() != Matrix().NumGlobalCols64())
    IFPACK_CHK_ERR(-2); // only square matrices

  NumMyRows_ = Matrix_->NumMyRows();
  NumMyNonzeros_ = Matrix_->NumMyNonzeros();
  NumGlobalRows_ = Matrix_->NumGlobalRows64();
  NumGlobalNonzeros_ = Matrix_->NumGlobalNonzeros64();

  if (Comm().NumProc() != 1)
    IsParallel_ = true;
  else
    IsParallel_ = false;

  ++NumInitialize_;
  InitializeTime_ += Time_->ElapsedTime();
  IsInitialized_ = true;
  return(0);
}

//==============================================================================
int Ifpack_PointRelaxation::Compute()
{
  int ierr = 0;
  if (!IsInitialized())
    IFPACK_CHK_ERR(Initialize());

  Time_->ResetStartTime();

  // reset values
  IsComputed_ = false;
  Condest_ = -1.0;

  if (NumSweeps_ == 0) ierr = 1; // Warning: no sweeps performed.
  if (NumSweeps_ < 0)
    IFPACK_CHK_ERR(-2); // at least one application

  // NOTE: RowMatrixRowMap doesn't work correctly for Epetra_VbrMatrix
  const Epetra_VbrMatrix * VbrMat = dynamic_cast<const Epetra_VbrMatrix*>(&Matrix());
  if(VbrMat)  Diagonal_ = Teuchos::rcp( new Epetra_Vector(VbrMat->RowMap()) );
  else Diagonal_ = Teuchos::rcp( new Epetra_Vector(Matrix().RowMatrixRowMap()) );

  if (Diagonal_ == Teuchos::null)
    IFPACK_CHK_ERR(-5);

  IFPACK_CHK_ERR(Matrix().ExtractDiagonalCopy(*Diagonal_));

  // Setup for L1 Methods.
  // Here we add half the value of the off-processor entries in the row,
  // but only if diagonal isn't sufficiently large.
  //
  // Note: This is only done in the slower-but-more-general "RowMatrix" mode.
  //
  // This follows from Equation (6.5) in:
  // Baker, Falgout, Kolev and Yang.  Multigrid Smoothers for Ultraparallel Computing.
  // SIAM J. Sci. Comput., Vol. 33, No. 5. (2011), pp. 2864--2887.
  if(DoL1Method_ && IsParallel_) {
    int maxLength = Matrix().MaxNumEntries();
    std::vector<int> Indices(maxLength);
    std::vector<double> Values(maxLength);
    int NumEntries;

    for (int i = 0 ; i < NumMyRows_ ; ++i) {
      IFPACK_CHK_ERR(Matrix_->ExtractMyRowCopy(i, maxLength,NumEntries,
                                               &Values[0], &Indices[0]));
      double diagonal_boost=0.0;
      for (int k = 0 ; k < NumEntries ; ++k)
        if(Indices[k] > NumMyRows_)
          diagonal_boost+=std::abs(Values[k]/2.0);
      if ((*Diagonal_)[i] < L1Eta_*diagonal_boost)
        (*Diagonal_)[i]+=diagonal_boost;
    }
  }

  // check diagonal elements, store the inverses, and verify that
  // no zeros are around. If an element is zero, then by default
  // its inverse is zero as well (that is, the row is ignored).
  for (int i = 0 ; i < NumMyRows_ ; ++i) {
    double& diag = (*Diagonal_)[i];
    if (IFPACK_ABS(diag) < MinDiagonalValue_)
      diag = MinDiagonalValue_;
    if (diag != 0.0)
      diag = 1.0 / diag;
  }
#ifdef IFPACK_FLOPCOUNTERS
  ComputeFlops_ += NumMyRows_;
#endif

#if 0
  // some methods require the inverse of the diagonal, compute it
  // now and store it in `Diagonal_'
  if ((PrecType_ == IFPACK_JACOBI) || (PrecType_ == IFPACK_GS)) {
    Diagonal_->Reciprocal(*Diagonal_);
    // update flops
    ComputeFlops_ += NumMyRows_;
  }
#endif

  // We need to import data from external processors. Here I create an
  // Epetra_Import object because I cannot assume that Matrix_ has one.
  // This is a bit of waste of resources (but the code is more robust).
  // Note that I am doing some strange stuff to set the components of Y
  // from Y2 (to save some time).
  //
  if (IsParallel_ && ((PrecType_ == IFPACK_GS) || (PrecType_ == IFPACK_SGS))) {
    Importer_ = Teuchos::rcp( new Epetra_Import(Matrix().RowMatrixColMap(),
                                  Matrix().RowMatrixRowMap()) );
    if (Importer_ == Teuchos::null) IFPACK_CHK_ERR(-5);
  }

  ++NumCompute_;
  ComputeTime_ += Time_->ElapsedTime();
  IsComputed_ = true;

  IFPACK_CHK_ERR(ierr);
  return(0);
}

//==============================================================================
std::ostream& Ifpack_PointRelaxation::Print(std::ostream & os) const
{
  using std::endl;

  double MyMinVal, MyMaxVal;
  double MinVal, MaxVal;

  if (IsComputed_) {
    Diagonal_->MinValue(&MyMinVal);
    Diagonal_->MaxValue(&MyMaxVal);
    Comm().MinAll(&MyMinVal,&MinVal,1);
    Comm().MinAll(&MyMaxVal,&MaxVal,1);
  }

  if (!Comm().MyPID()) {
    os << endl;
    os << "================================================================================" << endl;
    os << "Ifpack_PointRelaxation" << endl;
    os << "Sweeps         = " << NumSweeps_ << endl;
    os << "damping factor = " << DampingFactor_ << endl;
    if (PrecType_ == IFPACK_JACOBI)
      os << "Type           = Jacobi" << endl;
    else if (PrecType_ == IFPACK_GS)
      os << "Type           = Gauss-Seidel" << endl;
    else if (PrecType_ == IFPACK_SGS)
      os << "Type           = symmetric Gauss-Seidel" << endl;
    if (DoBackwardGS_)
      os << "Using backward mode (GS only)" << endl;
    if (ZeroStartingSolution_)
      os << "Using zero starting solution" << endl;
    else
      os << "Using input starting solution" << endl;
    os << "Condition number estimate = " << Condest() << endl;
    os << "Global number of rows            = " << Matrix_->NumGlobalRows64() << endl;
    if (IsComputed_) {
      os << "Minimum value on stored diagonal = " << MinVal << endl;
      os << "Maximum value on stored diagonal = " << MaxVal << endl;
    }
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
double Ifpack_PointRelaxation::
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
void Ifpack_PointRelaxation::SetLabel()
{
  std::string PT;
  if (PrecType_ == IFPACK_JACOBI)
    PT = "Jacobi";
  else if (PrecType_ == IFPACK_GS){
    PT = "GS";
    if(DoBackwardGS_)
      PT = "Backward " + PT;
  }
  else if (PrecType_ == IFPACK_SGS)
    PT = "SGS";

  if(NumLocalSmoothingIndices_!=NumMyRows_) PT="Local " + PT;
  else if(LocalSmoothingIndices_) PT="Reordered " + PT;

  Label_ = "IFPACK (" + PT + ", sweeps=" + Ifpack_toString(NumSweeps_)
    + ", damping=" + Ifpack_toString(DampingFactor_) + ")";
}

//==============================================================================
// Note that Ifpack_PointRelaxation and Jacobi is much faster than
// Ifpack_AdditiveSchwarz<Ifpack_PointRelaxation> (because of the
// way the matrix-vector produce is performed).
//
// Another ML-related observation is that the starting solution (in Y)
// is NOT supposed to be zero. This may slow down the application of just
// one sweep of Jacobi.
//
int Ifpack_PointRelaxation::
ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  if (!IsComputed())
    IFPACK_CHK_ERR(-3);

  if (X.NumVectors() != Y.NumVectors())
    IFPACK_CHK_ERR(-2);

  Time_->ResetStartTime();

  // AztecOO gives X and Y pointing to the same memory location,
  // need to create an auxiliary vector, Xcopy
  Teuchos::RefCountPtr< const Epetra_MultiVector > Xcopy;
  if (X.Pointers()[0] == Y.Pointers()[0])
    Xcopy = Teuchos::rcp( new Epetra_MultiVector(X) );
  else
    Xcopy = Teuchos::rcp( &X, false );

  // Flops are updated in each of the following.
  switch (PrecType_) {
  case IFPACK_JACOBI:
    IFPACK_CHK_ERR(ApplyInverseJacobi(*Xcopy,Y));
    break;
  case IFPACK_GS:
    IFPACK_CHK_ERR(ApplyInverseGS(*Xcopy,Y));
    break;
  case IFPACK_SGS:
    IFPACK_CHK_ERR(ApplyInverseSGS(*Xcopy,Y));
    break;
  default:
    IFPACK_CHK_ERR(-1); // something wrong
  }

  ++NumApplyInverse_;
  ApplyInverseTime_ += Time_->ElapsedTime();
  return(0);
}

//==============================================================================
// This preconditioner can be much slower than AztecOO and ML versions
// if the matrix-vector product requires all ExtractMyRowCopy()
// (as done through Ifpack_AdditiveSchwarz).
int Ifpack_PointRelaxation::
ApplyInverseJacobi(const Epetra_MultiVector& RHS, Epetra_MultiVector& LHS) const
{
  int NumVectors = LHS.NumVectors();

  int startIter = 0;
  if (NumSweeps_ > 0 && ZeroStartingSolution_) {
    // When we have a zero initial guess, we can skip the first
    // matrix apply call and zero initialization
    for (int v = 0; v < NumVectors; v++)
      IFPACK_CHK_ERR(LHS(v)->Multiply(DampingFactor_, *(RHS(v)), *Diagonal_, 0.0));

    startIter = 1;
  }

  bool zeroOut = false;
  Epetra_MultiVector A_times_LHS(LHS.Map(), NumVectors, zeroOut);
  for (int j = startIter; j < NumSweeps_ ; j++) {
    IFPACK_CHK_ERR(Apply(LHS, A_times_LHS));
    IFPACK_CHK_ERR(A_times_LHS.Update(1.0,RHS,-1.0));
    for (int v = 0 ; v < NumVectors ; ++v)
      IFPACK_CHK_ERR(LHS(v)->Multiply(DampingFactor_, *(A_times_LHS(v)),
                                     *Diagonal_, 1.0));

  }

  // Flops:
  // - matrix vector              (2 * NumGlobalNonzeros_)
  // - update                     (2 * NumGlobalRows_)
  // - Multiply:
  //   - DampingFactor            (NumGlobalRows_)
  //   - Diagonal                 (NumGlobalRows_)
  //   - A + B                    (NumGlobalRows_)
  //   - 1.0                      (NumGlobalRows_)
#ifdef IFPACK_FLOPCOUNTERS
  ApplyInverseFlops_ += NumVectors * (6 * NumGlobalRows_ + 2 * NumGlobalNonzeros_);
#endif

  return(0);
}

//==============================================================================
int Ifpack_PointRelaxation::
ApplyInverseGS(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  if (ZeroStartingSolution_)
      Y.PutScalar(0.0);

  const Epetra_CrsMatrix* CrsMatrix = dynamic_cast<const Epetra_CrsMatrix*>(&*Matrix_);
  // try to pick the best option; performances may be improved
  // if several sweeps are used.
  if (CrsMatrix != 0)
  {
    if (CrsMatrix->StorageOptimized() && LocalSmoothingIndices_)
      return(ApplyInverseGS_LocalFastCrsMatrix(CrsMatrix, X, Y));
    else if (CrsMatrix->StorageOptimized())
      return(ApplyInverseGS_FastCrsMatrix(CrsMatrix, X, Y));
    else
      return(ApplyInverseGS_CrsMatrix(CrsMatrix, X, Y));
  }
  else
    return(ApplyInverseGS_RowMatrix(X, Y));
} //ApplyInverseGS()



//==============================================================================
int Ifpack_PointRelaxation::
ApplyInverseGS_RowMatrix(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  int NumVectors = X.NumVectors();

  int Length = Matrix().MaxNumEntries();
  std::vector<int> Indices(Length);
  std::vector<double> Values(Length);

  Teuchos::RefCountPtr< Epetra_MultiVector > Y2;
  if (IsParallel_)
    Y2 = Teuchos::rcp( new Epetra_MultiVector(Importer_->TargetMap(), NumVectors) );
  else
    Y2 = Teuchos::rcp( &Y, false );

  // extract views (for nicer and faster code)
  double** y_ptr, ** y2_ptr, ** x_ptr, *d_ptr;
  X.ExtractView(&x_ptr);
  Y.ExtractView(&y_ptr);
  Y2->ExtractView(&y2_ptr);
  Diagonal_->ExtractView(&d_ptr);

  for (int j = 0; j < NumSweeps_ ; j++) {

    // data exchange is here, once per sweep
    if (IsParallel_)
      IFPACK_CHK_ERR(Y2->Import(Y,*Importer_,Insert));

    // FIXME: do I really need this code below?
    if (NumVectors == 1) {

      double* y0_ptr = y_ptr[0];
      double* y20_ptr = y2_ptr[0];
      double* x0_ptr = x_ptr[0];

      if(!DoBackwardGS_){
        /* Forward Mode */
        for (int ii = 0 ; ii < NumLocalSmoothingIndices_ ; ++ii) {
          int i = (!LocalSmoothingIndices_)? ii : LocalSmoothingIndices_[ii];

          int NumEntries;
          int col;
          IFPACK_CHK_ERR(Matrix_->ExtractMyRowCopy(i, Length,NumEntries,
                                                   &Values[0], &Indices[0]));

          double dtemp = 0.0;
          for (int k = 0 ; k < NumEntries ; ++k) {

            col = Indices[k];
            dtemp += Values[k] * y20_ptr[col];
          }

          y20_ptr[i] += DampingFactor_ * d_ptr[i] * (x0_ptr[i] - dtemp);
        }
      }
      else {
        /* Backward Mode */
        for (int ii = NumLocalSmoothingIndices_  - 1 ; ii > -1 ; --ii) {
          int i = (!LocalSmoothingIndices_)? ii : LocalSmoothingIndices_[ii];

          int NumEntries;
          // int col; // unused
          // (void) col; // Forestall compiler warning for unused variable.
          IFPACK_CHK_ERR(Matrix_->ExtractMyRowCopy(i, Length,NumEntries,
                                                   &Values[0], &Indices[0]));
          double dtemp = 0.0;
          for (int k = 0 ; k < NumEntries ; ++k) {
            //col = Indices[k]; // unused
            dtemp += Values[k] * y20_ptr[i];
          }

          y20_ptr[i] += DampingFactor_ * d_ptr[i] * (x0_ptr[i] - dtemp);
        }
      }

      // using Export() sounded quite expensive
      if (IsParallel_)
        for (int i = 0 ; i < NumMyRows_ ; ++i)
          y0_ptr[i] = y20_ptr[i];

    }
    else {
      if(!DoBackwardGS_){
        /* Forward Mode */
        for (int i = 0 ; i < NumMyRows_ ; ++i) {

          int NumEntries;
          int col;
          IFPACK_CHK_ERR(Matrix_->ExtractMyRowCopy(i, Length,NumEntries,
                                                   &Values[0], &Indices[0]));

          for (int m = 0 ; m < NumVectors ; ++m) {

            double dtemp = 0.0;
            for (int k = 0 ; k < NumEntries ; ++k) {

              col = Indices[k];
              dtemp += Values[k] * y2_ptr[m][col];
            }

            y2_ptr[m][i] += DampingFactor_ * d_ptr[i] * (x_ptr[m][i] - dtemp);
          }
        }
      }
      else {
        /* Backward Mode */
        for (int i = NumMyRows_  - 1 ; i > -1 ; --i) {
          int NumEntries;
          int col;
          IFPACK_CHK_ERR(Matrix_->ExtractMyRowCopy(i, Length,NumEntries,
                                                   &Values[0], &Indices[0]));

          for (int m = 0 ; m < NumVectors ; ++m) {

            double dtemp = 0.0;
            for (int k = 0 ; k < NumEntries ; ++k) {

              col = Indices[k];
              dtemp += Values[k] * y2_ptr[m][col];
            }

            y2_ptr[m][i] += DampingFactor_ * d_ptr[i] * (x_ptr[m][i] - dtemp);

          }
        }
      }

      // using Export() sounded quite expensive
      if (IsParallel_)
        for (int m = 0 ; m < NumVectors ; ++m)
          for (int ii = 0 ; ii < NumLocalSmoothingIndices_ ; ++ii) {
            int i = (!LocalSmoothingIndices_)? ii : LocalSmoothingIndices_[ii];
            y_ptr[m][i] = y2_ptr[m][i];
          }
    }
  }

#ifdef IFPACK_FLOPCOUNTERS
  ApplyInverseFlops_ += NumVectors * (4 * NumGlobalRows_ + 2 * NumGlobalNonzeros_);
#endif

  return(0);
} //ApplyInverseGS_RowMatrix()

//==============================================================================
int Ifpack_PointRelaxation::
ApplyInverseGS_CrsMatrix(const Epetra_CrsMatrix* A, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  int NumVectors = X.NumVectors();

  int* Indices;
  double* Values;

  Teuchos::RefCountPtr< Epetra_MultiVector > Y2;
  if (IsParallel_) {
    Y2 = Teuchos::rcp( new Epetra_MultiVector(Importer_->TargetMap(), NumVectors) );
  }
  else
    Y2 = Teuchos::rcp( &Y, false );

  double** y_ptr, ** y2_ptr, ** x_ptr, *d_ptr;
  X.ExtractView(&x_ptr);
  Y.ExtractView(&y_ptr);
  Y2->ExtractView(&y2_ptr);
  Diagonal_->ExtractView(&d_ptr);

  for (int iter = 0 ; iter < NumSweeps_ ; ++iter) {

    // only one data exchange per sweep
    if (IsParallel_)
      IFPACK_CHK_ERR(Y2->Import(Y,*Importer_,Insert));

    if(!DoBackwardGS_){
      /* Forward Mode */
      for (int ii = 0 ; ii < NumLocalSmoothingIndices_ ; ++ii) {
        int i = (!LocalSmoothingIndices_)? ii : LocalSmoothingIndices_[ii];

        int NumEntries;
        int col;
        double diag = d_ptr[i];

        IFPACK_CHK_ERR(A->ExtractMyRowView(i, NumEntries, Values, Indices));

        for (int m = 0 ; m < NumVectors ; ++m) {

          double dtemp = 0.0;

          for (int k = 0; k < NumEntries; ++k) {

            col = Indices[k];
            dtemp += Values[k] * y2_ptr[m][col];
          }

          y2_ptr[m][i] += DampingFactor_ * (x_ptr[m][i] - dtemp) * diag;
        }
      }
    }
    else {
      /* Backward Mode */
      for (int ii = NumLocalSmoothingIndices_  - 1 ; ii > -1 ; --ii) {
        int i = (!LocalSmoothingIndices_)? ii : LocalSmoothingIndices_[ii];

        int NumEntries;
        int col;
        double diag = d_ptr[i];

        IFPACK_CHK_ERR(A->ExtractMyRowView(i, NumEntries, Values, Indices));

        for (int m = 0 ; m < NumVectors ; ++m) {

          double dtemp = 0.0;
          for (int k = 0; k < NumEntries; ++k) {

            col = Indices[k];
            dtemp += Values[k] * y2_ptr[m][col];
          }

          y2_ptr[m][i] += DampingFactor_ * (x_ptr[m][i] - dtemp) * diag;

        }
      }
    }

    if (IsParallel_)
      for (int m = 0 ; m < NumVectors ; ++m)
        for (int ii = 0 ; ii < NumLocalSmoothingIndices_ ; ++ii) {
          int i = (!LocalSmoothingIndices_)? ii : LocalSmoothingIndices_[ii];
          y_ptr[m][i] = y2_ptr[m][i];
        }
  }

#ifdef IFPACK_FLOPCOUNTERS
  ApplyInverseFlops_ += NumVectors * (8 * NumGlobalRows_ + 4 * NumGlobalNonzeros_);
#endif
  return(0);
} //ApplyInverseGS_CrsMatrix()

//
//==============================================================================
// ApplyInverseGS_FastCrsMatrix() requires Epetra_CrsMatrix + OptimizeStorage()

int Ifpack_PointRelaxation::
ApplyInverseGS_FastCrsMatrix(const Epetra_CrsMatrix* A, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  int* IndexOffset;
  int* Indices;
  double* Values;
  IFPACK_CHK_ERR(A->ExtractCrsDataPointers(IndexOffset, Indices, Values));

  int NumVectors = X.NumVectors();

  Teuchos::RefCountPtr< Epetra_MultiVector > Y2;
  if (IsParallel_) {
    Y2 = Teuchos::rcp( new Epetra_MultiVector(Importer_->TargetMap(), NumVectors) );
  }
  else
    Y2 = Teuchos::rcp( &Y, false );

  double** y_ptr, ** y2_ptr, ** x_ptr, *d_ptr;
  X.ExtractView(&x_ptr);
  Y.ExtractView(&y_ptr);
  Y2->ExtractView(&y2_ptr);
  Diagonal_->ExtractView(&d_ptr);

  for (int iter = 0 ; iter < NumSweeps_ ; ++iter) {

    // only one data exchange per sweep
    if (IsParallel_)
      IFPACK_CHK_ERR(Y2->Import(Y,*Importer_,Insert));

    if(!DoBackwardGS_){
      /* Forward Mode */
      for (int i = 0 ; i < NumMyRows_ ; ++i) {

        int col;
        double diag = d_ptr[i];

        for (int m = 0 ; m < NumVectors ; ++m) {

          double dtemp = 0.0;

          for (int k = IndexOffset[i] ; k < IndexOffset[i + 1] ; ++k) {

            col = Indices[k];
            dtemp += Values[k] * y2_ptr[m][col];
          }

          y2_ptr[m][i] += DampingFactor_ * (x_ptr[m][i] - dtemp) * diag;
        }
      }
    }
    else {
      /* Backward Mode */
      for (int i = NumMyRows_  - 1 ; i > -1 ; --i) {

        int col;
        double diag = d_ptr[i];

        for (int m = 0 ; m < NumVectors ; ++m) {

          double dtemp = 0.0;
          for (int k = IndexOffset[i] ; k < IndexOffset[i + 1] ; ++k) {

            col = Indices[k];
            dtemp += Values[k] * y2_ptr[m][col];
          }

          y2_ptr[m][i] += DampingFactor_ * (x_ptr[m][i] - dtemp) * diag;

        }
      }
    }


    if (IsParallel_)
      for (int m = 0 ; m < NumVectors ; ++m)
        for (int i = 0 ; i < NumMyRows_ ; ++i)
          y_ptr[m][i] = y2_ptr[m][i];
  }

#ifdef IFPACK_FLOPCOUNTERS
  ApplyInverseFlops_ += NumVectors * (8 * NumGlobalRows_ + 4 * NumGlobalNonzeros_);
#endif
  return(0);
} //ApplyInverseGS_FastCrsMatrix()


//
//==============================================================================
// ApplyInverseGS_LocalFastCrsMatrix() requires Epetra_CrsMatrix + OptimizeStorage() + LocalSmoothingIndices

int Ifpack_PointRelaxation::
ApplyInverseGS_LocalFastCrsMatrix(const Epetra_CrsMatrix* A, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  int* IndexOffset;
  int* Indices;
  double* Values;
  IFPACK_CHK_ERR(A->ExtractCrsDataPointers(IndexOffset, Indices, Values));

  int NumVectors = X.NumVectors();

  Teuchos::RefCountPtr< Epetra_MultiVector > Y2;
  if (IsParallel_) {
    Y2 = Teuchos::rcp( new Epetra_MultiVector(Importer_->TargetMap(), NumVectors) );
  }
  else
    Y2 = Teuchos::rcp( &Y, false );

  double** y_ptr, ** y2_ptr, ** x_ptr, *d_ptr;
  X.ExtractView(&x_ptr);
  Y.ExtractView(&y_ptr);
  Y2->ExtractView(&y2_ptr);
  Diagonal_->ExtractView(&d_ptr);

  for (int iter = 0 ; iter < NumSweeps_ ; ++iter) {

    // only one data exchange per sweep
    if (IsParallel_)
      IFPACK_CHK_ERR(Y2->Import(Y,*Importer_,Insert));

    if(!DoBackwardGS_){
      /* Forward Mode */
      for (int ii = 0 ; ii < NumLocalSmoothingIndices_ ; ++ii) {
        int i=LocalSmoothingIndices_[ii];

        int col;
        double diag = d_ptr[i];

        for (int m = 0 ; m < NumVectors ; ++m) {

          double dtemp = 0.0;

          for (int k = IndexOffset[i] ; k < IndexOffset[i + 1] ; ++k) {

            col = Indices[k];
            dtemp += Values[k] * y2_ptr[m][col];
          }

          y2_ptr[m][i] += DampingFactor_ * (x_ptr[m][i] - dtemp) * diag;
        }
      }
    }
    else {
      /* Backward Mode */
      for (int ii = NumLocalSmoothingIndices_  - 1 ; ii > -1 ; --ii) {
        int i=LocalSmoothingIndices_[ii];

        int col;
        double diag = d_ptr[i];

        for (int m = 0 ; m < NumVectors ; ++m) {

          double dtemp = 0.0;
          for (int k = IndexOffset[i] ; k < IndexOffset[i + 1] ; ++k) {

            col = Indices[k];
            dtemp += Values[k] * y2_ptr[m][col];
          }

          y2_ptr[m][i] += DampingFactor_ * (x_ptr[m][i] - dtemp) * diag;

        }
      }
    }


    if (IsParallel_)
      for (int m = 0 ; m < NumVectors ; ++m)
        for (int ii = 0 ; ii < NumLocalSmoothingIndices_ ; ++ii) {
          int i=LocalSmoothingIndices_[ii];
          y_ptr[m][i] = y2_ptr[m][i];
        }
  }

#ifdef IFPACK_FLOPCOUNTERS
  ApplyInverseFlops_ += NumVectors * (8 * NumGlobalRows_ + 4 * NumGlobalNonzeros_);
#endif
  return(0);
} //ApplyInverseGS_LocalFastCrsMatrix()


//==============================================================================
int Ifpack_PointRelaxation::
ApplyInverseSGS(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  if (ZeroStartingSolution_)
    Y.PutScalar(0.0);

  const Epetra_CrsMatrix* CrsMatrix = dynamic_cast<const Epetra_CrsMatrix*>(&*Matrix_);
  // try to pick the best option; performances may be improved
  // if several sweeps are used.
  if (CrsMatrix != 0)
  {
    if (CrsMatrix->StorageOptimized() && LocalSmoothingIndices_)
      return(ApplyInverseSGS_LocalFastCrsMatrix(CrsMatrix, X, Y));
    else if (CrsMatrix->StorageOptimized())
      return(ApplyInverseSGS_FastCrsMatrix(CrsMatrix, X, Y));
    else
      return(ApplyInverseSGS_CrsMatrix(CrsMatrix, X, Y));
  }
  else
    return(ApplyInverseSGS_RowMatrix(X, Y));
}

//==============================================================================
int Ifpack_PointRelaxation::
ApplyInverseSGS_RowMatrix(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  int NumVectors = X.NumVectors();
  int Length = Matrix().MaxNumEntries();
  std::vector<int> Indices(Length);
  std::vector<double> Values(Length);

  Teuchos::RefCountPtr< Epetra_MultiVector > Y2;
  if (IsParallel_) {
    Y2 = Teuchos::rcp( new Epetra_MultiVector(Importer_->TargetMap(), NumVectors) );
  }
  else
    Y2 = Teuchos::rcp( &Y, false );

  double** y_ptr, ** y2_ptr, ** x_ptr, *d_ptr;
  X.ExtractView(&x_ptr);
  Y.ExtractView(&y_ptr);
  Y2->ExtractView(&y2_ptr);
  Diagonal_->ExtractView(&d_ptr);

  for (int iter = 0 ; iter < NumSweeps_ ; ++iter) {

    // only one data exchange per sweep
    if (IsParallel_)
      IFPACK_CHK_ERR(Y2->Import(Y,*Importer_,Insert));

    for (int ii = 0 ; ii < NumLocalSmoothingIndices_ ; ++ii) {
      int i = (!LocalSmoothingIndices_)? ii : LocalSmoothingIndices_[ii];

      int NumEntries;
      int col;
      double diag = d_ptr[i];

      IFPACK_CHK_ERR(Matrix_->ExtractMyRowCopy(i, Length,NumEntries,
                                               &Values[0], &Indices[0]));

      for (int m = 0 ; m < NumVectors ; ++m) {

        double dtemp = 0.0;

        for (int k = 0 ; k < NumEntries ; ++k) {

          col = Indices[k];
          dtemp += Values[k] * y2_ptr[m][col];
        }

        y2_ptr[m][i] += DampingFactor_ * (x_ptr[m][i] - dtemp) * diag;
      }
    }
    for (int ii = NumLocalSmoothingIndices_  - 1 ; ii > -1 ; --ii) {
      int i = (!LocalSmoothingIndices_)? ii : LocalSmoothingIndices_[ii];

      int NumEntries;
      int col;
      double diag = d_ptr[i];

      IFPACK_CHK_ERR(Matrix_->ExtractMyRowCopy(i, Length,NumEntries,
                                               &Values[0], &Indices[0]));

      for (int m = 0 ; m < NumVectors ; ++m) {

        double dtemp = 0.0;
        for (int k = 0 ; k < NumEntries ; ++k) {

          col = Indices[k];
          dtemp += Values[k] * y2_ptr[m][col];
        }

        y2_ptr[m][i] += DampingFactor_ * (x_ptr[m][i] - dtemp) * diag;

      }
    }

    if (IsParallel_)
      for (int m = 0 ; m < NumVectors ; ++m)
        for (int ii = 0 ; ii < NumLocalSmoothingIndices_ ; ++ii) {
          int i = (!LocalSmoothingIndices_)? ii : LocalSmoothingIndices_[ii];
          y_ptr[m][i] = y2_ptr[m][i];
        }
  }

#ifdef IFPACK_FLOPCOUNTERS
  ApplyInverseFlops_ += NumVectors * (8 * NumGlobalRows_ + 4 * NumGlobalNonzeros_);
#endif
  return(0);
}

//==============================================================================
int Ifpack_PointRelaxation::
ApplyInverseSGS_CrsMatrix(const Epetra_CrsMatrix* A, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  int NumVectors = X.NumVectors();

  int* Indices;
  double* Values;

  Teuchos::RefCountPtr< Epetra_MultiVector > Y2;
  if (IsParallel_) {
    Y2 = Teuchos::rcp( new Epetra_MultiVector(Importer_->TargetMap(), NumVectors) );
  }
  else
    Y2 = Teuchos::rcp( &Y, false );

  double** y_ptr, ** y2_ptr, ** x_ptr, *d_ptr;
  X.ExtractView(&x_ptr);
  Y.ExtractView(&y_ptr);
  Y2->ExtractView(&y2_ptr);
  Diagonal_->ExtractView(&d_ptr);

  for (int iter = 0 ; iter < NumSweeps_ ; ++iter) {

    // only one data exchange per sweep
    if (IsParallel_)
      IFPACK_CHK_ERR(Y2->Import(Y,*Importer_,Insert));

    for (int ii = 0 ; ii < NumLocalSmoothingIndices_ ; ++ii) {
      int i = (!LocalSmoothingIndices_)? ii : LocalSmoothingIndices_[ii];

      int NumEntries;
      int col;
      double diag = d_ptr[i];

      IFPACK_CHK_ERR(A->ExtractMyRowView(i, NumEntries, Values, Indices));

      for (int m = 0 ; m < NumVectors ; ++m) {

        double dtemp = 0.0;

        for (int k = 0; k < NumEntries; ++k) {

          col = Indices[k];
          dtemp += Values[k] * y2_ptr[m][col];
        }

        y2_ptr[m][i] += DampingFactor_ * (x_ptr[m][i] - dtemp) * diag;
      }
    }
    for (int ii = NumLocalSmoothingIndices_  - 1 ; ii > -1 ; --ii) {
      int i = (!LocalSmoothingIndices_)? ii : LocalSmoothingIndices_[ii];

      int NumEntries;
      int col;
      double diag = d_ptr[i];

      IFPACK_CHK_ERR(A->ExtractMyRowView(i, NumEntries, Values, Indices));

      for (int m = 0 ; m < NumVectors ; ++m) {

        double dtemp = 0.0;
        for (int k = 0; k < NumEntries; ++k) {

          col = Indices[k];
          dtemp += Values[k] * y2_ptr[m][col];
        }

        y2_ptr[m][i] += DampingFactor_ * (x_ptr[m][i] - dtemp) * diag;

      }
    }

    if (IsParallel_)
      for (int m = 0 ; m < NumVectors ; ++m)
        for (int ii = 0 ; ii < NumLocalSmoothingIndices_ ; ++ii) {
          int i = (!LocalSmoothingIndices_)? ii : LocalSmoothingIndices_[ii];
          y_ptr[m][i] = y2_ptr[m][i];
        }
  }

#ifdef IFPACK_FLOPCOUNTERS
  ApplyInverseFlops_ += NumVectors * (8 * NumGlobalRows_ + 4 * NumGlobalNonzeros_);
#endif
  return(0);
}

//==============================================================================
// this requires Epetra_CrsMatrix + OptimizeStorage()
int Ifpack_PointRelaxation::
ApplyInverseSGS_FastCrsMatrix(const Epetra_CrsMatrix* A, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  int* IndexOffset;
  int* Indices;
  double* Values;
  IFPACK_CHK_ERR(A->ExtractCrsDataPointers(IndexOffset, Indices, Values));

  int NumVectors = X.NumVectors();

  Teuchos::RefCountPtr< Epetra_MultiVector > Y2;
  if (IsParallel_) {
    Y2 = Teuchos::rcp( new Epetra_MultiVector(Importer_->TargetMap(), NumVectors) );
  }
  else
    Y2 = Teuchos::rcp( &Y, false );

  double** y_ptr, ** y2_ptr, ** x_ptr, *d_ptr;
  X.ExtractView(&x_ptr);
  Y.ExtractView(&y_ptr);
  Y2->ExtractView(&y2_ptr);
  Diagonal_->ExtractView(&d_ptr);

  for (int iter = 0 ; iter < NumSweeps_ ; ++iter) {

    // only one data exchange per sweep
    if (IsParallel_)
      IFPACK_CHK_ERR(Y2->Import(Y,*Importer_,Insert));

    for (int i = 0 ; i < NumMyRows_ ; ++i) {

      int col;
      double diag = d_ptr[i];

      // no need to extract row i
      //IFPACK_CHK_ERR(A->ExtractMyRowView(i, NumEntries, Values, Indices));

      for (int m = 0 ; m < NumVectors ; ++m) {

        double dtemp = 0.0;

        for (int k = IndexOffset[i] ; k < IndexOffset[i + 1] ; ++k) {

          col = Indices[k];
          dtemp += Values[k] * y2_ptr[m][col];
        }

        y2_ptr[m][i] += DampingFactor_ * (x_ptr[m][i] - dtemp) * diag;
      }
    }

    for (int i = NumMyRows_  - 1 ; i > -1 ; --i) {

      int col;
      double diag = d_ptr[i];

      // no need to extract row i
      //IFPACK_CHK_ERR(A->ExtractMyRowView(i, NumEntries, Values, Indices));

      for (int m = 0 ; m < NumVectors ; ++m) {

        double dtemp = 0.0;
        for (int k = IndexOffset[i] ; k < IndexOffset[i + 1] ; ++k) {

          col = Indices[k];
          dtemp += Values[k] * y2_ptr[m][col];
        }

        y2_ptr[m][i] += DampingFactor_ * (x_ptr[m][i] - dtemp) * diag;

      }
    }

    if (IsParallel_)
      for (int m = 0 ; m < NumVectors ; ++m)
        for (int i = 0 ; i < NumMyRows_ ; ++i)
          y_ptr[m][i] = y2_ptr[m][i];
  }

#ifdef IFPACK_FLOPCOUNTERS
  ApplyInverseFlops_ += NumVectors * (8 * NumGlobalRows_ + 4 * NumGlobalNonzeros_);
#endif
  return(0);
}


//==============================================================================
// this requires Epetra_CrsMatrix + OptimizeStorage() + LocalSmoothingIndices_
int Ifpack_PointRelaxation::
ApplyInverseSGS_LocalFastCrsMatrix(const Epetra_CrsMatrix* A, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  int* IndexOffset;
  int* Indices;
  double* Values;
  IFPACK_CHK_ERR(A->ExtractCrsDataPointers(IndexOffset, Indices, Values));

  int NumVectors = X.NumVectors();

  Teuchos::RefCountPtr< Epetra_MultiVector > Y2;
  if (IsParallel_) {
    Y2 = Teuchos::rcp( new Epetra_MultiVector(Importer_->TargetMap(), NumVectors) );
  }
  else
    Y2 = Teuchos::rcp( &Y, false );

  double** y_ptr, ** y2_ptr, ** x_ptr, *d_ptr;
  X.ExtractView(&x_ptr);
  Y.ExtractView(&y_ptr);
  Y2->ExtractView(&y2_ptr);
  Diagonal_->ExtractView(&d_ptr);

  for (int iter = 0 ; iter < NumSweeps_ ; ++iter) {

    // only one data exchange per sweep
    if (IsParallel_)
      IFPACK_CHK_ERR(Y2->Import(Y,*Importer_,Insert));

    for (int ii = 0 ; ii < NumLocalSmoothingIndices_ ; ++ii) {
      int i = (!LocalSmoothingIndices_)? ii : LocalSmoothingIndices_[ii];

      int col;
      double diag = d_ptr[i];

      // no need to extract row i
      //IFPACK_CHK_ERR(A->ExtractMyRowView(i, NumEntries, Values, Indices));

      for (int m = 0 ; m < NumVectors ; ++m) {

        double dtemp = 0.0;

        for (int k = IndexOffset[i] ; k < IndexOffset[i + 1] ; ++k) {

          col = Indices[k];
          dtemp += Values[k] * y2_ptr[m][col];
        }

        y2_ptr[m][i] += DampingFactor_ * (x_ptr[m][i] - dtemp) * diag;
      }
    }

    for (int ii = NumLocalSmoothingIndices_  - 1 ; ii > -1 ; --ii) {
      int i = (!LocalSmoothingIndices_)? ii : LocalSmoothingIndices_[ii];

      int col;
      double diag = d_ptr[i];

      // no need to extract row i
      //IFPACK_CHK_ERR(A->ExtractMyRowView(i, NumEntries, Values, Indices));

      for (int m = 0 ; m < NumVectors ; ++m) {

        double dtemp = 0.0;
        for (int k = IndexOffset[i] ; k < IndexOffset[i + 1] ; ++k) {

          col = Indices[k];
          dtemp += Values[k] * y2_ptr[m][col];
        }

        y2_ptr[m][i] += DampingFactor_ * (x_ptr[m][i] - dtemp) * diag;

      }
    }

    if (IsParallel_)
      for (int m = 0 ; m < NumVectors ; ++m)
        for (int ii = 0 ; ii < NumLocalSmoothingIndices_ ; ++ii) {
          int i = (!LocalSmoothingIndices_)? ii : LocalSmoothingIndices_[ii];
          y_ptr[m][i] = y2_ptr[m][i];
        }
  }

#ifdef IFPACK_FLOPCOUNTERS
  ApplyInverseFlops_ += NumVectors * (8 * NumGlobalRows_ + 4 * NumGlobalNonzeros_);
#endif
  return(0);
}
