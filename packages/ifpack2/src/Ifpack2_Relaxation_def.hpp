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

#ifndef IFPACK2_RELAXATION_DEF_HPP
#define IFPACK2_RELAXATION_DEF_HPP

#include "Ifpack2_Relaxation_decl.hpp"

namespace Ifpack2 {

//==========================================================================
template<class MatrixType>
Relaxation<MatrixType>::Relaxation(const Teuchos::RCP<const Tpetra::RowMatrix<scalar_type,local_ordinal_type,global_ordinal_type,node_type> >& A)
: A_(A),
  Comm_ (A->getRowMap ()->getComm ()),
  Time_ (Teuchos::rcp (new Teuchos::Time("Ifpack2::Relaxation"))),
  NumSweeps_ (1),
  PrecType_ (Ifpack2::JACOBI),
  MinDiagonalValue_ (Teuchos::as<scalar_type> (0.0)),
  DampingFactor_ (Teuchos::as<scalar_type> (1.0)),
  IsParallel_ (A->getRowMap ()->getComm ()->getSize () > 1),
  ZeroStartingSolution_ (true),
  DoBackwardGS_ (false),
  DoL1Method_ (false),
  L1Eta_ (Teuchos::as<magnitude_type> (1.5)),
  Condest_ (Teuchos::as<magnitude_type> (-1)),
  IsInitialized_ (false),
  IsComputed_ (false),
  NumInitialize_ (0),
  NumCompute_ (0),
  NumApply_ (0),
  InitializeTime_ (0.0), // Times are double anyway, so no need for ScalarTraits.
  ComputeTime_ (0.0),
  ApplyTime_ (0.0),
  ComputeFlops_ (0.0),
  ApplyFlops_ (0.0),
  NumMyRows_ (0),
  NumGlobalRows_ (0),
  NumGlobalNonzeros_ (0)
{
  this->setObjectLabel("Ifpack2::Relaxation");
  TEUCHOS_TEST_FOR_EXCEPTION(
    A_.is_null (),
    std::runtime_error,
    "Ifpack2::Relaxation(): The constructor needs a non-null input matrix.");
}

//==========================================================================
template<class MatrixType>
Relaxation<MatrixType>::~Relaxation() {
}

//==========================================================================
template<class MatrixType>
void Relaxation<MatrixType>::setParameters(const Teuchos::ParameterList& List)
{
  Teuchos::ParameterList validparams;
  Ifpack2::getValidParameters(validparams);
  List.validateParameters(validparams);

  std::string PT;
  if (PrecType_ == Ifpack2::JACOBI)
    PT = "Jacobi";
  else if (PrecType_ == Ifpack2::GS)
    PT = "Gauss-Seidel";
  else if (PrecType_ == Ifpack2::SGS)
    PT = "Symmetric Gauss-Seidel";

  Ifpack2::getParameter(List, "relaxation: type", PT);

  if (PT == "Jacobi") {
    PrecType_ = Ifpack2::JACOBI;
  }
  else if (PT == "Gauss-Seidel") {
    PrecType_ = Ifpack2::GS;
  }
  else if (PT == "Symmetric Gauss-Seidel") {
    PrecType_ = Ifpack2::SGS;
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(
      true,
      std::invalid_argument,
      "Ifpack2::Relaxation::setParameters: unsupported value for 'relaxation: type' parameter (\"" << PT << "\")");
  }

  Ifpack2::getParameter(List, "relaxation: sweeps",NumSweeps_);
  Ifpack2::getParameter(List, "relaxation: damping factor", DampingFactor_);
  Ifpack2::getParameter(List, "relaxation: min diagonal value", MinDiagonalValue_);
  Ifpack2::getParameter(List, "relaxation: zero starting solution", ZeroStartingSolution_);
  Ifpack2::getParameter(List, "relaxation: backward mode",DoBackwardGS_);
  Ifpack2::getParameter(List, "relaxation: use l1",DoL1Method_);
  Ifpack2::getParameter(List, "relaxation: l1 eta",L1Eta_);
}

//==========================================================================
template<class MatrixType>
const Teuchos::RCP<const Teuchos::Comm<int> > &
Relaxation<MatrixType>::getComm() const{
  return(Comm_);
}

//==========================================================================
template<class MatrixType>
Teuchos::RCP<const Tpetra::RowMatrix<typename MatrixType::scalar_type,
                                     typename MatrixType::local_ordinal_type,
                                     typename MatrixType::global_ordinal_type,
                                     typename MatrixType::node_type> >
Relaxation<MatrixType>::getMatrix() const {
  return(A_);
}

//==========================================================================
template<class MatrixType>
const Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,
                                     typename MatrixType::global_ordinal_type,
                                     typename MatrixType::node_type> >&
Relaxation<MatrixType>::getDomainMap() const {
  return A_->getDomainMap();
}

//==========================================================================
template<class MatrixType>
const Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,
                                     typename MatrixType::global_ordinal_type,
                                     typename MatrixType::node_type> >&
Relaxation<MatrixType>::getRangeMap() const {
  return A_->getRangeMap();
}

//==========================================================================
template<class MatrixType>
bool Relaxation<MatrixType>::hasTransposeApply() const {
  return true;
}

//==========================================================================
template<class MatrixType>
int Relaxation<MatrixType>::getNumInitialize() const {
  return(NumInitialize_);
}

//==========================================================================
template<class MatrixType>
int Relaxation<MatrixType>::getNumCompute() const {
  return(NumCompute_);
}

//==========================================================================
template<class MatrixType>
int Relaxation<MatrixType>::getNumApply() const {
  return(NumApply_);
}

//==========================================================================
template<class MatrixType>
double Relaxation<MatrixType>::getInitializeTime() const {
  return(InitializeTime_);
}

//==========================================================================
template<class MatrixType>
double Relaxation<MatrixType>::getComputeTime() const {
  return(ComputeTime_);
}

//==========================================================================
template<class MatrixType>
double Relaxation<MatrixType>::getApplyTime() const {
  return(ApplyTime_);
}

//==========================================================================
template<class MatrixType>
double Relaxation<MatrixType>::getComputeFlops() const {
  return(ComputeFlops_);
}

//==========================================================================
template<class MatrixType>
double Relaxation<MatrixType>::getApplyFlops() const {
  return(ApplyFlops_);
}

//==========================================================================
template<class MatrixType>
typename Teuchos::ScalarTraits<typename MatrixType::scalar_type>::magnitudeType
Relaxation<MatrixType>::getCondEst() const
{
  return(Condest_);
}

//==========================================================================
template<class MatrixType>
typename Teuchos::ScalarTraits<typename MatrixType::scalar_type>::magnitudeType
Relaxation<MatrixType>::computeCondEst(
                     CondestType CT,
                     typename MatrixType::local_ordinal_type MaxIters,
                     magnitude_type Tol,
     const Teuchos::Ptr<const Tpetra::RowMatrix<typename MatrixType::scalar_type,
                                                typename MatrixType::local_ordinal_type,
                                                typename MatrixType::global_ordinal_type,
                                                typename MatrixType::node_type> > &matrix)
{
  if (!isComputed()) // cannot compute right now
    return(-1.0);

  // always compute it. Call Condest() with no parameters to get
  // the previous estimate.
  Condest_ = Ifpack2::Condest(*this, CT, MaxIters, Tol, matrix);

  return(Condest_);
}

//==========================================================================
template<class MatrixType>
void Relaxation<MatrixType>::apply(
          const Tpetra::MultiVector<typename MatrixType::scalar_type,
                                    typename MatrixType::local_ordinal_type,
                                    typename MatrixType::global_ordinal_type,
                                    typename MatrixType::node_type>& X,
                Tpetra::MultiVector<typename MatrixType::scalar_type,
                                    typename MatrixType::local_ordinal_type,
                                    typename MatrixType::global_ordinal_type,
                                    typename MatrixType::node_type>& Y,
                Teuchos::ETransp mode,
                 scalar_type alpha,
                 scalar_type beta) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  typedef Teuchos::ScalarTraits<scalar_type> STS;
  typedef Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> MV;

  TEUCHOS_TEST_FOR_EXCEPTION(
    ! isComputed (),
    std::runtime_error,
    "Ifpack2::Relaxation::apply: You must call compute() on this Ifpack2 "
    "preconditioner instance before you may call apply().  You may call "
    "isComputed() to find out if compute() has been called already.");

  TEUCHOS_TEST_FOR_EXCEPTION(
    X.getNumVectors() != Y.getNumVectors(),
    std::runtime_error,
    "Ifpack2::Relaxation::apply: X and Y have different numbers of columns.  "
    "X has " << X.getNumVectors() << " columns, but Y has "
    << Y.getNumVectors() << " columns.");

  Time_->start(true);

  // If X and Y are pointing to the same memory location,
  // we need to create an auxiliary vector, Xcopy
  RCP<const MV> Xcopy;
  if (X.getLocalMV().getValues() == Y.getLocalMV().getValues()) {
    Xcopy = rcp (new MV (X));
  }
  else {
    Xcopy = rcpFromRef (X);
  }

  if (ZeroStartingSolution_) {
    Y.putScalar (STS::zero ());
  }

  // Flops are updated in each of the following.
  switch (PrecType_) {
  case Ifpack2::JACOBI:
    ApplyInverseJacobi(*Xcopy,Y);
    break;
  case Ifpack2::GS:
    ApplyInverseGS(*Xcopy,Y);
    break;
  case Ifpack2::SGS:
    ApplyInverseSGS(*Xcopy,Y);
    break;
  default:
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
      "Ifpack2::Relaxation::apply: Invalid preconditioner type enum value "
      << PrecType_ << ".  Please report this bug to the Ifpack2 developers.");
  }

  ++NumApply_;
  Time_->stop();
  ApplyTime_ += Time_->totalElapsedTime();
}

//==========================================================================
template<class MatrixType>
void Relaxation<MatrixType>::applyMat(
          const Tpetra::MultiVector<typename MatrixType::scalar_type,
                                    typename MatrixType::local_ordinal_type,
                                    typename MatrixType::global_ordinal_type,
                                    typename MatrixType::node_type>& X,
                Tpetra::MultiVector<typename MatrixType::scalar_type,
                                    typename MatrixType::local_ordinal_type,
                                    typename MatrixType::global_ordinal_type,
                                    typename MatrixType::node_type>& Y,
             Teuchos::ETransp mode) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(isComputed() == false, std::runtime_error,
     "Ifpack2::Relaxation::applyMat() ERROR: isComputed() must be true prior to calling applyMat().");
  TEUCHOS_TEST_FOR_EXCEPTION(X.getNumVectors() != Y.getNumVectors(), std::runtime_error,
     "Ifpack2::Relaxation::applyMat() ERROR: X.getNumVectors() != Y.getNumVectors().");
  A_->apply(X, Y, mode);
}

//==========================================================================
template<class MatrixType>
void Relaxation<MatrixType>::initialize() {
  IsInitialized_ = false;

  TEUCHOS_TEST_FOR_EXCEPTION(A_ == Teuchos::null, std::runtime_error,
    "Ifpack2::Relaxation::Initialize ERROR, Matrix is NULL");

  Time_->start(true);

  NumMyRows_ = A_->getNodeNumRows();
  NumGlobalRows_ = A_->getGlobalNumRows();
  NumGlobalNonzeros_ = A_->getGlobalNumEntries();

  if (Comm_->getSize() != 1) {
    IsParallel_ = true;
  }
  else {
    IsParallel_ = false;
  }

  ++NumInitialize_;
  Time_->stop();
  InitializeTime_ += Time_->totalElapsedTime();
  IsInitialized_ = true;
}

//==========================================================================
template<class MatrixType>
void Relaxation<MatrixType>::compute()
{
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::as;
  using Teuchos::rcp;
  typedef Teuchos::ScalarTraits<scalar_type> STS;
  typedef Teuchos::ScalarTraits<magnitude_type> STM;
  typedef Tpetra::Import<local_ordinal_type,global_ordinal_type,node_type> import_type;
  typedef Tpetra::Vector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> vector_type;

  if (! isInitialized ()) {
    initialize ();
  }

  Time_->start(true);

  // reset values
  IsComputed_ = false;
  Condest_ = -1.0;

  TEUCHOS_TEST_FOR_EXCEPTION(
    NumSweeps_ < 0,
    std::logic_error,
    "Ifpack2::Relaxation::compute: NumSweeps_ = " << NumSweeps_ << " < 0.  "
    "Please report this bug to the Ifpack2 developers.");

  Diagonal_ = rcp (new vector_type (A_->getRowMap ()));

  TEUCHOS_TEST_FOR_EXCEPTION(
    Diagonal_.is_null (),
    std::logic_error,
    "Ifpack2::Relaxation::compute: Vector of diagonal entries has not been "
    "created yet.  Please report this bug to the Ifpack2 developers.");

  A_->getLocalDiagCopy (*Diagonal_);
  ArrayRCP<scalar_type> DiagView = Diagonal_->get1dViewNonConst ();

  // Setup for L1 Methods.
  // Here we add half the value of the off-processor entries in the row,
  // but only if diagonal isn't sufficiently large.
  //
  // Note: This is only done in the slower-but-more-general "RowMatrix" mode.
  //
  // This follows from Equation (6.5) in:
  // Baker, Falgout, Kolev and Yang.  Multigrid Smoothers for Ultraparallel Computing.
  // SIAM J. Sci. Comput., Vol. 33, No. 5. (2011), pp. 2864--2887.
  if (DoL1Method_ && IsParallel_) {
    size_t maxLength = A_->getNodeMaxNumRowEntries ();
    Array<local_ordinal_type> Indices(maxLength);
    Array<scalar_type> Values(maxLength);
    size_t NumEntries;

    scalar_type two = STS::one () + STS::one ();

    for (size_t i = 0; i < NumMyRows_; ++i) {
      A_->getLocalRowCopy (i, Indices (), Values (), NumEntries);
      magnitude_type diagonal_boost = STM::zero ();
      for (size_t k = 0 ; k < NumEntries ; ++k) {
        if (as<size_t> (Indices[k]) > i) {
          diagonal_boost += STS::magnitude (Values[k] / two);
        }
      }
      if (STS::magnitude (DiagView[i]) < L1Eta_ * diagonal_boost) {
        DiagView[i] += diagonal_boost;
      }
    }
  }

  // check diagonal elements, store the inverses, and verify that
  // no zeros are around. If an element is zero, then by default
  // its inverse is zero as well (that is, the row is ignored).
  for (size_t i = 0 ; i < NumMyRows_ ; ++i) {
    scalar_type& diag = DiagView[i];
    if (STS::magnitude (diag) < STS::magnitude (MinDiagonalValue_)) {
      diag = MinDiagonalValue_;
    }
    if (diag != STS::zero ()) {
      diag = STS::one () / diag;
    }
  }
  ComputeFlops_ += NumMyRows_;


  // We need to import data from external processors. Here I create a
  // Tpetra::Import object if needed (stealing from A_ if possible)
  // Marzio's comment:
  // Note that I am doing some strange stuff to set the components of Y
  // from Y2 (to save some time).
  //
  if (IsParallel_ && ((PrecType_ == Ifpack2::GS) || (PrecType_ == Ifpack2::SGS))) {
    Importer_=A_->getGraph()->getImporter();
    if (Importer_.is_null ()) {
      Importer_ = rcp (new import_type (A_->getDomainMap (),
                                        A_->getColMap ()));
    }
  }

  ++NumCompute_;
  Time_->stop();
  ComputeTime_ += Time_->totalElapsedTime();
  IsComputed_ = true;
}

//==========================================================================
template<class MatrixType>
void Relaxation<MatrixType>::ApplyInverseJacobi(
        const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
              Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y) const
{
  typedef Teuchos::ScalarTraits<scalar_type> STS;
  typedef Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> MV;

  int NumVectors = X.getNumVectors();
  MV A_times_Y (Y.getMap (), NumVectors);

  for (int j = 0; j < NumSweeps_; ++j) {
    applyMat (Y, A_times_Y);
    A_times_Y.update (STS::one (), X, -STS::one ());
    Y.elementWiseMultiply (DampingFactor_, *Diagonal_, A_times_Y, STS::one ());
  }

  // For each column of output, for each pass over the matrix:
  //
  // - One + and one * for each matrix entry
  // - One / and one + for each row of the matrix
  // - If the damping factor is not one: one * for each row of the
  //   matrix.  (It's not fair to count this if the damping factor is
  //   one, since the implementation could skip it.  Whether it does
  //   or not is the implementation's choice.)

  // Floating-point operations due to the damping factor, per matrix
  // row, per direction, per columm of output.
  const size_t dampingFlops = (DampingFactor_ == STS::one()) ? 0 : 1;
  const size_t numVectors = X.getNumVectors ();

  ApplyFlops_ += NumSweeps_ * numVectors *
    (2 * NumGlobalRows_ + 2 * NumGlobalNonzeros_ + dampingFlops);
}

//==========================================================================
template<class MatrixType>
void Relaxation<MatrixType>::ApplyInverseGS(
        const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
              Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y) const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  // FIXME (mfh 02 Jan 2013) This assumes that MatrixType is a
  // CrsMatrix specialization.
  RCP<const MatrixType> crsMat = rcp_dynamic_cast<const MatrixType> (A_);

  // The CrsMatrix version is faster, because it can access the sparse
  // matrix data directly, rather than by copying out each row's data
  // in turn.  Thus, we check whether the RowMatrix is really a
  // CrsMatrix.
  if (! crsMat.is_null ()) {
    ApplyInverseGS_CrsMatrix (*crsMat, X, Y);
  }
  else {
    ApplyInverseGS_RowMatrix (X, Y);
  }
}

//==========================================================================
template<class MatrixType>
void Relaxation<MatrixType>::ApplyInverseGS_RowMatrix(
        const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
              Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y) const
{
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  typedef Teuchos::ScalarTraits<scalar_type> STS;
  typedef Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> MV;

  size_t NumVectors = X.getNumVectors();

  size_t maxLength = A_->getNodeMaxNumRowEntries();
  Array<local_ordinal_type> Indices(maxLength);
  Array<scalar_type> Values(maxLength);

  RCP<MV> Y2;
  if (IsParallel_) {
    Y2 = rcp (new MV (Importer_->getTargetMap (), NumVectors));
  }
  else {
    Y2 = rcpFromRef (Y);
  }

  // extract views (for nicer and faster code)
  ArrayRCP<ArrayRCP<scalar_type> > y_ptr = Y.get2dViewNonConst();
  ArrayRCP<ArrayRCP<scalar_type> > y2_ptr = Y2->get2dViewNonConst();
  ArrayRCP<ArrayRCP<const scalar_type> > x_ptr =  X.get2dView();
  ArrayRCP<const scalar_type> d_ptr = Diagonal_->get1dView();

  for (int j = 0; j < NumSweeps_; j++) {
    // data exchange is here, once per sweep
    if (IsParallel_) {
      Y2->doImport (Y, *Importer_, Tpetra::INSERT);
    }

    if (! DoBackwardGS_) { // Forward sweep
      for (size_t i = 0; i < NumMyRows_; ++i) {
        size_t NumEntries;
        A_->getLocalRowCopy (as<local_ordinal_type> (i), Indices (), Values (), NumEntries);

        for (size_t m = 0; m < NumVectors; ++m) {
          scalar_type dtemp = STS::zero ();
          for (size_t k = 0; k < NumEntries; ++k) {
            const local_ordinal_type col = Indices[k];
            dtemp += Values[k] * y2_ptr[m][col];
          }
          y2_ptr[m][i] += DampingFactor_ * d_ptr[i] * (x_ptr[m][i] - dtemp);
        }
      }
    }
    else { // Backward sweep
      // ptrdiff_t is the same size as size_t, but is signed.  Being
      // signed is important so that i >= 0 is not trivially true.
      for (ptrdiff_t i = NumMyRows_ - 1; i >= 0; --i) {
        size_t NumEntries;
        A_->getLocalRowCopy (as<local_ordinal_type> (i), Indices (), Values (), NumEntries);

        for (size_t m = 0; m < NumVectors; ++m) {
          scalar_type dtemp = STS::zero ();
          for (size_t k = 0; k < NumEntries; ++k) {
            const local_ordinal_type col = Indices[k];
            dtemp += Values[k] * y2_ptr[m][col];
          }
          y2_ptr[m][i] += DampingFactor_ * d_ptr[i] * (x_ptr[m][i] - dtemp);
        }
      }
    }

    // FIXME (mfh 02 Jan 2013) This is only correct if row Map == range Map.
    if (IsParallel_) {
      Y = *Y2;
    }
  }

  // See flop count discussion in implementation of ApplyInverseGS_CrsMatrix().
  const size_t dampingFlops = (DampingFactor_ == STS::one()) ? 0 : 1;
  const size_t numVectors = X.getNumVectors ();

  ApplyFlops_ += 2 * NumSweeps_ * numVectors *
    (2 * NumGlobalRows_ + 2 * NumGlobalNonzeros_ + dampingFlops);
}

//==========================================================================
template<class MatrixType>
void
Relaxation<MatrixType>::
ApplyInverseGS_CrsMatrix (const MatrixType& A,
                          const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
                          Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y) const
{
  typedef Teuchos::ScalarTraits<scalar_type> STS;

  const Tpetra::ESweepDirection direction =
    DoBackwardGS_ ? Tpetra::Backward : Tpetra::Forward;
  A.gaussSeidelCopy (Y, X, *Diagonal_, DampingFactor_, direction, NumSweeps_);

  // For each column of output, for each sweep over the matrix:
  //
  // - One + and one * for each matrix entry
  // - One / and one + for each row of the matrix
  // - If the damping factor is not one: one * for each row of the
  //   matrix.  (It's not fair to count this if the damping factor is
  //   one, since the implementation could skip it.  Whether it does
  //   or not is the implementation's choice.)

  // Floating-point operations due to the damping factor, per matrix
  // row, per direction, per columm of output.
  const size_t dampingFlops = (DampingFactor_ == STS::one()) ? 0 : 1;
  const size_t numVectors = X.getNumVectors ();

  ApplyFlops_ += NumSweeps_ * numVectors *
    (2 * NumGlobalRows_ + 2 * NumGlobalNonzeros_ + dampingFlops);
}

//==========================================================================
template<class MatrixType>
void Relaxation<MatrixType>::ApplyInverseSGS(
        const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
              Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y) const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  // FIXME (mfh 02 Jan 2013) This assumes that MatrixType is a
  // CrsMatrix specialization.  Alas, C++ doesn't let me do pattern
  // matching on template types.  I could just cast to the four
  // template argument specialization of CrsMatrix, but that would
  // miss nondefault values of the fifth template parameter of
  // CrsMatrix.
  //
  // Another way to solve this problem would be to move
  // implementations of relaxations into Tpetra.  Tpetra::RowMatrix
  // could provide a gaussSeidel, etc. interface, with a default
  // implementation (same as the RowMatrix version here in Ifpack2
  // now), and Tpetra::CrsMatrix specializations could reimplement
  // this however they wish.
  RCP<const MatrixType> crsMat = rcp_dynamic_cast<const MatrixType> (A_);

  // The CrsMatrix version is faster, because it can access the sparse
  // matrix data directly, rather than by copying out each row's data
  // in turn.  Thus, we check whether the RowMatrix is really a
  // CrsMatrix.
  if (! crsMat.is_null ()) {
    ApplyInverseSGS_CrsMatrix (*crsMat, X, Y);
  }
  else {
    ApplyInverseSGS_RowMatrix (X, Y);
  }
}

//==========================================================================
template<class MatrixType>
void Relaxation<MatrixType>::ApplyInverseSGS_RowMatrix(
        const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
              Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y) const
{
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  typedef Teuchos::ScalarTraits<scalar_type> STS;
  typedef Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> MV;

  size_t NumVectors = X.getNumVectors ();
  size_t maxLength = A_->getNodeMaxNumRowEntries ();
  Array<local_ordinal_type> Indices (maxLength);
  Array<scalar_type> Values (maxLength);

  RCP<MV> Y2;
  if (IsParallel_) {
    Y2 = rcp (new MV (Importer_->getTargetMap (), NumVectors));
  }
  else {
    Y2 = rcpFromRef (Y);
  }

  ArrayRCP<ArrayRCP<scalar_type> > y_ptr = Y.get2dViewNonConst ();
  ArrayRCP<ArrayRCP<scalar_type> > y2_ptr = Y2->get2dViewNonConst ();
  ArrayRCP<ArrayRCP<const scalar_type> > x_ptr =  X.get2dView ();
  ArrayRCP<const scalar_type> d_ptr = Diagonal_->get1dView ();

  for (int iter = 0; iter < NumSweeps_; ++iter) {
    // only one data exchange per sweep
    if (IsParallel_) {
      Y2->doImport (Y, *Importer_, Tpetra::INSERT);
    }

    for (size_t i = 0; i < NumMyRows_; ++i) {
      const scalar_type diag = d_ptr[i];
      size_t NumEntries;
      A_->getLocalRowCopy (as<local_ordinal_type> (i), Indices (), Values (), NumEntries);

      for (size_t m = 0; m < NumVectors; ++m) {
        scalar_type dtemp = STS::zero ();
        for (size_t k = 0; k < NumEntries; ++k) {
          const local_ordinal_type col = Indices[k];
          dtemp += Values[k] * y2_ptr[m][col];
        }
        y2_ptr[m][i] += DampingFactor_ * (x_ptr[m][i] - dtemp) * diag;
      }
    }

    // ptrdiff_t is the same size as size_t, but is signed.  Being
    // signed is important so that i >= 0 is not trivially true.
    for (ptrdiff_t i = NumMyRows_  - 1; i >= 0; --i) {
      const scalar_type diag = d_ptr[i];
      size_t NumEntries;
      A_->getLocalRowCopy (as<local_ordinal_type> (i), Indices (), Values (), NumEntries);

      for (size_t m = 0; m < NumVectors; ++m) {
        scalar_type dtemp = STS::zero ();
        for (size_t k = 0; k < NumEntries; ++k) {
          const local_ordinal_type col = Indices[k];
          dtemp += Values[k] * y2_ptr[m][col];
        }
        y2_ptr[m][i] += DampingFactor_ * (x_ptr[m][i] - dtemp) * diag;
      }
    }

    // FIXME (mfh 02 Jan 2013) This is only correct if row Map == range Map.
    if (IsParallel_) {
      Y = *Y2;
    }
  }

  // See flop count discussion in implementation of ApplyInverseSGS_CrsMatrix().
  const size_t dampingFlops = (DampingFactor_ == STS::one()) ? 0 : 1;
  const size_t numVectors = X.getNumVectors ();

  ApplyFlops_ += 2 * NumSweeps_ * numVectors *
    (2 * NumGlobalRows_ + 2 * NumGlobalNonzeros_ + dampingFlops);
}

//==========================================================================
template<class MatrixType>
void
Relaxation<MatrixType>::
ApplyInverseSGS_CrsMatrix (const MatrixType& A,
                           const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
                           Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y) const
{
  typedef Teuchos::ScalarTraits<scalar_type> STS;

  const Tpetra::ESweepDirection direction = Tpetra::Symmetric;
  A.gaussSeidelCopy (Y, X, *Diagonal_, DampingFactor_, direction, NumSweeps_);

  // For each column of output, for each sweep over the matrix:
  //
  // - One + and one * for each matrix entry
  // - One / and one + for each row of the matrix
  // - If the damping factor is not one: one * for each row of the
  //   matrix.  (It's not fair to count this if the damping factor is
  //   one, since the implementation could skip it.  Whether it does
  //   or not is the implementation's choice.)
  //
  // Each sweep of symmetric Gauss-Seidel / SOR counts as two sweeps,
  // one forward and one backward.

  // Floating-point operations due to the damping factor, per matrix
  // row, per direction, per columm of output.
  const size_t dampingFlops = (DampingFactor_ == STS::one()) ? 0 : 1;
  const size_t numVectors = X.getNumVectors ();

  ApplyFlops_ += 2 * NumSweeps_ * numVectors *
    (2 * NumGlobalRows_ + 2 * NumGlobalNonzeros_ + dampingFlops);
}

//==========================================================================
template<class MatrixType>
std::string Relaxation<MatrixType>::description() const {
  using Teuchos::TypeNameTraits;
  std::ostringstream os;

  std::string status;
  if (isInitialized ()) {
    status = "initialized";
    if (isComputed ()) {
      status += ", computed";
    } else {
      status += ", not computed";
    }
  } else {
    status = "not initialized";
  }

  std::string type;
  if (PrecType_ == Ifpack2::JACOBI) {
    type = "Jacobi";
  } else if (PrecType_ == Ifpack2::GS) {
    type = "Gauss-Seidel";
  } else if (PrecType_ == Ifpack2::SGS) {
    type = "Symmetric Gauss-Seidel";
  } else {
    type = "INVALID";
  }

  // Output is a valid YAML dictionary in flow style.  If you don't
  // like everything on a single line, you should call describe()
  // instead.
  os << "\"Ifpack2::Relaxation\": { "
     << "MatrixType: \"" << TypeNameTraits<MatrixType>::name () << "\", "
     << "Status: " << status << ", "
     << "\"relaxation: type\": " << type << ", "
     << "\"relaxation: sweeps\": " << NumSweeps_ << ", "
     << "\"relaxation: damping factor\": " << DampingFactor_ << ", ";
  if (DoL1Method_) {
    os << "\"relaxation: use l1\": " << DoL1Method_ << ", "
       << "\"relaxation: l1 eta\": " << L1Eta_ << ", ";
  }
  os << "\"Global number of rows\": " << A_->getGlobalNumRows () << ", "
     << "\"Global number of columns\": " << A_->getGlobalNumCols ()
     << " }";
  return os.str ();
}

//==========================================================================
template<class MatrixType>
void
Relaxation<MatrixType>::
describe (Teuchos::FancyOStream &out,
          const Teuchos::EVerbosityLevel verbLevel) const
{
  using Teuchos::OSTab;
  using Teuchos::TypeNameTraits;
  using Teuchos::VERB_DEFAULT;
  using Teuchos::VERB_NONE;
  using Teuchos::VERB_LOW;
  using Teuchos::VERB_MEDIUM;
  using Teuchos::VERB_HIGH;
  using Teuchos::VERB_EXTREME;
  using std::endl;
  typedef Teuchos::ScalarTraits<scalar_type> STS;

  const Teuchos::EVerbosityLevel vl =
    (verbLevel == VERB_DEFAULT) ? VERB_LOW : verbLevel;

  const int myRank = Comm_->getRank ();

  //    none: print nothing
  //     low: print O(1) info from Proc 0
  //  medium:
  //    high:
  // extreme:

  if (vl != VERB_NONE && myRank == 0) {
    // Describable interface asks each implementation to start with a tab.
    OSTab tab1 (out);
    // Output is valid YAML; hence the quotes, to protect the colons.
    out << "\"Ifpack2::Relaxation\":" << endl;
    OSTab tab2 (out);
    out << "MatrixType: \"" << TypeNameTraits<MatrixType>::name () << "\"" << endl
        << "Label: " << this->getObjectLabel () << endl
        << "Parameters: " << endl;
    {
      OSTab tab3 (out);
      out << "\"relaxation: type\": ";
      if (PrecType_ == Ifpack2::JACOBI) {
        out << "Jacobi";
      } else if (PrecType_ == Ifpack2::GS) {
        out << "Gauss-Seidel";
      } else if (PrecType_ == Ifpack2::SGS) {
        out << "Symmetric Gauss-Seidel";
      } else {
        out << "INVALID";
      }
      out << endl
          << "\"relaxation: sweeps\": " << NumSweeps_ << endl
          << "\"relaxation: damping factor\": " << DampingFactor_ << endl
          << "\"relaxation: min diagonal value\": " << MinDiagonalValue_ << endl
          << "\"relaxation: zero starting solution\": " << ZeroStartingSolution_ << endl
          << "\"relaxation: backward mode\": " << DoBackwardGS_ << endl
          << "\"relaxation: use l1\": " << DoL1Method_ << endl
          << "\"relaxation: l1 eta\": " << L1Eta_ << endl;
    }
    out << "Computed quantities: " << endl;
    {
      OSTab tab3 (out);
      out << "initialized: " << (isInitialized () ? "true" : "false") << endl
          << "computed: " << (isComputed () ? "true" : "false") << endl
          << "Condition number estimate: " << Condest_ << endl
          << "Global number of rows: " << A_->getGlobalNumRows () << endl
          << "Global number of columns: " << A_->getGlobalNumCols () << endl;
    }
    out << "Call counts and total times (in seconds): " << endl;
    {
      OSTab tab3 (out);
      out << "initialize: " << endl;
      {
        OSTab tab4 (out);
        out << "Call count: " << NumInitialize_ << endl;
        out << "Total time: " << InitializeTime_ << endl;
      }
      out << "compute: " << endl;
      {
        OSTab tab4 (out);
        out << "Call count: " << NumCompute_ << endl;
        out << "Total time: " << ComputeTime_ << endl;
      }
      out << "apply: " << endl;
      {
        OSTab tab4 (out);
        out << "Call count: " << NumApply_ << endl;
        out << "Total time: " << ApplyTime_ << endl;
      }
    }
  }
}

} // namespace Ifpack2

#endif // IFPACK2_RELAXATION_DEF_HPP

