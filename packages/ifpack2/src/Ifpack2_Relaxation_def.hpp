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
  Comm_(A->getRowMap()->getComm()),
  Time_( Teuchos::rcp( new Teuchos::Time("Ifpack2::Relaxation") ) ),
  NumSweeps_(1),
  PrecType_(Ifpack2::JACOBI),
  MinDiagonalValue_(0.0),
  DampingFactor_(1.0),
  IsParallel_(false),
  ZeroStartingSolution_(true),
  DoBackwardGS_(false),
  DoL1Method_(false),
  L1Eta_(1.5),
  Condest_(-1.0),
  IsInitialized_(false),
  IsComputed_(false),
  NumInitialize_(0),
  NumCompute_(0),
  NumApply_(0),
  InitializeTime_(0.0),
  ComputeTime_(0.0),
  ApplyTime_(0.0),
  ComputeFlops_(0.0),
  ApplyFlops_(0.0),
  NumMyRows_(0),
  NumGlobalRows_(0),
  NumGlobalNonzeros_(0)
{
  TEUCHOS_TEST_FOR_EXCEPTION(A_ == Teuchos::null, std::runtime_error,
      Teuchos::typeName(*this) << "::Relaxation(): input matrix reference was null.");
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

  if (PT == "Jacobi")
    PrecType_ = Ifpack2::JACOBI;
  else if (PT == "Gauss-Seidel")
    PrecType_ = Ifpack2::GS;
  else if (PT == "Symmetric Gauss-Seidel")
    PrecType_ = Ifpack2::SGS;
  else {
    std::ostringstream osstr;
    osstr << "Ifpack2::Relaxation::setParameters: unsupported parameter-value for 'relaxation: type' (" << PT << ")";
    std::string str = osstr.str();
    throw std::runtime_error(str);
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

  // Flops:
  // - matrix vector              (2 * NumGlobalNonzeros_)
  // - update                     (2 * NumGlobalRows_)
  // - Multiply:
  //   - DampingFactor            (NumGlobalRows_)
  //   - Diagonal                 (NumGlobalRows_)
  //   - A + B                    (NumGlobalRows_)
  //   - 1.0                      (NumGlobalRows_)
  ApplyFlops_ += NumVectors * (6 * NumGlobalRows_ + 2 * NumGlobalNonzeros_);
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
  size_t NumVectors = X.getNumVectors();

  size_t maxLength = A_->getNodeMaxNumRowEntries();
  Teuchos::Array<local_ordinal_type> Indices(maxLength);
  Teuchos::Array<scalar_type> Values(maxLength);

  Teuchos::RCP< Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> > Y2;
  if (IsParallel_)
    Y2 = Teuchos::rcp( new Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>(Importer_->getTargetMap(), NumVectors) );
  else
    Y2 = Teuchos::rcp( &Y, false );

  // extract views (for nicer and faster code)
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<scalar_type> > y_ptr = Y.get2dViewNonConst();
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<scalar_type> > y2_ptr = Y2->get2dViewNonConst();
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<const scalar_type> > x_ptr =  X.get2dView();
  Teuchos::ArrayRCP<const scalar_type> d_ptr = Diagonal_->get1dView();

  for (int j = 0; j < NumSweeps_ ; j++) {

    // data exchange is here, once per sweep
    if (IsParallel_)
      Y2->doImport(Y,*Importer_,Tpetra::INSERT);

    if(!DoBackwardGS_){
      /* Forward Mode */
      for (size_t i = 0 ; i < NumMyRows_ ; ++i) {

        size_t NumEntries;
        A_->getLocalRowCopy(i, Indices(), Values(), NumEntries);

        for (size_t m = 0 ; m < NumVectors ; ++m) {

          scalar_type dtemp = 0.0;
          for (size_t k = 0 ; k < NumEntries ; ++k) {
            local_ordinal_type col = Indices[k];
            dtemp += Values[k] * y2_ptr[m][col];
          }

          y2_ptr[m][i] += DampingFactor_ * d_ptr[i] * (x_ptr[m][i] - dtemp);
        }
      }
    }
    else {
      /* Backward Mode */
      for (int i = NumMyRows_  - 1 ; i > -1 ; --i) {
        size_t NumEntries;
        A_->getLocalRowCopy(i, Indices(), Values(), NumEntries);

        for (size_t m = 0 ; m < NumVectors ; ++m) {

          scalar_type dtemp = 0.0;
          for (size_t k = 0 ; k < NumEntries ; ++k) {
            local_ordinal_type col = Indices[k];
            dtemp += Values[k] * y2_ptr[m][col];
          }

          y2_ptr[m][i] += DampingFactor_ * d_ptr[i] * (x_ptr[m][i] - dtemp);
        }
      }
    }

    // using Export() sounded quite expensive
    if (IsParallel_) {
      for (size_t m = 0 ; m < NumVectors ; ++m)
        for (size_t i = 0 ; i < NumMyRows_ ; ++i)
          y_ptr[m][i] = y2_ptr[m][i];
    }
  }

  ApplyFlops_ += NumVectors * (4 * NumGlobalRows_ + 2 * NumGlobalNonzeros_);
}

//==========================================================================
template<class MatrixType>
void Relaxation<MatrixType>::ApplyInverseGS_CrsMatrix(
        const MatrixType& A,
        const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
              Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y) const
{
  size_t NumVectors = X.getNumVectors();

  Teuchos::ArrayView<const local_ordinal_type> Indices;
  Teuchos::ArrayView<const scalar_type> Values;

  Teuchos::RCP< Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> > Y2;
  if (IsParallel_) {
    Y2 = Teuchos::rcp( new Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>(Importer_->getTargetMap(), NumVectors) );
  }
  else
    Y2 = Teuchos::rcp( &Y, false );

  Teuchos::ArrayRCP<Teuchos::ArrayRCP<scalar_type> > y_ptr = Y.get2dViewNonConst();
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<scalar_type> > y2_ptr = Y2->get2dViewNonConst();
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<const scalar_type> > x_ptr =  X.get2dView();
  Teuchos::ArrayRCP<const scalar_type> d_ptr = Diagonal_->get1dView();

  for (int iter = 0 ; iter < NumSweeps_ ; ++iter) {

    // only one data exchange per sweep
    if (IsParallel_)
      Y2->doImport(Y,*Importer_,Tpetra::INSERT);

    if(!DoBackwardGS_){
      /* Forward Mode */
      for (size_t i = 0 ; i < NumMyRows_ ; ++i) {

        local_ordinal_type col;
        scalar_type diag = d_ptr[i];

        A.getLocalRowView(i, Indices, Values);
        size_t NumEntries = Indices.size();

        for (size_t m = 0 ; m < NumVectors ; ++m) {

          scalar_type dtemp = 0.0;

          for (size_t k = 0; k < NumEntries; ++k) {
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

        local_ordinal_type col;
        scalar_type diag = d_ptr[i];

        A.getLocalRowView(i, Indices, Values);
        size_t NumEntries = Indices.size();

        for (size_t m = 0 ; m < NumVectors ; ++m) {

          scalar_type dtemp = 0.0;
          for (size_t k = 0; k < NumEntries; ++k) {
            col = Indices[k];
            dtemp += Values[k] * y2_ptr[m][col];
          }

          y2_ptr[m][i] += DampingFactor_ * (x_ptr[m][i] - dtemp) * diag;

        }
      }
    }

    if (IsParallel_) {
      for (size_t m = 0 ; m < NumVectors ; ++m)
        for (size_t i = 0 ; i < NumMyRows_ ; ++i)
          y_ptr[m][i] = y2_ptr[m][i];
    }
  }

  ApplyFlops_ += NumVectors * (8 * NumGlobalRows_ + 4 * NumGlobalNonzeros_);
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
  // CrsMatrix specialization.
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
  size_t NumVectors = X.getNumVectors();
  size_t maxLength = A_->getNodeMaxNumRowEntries();
  Teuchos::Array<local_ordinal_type> Indices(maxLength);
  Teuchos::Array<scalar_type> Values(maxLength);

  Teuchos::RCP< Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> > Y2;
  if (IsParallel_) {
    Y2 = Teuchos::rcp( new Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>(Importer_->getTargetMap(), NumVectors) );
  }
  else
    Y2 = Teuchos::rcp( &Y, false );

  Teuchos::ArrayRCP<Teuchos::ArrayRCP<scalar_type> > y_ptr = Y.get2dViewNonConst();
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<scalar_type> > y2_ptr = Y2->get2dViewNonConst();
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<const scalar_type> > x_ptr =  X.get2dView();
  Teuchos::ArrayRCP<const scalar_type> d_ptr = Diagonal_->get1dView();

  for (int iter = 0 ; iter < NumSweeps_ ; ++iter) {

    // only one data exchange per sweep
    if (IsParallel_)
      Y2->doImport(Y,*Importer_,Tpetra::INSERT);

    for (size_t i = 0 ; i < NumMyRows_ ; ++i) {

      size_t NumEntries;
      scalar_type diag = d_ptr[i];

      A_->getLocalRowCopy(i, Indices(), Values(), NumEntries);

      for (size_t m = 0 ; m < NumVectors ; ++m) {

        scalar_type dtemp = 0.0;

        for (size_t k = 0 ; k < NumEntries ; ++k) {
          local_ordinal_type col = Indices[k];
          dtemp += Values[k] * y2_ptr[m][col];
        }

        y2_ptr[m][i] += DampingFactor_ * (x_ptr[m][i] - dtemp) * diag;
      }
    }

    for (int i = NumMyRows_  - 1 ; i > -1 ; --i) {

      size_t NumEntries;
      scalar_type diag = d_ptr[i];

      A_->getLocalRowCopy(i, Indices(), Values(), NumEntries);

      for (size_t m = 0 ; m < NumVectors ; ++m) {

        scalar_type dtemp = Teuchos::ScalarTraits<scalar_type>::zero();
        for (size_t k = 0 ; k < NumEntries ; ++k) {
          local_ordinal_type col = Indices[k];
          dtemp += Values[k] * y2_ptr[m][col];
        }

        y2_ptr[m][i] += DampingFactor_ * (x_ptr[m][i] - dtemp) * diag;
      }
    }

    if (IsParallel_) {
      for (size_t m = 0 ; m < NumVectors ; ++m)
        for (size_t i = 0 ; i < NumMyRows_ ; ++i)
          y_ptr[m][i] = y2_ptr[m][i];
    }
  }

  ApplyFlops_ += NumVectors * (8 * NumGlobalRows_ + 4 * NumGlobalNonzeros_);
}

//==========================================================================
template<class MatrixType>
void Relaxation<MatrixType>::ApplyInverseSGS_CrsMatrix(
        const MatrixType& A,
        const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
              Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y) const
{
  size_t NumVectors = X.getNumVectors();

  Teuchos::ArrayView<const local_ordinal_type> Indices;
  Teuchos::ArrayView<const scalar_type> Values;

  Teuchos::RCP< Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> > Y2;
  if (IsParallel_) {
    Y2 = Teuchos::rcp( new Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>(Importer_->getTargetMap(), NumVectors) );
  }
  else
    Y2 = Teuchos::rcp( &Y, false );

  Teuchos::ArrayRCP<Teuchos::ArrayRCP<scalar_type> > y_ptr = Y.get2dViewNonConst();
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<scalar_type> > y2_ptr = Y2->get2dViewNonConst();
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<const scalar_type> > x_ptr =  X.get2dView();
  Teuchos::ArrayRCP<const scalar_type> d_ptr = Diagonal_->get1dView();

  for (int iter = 0 ; iter < NumSweeps_ ; ++iter) {

    // only one data exchange per sweep
    if (IsParallel_)
      Y2->doImport(Y,*Importer_,Tpetra::INSERT);


    for (size_t i = 0 ; i < NumMyRows_ ; ++i) {

      scalar_type diag = d_ptr[i];

      A.getLocalRowView(i, Indices, Values);
      size_t NumEntries = Indices.size();

      for (size_t m = 0 ; m < NumVectors ; ++m) {

        scalar_type dtemp = Teuchos::ScalarTraits<scalar_type>::zero();

        for (size_t k = 0; k < NumEntries; ++k) {

          local_ordinal_type col = Indices[k];
          dtemp += Values[k] * y2_ptr[m][col];
        }

        y2_ptr[m][i] += DampingFactor_ * (x_ptr[m][i] - dtemp) * diag;
      }
    }

    for (int i = NumMyRows_  - 1 ; i > -1 ; --i) {

      scalar_type diag = d_ptr[i];

      A.getLocalRowView(i, Indices, Values);
      size_t NumEntries = Indices.size();

      for (size_t m = 0 ; m < NumVectors ; ++m) {

        scalar_type dtemp = Teuchos::ScalarTraits<scalar_type>::zero();
        for (size_t k = 0; k < NumEntries; ++k) {

          local_ordinal_type col = Indices[k];
          dtemp += Values[k] * y2_ptr[m][col];
        }

        y2_ptr[m][i] += DampingFactor_ * (x_ptr[m][i] - dtemp) * diag;

      }
    }

    if (IsParallel_) {
      for (size_t m = 0 ; m < NumVectors ; ++m)
        for (size_t i = 0 ; i < NumMyRows_ ; ++i)
          y_ptr[m][i] = y2_ptr[m][i];
    }
  }

  ApplyFlops_ += NumVectors * (8 * NumGlobalRows_ + 4 * NumGlobalNonzeros_);
}

//==========================================================================
template<class MatrixType>
std::string Relaxation<MatrixType>::description() const {
  std::ostringstream oss;
  oss << Teuchos::Describable::description();
  if (isInitialized()) {
    if (isComputed()) {
      oss << "{status = initialized, computed";
    }
    else {
      oss << "{status = initialized, not computed";
    }
  }
  else {
    oss << "{status = not initialized, not computed";
  }
  //
  if (PrecType_ == Ifpack2::JACOBI)   oss << "Type = Jacobi, " << std::endl;
  else if (PrecType_ == Ifpack2::GS)  oss << "Type = Gauss-Seidel, " << std::endl;
  else if (PrecType_ == Ifpack2::SGS) oss << "Type = Sym. Gauss-Seidel, " << std::endl;
  //
  oss << ", global rows = " << A_->getGlobalNumRows()
      << ", global cols = " << A_->getGlobalNumCols();

  if (DoL1Method_)
    oss<<", using L1 correction with eta "<<L1Eta_;

  oss << "}";
  return oss.str();
}

//==========================================================================
template<class MatrixType>
void Relaxation<MatrixType>::describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const {
  using std::endl;
  using std::setw;
  using Teuchos::VERB_DEFAULT;
  using Teuchos::VERB_NONE;
  using Teuchos::VERB_LOW;
  using Teuchos::VERB_MEDIUM;
  using Teuchos::VERB_HIGH;
  using Teuchos::VERB_EXTREME;
  typedef Teuchos::ScalarTraits<scalar_type> STS;

  Teuchos::EVerbosityLevel vl = verbLevel;
  if (vl == VERB_DEFAULT) vl = VERB_LOW;
  const int myImageID = Comm_->getRank();
  Teuchos::OSTab tab(out);

  scalar_type MinVal = STS::zero();
  scalar_type MaxVal = STS::zero();

  if (IsComputed_) {
    Teuchos::ArrayRCP<scalar_type> DiagView = Diagonal_->get1dViewNonConst();
    scalar_type myMinVal = DiagView[0];
    scalar_type myMaxVal = DiagView[0];
    for(typename Teuchos::ArrayRCP<scalar_type>::size_type i=0; i<DiagView.size(); ++i) {
      if (STS::magnitude(myMinVal) > STS::magnitude(DiagView[i])) myMinVal = DiagView[i];
      if (STS::magnitude(myMaxVal) < STS::magnitude(DiagView[i])) myMaxVal = DiagView[i];
    }

    Teuchos::reduceAll(*Comm_, Teuchos::REDUCE_MIN, 1, &myMinVal, &MinVal);
    Teuchos::reduceAll(*Comm_, Teuchos::REDUCE_MAX, 1, &myMaxVal, &MaxVal);
  }

  //    none: print nothing
  //     low: print O(1) info from node 0
  //  medium:
  //    high:
  // extreme:
  if (vl != VERB_NONE && myImageID == 0) {
    out << this->description() << endl;
    out << endl;
    out << "===============================================================================" << endl;
    out << "Sweeps         = " << NumSweeps_ << endl;
    out << "damping factor = " << DampingFactor_ << endl;
    if (PrecType_ == Ifpack2::GS && DoBackwardGS_) {
      out << "Using backward mode (GS only)" << endl;
    }
    if   (ZeroStartingSolution_) { out << "Using zero starting solution" << endl; }
    else                         { out << "Using input starting solution" << endl; }
    if   (Condest_ == -1.0) { out << "Condition number estimate       = N/A" << endl; }
    else                    { out << "Condition number estimate       = " << Condest_ << endl; }
    if (IsComputed_) {
      out << "Minimum value on stored diagonal = " << MinVal << endl;
      out << "Maximum value on stored diagonal = " << MaxVal << endl;
    }
    out << endl;
    out << "Phase           # calls    Total Time (s)     Total MFlops      MFlops/s       " << endl;
    out << "------------    -------    ---------------    ---------------   ---------------" << endl;
    out << setw(12) << "initialize()" << setw(5) << getNumInitialize() << "    " << setw(15) << getInitializeTime() << endl;
    out << setw(12) << "compute()" << setw(5) << getNumCompute()    << "    " << setw(15) << getComputeTime() << "    "
        << setw(15) << getComputeFlops() << "    "
        << setw(15) << (getComputeTime() != 0.0 ? getComputeFlops() / getComputeTime() * 1.0e-6 : 0.0) << endl;
    out << setw(12) << "apply()" << setw(5) << getNumApply()    << "    " << setw(15) << getApplyTime() << "    "
        << setw(15) << getApplyFlops() << "    "
        << setw(15) << (getApplyTime() != 0.0 ? getApplyFlops() / getApplyTime() * 1.0e-6 : 0.0) << endl;
    out << "===============================================================================" << endl;
    out << endl;
  }
}

}//namespace Ifpack2

#endif // IFPACK2_RELAXATION_DEF_HPP

