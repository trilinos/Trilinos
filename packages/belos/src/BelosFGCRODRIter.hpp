// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef BELOS_FGCRODR_ITER_HPP
#define BELOS_FGCRODR_ITER_HPP

/*! \file BelosFGCRODRIter.hpp
    \brief Belos concrete class for performing the flexible GCRO-DR iteration.
*/

#include "BelosGCRODRIter.hpp"

namespace Belos {

class FGCRODRIterInitFailure : public GCRODRIterInitFailure {
public:
  FGCRODRIterInitFailure(const std::string& what_arg) :
    GCRODRIterInitFailure(what_arg) {}
};

class FGCRODRIterOrthoFailure : public GCRODRIterOrthoFailure {
public:
  FGCRODRIterOrthoFailure(const std::string& what_arg) :
    GCRODRIterOrthoFailure(what_arg) {}
};

template<class ScalarType, class MV, class OP, class DM = DefaultDenseMatrix<int, ScalarType> >
class FGCRODRIter : virtual public GCRODRIteration<ScalarType,MV,OP,DM> {
public:
  typedef MultiVecTraits<ScalarType,MV,DM> MVT;
  typedef OperatorTraits<ScalarType,MV,OP> OPT;
  typedef DenseMatTraits<ScalarType,DM>    DMT;
  typedef Teuchos::ScalarTraits<ScalarType> SCT;
  typedef typename SCT::magnitudeType MagnitudeType;

  FGCRODRIter(const Teuchos::RCP<LinearProblem<ScalarType,MV,OP,DM> > &problem,
              const Teuchos::RCP<OutputManager<ScalarType> > &printer,
              const Teuchos::RCP<StatusTest<ScalarType,MV,OP,DM> > &tester,
              const Teuchos::RCP<MatOrthoManager<ScalarType,MV,OP,DM> > &ortho,
              Teuchos::ParameterList &params);

  virtual ~FGCRODRIter() {}

  void iterate();

  void initialize(GCRODRIterState<ScalarType,MV,DM>& newstate);

  void initialize() {
    GCRODRIterState<ScalarType,MV,DM> empty;
    initialize(empty);
  }

  GCRODRIterState<ScalarType,MV,DM> getState() const {
    GCRODRIterState<ScalarType,MV,DM> state;
    state.curDim = curDim_;
    state.V = V_;
    state.Z = Z_;
    state.U = U_;
    state.C = C_;
    state.H2 = H2_;
    state.H = H_;
    state.B = B_;
    return state;
  }

  int getNumIters() const { return iter_; }

  void resetNumIters(int iter = 0) { iter_ = iter; }

  Teuchos::RCP<const MV>
  getNativeResiduals(std::vector<MagnitudeType> *norms) const;

  Teuchos::RCP<MV> getCurrentUpdate() const;

  void updateLSQR(int dim = -1);

  int getCurSubspaceDim() const {
    if (!initialized_) return 0;
    return curDim_;
  }

  int getMaxSubspaceDim() const { return numBlocks_; }

  const LinearProblem<ScalarType,MV,OP,DM>& getProblem() const {
    return *lp_;
  }

  int getNumBlocks() const { return numBlocks_; }

  void setNumBlocks(int numBlocks) {
    setSize(recycledBlocks_, numBlocks);
  }

  int getBlockSize() const { return 1; }

  void setBlockSize(int blockSize) {
    TEUCHOS_TEST_FOR_EXCEPTION(
      blockSize != 1,
      std::invalid_argument,
      "Belos::FGCRODRIter::setBlockSize(): Cannot use a block size that is not one.");
  }

  void setSize(int recycledBlocks, int numBlocks) {
    if (recycledBlocks_ != recycledBlocks)
      recycledBlocks_ = recycledBlocks;

    if (numBlocks_ != numBlocks) {
      numBlocks_ = numBlocks;
      cs_.resize(numBlocks_ + 1);
      sn_.resize(numBlocks_ + 1);
      z_ = DMT::Create(numBlocks_ + 1, 1, false);
      R_ = DMT::Create(numBlocks_ + 1, numBlocks_, false);
    }
  }

  bool isInitialized() { return initialized_; }

private:
  const Teuchos::RCP<LinearProblem<ScalarType,MV,OP,DM> > lp_;
  const Teuchos::RCP<OutputManager<ScalarType> > om_;
  const Teuchos::RCP<StatusTest<ScalarType,MV,OP,DM> > stest_;
  const Teuchos::RCP<OrthoManager<ScalarType,MV,DM> > ortho_;

  int numBlocks_;
  int recycledBlocks_;

  std::vector<ScalarType> sn_;
  std::vector<MagnitudeType> cs_;

  bool initialized_;
  int curDim_, iter_;
  int ptrH00_;

  Teuchos::RCP<MV> V_;
  Teuchos::RCP<MV> Z_;
  Teuchos::RCP<MV> U_, C_;

  Teuchos::RCP<DM> H2_;
  Teuchos::RCP<DM> H_;
  Teuchos::RCP<DM> B_;

  Teuchos::RCP<DM> R_;
  Teuchos::RCP<DM> z_;
};


// Constructor
template<class ScalarType, class MV, class OP, class DM>
FGCRODRIter<ScalarType,MV,OP,DM>::
FGCRODRIter(const Teuchos::RCP<LinearProblem<ScalarType,MV,OP,DM> > &problem,
            const Teuchos::RCP<OutputManager<ScalarType> > &printer,
            const Teuchos::RCP<StatusTest<ScalarType,MV,OP,DM> > &tester,
            const Teuchos::RCP<MatOrthoManager<ScalarType,MV,OP,DM> > &ortho,
            Teuchos::ParameterList &params) :
  lp_(problem),
  om_(printer),
  stest_(tester),
  ortho_(ortho),
  numBlocks_(0),
  recycledBlocks_(0),
  initialized_(false),
  curDim_(0),
  iter_(0),
  ptrH00_(0),
  V_(Teuchos::null),
  Z_(Teuchos::null),
  U_(Teuchos::null),
  C_(Teuchos::null),
  H2_(Teuchos::null),
  H_(Teuchos::null),
  B_(Teuchos::null)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    !params.isParameter("Num Blocks"),
    std::invalid_argument,
    "Belos::FGCRODRIter::constructor: mandatory parameter \"Num Blocks\" is not specified.");
  int nb = Teuchos::getParameter<int>(params, "Num Blocks");

  TEUCHOS_TEST_FOR_EXCEPTION(
    !params.isParameter("Recycled Blocks"),
    std::invalid_argument,
    "Belos::FGCRODRIter::constructor: mandatory parameter \"Recycled Blocks\" is not specified.");
  int rb = Teuchos::getParameter<int>(params, "Recycled Blocks");

  TEUCHOS_TEST_FOR_EXCEPTION(
    nb <= 0,
    std::invalid_argument,
    "Belos::FGCRODRIter() was passed a non-positive argument for \"Num Blocks\".");

  TEUCHOS_TEST_FOR_EXCEPTION(
    rb >= nb,
    std::invalid_argument,
    "Belos::FGCRODRIter() the number of recycled blocks is larger than the allowable subspace.");

  numBlocks_ = nb;
  recycledBlocks_ = rb;

  cs_.resize(numBlocks_ + 1);
  sn_.resize(numBlocks_ + 1);
  z_ = DMT::Create(numBlocks_ + 1, 1, false);
  R_ = DMT::Create(numBlocks_ + 1, numBlocks_, false);
}


// Get current flexible update.
template<class ScalarType, class MV, class OP, class DM>
Teuchos::RCP<MV>
FGCRODRIter<ScalarType,MV,OP,DM>::getCurrentUpdate() const
{
  Teuchos::RCP<MV> currentUpdate = Teuchos::null;

  if (curDim_ == 0) {
    return currentUpdate;
  }

  const ScalarType one = SCT::one();
  const ScalarType zero = SCT::zero();

  Teuchos::BLAS<int,ScalarType> blas;

  currentUpdate = MVT::Clone(*Z_, 1);

  Teuchos::RCP<DM> y = DMT::SubviewCopy(*z_, curDim_, 1);

  DMT::SyncDeviceToHost(*y);
  DMT::SyncDeviceToHost(*R_);

  blas.TRSM(Teuchos::LEFT_SIDE,
            Teuchos::UPPER_TRI,
            Teuchos::NO_TRANS,
            Teuchos::NON_UNIT_DIAG,
            curDim_,
            1,
            one,
            DMT::GetConstRawHostPtr(*R_),
            DMT::GetStride(*R_),
            DMT::GetRawHostPtr(*y),
            DMT::GetStride(*y));

  DMT::SyncHostToDevice(*y);

  // Flexible part: update = Z(:,1:curDim) * y.
  std::vector<int> index(curDim_);
  for (int i = 0; i < curDim_; ++i) {
    index[i] = i;
  }

  Teuchos::RCP<const MV> Zjp1 = MVT::CloneView(*Z_, index);
  MVT::MvTimesMatAddMv(one, *Zjp1, *y, zero, *currentUpdate);

  // Recycle correction: update -= U * B * y.
  if (U_ != Teuchos::null) {
    Teuchos::RCP<DM> z = DMT::Create(recycledBlocks_, 1);

    DMT::SyncDeviceToHost(*H2_);

    blas.GEMM(Teuchos::NO_TRANS,
              Teuchos::NO_TRANS,
              recycledBlocks_,
              1,
              curDim_,
              one,
              DMT::GetConstRawHostPtr(*B_),
              DMT::GetStride(*B_),
              DMT::GetConstRawHostPtr(*y),
              DMT::GetStride(*y),
              zero,
              DMT::GetRawHostPtr(*z),
              DMT::GetStride(*z));

    DMT::SyncHostToDevice(*z);

    MVT::MvTimesMatAddMv(-one, *U_, *z, one, *currentUpdate);
  }

  return currentUpdate;
}


// Native residual norms.
template<class ScalarType, class MV, class OP, class DM>
Teuchos::RCP<const MV>
FGCRODRIter<ScalarType,MV,OP,DM>::
getNativeResiduals(std::vector<MagnitudeType> *norms) const
{
  if (norms && static_cast<int>(norms->size()) == 0) {
    norms->resize(1);
  }

  if (norms) {
    DMT::SyncDeviceToHost(*z_);
    const ScalarType curNativeResid = DMT::ValueConst(*z_, curDim_, 0);
    (*norms)[0] = SCT::magnitude(curNativeResid);
  }

  return Teuchos::null;
}


// Initialize.
template<class ScalarType, class MV, class OP, class DM>
void
FGCRODRIter<ScalarType,MV,OP,DM>::
initialize(GCRODRIterState<ScalarType,MV,DM>& newstate)
{
  if (newstate.V != Teuchos::null &&
      newstate.Z != Teuchos::null &&
      newstate.H2 != Teuchos::null) {
    curDim_ = newstate.curDim;
    V_ = newstate.V;
    Z_ = newstate.Z;
    U_ = newstate.U;
    C_ = newstate.C;
    H2_ = newstate.H2;

    // No recycled space; this cycle primes the recycle space.
    if (newstate.U == Teuchos::null) {
      ptrH00_ = recycledBlocks_ + 1;
      H_ = DMT::Subview(*H2_, numBlocks_ + 1, numBlocks_, ptrH00_, ptrH00_);
      B_ = Teuchos::null;
    }
    else {
      ptrH00_ = recycledBlocks_;
      H_ = DMT::Subview(*H2_, numBlocks_ + 1, numBlocks_, ptrH00_, ptrH00_);
      B_ = DMT::Subview(*H2_, recycledBlocks_, numBlocks_, 0, ptrH00_);
    }
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(
      newstate.V == Teuchos::null,
      std::invalid_argument,
      "Belos::FGCRODRIter::initialize(): GCRODRIterState does not have V initialized.");
    TEUCHOS_TEST_FOR_EXCEPTION(
      newstate.Z == Teuchos::null,
      std::invalid_argument,
      "Belos::FGCRODRIter::initialize(): GCRODRIterState does not have Z initialized.");
    TEUCHOS_TEST_FOR_EXCEPTION(
      newstate.H2 == Teuchos::null,
      std::invalid_argument,
      "Belos::FGCRODRIter::initialize(): GCRODRIterState does not have H2 initialized.");
  }

  initialized_ = true;
}


// Iterate.
template<class ScalarType, class MV, class OP, class DM>
void
FGCRODRIter<ScalarType,MV,OP,DM>::iterate()
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    initialized_ == false,
    FGCRODRIterInitFailure,
    "Belos::FGCRODRIter::iterate(): FGCRODRIter class not initialized.");

  setSize(recycledBlocks_, numBlocks_);

  Teuchos::RCP<MV> Vnext;
  Teuchos::RCP<MV> Znext;
  Teuchos::RCP<const MV> Vprev;

  std::vector<int> curind(1);

  DMT::PutScalar(*z_);

  // Orthonormalize the initial residual vector in V(:,0).
  curind[0] = 0;
  Vnext = MVT::CloneViewNonConst(*V_, curind);

  Teuchos::RCP<DM> z0 = DMT::Subview(*z_, 1, 1);
  int rank = ortho_->normalize(*Vnext, z0);

  TEUCHOS_TEST_FOR_EXCEPTION(
    rank != 1,
    FGCRODRIterOrthoFailure,
    "Belos::FGCRODRIter::iterate(): couldn't generate initial basis of full rank.");

  std::vector<int> prevind(numBlocks_ + 1);

  while (stest_->checkStatus(this) != Passed && curDim_ + 1 <= numBlocks_) {
    iter_++;

    const int lclDim = curDim_ + 1;

    // Next V basis vector storage.
    curind[0] = lclDim;
    Vnext = MVT::CloneViewNonConst(*V_, curind);

    // Current V vector.
    curind[0] = curDim_;
    Vprev = MVT::CloneView(*V_, curind);

    // Current flexible correction vector.
    Znext = MVT::CloneViewNonConst(*Z_, curind);

    // z_j = M_j(v_j).  If no right preconditioner exists,
    // LinearProblem::applyRightPrec copies v_j into z_j.
    lp_->applyRightPrec(*Vprev, *Znext);
    Vprev = Teuchos::null;

    // w = A z_j.
    lp_->applyOp(*Znext, *Vnext);
    Znext = Teuchos::null;

    if (U_ != Teuchos::null) {
      DMT::SyncHostToDevice(*H2_);

      // Project out recycled image space C and store coefficients in B.
      Teuchos::Array<Teuchos::RCP<const MV> > C(1, C_);
      Teuchos::RCP<DM> subB =
        DMT::Subview(*H2_, recycledBlocks_, 1, 0, ptrH00_ + curDim_);
      Teuchos::Array<Teuchos::RCP<DM> > AsubB(1, subB);

      ortho_->project(*Vnext, AsubB, C);
    }

    // Orthogonalize against previous Krylov basis vectors.
    prevind.resize(lclDim);
    for (int i = 0; i < lclDim; ++i) {
      prevind[i] = i;
    }

    Vprev = MVT::CloneView(*V_, prevind);
    Teuchos::Array<Teuchos::RCP<const MV> > AVprev(1, Vprev);

    Teuchos::RCP<DM> subH =
      DMT::Subview(*H2_, lclDim, 1, ptrH00_, ptrH00_ + curDim_);
    Teuchos::Array<Teuchos::RCP<DM> > AsubH(1, subH);

    Teuchos::RCP<DM> subR =
      DMT::Subview(*H2_, 1, 1, ptrH00_ + lclDim, ptrH00_ + curDim_);

    rank = ortho_->projectAndNormalize(*Vnext, AsubH, subR, AVprev);

    Teuchos::RCP<DM> subR2 =
      DMT::Subview(*R_, lclDim + 1, 1, 0, curDim_);
    Teuchos::RCP<const DM> subH2 =
      DMT::SubviewConst(*H2_, lclDim + 1, 1, ptrH00_, ptrH00_ + curDim_);

    DMT::Assign(*subR2, *subH2);

    TEUCHOS_TEST_FOR_EXCEPTION(
      rank != 1,
      FGCRODRIterOrthoFailure,
      "Belos::FGCRODRIter::iterate(): couldn't generate basis of full rank.");

    updateLSQR();

    curDim_++;
  }
}


// Update QR factorization.
template<class ScalarType, class MV, class OP, class DM>
void
FGCRODRIter<ScalarType,MV,OP,DM>::updateLSQR(int dim)
{
  int i;
  const ScalarType zero = SCT::zero();

  int curDim = curDim_;
  if ((dim >= curDim_) && (dim < getMaxSubspaceDim())) {
    curDim = dim;
  }

  Teuchos::BLAS<int, ScalarType> blas;

  DMT::SyncDeviceToHost(*R_);
  DMT::SyncDeviceToHost(*z_);

  for (i = 0; i < curDim; ++i) {
    blas.ROT(1,
             &(DMT::Value(*R_, i, curDim)),
             1,
             &(DMT::Value(*R_, i+1, curDim)),
             1,
             &cs_[i],
             &sn_[i]);
  }

  blas.ROTG(&(DMT::Value(*R_, curDim, curDim)),
            &(DMT::Value(*R_, curDim+1, curDim)),
            &cs_[curDim],
            &sn_[curDim]);

  DMT::Value(*R_, curDim+1, curDim) = zero;

  blas.ROT(1,
           &(DMT::Value(*z_, curDim, 0)),
           1,
           &(DMT::Value(*z_, curDim+1, 0)),
           1,
           &cs_[curDim],
           &sn_[curDim]);

  DMT::SyncHostToDevice(*R_);
  DMT::SyncHostToDevice(*z_);
}

} // namespace Belos

#endif // BELOS_FGCRODR_ITER_HPP
