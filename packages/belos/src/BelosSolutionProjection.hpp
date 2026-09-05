// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef BELOS_SOLUTION_PROJECTION_HPP
#define BELOS_SOLUTION_PROJECTION_HPP

/*! \file BelosSolutionProjection.hpp
    \brief Utility implementing Fischer's RHS / solution projection scheme.

    This utility stores paired bases U and C satisfying

      A U = C,   C^H C = I,

    and applies the projection

      x <- x + U C^H r,
      r <- r - C C^H r,

    where r = b - A x.

    This is Fischer's Method 1 when U stores previous solution-like
    vectors and C stores their corresponding operator images.

    The projection tolerance is relative: a candidate image vector c is
    rejected if, after orthogonalization against the current C basis,

      ||c_orth|| <= projectionTol * ||c_original||.

    The default replacement strategy is "Restart": if the basis is full
    and a new vector is accepted, the previous basis is discarded and the
    new vector starts the next basis.  Thus, when the history space is
    full, the next accepted vector "kicks off" the new space using the
    most recent information.
*/

#include "BelosConfigDefs.hpp"
#include "BelosTypes.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosMultiVecTraits.hpp"
#include "BelosOperatorTraits.hpp"
#include "BelosTeuchosDenseAdapter.hpp"
#include "BelosKokkosDenseAdapter.hpp"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ScalarTraits.hpp"

#include <string>
#include <vector>

namespace Belos {

class SolutionProjectionFailure : public BelosError {
public:
  SolutionProjectionFailure(const std::string& what_arg) :
    BelosError(what_arg) {}
};

/*! \class Belos::SolutionProjection
    \brief Fischer Method 1 style projection utility.

    The class stores a correction basis U and an orthonormal image basis C
    satisfying A U = C.  The projection is minimum-residual over range(U).

    This utility is intended to be used as a black-box initial-guess
    generator for sequences of related linear systems.
*/
template<class ScalarType, class MV, class OP,
         class DM = DefaultDenseMatrix<int, ScalarType> >
class SolutionProjection {
public:
  typedef MultiVecTraits<ScalarType,MV,DM> MVT;
  typedef OperatorTraits<ScalarType,MV,OP> OPT;
  typedef DenseMatTraits<ScalarType,DM> DMT;
  typedef Teuchos::ScalarTraits<ScalarType> SCT;
  typedef typename SCT::magnitudeType MagnitudeType;
  typedef Teuchos::ScalarTraits<MagnitudeType> MT;

  //! Constructor with explicit maximum basis size.
  SolutionProjection(const Teuchos::RCP<const OP>& A,
                     const int maxBasisSize) :
    A_(A),
    maxBasisSize_(maxBasisSize),
    curBasisSize_(0),
    replacementStrategy_("Restart"),
    projectionTol_(static_cast<MagnitudeType>(100) * MT::eps()),
    reorthogonalizationSteps_(2),
    U_(Teuchos::null),
    C_(Teuchos::null)
  {
    validateParameters();
  }

  //! Constructor with parameter list.
  SolutionProjection(const Teuchos::RCP<const OP>& A,
                     const Teuchos::RCP<Teuchos::ParameterList>& params) :
    A_(A),
    maxBasisSize_(20),
    curBasisSize_(0),
    replacementStrategy_("Restart"),
    projectionTol_(static_cast<MagnitudeType>(100) * MT::eps()),
    reorthogonalizationSteps_(2),
    U_(Teuchos::null),
    C_(Teuchos::null)
  {
    if (!params.is_null()) {
      maxBasisSize_ =
        params->get("Maximum Basis Size", maxBasisSize_);
      replacementStrategy_ =
        params->get("Replacement Strategy", replacementStrategy_);
      projectionTol_ =
        params->get("Projection Tolerance", projectionTol_);
      reorthogonalizationSteps_ =
        params->get("Reorthogonalization Steps", reorthogonalizationSteps_);
    }

    validateParameters();
  }

  //! Set or replace the operator.
  /*!
      Changing the operator invalidates the stored relation A U = C, so
      the basis is reset.
  */
  void setOperator(const Teuchos::RCP<const OP>& A) {
    A_ = A;
    reset();
  }

  //! Clear the stored projection basis.
  void reset() {
    curBasisSize_ = 0;
  }

  //! Current number of stored basis vectors.
  int getBasisSize() const {
    return curBasisSize_;
  }

  //! Maximum number of stored basis vectors.
  int getMaxBasisSize() const {
    return maxBasisSize_;
  }

  //! Set replacement strategy.
  void setReplacementStrategy(const std::string& strategy) {
    replacementStrategy_ = strategy;
    validateParameters();
  }

  //! Get replacement strategy.
  std::string getReplacementStrategy() const {
    return replacementStrategy_;
  }

  //! Set relative projection tolerance.
  void setProjectionTolerance(const MagnitudeType tol) {
    projectionTol_ = tol;
    validateParameters();
  }

  //! Get relative projection tolerance.
  MagnitudeType getProjectionTolerance() const {
    return projectionTol_;
  }

  //! Set number of reorthogonalization passes.
  void setReorthogonalizationSteps(const int steps) {
    reorthogonalizationSteps_ = steps;
    validateParameters();
  }

  //! Get number of reorthogonalization passes.
  int getReorthogonalizationSteps() const {
    return reorthogonalizationSteps_;
  }

  //! Get correction basis U.
  Teuchos::RCP<const MV> getCorrectionBasis() const {
    if (curBasisSize_ == 0 || U_.is_null()) {
      return Teuchos::null;
    }

    std::vector<int> ind(curBasisSize_);
    for (int i = 0; i < curBasisSize_; ++i) {
      ind[i] = i;
    }

    return MVT::CloneView(*U_, ind);
  }

  //! Get image basis C.
  Teuchos::RCP<const MV> getImageBasis() const {
    if (curBasisSize_ == 0 || C_.is_null()) {
      return Teuchos::null;
    }

    std::vector<int> ind(curBasisSize_);
    for (int i = 0; i < curBasisSize_; ++i) {
      ind[i] = i;
    }

    return MVT::CloneView(*C_, ind);
  }

  //! Add solution vector(s) u.  The class computes c = A u.
  int addSolution(const MV& u) {
    TEUCHOS_TEST_FOR_EXCEPTION(
      A_.is_null(),
      SolutionProjectionFailure,
      "Belos::SolutionProjection::addSolution(): operator is null.");

    const int numVecs = MVT::GetNumberVecs(u);

    TEUCHOS_TEST_FOR_EXCEPTION(
      numVecs <= 0,
      std::invalid_argument,
      "Belos::SolutionProjection::addSolution(): input multivector must contain at least one vector.");

    Teuchos::RCP<MV> c = MVT::Clone(u, numVecs);

    OPT::Apply(*A_, u, *c);

    return addPair(u, *c);
  }

  //! Add pair(s) u and c, where c should equal A u.
  int addPair(const MV& u, const MV& c) {
    const int numVecs = MVT::GetNumberVecs(u);

    TEUCHOS_TEST_FOR_EXCEPTION(
      numVecs <= 0,
      std::invalid_argument,
      "Belos::SolutionProjection::addPair(): input multivector must contain at least one vector.");

    TEUCHOS_TEST_FOR_EXCEPTION(
      MVT::GetNumberVecs(c) != numVecs,
      std::invalid_argument,
      "Belos::SolutionProjection::addPair(): u and c must have the same number of vectors.");

    TEUCHOS_TEST_FOR_EXCEPTION(
      MVT::GetGlobalLength(u) != MVT::GetGlobalLength(c),
      std::invalid_argument,
      "Belos::SolutionProjection::addPair(): u and c must have the same global length.");

    if (maxBasisSize_ == 0) {
      return 0;
    }

    ensureStorage(u, c);

    int accepted = 0;

    for (int j = 0; j < numVecs; ++j) {
      std::vector<int> srcInd(1);
      srcInd[0] = j;

      Teuchos::RCP<const MV> uSrc = MVT::CloneView(u, srcInd);
      Teuchos::RCP<const MV> cSrc = MVT::CloneView(c, srcInd);

      Teuchos::RCP<MV> uNew = MVT::CloneCopy(*uSrc);
      Teuchos::RCP<MV> cNew = MVT::CloneCopy(*cSrc);

      const bool ok = addOnePair(*uNew, *cNew);
      if (ok) {
        ++accepted;
      }
    }

    return accepted;
  }

  //! Apply projection to x and r, where r is assumed to be b - A x.
  void project(MV& x, MV& r) const {
    if (curBasisSize_ == 0) {
      return;
    }

    validateProjectionInputs(x, r);

    const ScalarType one = SCT::one();

    const int numRhs = MVT::GetNumberVecs(r);

    std::vector<int> ind(curBasisSize_);
    for (int i = 0; i < curBasisSize_; ++i) {
      ind[i] = i;
    }

    Teuchos::RCP<const MV> Ucur = MVT::CloneView(*U_, ind);
    Teuchos::RCP<const MV> Ccur = MVT::CloneView(*C_, ind);

    Teuchos::RCP<DM> gamma = DMT::Create(curBasisSize_, numRhs);

    // gamma = C^H r.
    MVT::MvTransMv(one, *Ccur, r, *gamma);

    // x <- x + U gamma.
    MVT::MvTimesMatAddMv(one, *Ucur, *gamma, one, x);

    // r <- r - C gamma.
    MVT::MvTimesMatAddMv(-one, *Ccur, *gamma, one, r);
  }

  //! Compute r = b - A x, then apply projection to x and r.
  void project(MV& x,
               const MV& b,
               Teuchos::RCP<MV> projectedResidual = Teuchos::null) const {
    TEUCHOS_TEST_FOR_EXCEPTION(
      A_.is_null(),
      SolutionProjectionFailure,
      "Belos::SolutionProjection::project(): operator is null.");

    validateProjectionInputs(x, b);

    if (!projectedResidual.is_null()) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        MVT::GetNumberVecs(*projectedResidual) != MVT::GetNumberVecs(b),
        std::invalid_argument,
        "Belos::SolutionProjection::project(): projectedResidual must have the same number of vectors as b.");

      TEUCHOS_TEST_FOR_EXCEPTION(
        MVT::GetGlobalLength(*projectedResidual) != MVT::GetGlobalLength(b),
        std::invalid_argument,
        "Belos::SolutionProjection::project(): projectedResidual must have the same global length as b.");
    }

    Teuchos::RCP<MV> r = MVT::Clone(b, MVT::GetNumberVecs(b));

    OPT::Apply(*A_, x, *r);
    MVT::MvAddMv(-SCT::one(), *r, SCT::one(), b, *r);

    project(x, *r);

    if (!projectedResidual.is_null()) {
      MVT::Assign(*r, *projectedResidual);
    }
  }

  //! Project the current solution in a LinearProblem.
  void project(LinearProblem<ScalarType,MV,OP,DM>& problem,
               Teuchos::RCP<MV> projectedResidual = Teuchos::null) const {
    Teuchos::RCP<MV> x = problem.getCurrLHSVec();

    TEUCHOS_TEST_FOR_EXCEPTION(
      x.is_null(),
      SolutionProjectionFailure,
      "Belos::SolutionProjection::project(LinearProblem): current LHS is null. "
      "Did you call setLSIndex() or otherwise set the current linear system?");

    Teuchos::RCP<const MV> b = problem.getCurrRHSVec();

    TEUCHOS_TEST_FOR_EXCEPTION(
      b.is_null(),
      SolutionProjectionFailure,
      "Belos::SolutionProjection::project(LinearProblem): current RHS is null. "
      "Did you call setLSIndex() or otherwise set the current linear system?");

    if (!projectedResidual.is_null()) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        MVT::GetNumberVecs(*projectedResidual) != MVT::GetNumberVecs(*b),
        std::invalid_argument,
        "Belos::SolutionProjection::project(LinearProblem): projectedResidual has wrong number of vectors.");

      TEUCHOS_TEST_FOR_EXCEPTION(
        MVT::GetGlobalLength(*projectedResidual) != MVT::GetGlobalLength(*b),
        std::invalid_argument,
        "Belos::SolutionProjection::project(LinearProblem): projectedResidual has wrong global length.");
    }

    Teuchos::RCP<MV> r = MVT::Clone(*b, MVT::GetNumberVecs(*b));

    // True residual, not preconditioned residual.
    problem.computeCurrResVec(&*r);

    project(*x, *r);

    if (!projectedResidual.is_null()) {
      MVT::Assign(*r, *projectedResidual);
    }
  }

private:
  void validateParameters() const {
    TEUCHOS_TEST_FOR_EXCEPTION(
      maxBasisSize_ < 0,
      std::invalid_argument,
      "Belos::SolutionProjection: maximum basis size must be nonnegative.");

    TEUCHOS_TEST_FOR_EXCEPTION(
      projectionTol_ < MT::zero(),
      std::invalid_argument,
      "Belos::SolutionProjection: projection tolerance must be nonnegative.");

    TEUCHOS_TEST_FOR_EXCEPTION(
      reorthogonalizationSteps_ < 1,
      std::invalid_argument,
      "Belos::SolutionProjection: reorthogonalization steps must be at least one.");

    TEUCHOS_TEST_FOR_EXCEPTION(
      replacementStrategy_ != "Restart",
      std::invalid_argument,
      "Belos::SolutionProjection: currently only \"Restart\" replacement strategy is implemented.");
  }

  void validateProjectionInputs(const MV& x, const MV& r) const {
    TEUCHOS_TEST_FOR_EXCEPTION(
      MVT::GetNumberVecs(x) != MVT::GetNumberVecs(r),
      std::invalid_argument,
      "Belos::SolutionProjection::project(): x and residual/right-hand side must have the same number of vectors.");

    TEUCHOS_TEST_FOR_EXCEPTION(
      MVT::GetGlobalLength(x) != MVT::GetGlobalLength(r),
      std::invalid_argument,
      "Belos::SolutionProjection::project(): x and residual/right-hand side must have the same global length.");

    TEUCHOS_TEST_FOR_EXCEPTION(
      !U_.is_null() && MVT::GetGlobalLength(*U_) != MVT::GetGlobalLength(x),
      std::invalid_argument,
      "Belos::SolutionProjection::project(): x has incompatible global length with stored basis.");

    TEUCHOS_TEST_FOR_EXCEPTION(
      !C_.is_null() && MVT::GetGlobalLength(*C_) != MVT::GetGlobalLength(r),
      std::invalid_argument,
      "Belos::SolutionProjection::project(): residual/right-hand side has incompatible global length with stored basis.");
  }

  //! Allocate U_ and C_ if necessary.
  void ensureStorage(const MV& uPrototype, const MV& cPrototype) {
    if (U_.is_null()) {
      U_ = MVT::Clone(uPrototype, maxBasisSize_);
    }
    else if (MVT::GetNumberVecs(*U_) < maxBasisSize_) {
      Teuchos::RCP<const MV> tmp = U_;
      U_ = MVT::Clone(*tmp, maxBasisSize_);
    }

    // Important: clone C_ from the image vector cPrototype, not uPrototype.
    if (C_.is_null()) {
      C_ = MVT::Clone(cPrototype, maxBasisSize_);
    }
    else if (MVT::GetNumberVecs(*C_) < maxBasisSize_) {
      Teuchos::RCP<const MV> tmp = C_;
      C_ = MVT::Clone(*tmp, maxBasisSize_);
    }
  }

  //! Add one pair uNew, cNew, maintaining A U = C and C^H C = I.
  bool addOnePair(MV& uNew, MV& cNew) {
    const ScalarType one = SCT::one();

    std::vector<MagnitudeType> initialNorm(1);
    MVT::MvNorm(cNew, initialNorm);

    if (initialNorm[0] == MT::zero()) {
      return false;
    }

    if (curBasisSize_ == maxBasisSize_) {
      if (replacementStrategy_ == "Restart") {
        // The next accepted vector starts a new basis, using the most
        // recent information available.
        curBasisSize_ = 0;
      }
      else {
        TEUCHOS_TEST_FOR_EXCEPTION(
          true,
          SolutionProjectionFailure,
          "Belos::SolutionProjection: unsupported replacement strategy.");
      }
    }

    for (int pass = 0; pass < reorthogonalizationSteps_; ++pass) {
      if (curBasisSize_ == 0) {
        break;
      }

      std::vector<int> ind(curBasisSize_);
      for (int i = 0; i < curBasisSize_; ++i) {
        ind[i] = i;
      }

      Teuchos::RCP<const MV> Ucur = MVT::CloneView(*U_, ind);
      Teuchos::RCP<const MV> Ccur = MVT::CloneView(*C_, ind);

      Teuchos::RCP<DM> h = DMT::Create(curBasisSize_, 1);

      // h = C^H cNew.
      MVT::MvTransMv(one, *Ccur, cNew, *h);

      // cNew <- cNew - C h.
      MVT::MvTimesMatAddMv(-one, *Ccur, *h, one, cNew);

      // uNew <- uNew - U h.
      //
      // This keeps A uNew = cNew, assuming A U = C before this step.
      MVT::MvTimesMatAddMv(-one, *Ucur, *h, one, uNew);
    }

    std::vector<MagnitudeType> finalNorm(1);
    MVT::MvNorm(cNew, finalNorm);

    if (finalNorm[0] <= projectionTol_ * initialNorm[0]) {
      return false;
    }

    const ScalarType scale =
      static_cast<ScalarType>(MT::one() / finalNorm[0]);

    MVT::MvScale(cNew, scale);
    MVT::MvScale(uNew, scale);

    std::vector<int> dstInd(1);
    dstInd[0] = curBasisSize_;

    MVT::SetBlock(uNew, dstInd, *U_);
    MVT::SetBlock(cNew, dstInd, *C_);

    ++curBasisSize_;

    return true;
  }

private:
  Teuchos::RCP<const OP> A_;

  int maxBasisSize_;
  int curBasisSize_;

  std::string replacementStrategy_;

  // Relative tolerance.  Candidate pair is rejected if
  // ||c_orth|| <= projectionTol_ * ||c_original||.
  MagnitudeType projectionTol_;

  int reorthogonalizationSteps_;

  // Correction basis and image basis.
  //
  // Invariant:
  //
  //   A U = C,
  //   C^H C = I.
  Teuchos::RCP<MV> U_;
  Teuchos::RCP<MV> C_;
};

} // namespace Belos

#endif // BELOS_SOLUTION_PROJECTION_HPP
