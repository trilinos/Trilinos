// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//                    Tobias Wiesner    (tawiesn@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef BELOS_XPETRA_STATUS_TEST_GEN_RES_SUB_NORM_HPP
#define BELOS_XPETRA_STATUS_TEST_GEN_RES_SUB_NORM_HPP

#include "Xpetra_ConfigDefs.hpp"

#include "Xpetra_BlockedCrsMatrix.hpp"

#include "MueLu_Exceptions.hpp"

#include <BelosConfigDefs.hpp>
#include <BelosTypes.hpp>
#include <BelosOperatorT.hpp>
#include <BelosXpetraAdapterOperator.hpp>
#include <BelosStatusTestGenResSubNorm.hpp>

namespace Belos {

/*! \brief Template specialization of Belos::StatusTestGenResSubNorm class using the
 * Xpetra::MultiVector and Belos::OperatorT MueLu adapter class.
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class StatusTestGenResSubNorm<Scalar, Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>, Belos::OperatorT<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > >
  : public StatusTestResNorm<Scalar, Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>, Belos::OperatorT<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > > {
 public:
  // Convenience typedefs
  typedef Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> MV;
  typedef Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> BCRS;
  typedef Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> ME;
  typedef Belos::OperatorT<MV> OP;

  typedef Teuchos::ScalarTraits<Scalar> SCT;
  typedef typename SCT::magnitudeType MagnitudeType;
  typedef MultiVecTraits<Scalar, MV> MVT;
  typedef OperatorTraits<Scalar, MV, OP> OT;

  //! @name Constructors/destructors.
  //@{
  //! Constructor
  /*! The constructor takes a single argument specifying the tolerance (\f$\tau\f$).
    If none of the form definition methods are called, we use \f$\|r\|/\|r^{(0)}\| \le \tau\f$
    as the stopping criterion, where \f$\|r\|\f$ always uses the true residual and
    \f$\|r^{(0)}\|\f$ is the corresponding norm of the initial residual.
    The used norm can be specified by defineResForm and defineScaleForm.

    @param Tolerance: Specifies tolerance \f$\tau\f
    @param subIdx: index of block row in the n x n block system we want to check the residual of
    @param quorum: Number of residual (sub-)vectors which are needed to be within the tolerance before check is considered to be passed
    @param showMaxResNormOnly: for output only

  */
  StatusTestGenResSubNorm(MagnitudeType Tolerance, size_t subIdx, int quorum = -1, bool showMaxResNormOnly = false)
    : tolerance_(Tolerance)
    , subIdx_(subIdx)
    , quorum_(quorum)
    , showMaxResNormOnly_(showMaxResNormOnly)
    , resnormtype_(TwoNorm)
    , scaletype_(NormOfInitRes)
    , scalenormtype_(TwoNorm)
    , scalevalue_(Teuchos::ScalarTraits<MagnitudeType>::one())
    , status_(Undefined)
    , curBlksz_(0)
    , curNumRHS_(0)
    , curLSNum_(0)
    , numrhs_(0)
    , firstcallCheckStatus_(true)
    , firstcallDefineResForm_(true)
    , firstcallDefineScaleForm_(true)
    , mapExtractor_(Teuchos::null) {}

  //! Destructor
  virtual ~StatusTestGenResSubNorm(){};
  //@}

  //! @name Form and parameter definition methods.
  //@{

  //! Define form of the residual, its norm and optional weighting std::vector.
  /*! This method defines the form of \f$\|r\|\f$.  We specify:
    <ul>
    <li> The norm to be used on the residual (this may be different than the norm used in
    DefineScaleForm()).
    </ul>
  */
  int defineResForm(NormType TypeOfNorm) {
    TEUCHOS_TEST_FOR_EXCEPTION(firstcallDefineResForm_ == false, StatusTestError,
                               "StatusTestGenResSubNorm::defineResForm(): The residual form has already been defined.");
    firstcallDefineResForm_ = false;

    resnormtype_ = TypeOfNorm;

    return (0);
  }

  //! Define form of the scaling, its norm, its optional weighting std::vector, or, alternatively, define an explicit value.
  /*! This method defines the form of how the residual is scaled (if at all).  It operates in two modes:
    <ol>
    <li> User-provided scaling value:
    <ul>
    <li> Set argument TypeOfScaling to UserProvided.
    <li> Set ScaleValue to a non-zero value that the residual norm will be divided by.
    <li> TypeOfNorm argument will be ignored.
    <li> Sample use:  Define ScaleValue = \f$\|A\|_{\infty}\f$ where \f$ A \f$ is the matrix
    of the linear problem.
    </ul>

    <li> Use a supported Scaling Form:
    <ul>
    <li> Define TypeOfScaling to be the norm of the right hand side, the initial residual std::vector,
    or to none.
    <li> Define norm to be used on the scaling std::vector (this may be different than the norm used
    in DefineResForm()).
    </ul>
    </ol>
  */
  int defineScaleForm(ScaleType TypeOfScaling, NormType TypeOfNorm, MagnitudeType ScaleValue = Teuchos::ScalarTraits<MagnitudeType>::one()) {
    TEUCHOS_TEST_FOR_EXCEPTION(firstcallDefineScaleForm_ == false, StatusTestError,
                               "StatusTestGenResSubNorm::defineScaleForm(): The scaling type has already been defined.");
    firstcallDefineScaleForm_ = false;

    scaletype_     = TypeOfScaling;
    scalenormtype_ = TypeOfNorm;
    scalevalue_    = ScaleValue;

    return (0);
  }

  //! Set the value of the tolerance
  /*! We allow the tolerance to be reset for cases where, in the process of testing the residual,
    we find that the initial tolerance was too tight or too lax.
  */
  int setTolerance(MagnitudeType tolerance) {
    tolerance_ = tolerance;
    return (0);
  }

  //! Set the block index of which we want to check the norm of the sub-residuals
  /*! It does not really make sense to change/reset the index during the solution process
   */
  int setSubIdx(size_t subIdx) {
    subIdx_ = subIdx;
    return (0);
  }

  //! Sets the number of residuals that must pass the convergence test before Passed is returned.
  //! \note If \c quorum=-1 then all residuals must pass the convergence test before Passed is returned.
  int setQuorum(int quorum) {
    quorum_ = quorum;
    return (0);
  }

  //! Set whether the only maximum residual norm is displayed when the print() method is called
  int setShowMaxResNormOnly(bool showMaxResNormOnly) {
    showMaxResNormOnly_ = showMaxResNormOnly;
    return (0);
  }

  //@}

  //! @name Status methods
  //@{
  //! Check convergence status: Passed, Failed, or Undefined.
  /*! This method checks to see if the convergence criteria are met.
    Depending on how the residual test is constructed this method will return
    the appropriate status type.

    \return StatusType: Passed, Failed, or Undefined.
  */
  StatusType checkStatus(Iteration<Scalar, MV, OP>* iSolver) {
    MagnitudeType zero                      = Teuchos::ScalarTraits<MagnitudeType>::zero();
    const LinearProblem<Scalar, MV, OP>& lp = iSolver->getProblem();
    // Compute scaling term (done once for each block that's being solved)
    if (firstcallCheckStatus_) {
      StatusType status = firstCallCheckStatusSetup(iSolver);
      if (status == Failed) {
        status_ = Failed;
        return (status_);
      }
    }

    //
    // This section computes the norm of the residual std::vector
    //
    if (curLSNum_ != lp.getLSNumber()) {
      //
      // We have moved on to the next rhs block
      //
      curLSNum_   = lp.getLSNumber();
      curLSIdx_   = lp.getLSIndex();
      curBlksz_   = (int)curLSIdx_.size();
      int validLS = 0;
      for (int i = 0; i < curBlksz_; ++i) {
        if (curLSIdx_[i] > -1 && curLSIdx_[i] < numrhs_)
          validLS++;
      }
      curNumRHS_ = validLS;
      curSoln_   = Teuchos::null;
      //
    } else {
      //
      // We are in the same rhs block, return if we are converged
      //
      if (status_ == Passed) {
        return status_;
      }
    }

    //
    // Request the true residual for this block of right-hand sides.
    //
    Teuchos::RCP<MV> cur_update = iSolver->getCurrentUpdate();
    curSoln_                    = lp.updateSolution(cur_update);
    Teuchos::RCP<MV> cur_res    = MVT::Clone(*curSoln_, MVT::GetNumberVecs(*curSoln_));
    lp.computeCurrResVec(&*cur_res, &*curSoln_);
    std::vector<MagnitudeType> tmp_resvector(MVT::GetNumberVecs(*cur_res));
    MvSubNorm(*cur_res, subIdx_, tmp_resvector, resnormtype_);

    typename std::vector<int>::iterator p = curLSIdx_.begin();
    for (int i = 0; p < curLSIdx_.end(); ++p, ++i) {
      // Check if this index is valid
      if (*p != -1)
        resvector_[*p] = tmp_resvector[i];
    }

    //
    // Compute the new linear system residuals for testing.
    // (if any of them don't meet the tolerance or are NaN, then we exit with that status)
    //
    if (scalevector_.size() > 0) {
      typename std::vector<int>::iterator pp = curLSIdx_.begin();
      for (; pp < curLSIdx_.end(); ++pp) {
        // Check if this index is valid
        if (*pp != -1) {
          // Scale the std::vector accordingly
          if (scalevector_[*pp] != zero) {
            // Don't intentionally divide by zero.
            testvector_[*pp] = resvector_[*pp] / scalevector_[*pp] / scalevalue_;
          } else {
            testvector_[*pp] = resvector_[*pp] / scalevalue_;
          }
        }
      }
    } else {
      typename std::vector<int>::iterator pp = curLSIdx_.begin();
      for (; pp < curLSIdx_.end(); ++pp) {
        // Check if this index is valid
        if (*pp != -1)
          testvector_[*pp] = resvector_[*pp] / scalevalue_;
      }
    }
    // Check status of new linear system residuals and see if we have the quorum.
    int have = 0;
    ind_.resize(curLSIdx_.size());
    typename std::vector<int>::iterator p2 = curLSIdx_.begin();
    for (; p2 < curLSIdx_.end(); ++p2) {
      // Check if this index is valid
      if (*p2 != -1) {
        // Check if any of the residuals are larger than the tolerance.
        if (testvector_[*p2] > tolerance_) {
          // do nothing.
        } else if (testvector_[*p2] == Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::zero())) {
          reset();
        } else if (testvector_[*p2] <= tolerance_) {
          ind_[have] = *p2;
          have++;
        } else {
          // Throw an std::exception if a NaN is found.
          status_ = Failed;
          TEUCHOS_TEST_FOR_EXCEPTION(true, StatusTestError, "StatusTestGenResSubNorm::checkStatus(): NaN has been detected.");
        }
      }
    }
    ind_.resize(have);
    int need = (quorum_ == -1) ? curNumRHS_ : quorum_;
    status_  = (have >= need) ? Passed : Failed;
    // Return the current status
    return status_;
  }

  //! Return the result of the most recent CheckStatus call.
  StatusType getStatus() const { return (status_); };
  //@}

  //! @name Reset methods
  //@{

  //! Resets the internal configuration to the initial state.
  void reset() {
    status_   = Undefined;
    curBlksz_ = 0;
    curLSNum_ = 0;
    curLSIdx_.resize(0);
    numrhs_ = 0;
    ind_.resize(0);
    firstcallCheckStatus_ = true;
    curSoln_              = Teuchos::null;
  }

  //@}

  //! @name Print methods
  //@{

  //! Output formatted description of stopping test to output stream.
  void print(std::ostream& os, int indent = 0) const {
    os.setf(std::ios_base::scientific);
    for (int j = 0; j < indent; j++)
      os << ' ';
    printStatus(os, status_);
    os << resFormStr();
    if (status_ == Undefined)
      os << ", tol = " << tolerance_ << std::endl;
    else {
      os << std::endl;
      if (showMaxResNormOnly_ && curBlksz_ > 1) {
        const MagnitudeType maxRelRes = *std::max_element(
            testvector_.begin() + curLSIdx_[0], testvector_.begin() + curLSIdx_[curBlksz_ - 1]);
        for (int j = 0; j < indent + 13; j++)
          os << ' ';
        os << "max{residual[" << curLSIdx_[0] << "..." << curLSIdx_[curBlksz_ - 1] << "]} = " << maxRelRes
           << (maxRelRes <= tolerance_ ? " <= " : " > ") << tolerance_ << std::endl;
      } else {
        for (int i = 0; i < numrhs_; i++) {
          for (int j = 0; j < indent + 13; j++)
            os << ' ';
          os << "residual [ " << i << " ] = " << testvector_[i];
          os << ((testvector_[i] < tolerance_) ? " < " : (testvector_[i] == tolerance_) ? " == "
                                                     : (testvector_[i] > tolerance_)    ? " > "
                                                                                        : " ")
             << tolerance_ << std::endl;
        }
      }
    }
    os << std::endl;
  }

  //! Print message for each status specific to this stopping test.
  void printStatus(std::ostream& os, StatusType type) const {
    os << std::left << std::setw(13) << std::setfill('.');
    switch (type) {
      case Passed:
        os << "Converged";
        break;
      case Failed:
        os << "Unconverged";
        break;
      case Undefined:
      default:
        os << "**";
        break;
    }
    os << std::left << std::setfill(' ');
    return;
  }
  //@}

  //! @name Methods to access data members.
  //@{

  //! Returns the current solution estimate that was computed for the most recent residual test.
  Teuchos::RCP<MV> getSolution() { return curSoln_; }

  //! Returns the number of residuals that must pass the convergence test before Passed is returned.
  //! \note If \c quorum=-1 then all residuals must pass the convergence test before Passed is returned.
  int getQuorum() const { return quorum_; }

  //! Returns the index of the block row the norms are calculated for
  size_t getSubIdx() const { return subIdx_; }

  //! Returns whether the only maximum residual norm is displayed when the print() method is called
  bool getShowMaxResNormOnly() { return showMaxResNormOnly_; }

  //! Returns the std::vector containing the indices of the residuals that passed the test.
  std::vector<int> convIndices() { return ind_; }

  //! Returns the value of the tolerance, \f$ \tau \f$, set in the constructor.
  MagnitudeType getTolerance() const { return (tolerance_); };

  //! Returns the test value, \f$ \frac{\|r\|}{\sigma} \f$, computed in most recent call to CheckStatus.
  const std::vector<MagnitudeType>* getTestValue() const { return (&testvector_); };

  //! Returns the residual norm value, \f$ \|r\| \f$, computed in most recent call to CheckStatus.
  const std::vector<MagnitudeType>* getResNormValue() const { return (&resvector_); };

  //! Returns the scaled norm value, \f$ \sigma \f$.
  const std::vector<MagnitudeType>* getScaledNormValue() const { return (&scalevector_); };

  //! Returns a boolean indicating a loss of accuracy has been detected in computing the residual.
  //! \note This status test does not check for loss of accuracy, so this method will always return false.
  bool getLOADetected() const { return false; }

  //@}

  /** @name Misc. */
  //@{

  /** \brief Call to setup initial scaling std::vector.
   *
   * After this function is called <tt>getScaledNormValue()</tt> can be called
   * to get the scaling std::vector.
   */
  StatusType firstCallCheckStatusSetup(Iteration<Scalar, MV, OP>* iSolver) {
    int i;
    MagnitudeType zero                      = Teuchos::ScalarTraits<MagnitudeType>::zero();
    MagnitudeType one                       = Teuchos::ScalarTraits<MagnitudeType>::one();
    const LinearProblem<Scalar, MV, OP>& lp = iSolver->getProblem();
    // Compute scaling term (done once for each block that's being solved)
    if (firstcallCheckStatus_) {
      //
      // Get some current solver information.
      //
      firstcallCheckStatus_ = false;

      // try to access the underlying blocked operator
      Teuchos::RCP<const OP> Op = lp.getOperator();
      Teuchos::RCP<const Belos::XpetraOp<Scalar, LocalOrdinal, GlobalOrdinal, Node> > xOp =
          Teuchos::rcp_dynamic_cast<const Belos::XpetraOp<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(Op);
      TEUCHOS_TEST_FOR_EXCEPTION(xOp.is_null(), MueLu::Exceptions::BadCast, "Bad cast from \'const Belos::OperatorT\' to \'const Belos::XpetraOp\'. The origin type is " << typeid(const OP).name() << ".");
      Teuchos::RCP<const Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> > xIntOp =
          xOp->getOperator();
      TEUCHOS_TEST_FOR_EXCEPTION(xIntOp.is_null(), MueLu::Exceptions::BadCast, "Cannot access Xpetra::Operator stored in Belos::XpetraOperator.");
      Teuchos::RCP<const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > xMat =
          Teuchos::rcp_dynamic_cast<const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(xIntOp);
      TEUCHOS_TEST_FOR_EXCEPTION(xMat.is_null(), MueLu::Exceptions::RuntimeError, "Cannot access Xpetra::Matrix stored in Belos::XpetraOp. Error.");
      Teuchos::RCP<const Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > bMat = Teuchos::rcp_dynamic_cast<const Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(xMat);
      TEUCHOS_TEST_FOR_EXCEPTION(bMat.is_null(), MueLu::Exceptions::BadCast, "Bad cast from \'const Xpetra::Matrix\' to \'const Xpetra::BlockedCrsMatrix\'. The origin type is " << typeid(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>).name() << ". Note: you need a BlockedCrsMatrix object for the StatusTestGenResSubNorm to work!");
      mapExtractor_ = bMat->getRangeMapExtractor();
      TEUCHOS_TEST_FOR_EXCEPTION(mapExtractor_.is_null(), MueLu::Exceptions::RuntimeError, "Could not extract map extractor from BlockedCrsMatrix. Error.");
      TEUCHOS_TEST_FOR_EXCEPTION(mapExtractor_->NumMaps() <= subIdx_, MueLu::Exceptions::RuntimeError, "The multivector is only split into " << mapExtractor_->NumMaps() << " sub parts. Cannot access sub-block " << subIdx_ << ".");

      // calculate initial norms
      if (scaletype_ == NormOfRHS) {
        Teuchos::RCP<const MV> rhs = lp.getRHS();
        numrhs_                    = MVT::GetNumberVecs(*rhs);
        scalevector_.resize(numrhs_);
        MvSubNorm(*rhs, subIdx_, scalevector_, scalenormtype_);
      } else if (scaletype_ == NormOfInitRes) {
        Teuchos::RCP<const MV> init_res = lp.getInitResVec();
        numrhs_                         = MVT::GetNumberVecs(*init_res);
        scalevector_.resize(numrhs_);
        MvSubNorm(*init_res, subIdx_, scalevector_, scalenormtype_);
      } else if (scaletype_ == NormOfPrecInitRes) {
        Teuchos::RCP<const MV> init_res = lp.getInitPrecResVec();
        numrhs_                         = MVT::GetNumberVecs(*init_res);
        scalevector_.resize(numrhs_);
        MvSubNorm(*init_res, subIdx_, scalevector_, scalenormtype_);
      } else if (scaletype_ == NormOfFullInitRes) {
        Teuchos::RCP<const MV> init_res = lp.getInitResVec();
        numrhs_                         = MVT::GetNumberVecs(*init_res);
        scalevector_.resize(numrhs_);
        MVT::MvNorm(*init_res, scalevector_, scalenormtype_);
        scalevalue_ = one;
      } else if (scaletype_ == NormOfFullPrecInitRes) {
        Teuchos::RCP<const MV> init_res = lp.getInitPrecResVec();
        numrhs_                         = MVT::GetNumberVecs(*init_res);
        scalevector_.resize(numrhs_);
        MVT::MvNorm(*init_res, scalevector_, scalenormtype_);
        scalevalue_ = one;
      } else if (scaletype_ == NormOfFullScaledInitRes) {
        Teuchos::RCP<const MV> init_res = lp.getInitResVec();
        numrhs_                         = MVT::GetNumberVecs(*init_res);
        scalevector_.resize(numrhs_);
        MVT::MvNorm(*init_res, scalevector_, scalenormtype_);
        MvScalingRatio(*init_res, subIdx_, scalevalue_);
      } else if (scaletype_ == NormOfFullScaledPrecInitRes) {
        Teuchos::RCP<const MV> init_res = lp.getInitPrecResVec();
        numrhs_                         = MVT::GetNumberVecs(*init_res);
        scalevector_.resize(numrhs_);
        MVT::MvNorm(*init_res, scalevector_, scalenormtype_);
        MvScalingRatio(*init_res, subIdx_, scalevalue_);
      } else {
        numrhs_ = MVT::GetNumberVecs(*(lp.getRHS()));
      }

      resvector_.resize(numrhs_);
      testvector_.resize(numrhs_);

      curLSNum_   = lp.getLSNumber();
      curLSIdx_   = lp.getLSIndex();
      curBlksz_   = (int)curLSIdx_.size();
      int validLS = 0;
      for (i = 0; i < curBlksz_; ++i) {
        if (curLSIdx_[i] > -1 && curLSIdx_[i] < numrhs_)
          validLS++;
      }
      curNumRHS_ = validLS;
      //
      // Initialize the testvector.
      for (i = 0; i < numrhs_; i++) {
        testvector_[i] = one;
      }

      // Return an error if the scaling is zero.
      if (scalevalue_ == zero) {
        return Failed;
      }
    }
    return Undefined;
  }
  //@}

  /** \name Overridden from Teuchos::Describable */
  //@{

  /** \brief Method to return description of the maximum iteration status test  */
  std::string description() const {
    std::ostringstream oss;
    oss << "Belos::StatusTestGenResSubNorm<>: " << resFormStr();
    oss << ", tol = " << tolerance_;
    return oss.str();
  }
  //@}

 protected:
 private:
  //! @name Private methods.
  //@{
  /** \brief Description of current residual form */
  std::string resFormStr() const {
    std::ostringstream oss;
    oss << "(";
    oss << ((resnormtype_ == OneNorm) ? "1-Norm" : (resnormtype_ == TwoNorm) ? "2-Norm"
                                                                             : "Inf-Norm");
    oss << " Exp";
    oss << " Res Vec [" << subIdx_ << "]) ";

    // If there is no residual scaling, return current string.
    if (scaletype_ != None) {
      // Insert division sign.
      oss << "/ ";

      // Determine output string for scaling, if there is any.
      if (scaletype_ == UserProvided)
        oss << " (User Scale)";
      else {
        oss << "(";
        oss << ((scalenormtype_ == OneNorm) ? "1-Norm" : (resnormtype_ == TwoNorm) ? "2-Norm"
                                                                                   : "Inf-Norm");
        if (scaletype_ == NormOfInitRes)
          oss << " Res0 [" << subIdx_ << "]";
        else if (scaletype_ == NormOfPrecInitRes)
          oss << " Prec Res0 [" << subIdx_ << "]";
        else if (scaletype_ == NormOfFullInitRes)
          oss << " Full Res0 [" << subIdx_ << "]";
        else if (scaletype_ == NormOfFullPrecInitRes)
          oss << " Full Prec Res0 [" << subIdx_ << "]";
        else if (scaletype_ == NormOfFullScaledInitRes)
          oss << " scaled Full Res0 [" << subIdx_ << "]";
        else if (scaletype_ == NormOfFullScaledPrecInitRes)
          oss << " scaled Full Prec Res0 [" << subIdx_ << "]";
        else
          oss << " RHS [" << subIdx_ << "]";
        oss << ")";
      }
    }

    // TODO add a tagging name

    return oss.str();
  }

  //@}

  //! @name Private helper functions
  //@{

  // calculate norm of partial multivector
  void MvSubNorm(const MV& mv, size_t block, std::vector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>& normVec, NormType type = TwoNorm) {
    Teuchos::RCP<const MV> input = Teuchos::rcpFromRef(mv);

    Teuchos::RCP<const MV> SubVec = mapExtractor_->ExtractVector(input, block);
    MVT::MvNorm(*SubVec, normVec, type);
  }

  // calculate ration of sub-vector length to full vector length (for scalevalue_)
  void MvScalingRatio(const MV& mv, size_t block, MagnitudeType& lengthRatio) {
    Teuchos::RCP<const MV> input = Teuchos::rcpFromRef(mv);

    Teuchos::RCP<const MV> SubVec = mapExtractor_->ExtractVector(input, block);

    lengthRatio = Teuchos::as<MagnitudeType>(SubVec->getGlobalLength()) / Teuchos::as<MagnitudeType>(input->getGlobalLength());
  }
  //@}

  //! @name Private data members.
  //@{

  //! Tolerance used to determine convergence
  MagnitudeType tolerance_;

  //! Index of block row in n x n block system of which we want to check the sub-residuals
  size_t subIdx_;

  //! Number of residuals that must pass the convergence test before Passed is returned.
  int quorum_;

  //! Determines if the entries for all of the residuals are shown or just the max.
  bool showMaxResNormOnly_;

  //! Type of norm to use on residual (OneNorm, TwoNorm, or InfNorm).
  NormType resnormtype_;

  //! Type of scaling to use (Norm of RHS, Norm of Initial Residual, None or User provided)
  ScaleType scaletype_;

  //! Type of norm to use on the scaling (OneNorm, TwoNorm, or InfNorm)
  NormType scalenormtype_;

  //! Scaling value.
  MagnitudeType scalevalue_;

  //! Scaling std::vector.
  std::vector<MagnitudeType> scalevector_;

  //! Residual norm std::vector.
  std::vector<MagnitudeType> resvector_;

  //! Test std::vector = resvector_ / scalevector_
  std::vector<MagnitudeType> testvector_;

  //! Vector containing the indices for the vectors that passed the test.
  std::vector<int> ind_;

  //! Most recent solution vector used by this status test.
  Teuchos::RCP<MV> curSoln_;

  //! Status
  StatusType status_;

  //! The current blocksize of the linear system being solved.
  int curBlksz_;

  //! The current number of right-hand sides being solved for.
  int curNumRHS_;

  //! The indices of the current number of right-hand sides being solved for.
  std::vector<int> curLSIdx_;

  //! The current number of linear systems that have been loaded into the linear problem.
  int curLSNum_;

  //! The total number of right-hand sides being solved for.
  int numrhs_;

  //! Is this the first time CheckStatus is called?
  bool firstcallCheckStatus_;

  //! Is this the first time DefineResForm is called?
  bool firstcallDefineResForm_;

  //! Is this the first time DefineScaleForm is called?
  bool firstcallDefineScaleForm_;

  //! MapExtractor for range space
  Teuchos::RCP<const ME> mapExtractor_;
  //@}
};

}  // namespace Belos

#endif /* BELOS_XPETRA_STATUS_TEST_GEN_RES_SUB_NORM_HPP */
