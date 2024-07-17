// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef BELOS_STATUS_TEST_GEN_RESNORM_H
#define BELOS_STATUS_TEST_GEN_RESNORM_H

/*!
  \file BelosStatusTestGenResNorm.hpp
  \brief Belos::StatusTestResNorm for specifying general residual norm stopping criteria.
*/

#include "BelosStatusTestResNorm.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosMultiVecTraits.hpp"

/*!
  \class Belos::StatusTestGenResNorm
  \brief An implementation of StatusTestResNorm using a family of residual norms.

  StatusTestGenResNorm is an implementation of StatusTestResNorm that allows a user to construct
  one of a family of residual tests for use as a status/convergence test for Belos.
  The form of the test is
   \f[
   \frac{\|r_i\|}{\sigma_i} \le \tau
   \f]
   where
   <ul>
   <li> \f$r_i\f$ is the i-th residual std::vector, implicitly or explicitly computed (determined by enum ResType),
   <li> \f$\|r_i\|\f$ is the i-th residual norm determined by the enum NormType  (1-norm, 2-norm or inf-norm),
   <li> \f$\sigma_i\f$ is the i-th scale factor that can be passed in as a precomputed number of the templated type,
   or can be selected from by the enum ScaleType (norm of RHS, norm of initial residual).
   <li> \f$\tau\f$ is the tolerance that is passed in as a number of the templated type to the constructor.
   The value of \f$\tau\f$ can be reset using the setTolerance() method.
   </ul>

*/

namespace Belos {

template <class ScalarType, class MV, class OP>
class StatusTestGenResNorm: public StatusTestResNorm<ScalarType,MV,OP> {

 public:

  // Convenience typedefs
  typedef Teuchos::ScalarTraits<ScalarType> SCT;
  typedef typename SCT::magnitudeType MagnitudeType;
  typedef MultiVecTraits<ScalarType,MV>  MVT;

  //! @name Enums.
  //@{

  /*!
    \brief Select how the residual std::vector is produced.
  */
  enum ResType {Implicit, /*!< Use the residual std::vector produced by the iterative solver. */
                Explicit  /*!< Explicitly compute the residual std::vector r = b - A*x using the
                            linear problem. */
  };

  //@}

  //! @name Constructors/destructors.
  //@{
  //! Constructor
  /*! The constructor takes a single argument specifying the tolerance (\f$\tau\f$).
    If none of the form definition methods are called, we use \f$\|r\|_2/\|r^{(0)}\|_2 \le \tau\f$
    as the stopping criterion, where \f$\|r\|_2\f$ uses the least costly form of the 2-norm of
    residual available from the iterative method and \f$\|r^{(0)}\|_2\f$ is the corresponding norm
    of the initial residual.  The least costly form of the 2-norm depends on the chosen iterative
    method.  Most Krylov methods produce the preconditioned residual std::vector in a form that would be
    exact in infinite precision arithmetic.  This std::vector may be different from the true residual
    either because left preconditioning was used, or because round-off error has
    introduced significant error, or both.

    You can also state the number of vectors that must pass the convergence criteria before the
    status test passes by using the \c quorum argument.
  */
  StatusTestGenResNorm( MagnitudeType Tolerance, int quorum = -1, bool showMaxResNormOnly = false );

  //! Destructor
  virtual ~StatusTestGenResNorm();
  //@}

  //! @name Form and parameter definition methods.
  //@{

  //! Define form of the residual, its norm and optional weighting std::vector.
  /*! This method defines the form of \f$\|r\|\f$.  We specify:
    <ul>
    <li> Whether the residual std::vector should be explicitly computed, or taken from the iterative method.
    <li> The norm to be used on the residual (this may be different than the norm used in
    DefineScaleForm()).
    </ul>
  */
  int defineResForm( ResType TypeOfResidual, NormType TypeOfNorm);

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
  int defineScaleForm( ScaleType TypeOfScaling, NormType TypeOfNorm, MagnitudeType ScaleValue = Teuchos::ScalarTraits<MagnitudeType>::one());

  NormType getResNormType() {return resnormtype_;}

  //! Set the value of the tolerance
  /*! We allow the tolerance to be reset for cases where, in the process of testing the residual,
    we find that the initial tolerance was too tight or too lax.
  */
  int setTolerance(MagnitudeType tolerance) {tolerance_ = tolerance; return(0);}

  //! Sets the number of residuals that must pass the convergence test before Passed is returned.
  //! \note If \c quorum=-1 then all residuals must pass the convergence test before Passed is returned.
  int setQuorum(int quorum) {quorum_ = quorum; return(0);}

  //! Set whether the only maximum residual norm is displayed when the print() method is called
  int setShowMaxResNormOnly(bool showMaxResNormOnly) {showMaxResNormOnly_ = showMaxResNormOnly; return(0);}

  //@}

  //! @name Status methods
  //@{
  //! Check convergence status: Passed, Failed, or Undefined.
  /*! This method checks to see if the convergence criteria are met.
    Depending on how the residual test is constructed this method will return
    the appropriate status type.

    \return StatusType: Passed, Failed, or Undefined.
  */
  StatusType checkStatus(Iteration<ScalarType,MV,OP>* iSolver);

  //! Return the result of the most recent CheckStatus call.
  StatusType getStatus() const {return(status_);};
  //@}

  //! @name Reset methods
  //@{

  //! Resets the internal configuration to the initial state.
  void reset();

  //@}

  //! @name Print methods
  //@{

  //! Output formatted description of stopping test to output stream.
  void print(std::ostream& os, int indent = 0) const;

  //! Print message for each status specific to this stopping test.
  void printStatus(std::ostream& os, StatusType type) const;
  //@}

  //! @name Methods to access data members.
  //@{

  //! Returns the current solution estimate that was computed for the most recent residual test.
  //! \note This is useful for explicit residual tests, if this test is an implicit residual test
  //! a null pointer will be returned.
  Teuchos::RCP<MV> getSolution() { if (restype_==Implicit) { return Teuchos::null; } else { return curSoln_; } }

  //! Returns the number of residuals that must pass the convergence test before Passed is returned.
  //! \note If \c quorum=-1 then all residuals must pass the convergence test before Passed is returned.
  int getQuorum() const { return quorum_; }

  //! Returns whether the only maximum residual norm is displayed when the print() method is called
  bool getShowMaxResNormOnly() { return showMaxResNormOnly_; }

  //! Returns the std::vector containing the indices of the residuals that passed the test.
  std::vector<int> convIndices() { return ind_; }

  //! Returns the value of the tolerance, \f$ \tau \f$, set in the constructor.
  MagnitudeType getTolerance() const {return(tolerance_);};

  //! Returns the test value, \f$ \frac{\|r\|}{\sigma} \f$, computed in most recent call to CheckStatus.
  const std::vector<MagnitudeType>* getTestValue() const {return(&testvector_);};

  //! Returns the residual norm value, \f$ \|r\| \f$, computed in most recent call to CheckStatus.
  const std::vector<MagnitudeType>* getResNormValue() const {return(&resvector_);};

  //! Returns the scaled norm value, \f$ \sigma \f$.
  const std::vector<MagnitudeType>* getScaledNormValue() const {return(&scalevector_);};

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
  StatusType firstCallCheckStatusSetup(Iteration<ScalarType,MV,OP>* iSolver);
  //@}

  /** \name Overridden from Teuchos::Describable */
  //@{

  /** \brief Method to return description of the maximum iteration status test  */
  std::string description() const
  {
    std::ostringstream oss;
    oss << "Belos::StatusTestGenResNorm<>: " << resFormStr();
    oss << ", tol = " << tolerance_;
    return oss.str();
  }
  //@}

 protected:

 private:

  //! @name Private methods.
  //@{
  /** \brief Description of current residual form */
  std::string resFormStr() const
  {
    std::ostringstream oss;
    oss << "(";
    oss << ((resnormtype_==OneNorm) ? "1-Norm" : (resnormtype_==TwoNorm) ? "2-Norm" : "Inf-Norm");
    oss << ((restype_==Explicit) ? " Exp" : " Imp");
    oss << " Res Vec) ";

    // If there is no residual scaling, return current string.
    if (scaletype_!=None)
    {
      // Insert division sign.
      oss << "/ ";

      // Determine output string for scaling, if there is any.
      if (scaletype_==UserProvided)
        oss << " (User Scale)";
      else {
        oss << "(";
        oss << ((scalenormtype_==OneNorm) ? "1-Norm" : (resnormtype_==TwoNorm) ? "2-Norm" : "Inf-Norm");
        if (scaletype_==NormOfInitRes || scaletype_==NormOfFullInitRes || scaletype_==NormOfFullScaledInitRes)
          oss << " Res0";
        else if (scaletype_==NormOfPrecInitRes || scaletype_==NormOfFullPrecInitRes || scaletype_==NormOfFullScaledPrecInitRes)
          oss << " Prec Res0";
        else
          oss << " RHS ";
        oss << ")";
      }
    }

    return oss.str();
  }

  //@}

  //! @name Private data members.
  //@{

  //! Tolerance used to determine convergence
  MagnitudeType tolerance_;

  //! Number of residuals that must pass the convergence test before Passed is returned.
  int quorum_;

  //! Determines if the entries for all of the residuals are shown or just the max.
  bool showMaxResNormOnly_;

  //! Type of residual to use (explicit or implicit)
  ResType restype_;

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

  //@}

};

template <class ScalarType, class MV, class OP>
StatusTestGenResNorm<ScalarType,MV,OP>::
StatusTestGenResNorm (MagnitudeType Tolerance, int quorum, bool showMaxResNormOnly)
  : tolerance_(Tolerance),
    quorum_(quorum),
    showMaxResNormOnly_(showMaxResNormOnly),
    restype_(Implicit),
    resnormtype_(TwoNorm),
    scaletype_(NormOfInitRes),
    scalenormtype_(TwoNorm),
    scalevalue_(Teuchos::ScalarTraits<MagnitudeType>::one ()),
    status_(Undefined),
    curBlksz_(0),
    curNumRHS_(0),
    curLSNum_(0),
    numrhs_(0),
    firstcallCheckStatus_(true),
    firstcallDefineResForm_(true),
    firstcallDefineScaleForm_(true)
{
  // This constructor will compute the residual ||r_i||/||r0_i|| <= tolerance using the 2-norm of
  // the implicit residual std::vector.
}

template <class ScalarType, class MV, class OP>
StatusTestGenResNorm<ScalarType,MV,OP>::~StatusTestGenResNorm()
{}

template <class ScalarType, class MV, class OP>
void StatusTestGenResNorm<ScalarType,MV,OP>::reset()
{
  status_ = Undefined;
  curBlksz_ = 0;
  curLSNum_ = 0;
  curLSIdx_.resize(0);
  numrhs_ = 0;
  ind_.resize(0);
  firstcallCheckStatus_ = true;
  curSoln_ = Teuchos::null;
}

template <class ScalarType, class MV, class OP>
int StatusTestGenResNorm<ScalarType,MV,OP>::defineResForm( ResType TypeOfResidual, NormType TypeOfNorm )
{
  TEUCHOS_TEST_FOR_EXCEPTION(firstcallDefineResForm_==false,StatusTestError,
        "StatusTestGenResNorm::defineResForm(): The residual form has already been defined.");
  firstcallDefineResForm_ = false;

  restype_ = TypeOfResidual;
  resnormtype_ = TypeOfNorm;

  return(0);
}

template <class ScalarType, class MV, class OP>
int StatusTestGenResNorm<ScalarType,MV,OP>::defineScaleForm(ScaleType TypeOfScaling, NormType TypeOfNorm,
                                                         MagnitudeType ScaleValue )
{
  TEUCHOS_TEST_FOR_EXCEPTION(firstcallDefineScaleForm_==false,StatusTestError,
        "StatusTestGenResNorm::defineScaleForm(): The scaling type has already been defined.");
  firstcallDefineScaleForm_ = false;

  scaletype_ = TypeOfScaling;
  scalenormtype_ = TypeOfNorm;
  scalevalue_ = ScaleValue;

  return(0);
}

template <class ScalarType, class MV, class OP>
StatusType StatusTestGenResNorm<ScalarType,MV,OP>::checkStatus( Iteration<ScalarType,MV,OP>* iSolver )
{
  MagnitudeType zero = Teuchos::ScalarTraits<MagnitudeType>::zero();
  const LinearProblem<ScalarType,MV,OP>& lp = iSolver->getProblem();
  // Compute scaling term (done once for each block that's being solved)
  if (firstcallCheckStatus_) {
    StatusType status = firstCallCheckStatusSetup(iSolver);
    if(status==Failed) {
      status_ = Failed;
      return(status_);
    }
  }
  //
  // This section computes the norm of the residual std::vector
  //
  if ( curLSNum_ != lp.getLSNumber() ) {
    //
    // We have moved on to the next rhs block
    //
    curLSNum_ = lp.getLSNumber();
    curLSIdx_ = lp.getLSIndex();
    curBlksz_ = (int)curLSIdx_.size();
    int validLS = 0;
    for (int i=0; i<curBlksz_; ++i) {
      if (curLSIdx_[i] > -1 && curLSIdx_[i] < numrhs_)
        validLS++;
    }
    curNumRHS_ = validLS;
    curSoln_ = Teuchos::null;
    //
  } else {
    //
    // We are in the same rhs block, return if we are converged
    //
    if (status_==Passed) { return status_; }
  }
  if (restype_==Implicit) {
    //
    // get the native residual norms from the solver for this block of right-hand sides.
    // If the residual is returned in multivector form, use the resnormtype to compute the residual norms.
    // Otherwise the native residual is assumed to be stored in the resvector_.
    //
    std::vector<MagnitudeType> tmp_resvector( curBlksz_ );
    Teuchos::RCP<const MV> residMV = iSolver->getNativeResiduals( &tmp_resvector );
    if ( residMV != Teuchos::null ) {
      tmp_resvector.resize( MVT::GetNumberVecs( *residMV ) );
      MVT::MvNorm( *residMV, tmp_resvector, resnormtype_ );
      typename std::vector<int>::iterator p = curLSIdx_.begin();
      for (int i=0; p<curLSIdx_.end(); ++p, ++i) {
        // Check if this index is valid
        if (*p != -1)
          resvector_[*p] = tmp_resvector[i];
      }
    } else {
      typename std::vector<int>::iterator p = curLSIdx_.begin();
      for (int i=0; p<curLSIdx_.end(); ++p, ++i) {
        // Check if this index is valid
        if (*p != -1)
          resvector_[*p] = tmp_resvector[i];
      }
    }
  }
  else if (restype_==Explicit) {
    //
    // Request the true residual for this block of right-hand sides.
    //
    Teuchos::RCP<MV> cur_update = iSolver->getCurrentUpdate();
    curSoln_ = lp.updateSolution( cur_update );
    Teuchos::RCP<MV> cur_res = MVT::Clone( *curSoln_, MVT::GetNumberVecs( *curSoln_ ) );
    lp.computeCurrResVec( &*cur_res, &*curSoln_ );
    std::vector<MagnitudeType> tmp_resvector( MVT::GetNumberVecs( *cur_res ) );
    MVT::MvNorm( *cur_res, tmp_resvector, resnormtype_ );
    typename std::vector<int>::iterator p = curLSIdx_.begin();
    for (int i=0; p<curLSIdx_.end(); ++p, ++i) {
      // Check if this index is valid
      if (*p != -1)
        resvector_[*p] = tmp_resvector[i];
    }
  }
  //
  // Compute the new linear system residuals for testing.
  // (if any of them don't meet the tolerance or are NaN, then we exit with that status)
  //
  if ( scalevector_.size() > 0 ) {
    typename std::vector<int>::iterator p = curLSIdx_.begin();
    for (; p<curLSIdx_.end(); ++p) {
      // Check if this index is valid
      if (*p != -1) {
        // Scale the std::vector accordingly
        testvector_[ *p ] =
            scalevector_[ *p ] != zero? resvector_[ *p ] / (scalevector_[ *p ] * scalevalue_) : resvector_[ *p ] / scalevalue_;
      }
    }
  }
  else {
    typename std::vector<int>::iterator p = curLSIdx_.begin();
    for (; p<curLSIdx_.end(); ++p) {
      // Check if this index is valid
      if (*p != -1)
        testvector_[ *p ] = resvector_[ *p ] / scalevalue_;
    }
  }

  // Check status of new linear system residuals and see if we have the quorum.
  int have = 0;
  ind_.resize( curLSIdx_.size() );
  typename std::vector<int>::iterator p = curLSIdx_.begin();
  for (; p<curLSIdx_.end(); ++p) {
    // Check if this index is valid
    if (*p != -1) {
      // Check if any of the residuals are larger than the tolerance.
      if (testvector_[ *p ] > tolerance_) {
        // do nothing.
      } else if (testvector_[ *p ] <= tolerance_) {
        ind_[have] = *p;
        have++;
      } else {
        // Throw an std::exception if a NaN is found.
        status_ = Failed;
        TEUCHOS_TEST_FOR_EXCEPTION(true,StatusTestNaNError,"StatusTestGenResNorm::checkStatus(): NaN has been detected.");
      }
    }
  }
  ind_.resize(have);
  int need = (quorum_ == -1) ? curNumRHS_: quorum_;
  status_ = (have >= need) ? Passed : Failed;

  // Return the current status
  return status_;
}

template <class ScalarType, class MV, class OP>
void StatusTestGenResNorm<ScalarType,MV,OP>::print(std::ostream& os, int indent) const
{
  for (int j = 0; j < indent; j ++)
    os << ' ';
  printStatus(os, status_);
  os << resFormStr();
  if (status_==Undefined)
    os << ", tol = " << tolerance_ << std::endl;
  else {
    os << std::endl;
    if(showMaxResNormOnly_ && curBlksz_ > 1) {
      const MagnitudeType maxRelRes = *std::max_element(
        testvector_.begin()+curLSIdx_[0],testvector_.begin()+curLSIdx_[curBlksz_-1]
        );
      for (int j = 0; j < indent + 13; j ++)
        os << ' ';
      os << "max{residual["<<curLSIdx_[0]<<"..."<<curLSIdx_[curBlksz_-1]<<"]} = " << maxRelRes
         << ( maxRelRes <= tolerance_ ? " <= " : " > " ) << tolerance_ << std::endl;
    }
    else {
      for ( int i=0; i<numrhs_; i++ ) {
        for (int j = 0; j < indent + 13; j ++)
          os << ' ';
        os << "residual [ " << i << " ] = " << testvector_[ i ];
        os << ((testvector_[i]<tolerance_) ? " < " : (testvector_[i]==tolerance_) ? " == " : (testvector_[i]>tolerance_) ? " > " : " "  ) << tolerance_ << std::endl;
      }
    }
  }
  os << std::endl;
}

template <class ScalarType, class MV, class OP>
void StatusTestGenResNorm<ScalarType,MV,OP>::printStatus(std::ostream& os, StatusType type) const
{
  os << std::left << std::setw(13) << std::setfill('.');
  switch (type) {
  case  Passed:
    os << "Converged";
    break;
  case  Failed:
    os << "Unconverged";
    break;
  case  Undefined:
  default:
    os << "**";
    break;
  }
  os << std::left << std::setfill(' ');
    return;
}

template <class ScalarType, class MV, class OP>
StatusType StatusTestGenResNorm<ScalarType,MV,OP>::firstCallCheckStatusSetup( Iteration<ScalarType,MV,OP>* iSolver )
{
  int i;
  MagnitudeType zero = Teuchos::ScalarTraits<MagnitudeType>::zero();
  MagnitudeType one = Teuchos::ScalarTraits<MagnitudeType>::one();
  const LinearProblem<ScalarType,MV,OP>& lp = iSolver->getProblem();
  // Compute scaling term (done once for each block that's being solved)
  if (firstcallCheckStatus_) {
    //
    // Get some current solver information.
    //
    firstcallCheckStatus_ = false;

    if (scaletype_== NormOfRHS) {
      Teuchos::RCP<const MV> rhs = lp.getRHS();
      numrhs_ = MVT::GetNumberVecs( *rhs );
      scalevector_.resize( numrhs_ );
      MVT::MvNorm( *rhs, scalevector_, scalenormtype_ );
    }
    else if (scaletype_==NormOfInitRes || scaletype_==NormOfFullInitRes || scaletype_==NormOfFullScaledInitRes) {
      Teuchos::RCP<const MV> init_res = lp.getInitResVec();
      numrhs_ = MVT::GetNumberVecs( *init_res );
      scalevector_.resize( numrhs_ );
      MVT::MvNorm( *init_res, scalevector_, scalenormtype_ );
    }
    else if (scaletype_==NormOfPrecInitRes || scaletype_==NormOfFullPrecInitRes || scaletype_==NormOfFullScaledPrecInitRes) {
      Teuchos::RCP<const MV> init_res = lp.getInitPrecResVec();
      numrhs_ = MVT::GetNumberVecs( *init_res );
      scalevector_.resize( numrhs_ );
      MVT::MvNorm( *init_res, scalevector_, scalenormtype_ );
    }
    else {
      numrhs_ = MVT::GetNumberVecs( *(lp.getRHS()) );
    }

    resvector_.resize( numrhs_ );
    testvector_.resize( numrhs_ );

    curLSNum_ = lp.getLSNumber();
    curLSIdx_ = lp.getLSIndex();
    curBlksz_ = (int)curLSIdx_.size();
    int validLS = 0;
    for (i=0; i<curBlksz_; ++i) {
      if (curLSIdx_[i] > -1 && curLSIdx_[i] < numrhs_)
        validLS++;
    }
    curNumRHS_ = validLS;
    //
    // Initialize the testvector.
    for (i=0; i<numrhs_; i++) { testvector_[i] = one; }

    // Return an error if the scaling is zero.
    if (scalevalue_ == zero) {
      return Failed;
    }
  }
  return Undefined;
}

} // end namespace Belos

#endif /* BELOS_STATUS_TEST_RESNORM_H */
