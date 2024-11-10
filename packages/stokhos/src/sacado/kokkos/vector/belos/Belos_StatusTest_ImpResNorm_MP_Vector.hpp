// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef BELOS_STATUS_TEST_GEN_IMPRESNORM_MP_VECTOR_HPP
#define BELOS_STATUS_TEST_GEN_IMPRESNORM_MP_VECTOR_HPP

/*!
  \file Belos_StatusTest_ImpResNorm_MP_Vector.hpp
  \brief Belos::StatusTest for specifying an implicit residual norm stopping criteria that checks for loss of accuracy.

  Specialized for Sacado::MP::Vector scalar type to track how many iterations
  each component of ensemble takes.
*/

#include "Stokhos_Sacado_Kokkos_MP_Vector.hpp"
#include "BelosStatusTestImpResNorm.hpp"

/*!
  \class Belos::StatusTestImpResNorm
  \brief Convergence test using the implicit residual norm(s), with an
    explicit residual norm(s) check for loss of accuracy if necessary.

  This test passes when a quorum of residual vectors (for all
  right-hand sides) passes a residual norm test.  For residual vector
  \f$r_i\f$, the form of the test is \f$\|r_i\| / \sigma_i \le
  \tau\f$, where \f$\tau\f$ is the convergence tolerance set in the
  constructor or by setTolerance().  The default quorum consists of
  all the vectors, but you can set the number of vectors that must
  pass before the whole test passes, either by the constructor
  argument or by calling setQuorum().

  The default residual vector norm (the norm used in the numerator of
  the test, and applied to the current residual vectors) is the
  2-norm, but you can change the norm type (1-norm, 2-norm, or
  infinity norm) using defineResForm().  The default scaling factor
  \f$\sigma_i\f$ for residual vector i is the 2-norm of the initial
  residual vector, but you can change this using defineScaleForm().
  The norm used for the scaling factors may be different than the norm
  used for the current residual vectors.  Furthermore, each residual
  vector may use a different scaling factor \f$\sigma_i\f$.

  This test starts by using the "implicit" residual norm, otherwise
  called the "recursive" or "native" residual norm.  It comes from the
  iterative method's projected problem, unlike the "explicit" residual
  norm, which comes from explicitly computing \f$R = B - A * X\f$.
  Once the implicit residual norm reaches the convergence tolerance,
  this test then checks the explicit residual norm to make sure that
  it has reached the convergence tolerance as well.

  We say that there is a "potential loss of accuracy" when we first
  detect that the implicit residual norm(s) have met the desired
  original convergence tolerance, but the explicit residual norm(s)
  have not.  We don't actually call this a "loss of accuracy" (in the
  sense of getLOADetected() returning true) unless we tried to remedy
  this situation, but couldn't fix it.  Upon detecting a potential
  loss of accuracy, this test tells the solver to "iterate a few more
  steps" by making the "current" residual tolerance (the value of
  getCurrTolerance()) smaller, and forcing the implicit residual
  norm(s) to meet the current tolerance before moving on to test the
  explicit norm(s).  If that doesn't work, this test declares a loss
  of accuracy.

  This implementation is the same as the default, but is specialized for
  the Sacado::MP::Vector scalar type to keep track of when each component
  of the ensemble converges.

*/

namespace Belos {

template <class Storage, class MV, class OP>
class StatusTestImpResNorm<Sacado::MP::Vector<Storage>, MV, OP> :
    public StatusTestResNorm<Sacado::MP::Vector<Storage>,MV,OP> {

public:
  typedef Sacado::MP::Vector<Storage> ScalarType;
  //! The type of the magnitude (absolute value) of a ScalarType.
  typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;

private:
  //! @name Abbreviations for method implementations
  //@{
  typedef Teuchos::ScalarTraits<ScalarType> STS;
  typedef Teuchos::ScalarTraits<MagnitudeType> STM;
  typedef MultiVecTraits<ScalarType,MV> MVT;
  //@}

public:
  //! @name Constructors and destructor.
  //@{

  /// \brief Constructor
  ///
  /// \param Tolerance [in] Convergence tolerance \f$\tau\f$.
  ///
  /// \param quorum [in] The number of vectors in the problem that
  ///   must pass the convergence criteria before this StatusTest
  ///   passes.  -1 means that all residuals must pass.
  ///
  /// \param showMaxResNormOnly [in] Whether only the maximum residual
  ///   norm (of all vectors) is displayed when the print() method is
  ///   called.
  StatusTestImpResNorm (MagnitudeType Tolerance,
                        int quorum = -1,
                        bool showMaxResNormOnly = false);

  //! Destructor (virtual for memory safety).
  virtual ~StatusTestImpResNorm();

  //@}
  //! @name Form and parameter definition methods.
  //@{

  //! Define form of the residual, its norm and optional weighting vector.
  /*! This method defines the form of \f$\|r\|\f$.  We specify:
    <ul>
    <li> The norm to be used on the residual (this may be different than the norm used in
    defineScaleForm()).
    </ul>
  */
  int defineResForm( NormType TypeOfNorm);

  //! Define form of the scaling, its norm, its optional weighting vector, or, alternatively, define an explicit value.
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
    <li> Define TypeOfScaling to be the norm of the right hand side, the initial residual vector,
    or to none.
    <li> Define norm to be used on the scaling vector (this may be different than the norm used
    in DefineResForm()).
    </ul>
    </ol>
  */
  int defineScaleForm( ScaleType TypeOfScaling, NormType TypeOfNorm, MagnitudeType ScaleValue = Teuchos::ScalarTraits<MagnitudeType>::one());

  //! Set the value of the tolerance
  /*! We allow the tolerance to be reset for cases where, in the process of testing the residual,
    we find that the initial tolerance was too tight or too lax.
  */
  int setTolerance (MagnitudeType tolerance) {
    tolerance_ = tolerance;
    return 0;
  }

  //! Sets the number of residuals that must pass the convergence test before Passed is returned.
  //! \note If \c quorum=-1 then all residuals must pass the convergence test before Passed is returned.
  int setQuorum (int quorum) {
    quorum_ = quorum;
    return 0;
  }

  //! Set whether the only maximum residual norm is displayed when the print() method is called
  int setShowMaxResNormOnly (bool showMaxResNormOnly) {
    showMaxResNormOnly_ = showMaxResNormOnly;
    return 0;
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
  StatusType checkStatus(Iteration<ScalarType,MV,OP>* iSolver);

  //! Return the result of the most recent CheckStatus call.
  StatusType getStatus() const {return(status_);}
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
  Teuchos::RCP<MV> getSolution() { return curSoln_; }

  //! Returns the number of residuals that must pass the convergence test before Passed is returned.
  //! \note If \c quorum=-1 then all residuals must pass the convergence test before Passed is returned.
  int getQuorum() const { return quorum_; }

  //! Returns whether the only maximum residual norm is displayed when the print() method is called
  bool getShowMaxResNormOnly() { return showMaxResNormOnly_; }

  //! Returns the vector containing the indices of the residuals that passed the test.
  std::vector<int> convIndices() { return ind_; }

  /// \brief "Original" convergence tolerance \f$\tau\f$ as set by user.
  ///
  /// This value is the convergence tolerance as set by the user in
  /// the constructor, or by the setTolerance() method.  See this
  /// class' main documentation and the documentation of
  /// getCurrTolerance() for an explanation of the difference between
  /// the "original" and "current" tolerances.
  MagnitudeType getTolerance () const {
    return tolerance_;
  }

  /// \brief Current convergence tolerance; may be changed to prevent loss of accuracy.
  ///
  /// The difference between "original" tolerance (the value of
  /// getTolerance()) and "current" tolerance (the value of this
  /// method) relates to the idea of "loss of accuracy."  See this
  /// class' main documentation for details.  We do not allow users to
  /// set the "current" tolerance.
  MagnitudeType getCurrTolerance () const {
    return currTolerance_;
  }

  //! Returns the test value, \f$ \frac{\|r\|}{\sigma} \f$, computed in most recent call to CheckStatus.
  const std::vector<MagnitudeType>* getTestValue() const {return(&testvector_);}

  //! Returns the residual norm value, \f$ \|r\| \f$, computed in most recent call to CheckStatus.
  const std::vector<MagnitudeType>* getResNormValue() const {return(&resvector_);}

  //! Returns the scaled norm value, \f$ \sigma \f$.
  const std::vector<MagnitudeType>* getScaledNormValue() const {return(&scalevector_);}

  //! Returns a boolean indicating a loss of accuracy has been detected in computing the residual.
  bool getLOADetected() const { return lossDetected_; }

  //! Returns number of ensemble iterations
  const std::vector<int> getEnsembleIterations() const { return ensemble_iterations; }
  //@}


  /** @name Misc. */
  //@{

  /** \brief Call to setup initial scaling vector.
    *
    * After this function is called <tt>getScaledNormValue()</tt> can be called
    * to get the scaling vector.
    */
  StatusType firstCallCheckStatusSetup(Iteration<ScalarType,MV,OP>* iSolver);
  //@}

  /** \name Overridden from Teuchos::Describable */
  //@{

  /** \brief Method to return description of the maximum iteration status test  */
  std::string description() const
  {
    std::ostringstream oss;
    oss << "Belos::StatusTestImpResNorm<>: " << resFormStr();
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
        if (scaletype_==NormOfInitRes)
          oss << " Res0";
        else if (scaletype_==NormOfPrecInitRes)
          oss << " Prec Res0";
        else
          oss << " RHS ";
        oss << ")";
      }
    }

    return oss.str();
  }

  //! @name Private data members.
  //@{

  //! Current tolerance used to determine convergence and the default tolerance set by user.
  MagnitudeType tolerance_, currTolerance_;

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

  //! Scaling vector.
  std::vector<MagnitudeType> scalevector_;

  //! Residual norm vector.
  std::vector<MagnitudeType> resvector_;

  //! Test vector = resvector_ / scalevector_
  std::vector<MagnitudeType> testvector_;

  //! Current solution vector
  Teuchos::RCP<MV> curSoln_;

  //! Vector containing the indices for the vectors that passed the test.
  std::vector<int> ind_;

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

  //! Has a loss in accuracy been detected?
  bool lossDetected_;

  //! The number of iterations at which point each ensemble component converges
  std::vector<int> ensemble_iterations;

  //@}

};  

template <class StorageType, class MV, class OP>
StatusTestImpResNorm<Sacado::MP::Vector<StorageType>,MV,OP>::
StatusTestImpResNorm (MagnitudeType Tolerance, int quorum, bool showMaxResNormOnly)
  : tolerance_(Tolerance),
    currTolerance_(Tolerance),
    quorum_(quorum),
    showMaxResNormOnly_(showMaxResNormOnly),
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
    firstcallDefineScaleForm_(true),
    lossDetected_(false),
    ensemble_iterations(StorageType::static_size, 0)
{
  // This constructor will compute the residual ||r_i||/||r0_i|| <= tolerance using the 2-norm of
  // the implicit residual vector.
}

template <class StorageType, class MV, class OP>
StatusTestImpResNorm<Sacado::MP::Vector<StorageType>,MV,OP>::~StatusTestImpResNorm()
{}

template <class StorageType, class MV, class OP>
void StatusTestImpResNorm<Sacado::MP::Vector<StorageType>,MV,OP>::reset()
{
  status_ = Undefined;
  curBlksz_ = 0;
  curLSNum_ = 0;
  curLSIdx_.resize(0);
  numrhs_ = 0;
  ind_.resize(0);
  currTolerance_ = tolerance_;
  firstcallCheckStatus_ = true;
  lossDetected_ = false;
  curSoln_ = Teuchos::null;
  ensemble_iterations = std::vector<int>(StorageType::static_size, 0);
}

template <class StorageType, class MV, class OP>
int StatusTestImpResNorm<Sacado::MP::Vector<StorageType>,MV,OP>::defineResForm( NormType TypeOfNorm )
{
  TEUCHOS_TEST_FOR_EXCEPTION(firstcallDefineResForm_==false,StatusTestError,
        "StatusTestResNorm::defineResForm(): The residual form has already been defined.");
  firstcallDefineResForm_ = false;

  resnormtype_ = TypeOfNorm;

  return(0);
}

template <class StorageType, class MV, class OP>
int StatusTestImpResNorm<Sacado::MP::Vector<StorageType>,MV,OP>::defineScaleForm(ScaleType TypeOfScaling, NormType TypeOfNorm,
                                                          MagnitudeType ScaleValue )
{
  TEUCHOS_TEST_FOR_EXCEPTION(firstcallDefineScaleForm_==false,StatusTestError,
        "StatusTestResNorm::defineScaleForm(): The scaling type has already been defined.");
  firstcallDefineScaleForm_ = false;

  scaletype_ = TypeOfScaling;
  scalenormtype_ = TypeOfNorm;
  scalevalue_ = ScaleValue;

  return(0);
}

template <class StorageType, class MV, class OP>
StatusType StatusTestImpResNorm<Sacado::MP::Vector<StorageType>,MV,OP>::
checkStatus (Iteration<ScalarType,MV,OP>* iSolver)
{
  //std::cout << "check status start" << std::endl;
  using Teuchos::as;
  using Teuchos::RCP;

  const MagnitudeType zero = STM::zero ();
  const MagnitudeType one  = STM::one ();
  const LinearProblem<ScalarType,MV,OP>& lp = iSolver->getProblem ();

  // Compute scaling term (done once for each block that's being solved)
  if (firstcallCheckStatus_) {
    StatusType status = firstCallCheckStatusSetup (iSolver);
    if (status == Failed) {
      status_ = Failed;
      return status_;
    }
  }

  // mfh 23 Apr 2012: I don't know exactly what this code does.  It
  // has something to do with picking the block of right-hand sides
  // which we're currently checking.
  if (curLSNum_ != lp.getLSNumber ()) {
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
  } else {
    //
    // We are in the same rhs block, return if we are converged
    //
    if (status_ == Passed) {
      //std::cout << "check status end 2" << std::endl;
      return status_;
    }
  }

  //
  // Get the "native" residual norms from the solver for this block of
  // right-hand sides.  If the solver's getNativeResiduals() method
  // actually returns a multivector, compute the norms of the columns
  // of the multivector explicitly.  Otherwise, we assume that
  // resvector_ contains the norms.
  //
  // Note that "compute the norms explicitly" doesn't necessarily mean
  // the "explicit" residual norms (in the sense discussed in this
  // class' documentation).  These are just some vectors returned by
  // the solver.  Some Krylov methods, like CG, compute a residual
  // vector "recursively."  This is an "implicit residual" in the
  // sense of this class' documentation.  It equals the explicit
  // residual in exact arithmetic, but due to rounding error, it is
  // usually different than the explicit residual.
  //
  // FIXME (mfh 23 Apr 2012) This method does _not_ respect the
  // OrthoManager used by the solver.
  //
  std::vector<MagnitudeType> tmp_resvector( curBlksz_ );
  RCP<const MV> residMV = iSolver->getNativeResiduals (&tmp_resvector);
  if (! residMV.is_null ()) {
    // We got a multivector back.  Compute the norms explicitly.
    tmp_resvector.resize (MVT::GetNumberVecs (*residMV));
    MVT::MvNorm (*residMV, tmp_resvector, resnormtype_);
    typename std::vector<int>::iterator p = curLSIdx_.begin();
    for (int i=0; p<curLSIdx_.end(); ++p, ++i) {
      // Check if this index is valid
      if (*p != -1) {
        resvector_[*p] = tmp_resvector[i];
      }
    }
  } else {
    typename std::vector<int>::iterator p = curLSIdx_.begin();
    for (int i=0; p<curLSIdx_.end(); ++p, ++i) {
      // Check if this index is valid
      if (*p != -1) {
        resvector_[*p] = tmp_resvector[i];
      }
    }
  }
  //
  // Scale the unscaled residual norms we computed or obtained above.
  //
  if (scalevector_.size () > 0) {
    // There are per-vector scaling factors to apply.
    typename std::vector<int>::iterator p = curLSIdx_.begin();
    for (; p<curLSIdx_.end(); ++p) {
      // Check if this index is valid
      if (*p != -1) {
        // Scale the vector accordingly
        testvector_[ *p ] = resvector_[ *p ] / scalevalue_;
        mask_assign(scalevector_[ *p ]!= zero, testvector_[ *p ]) /= scalevector_[ *p ];
      }
    }
  }
  else { // There are no per-vector scaling factors.
    typename std::vector<int>::iterator p = curLSIdx_.begin();
    for (; p<curLSIdx_.end(); ++p) {
      // Check if this index is valid
      if (*p != -1) {
        testvector_[ *p ] = resvector_[ *p ] / scalevalue_;
      }
    }
  }

  // Count how many scaled residual norms (in testvector_) pass, using
  // the current tolerance (currTolerance_) rather than the original
  // tolerance (tolerance_).  If at least quorum_ of them pass, we
  // have a quorum for the whole test to pass.
  //
  // We also check here whether any of the scaled residual norms is
  // NaN, and throw an exception in that case.
  int have = 0;
  ind_.resize( curLSIdx_.size() );
  std::vector<int> lclInd( curLSIdx_.size() );
  typename std::vector<int>::iterator p = curLSIdx_.begin();
  for (int i=0; p<curLSIdx_.end(); ++p, ++i) {
    // Check if this index is valid
    if (*p != -1) {
      // Check if any of the residuals are larger than the tolerance.
      const int ensemble_size = StorageType::static_size;
      bool all_converged = true;
      for (int ii=0; ii<ensemble_size; ++ii) {
        if (testvector_[ *p ].coeff(ii) > currTolerance_.coeff(ii)) {
          ++ensemble_iterations[ii];
          all_converged = false;
        }
        else if (!(testvector_[ *p ].coeff(ii) <= currTolerance_.coeff(ii))) {
          // Throw an std::exception if the current residual norm is
          // NaN.  We know that it's NaN because it is not less than,
          // equal to, or greater than the current tolerance.  This is
          // only possible if either the residual norm or the current
          // tolerance is NaN; we assume the former.  We also mark the
          // test as failed, in case you want to catch the exception.
          status_ = Failed;
          TEUCHOS_TEST_FOR_EXCEPTION(true, StatusTestError, "Belos::"
            "StatusTestImpResNorm::checkStatus(): One or more of the current "
            "implicit residual norms is NaN.");
        }
      }
      if (all_converged) {
        ind_[have] = *p;
        have++;
      }
    }
  }
  // "have" is the number of residual norms that passed.
  ind_.resize(have);
  lclInd.resize(have);

  // Now check the exact residuals
  if (have) { // At least one residual norm has converged.
    //
    // Compute the explicit residual norm(s) from the current solution update.
    //
    RCP<MV> cur_update = iSolver->getCurrentUpdate ();
    curSoln_ = lp.updateSolution (cur_update);
    RCP<MV> cur_res = MVT::Clone (*curSoln_, MVT::GetNumberVecs (*curSoln_));
    lp.computeCurrResVec (&*cur_res, &*curSoln_);
    tmp_resvector.resize (MVT::GetNumberVecs (*cur_res));
    std::vector<MagnitudeType> tmp_testvector (have);

    MVT::MvNorm (*cur_res, tmp_resvector, resnormtype_);

    // Scale the explicit residual norm(s), just like the implicit norm(s).
    if ( scalevector_.size() > 0 ) {
      for (int i=0; i<have; ++i) {
        // Scale the vector accordingly
        tmp_testvector[ i ] = tmp_resvector[ lclInd[i] ] / scalevalue_;
        mask_assign(scalevector_[ ind_[i] ]!=zero,tmp_testvector[ i ]) /= scalevector_[ ind_[i] ];
      }
    }
    else {
      for (int i=0; i<have; ++i) {
        tmp_testvector[ i ] = tmp_resvector[ lclInd[i] ] / scalevalue_;
      }
    }

    //
    // Check whether the explicit residual norms also pass the
    // convergence test.  If not, check whether we want to try
    // iterating a little more to force both implicit and explicit
    // residual norms to pass.
    //
    int have2 = 0;
    for (int i = 0; i < have; ++i) {
      // testvector_ contains the implicit (i.e., recursive, computed
      // by the algorithm) (possibly scaled) residuals.  All of these
      // pass the convergence test.
      //
      // tmp_testvector contains the explicit (i.e., ||B-AX||)
      // (possibly scaled) residuals.  We're checking whether these
      // pass as well.  The explicit residual norms only have to meet
      // the _original_ tolerance (tolerance_), not the current
      // tolerance (currTolerance_).


      // Absolute difference between the current explicit and
      // implicit residual norm.
      const MagnitudeType diff = STM::magnitude (testvector_[ind_[i]] - tmp_testvector[i]);

      // Check if any of the residuals are larger than the tolerance.
      const int ensemble_size = StorageType::static_size;
      typedef typename StorageType::value_type value_type;
      bool all_converged = true;
      {
        for (int ii=0; ii<ensemble_size; ++ii) {
          if (!(tmp_testvector[i].coeff(ii) <= tolerance_.coeff(ii))) {
            if (diff.coeff(ii) > currTolerance_.coeff(ii)) {
              lossDetected_ = true;
            }
            else {
              const value_type onePointFive = as<value_type>(3) / as<value_type> (2);
              const value_type oneTenth = as<value_type>(1) / as<value_type> (10);

              currTolerance_.fastAccessCoeff(ii) = currTolerance_.coeff(ii) - onePointFive * diff.coeff(ii);
              while (currTolerance_.coeff(ii) < as<value_type>(0)) {
                currTolerance_.fastAccessCoeff(ii) += oneTenth * diff.coeff(ii);
              }
              all_converged = false;
            }
          }
        }
      }
      if (all_converged) {
        ind_[have2] = ind_[i];
        have2++; // This right-hand side has converged.
      }
    }
    have = have2;
    ind_.resize(have);
  }

  // Check whether we've met the quorum of vectors necessary for the
  // whole test to pass.
  int need = (quorum_ == -1) ? curNumRHS_: quorum_;
  status_ = (have >= need) ? Passed : Failed;


  //std::cout << "check status end default" << std::endl;
  // Return the current status
  return status_;
}

template <class StorageType, class MV, class OP>
void StatusTestImpResNorm<Sacado::MP::Vector<StorageType>,MV,OP>::print(std::ostream& os, int indent) const
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

template <class StorageType, class MV, class OP>
void StatusTestImpResNorm<Sacado::MP::Vector<StorageType>,MV,OP>::printStatus(std::ostream& os, StatusType type) const
{
  os << std::left << std::setw(13) << std::setfill('.');
  switch (type) {
  case  Passed:
    os << "Converged";
    break;
  case  Failed:
    if (lossDetected_)
      os << "Unconverged (LoA)";
    else
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

template <class StorageType, class MV, class OP>
StatusType StatusTestImpResNorm<Sacado::MP::Vector<StorageType>,MV,OP>::
firstCallCheckStatusSetup (Iteration<ScalarType,MV,OP>* iSolver)
{
  int i;
  const MagnitudeType zero = STM::zero ();
  const MagnitudeType one = STM::one ();
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
    else if (scaletype_==NormOfInitRes) {
      Teuchos::RCP<const MV> init_res = lp.getInitResVec();
      numrhs_ = MVT::GetNumberVecs( *init_res );
      scalevector_.resize( numrhs_ );
      MVT::MvNorm( *init_res, scalevector_, scalenormtype_ );
    }
    else if (scaletype_==NormOfPrecInitRes) {
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

#endif /* BELOS_STATUS_TEST_GEN_IMPRESNORM_MP_VECTOR_HPP */
