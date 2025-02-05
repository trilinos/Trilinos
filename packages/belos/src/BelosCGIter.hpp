// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef BELOS_CG_ITER_HPP
#define BELOS_CG_ITER_HPP

/*! \file BelosCGIter.hpp
    \brief Belos concrete class for performing the conjugate-gradient (CG) iteration.
*/

#include "BelosConfigDefs.hpp"
#include "BelosTypes.hpp"
#include "BelosCGIteration.hpp"

#include "BelosLinearProblem.hpp"
#include "BelosOutputManager.hpp"
#include "BelosStatusTest.hpp"
#include "BelosStatusTestGenResNorm.hpp"
#include "BelosOperatorTraits.hpp"
#include "BelosMultiVecTraits.hpp"

#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TimeMonitor.hpp"

/*!
  \class Belos::CGIter

  \brief This class implements the preconditioned Conjugate Gradient (CG) iteration.

  \ingroup belos_solver_framework

  \author Teri Barth and Heidi Thornquist
*/

namespace Belos {

  //! @name CGIteration Structures
  //@{

  /** \brief Structure to contain pointers to CGIteration state variables.
   *
   * This struct is utilized by CGIteration::initialize() and CGIteration::getState().
   */
  template <class ScalarType, class MV>
  class CGIterationState : public CGIterationStateBase<ScalarType, MV> {

  public:
    CGIterationState() = default;

    CGIterationState(Teuchos::RCP<const MV> tmp) {
      initialize(tmp);
    }

    virtual ~CGIterationState() = default;

    void initialize(Teuchos::RCP<const MV> tmp, int _numVectors) {
      using MVT = MultiVecTraits<ScalarType, MV>;
      TEUCHOS_ASSERT(_numVectors == 1);

      // S = (R, Z)
      // This allows to compute the inner products (R, S) = ((R, R), (R, Z)) using a single reduction.
      S = MVT::Clone( *tmp, 2 );
      std::vector<int> index(1,0);
      index[0] = 0;
      this->R = MVT::CloneViewNonConst( *S, index );
      index[0] = 1;
      this->Z = MVT::CloneViewNonConst( *S, index );

      this->P = MVT::Clone( *tmp, 1 );
      this->AP = MVT::Clone(*tmp, 1);

      CGIterationStateBase<ScalarType, MV>::initialize(tmp, _numVectors);
    }

    bool matches(Teuchos::RCP<const MV> tmp, int _numVectors=1) const {
      return (CGIterationStateBase<ScalarType, MV>::matches(tmp, _numVectors) &&
              !S.is_null());
    }

    Teuchos::RCP<MV> S;

};

template<class ScalarType, class MV, class OP>
class CGIter : virtual public CGIteration<ScalarType,MV,OP> {

  public:

  //
  // Convenience typedefs
  //
  using MVT = MultiVecTraits<ScalarType, MV>;
  using OPT = OperatorTraits<ScalarType, MV, OP>;
  using SCT = Teuchos::ScalarTraits<ScalarType>;
  using MagnitudeType = typename SCT::magnitudeType;

  //! @name Constructors/Destructor
  //@{

  /*! \brief %CGIter constructor with linear problem, solver utilities, and parameter list of solver options.
   *
   * This constructor takes pointers required by the linear solver iteration, in addition
   * to a parameter list of options for the linear solver.
   */
  CGIter( const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem,
		  const Teuchos::RCP<OutputManager<ScalarType> > &printer,
		  const Teuchos::RCP<StatusTest<ScalarType,MV,OP> > &tester,
                  const Teuchos::RCP<StatusTestGenResNorm<ScalarType,MV,OP> > &convTester,
		  Teuchos::ParameterList &params );

  //! Destructor.
  virtual ~CGIter() = default;
  //@}


  //! @name Solver methods
  //@{

  /*! \brief This method performs CG iterations until the status
   * test indicates the need to stop or an error occurs (in which case, an
   * std::exception is thrown).
   *
   * iterate() will first determine whether the solver is initialized; if
   * not, it will call initialize() using default arguments. After
   * initialization, the solver performs CG iterations until the
   * status test evaluates as ::Passed, at which point the method returns to
   * the caller.
   *
   * The status test is queried at the beginning of the iteration.
   */
  void iterate();

  /*! \brief Initialize the solver to an iterate, providing a complete state.
   *
   * The %CGIter contains a certain amount of state, consisting of the current
   * residual, preconditioned residual, and decent direction.
   *
   * initialize() gives the user the opportunity to manually set these,
   * although only the current unpreconditioned residual is required.
   *
   * \post
   * <li>isInitialized() == \c true (see post-conditions of isInitialize())
   *
   * \note For any pointer in \c newstate which directly points to the multivectors in
   * the solver, the data is not copied.
   */
  void initializeCG(Teuchos::RCP<CGIterationStateBase<ScalarType,MV> > newstate, Teuchos::RCP<MV> R_0);

  /*! \brief Initialize the solver with the initial vectors from the linear problem
   *  or random data.
   */
  void initialize()
  {
    initializeCG(Teuchos::null, Teuchos::null);
  }

  /*! \brief Get the current state of the linear solver.
   *
   * The data is only valid if isInitialized() == \c true.
   *
   * \returns A CGIterationState object containing const pointers to the current solver state.
   */
  Teuchos::RCP<CGIterationStateBase<ScalarType,MV> > getState() const {
    auto state = Teuchos::rcp(new CGIterationState<ScalarType,MV>());
    state->R = R_;
    state->P = P_;
    state->AP = AP_;
    state->Z = Z_;
    state->S = S_;
    return state;
  }

  void setState(Teuchos::RCP<CGIterationStateBase<ScalarType,MV> > state) {
    auto s = Teuchos::rcp_dynamic_cast<CGIterationState<ScalarType,MV> >(state, true);
    R_ = s->R;
    Z_ = s->Z;
    P_ = s->P;
    AP_ = s->AP;
    S_ = s->S;
  }

  //@}


  //! @name Status methods
  //@{

  //! \brief Get the current iteration count.
  int getNumIters() const { return iter_; }

  //! \brief Reset the iteration count.
  void resetNumIters( int iter = 0 ) { iter_ = iter; }

  //! Get the norms of the residuals native to the solver.
  //! \return A std::vector of length blockSize containing the native residuals.
  Teuchos::RCP<const MV> getNativeResiduals( std::vector<MagnitudeType> * norms ) const {
    if (foldConvergenceDetectionIntoAllreduce_ && convTest_->getResNormType() == Belos::TwoNorm) {
      (*norms)[0] = std::sqrt(Teuchos::ScalarTraits<ScalarType>::magnitude(rHr_));
      return Teuchos::null;
    } else
      return R_;
  }

  //! Get the current update to the linear system.
  /*! \note This method returns a null pointer because the linear problem is current.
  */
  Teuchos::RCP<MV> getCurrentUpdate() const { return Teuchos::null; }

  //@}

  //! @name Accessor methods
  //@{

  //! Get a constant reference to the linear problem.
  const LinearProblem<ScalarType,MV,OP>& getProblem() const { return *lp_; }

  //! Get the blocksize to be used by the iterative solver in solving this linear problem.
  int getBlockSize() const { return 1; }

  //! \brief Set the blocksize to be used by the iterative solver in solving this linear problem.
  void setBlockSize(int blockSize) {
    TEUCHOS_TEST_FOR_EXCEPTION(blockSize!=1,std::invalid_argument,
		       "Belos::CGIter::setBlockSize(): Cannot use a block size that is not one.");
  }

  //! States whether the solver has been initialized or not.
  bool isInitialized() { return initialized_; }


  //! Sets whether or not to store the diagonal for condition estimation
  void setDoCondEst(bool val) {
   if (numEntriesForCondEst_ != 0) doCondEst_=val;
  }

  //! Gets the diagonal for condition estimation
  Teuchos::ArrayView<MagnitudeType> getDiag() {
    // NOTE (mfh 30 Jul 2015) See note on getOffDiag() below.
    // getDiag() didn't actually throw for me in that case, but why
    // not be cautious?
    using size_type = typename Teuchos::ArrayView<MagnitudeType>::size_type;
    if (static_cast<size_type> (iter_) >= diag_.size ()) {
      return diag_ ();
    } else {
      return diag_ (0, iter_);
    }
    }

  //! Gets the off-diagonal for condition estimation
  Teuchos::ArrayView<MagnitudeType> getOffDiag() {
    // NOTE (mfh 30 Jul 2015) The implementation as I found it
    // returned "offdiag(0,iter_)".  This breaks (Teuchos throws in
    // debug mode) when the maximum number of iterations has been
    // reached, because iter_ == offdiag_.size() in that case.  The
    // new logic fixes this.
    using size_type = typename Teuchos::ArrayView<MagnitudeType>::size_type;
    if (static_cast<size_type> (iter_) >= offdiag_.size ()) {
      return offdiag_ ();
    } else {
      return offdiag_ (0, iter_);
    }
  }

  //@}

  private:

  //
  // Internal methods
  //
  // Classes inputed through constructor that define the linear problem to be solved.
  //
  const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> >    lp_;
  const Teuchos::RCP<OutputManager<ScalarType> >          om_;
  const Teuchos::RCP<StatusTest<ScalarType,MV,OP> >       stest_;
  const Teuchos::RCP<StatusTestGenResNorm<ScalarType,MV,OP> >       convTest_;

  //
  // Current solver state
  //
  // initialized_ specifies that the basis vectors have been initialized and the iterate() routine
  // is capable of running; _initialize is controlled  by the initialize() member method
  // For the implications of the state of initialized_, please see documentation for initialize()
  bool initialized_;

  // Current number of iterations performed.
  int iter_;

  // Should the allreduce for convergence detection be merged with one of the inner products?
  bool foldConvergenceDetectionIntoAllreduce_;

  // <r,r>
  ScalarType rHr_;

  // Assert that the matrix is positive definite
  bool assertPositiveDefiniteness_;

  // Tridiagonal system for condition estimation (if needed)
  Teuchos::ArrayRCP<MagnitudeType> diag_, offdiag_;
  int numEntriesForCondEst_;
  bool doCondEst_;



  //
  // State Storage
  //
  // Residual
  Teuchos::RCP<MV> R_;
  //
  // Preconditioned residual
  Teuchos::RCP<MV> Z_;
  //
  // Direction vector
  Teuchos::RCP<MV> P_;
  //
  // Operator applied to direction vector
  Teuchos::RCP<MV> AP_;

  Teuchos::RCP<MV> S_;

};

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor.
  template<class ScalarType, class MV, class OP>
  CGIter<ScalarType,MV,OP>::CGIter(const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem,
						   const Teuchos::RCP<OutputManager<ScalarType> > &printer,
						   const Teuchos::RCP<StatusTest<ScalarType,MV,OP> > &tester,
                                                   const Teuchos::RCP<StatusTestGenResNorm<ScalarType,MV,OP> > &convTester,
						   Teuchos::ParameterList &params ):
    lp_(problem),
    om_(printer),
    stest_(tester),
    convTest_(convTester),
    initialized_(false),
    iter_(0),
    assertPositiveDefiniteness_( params.get("Assert Positive Definiteness", true) ),
    numEntriesForCondEst_(params.get("Max Size For Condest",0) ),
    doCondEst_(false)
  {
    foldConvergenceDetectionIntoAllreduce_ = params.get<bool>("Fold Convergence Detection Into Allreduce",false);
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Initialize this iteration object
  template <class ScalarType, class MV, class OP>
  void CGIter<ScalarType,MV,OP>::initializeCG(Teuchos::RCP<CGIterationStateBase<ScalarType,MV> > newstate, Teuchos::RCP<MV> R_0)
  {
    // Initialize the state storage if it isn't already.
    Teuchos::RCP<const MV> lhsMV = lp_->getLHS();
    Teuchos::RCP<const MV> rhsMV = lp_->getRHS();
    Teuchos::RCP<const MV> tmp = ( (rhsMV!=Teuchos::null)? rhsMV: lhsMV );
    TEUCHOS_ASSERT(!newstate.is_null());
    if (!Teuchos::rcp_dynamic_cast<CGIterationState<ScalarType,MV> >(newstate, true)->matches(tmp, 1))
      newstate->initialize(tmp, 1);
    setState(newstate);

    // Tracking information for condition number estimation
    if(numEntriesForCondEst_ > 0) {
      diag_.resize(numEntriesForCondEst_);
      offdiag_.resize(numEntriesForCondEst_-1);
    }

    std::string errstr("Belos::CGIter::initialize(): Specified multivectors must have a consistent length and width.");
    {

      TEUCHOS_TEST_FOR_EXCEPTION( MVT::GetGlobalLength(*R_0) != MVT::GetGlobalLength(*R_),
                          std::invalid_argument, errstr );
      TEUCHOS_TEST_FOR_EXCEPTION( MVT::GetNumberVecs(*R_0) != 1,
                          std::invalid_argument, errstr );

      // Copy basis vectors from newstate into V
      if (R_0 != R_) {
        // copy over the initial residual (unpreconditioned).
	MVT::Assign( *R_0, *R_ );
      }

      // Compute initial direction vectors
      // Initially, they are set to the preconditioned residuals
      //
      if ( lp_->getLeftPrec() != Teuchos::null ) {
        lp_->applyLeftPrec( *R_, *Z_ );
        if ( lp_->getRightPrec() != Teuchos::null ) {
          Teuchos::RCP<MV> tmp2 = MVT::CloneCopy( *Z_ );
          lp_->applyRightPrec( *tmp2, *Z_ );
        }
      }
      else if ( lp_->getRightPrec() != Teuchos::null ) {
        lp_->applyRightPrec( *R_, *Z_ );
      }
      else {
        MVT::Assign( *R_, *Z_ );
      }
      MVT::Assign( *Z_, *P_ );
    }

    // The solver is initialized
    initialized_ = true;
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Iterate until the status test informs us we should stop.
  template <class ScalarType, class MV, class OP>
  void CGIter<ScalarType,MV,OP>::iterate()
  {
    //
    // Allocate/initialize data structures
    //
    if (!initialized_) {
      initialize();
    }

    // Allocate memory for scalars.
    std::vector<ScalarType> alpha(1);
    std::vector<ScalarType> beta(1);
    std::vector<ScalarType> rHz(1);
    std::vector<ScalarType> rHz_old(1);
    std::vector<ScalarType> pAp(1);
    Teuchos::SerialDenseMatrix<int,ScalarType> rHs( 1, 2 );

    // Create convenience variables for zero and one.
    const ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();
    const MagnitudeType zero = Teuchos::ScalarTraits<MagnitudeType>::zero();

    // Scalars for condition estimation (if needed) - These will always use entry zero, for convenience
    ScalarType pAp_old = one;
    ScalarType beta_old = one;
    ScalarType rHz_old2 = one;

    // Get the current solution vector.
    Teuchos::RCP<MV> cur_soln_vec = lp_->getCurrLHSVec();

    // Check that the current solution vector only has one column.
    TEUCHOS_TEST_FOR_EXCEPTION( MVT::GetNumberVecs(*cur_soln_vec) != 1, CGIterateFailure,
                        "Belos::CGIter::iterate(): current linear system has more than one vector!" );

    // Compute first <r,z> a.k.a. rHz
    if (foldConvergenceDetectionIntoAllreduce_ && convTest_->getResNormType() == Belos::TwoNorm) {
      MVT::MvTransMv( one, *R_, *S_, rHs );
      rHr_ = rHs(0,0);
      rHz[0] = rHs(0,1);
    } else
      MVT::MvDot( *R_, *Z_, rHz );

    ////////////////////////////////////////////////////////////////
    // Iterate until the status test tells us to stop.
    //
    while (stest_->checkStatus(this) != Passed) {

      // Increment the iteration
      iter_++;

      // Multiply the current direction vector by A and store in AP_
      lp_->applyOp( *P_, *AP_ );

      // Compute alpha := <R_,Z_> / <P_,AP_>
      MVT::MvDot( *P_, *AP_, pAp );
      alpha[0] = rHz[0] / pAp[0];

      // Check that alpha is a positive number!
      if(assertPositiveDefiniteness_) {
        TEUCHOS_TEST_FOR_EXCEPTION( SCT::real(alpha[0]) <= zero, CGPositiveDefiniteFailure,
                                    "Belos::CGIter::iterate(): non-positive value for p^H*A*p encountered!" );
      }

      //
      // Update the solution vector x := x + alpha * P_
      //
      MVT::MvAddMv( one, *cur_soln_vec, alpha[0], *P_, *cur_soln_vec );
      lp_->updateSolution();
      //
      // Save the denominator of beta before residual is updated [ old <R_, Z_> ]
      //
      rHz_old[0] = rHz[0];
      //
      // Compute the new residual R_ := R_ - alpha * AP_
      //
      MVT::MvAddMv( one, *R_, -alpha[0], *AP_, *R_ );
      //
      // Compute beta := [ new <R_, Z_> ] / [ old <R_, Z_> ],
      // and the new direction vector p.
      //
      if ( lp_->getLeftPrec() != Teuchos::null ) {
        lp_->applyLeftPrec( *R_, *Z_ );
        if ( lp_->getRightPrec() != Teuchos::null ) {
          Teuchos::RCP<MV> tmp = MVT::CloneCopy( *Z_);
          lp_->applyRightPrec( *tmp, *Z_ );
        }
      }
      else if ( lp_->getRightPrec() != Teuchos::null ) {
        lp_->applyRightPrec( *R_, *Z_ );
      }
      else {
        MVT::Assign( *R_, *Z_ );
      }
      //
      if (foldConvergenceDetectionIntoAllreduce_ && convTest_->getResNormType() == Belos::TwoNorm) {
        MVT::MvTransMv( one, *R_, *S_, rHs );
        rHr_ = rHs(0,0);
        rHz[0] = rHs(0,1);
      } else
        MVT::MvDot( *R_, *Z_, rHz );
      //
      beta[0] = rHz[0] / rHz_old[0];
      //
      MVT::MvAddMv( one, *Z_, beta[0], *P_, *P_ );

      // Condition estimate (if needed)
      if (doCondEst_) {
        if (iter_ > 1) {
          diag_[iter_-1]    = Teuchos::ScalarTraits<ScalarType>::real((beta_old * beta_old * pAp_old + pAp[0]) / rHz_old[0]);
          offdiag_[iter_-2] = -Teuchos::ScalarTraits<ScalarType>::real(beta_old * pAp_old / (sqrt( rHz_old[0] * rHz_old2)));
        }
        else {
          diag_[iter_-1]    = Teuchos::ScalarTraits<ScalarType>::real(pAp[0] / rHz_old[0]);
        }
        rHz_old2 = rHz_old[0];
        beta_old = beta[0];
        pAp_old = pAp[0];
      }

    } // end while (sTest_->checkStatus(this) != Passed)
  }

} // end Belos namespace

#endif /* BELOS_CG_ITER_HPP */
