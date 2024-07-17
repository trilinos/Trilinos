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
  
template<class ScalarType, class MV, class OP>
class CGIter : virtual public CGIteration<ScalarType,MV,OP> {

  public:
    
  //
  // Convenience typedefs
  //
  typedef MultiVecTraits<ScalarType,MV> MVT;
  typedef OperatorTraits<ScalarType,MV,OP> OPT;
  typedef Teuchos::ScalarTraits<ScalarType> SCT;
  typedef typename SCT::magnitudeType MagnitudeType;

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
  virtual ~CGIter() {};
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
  void initializeCG(CGIterationState<ScalarType,MV>& newstate);

  /*! \brief Initialize the solver with the initial vectors from the linear problem
   *  or random data.
   */
  void initialize()
  {
    CGIterationState<ScalarType,MV> empty;
    initializeCG(empty);
  }
  
  /*! \brief Get the current state of the linear solver.
   *
   * The data is only valid if isInitialized() == \c true.
   *
   * \returns A CGIterationState object containing const pointers to the current solver state.
   */
  CGIterationState<ScalarType,MV> getState() const {
    CGIterationState<ScalarType,MV> state;
    state.R = R_;
    state.P = P_;
    state.AP = AP_;
    state.Z = Z_;
    return state;
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
   if (numEntriesForCondEst_) doCondEst_=val;
  }
  
  //! Gets the diagonal for condition estimation
  Teuchos::ArrayView<MagnitudeType> getDiag() {
    // NOTE (mfh 30 Jul 2015) See note on getOffDiag() below.
    // getDiag() didn't actually throw for me in that case, but why
    // not be cautious?
    typedef typename Teuchos::ArrayView<MagnitudeType>::size_type size_type;
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
    typedef typename Teuchos::ArrayView<MagnitudeType>::size_type size_type;
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
  //! Method for initalizing the state storage needed by CG.
  void setStateSize();
  
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

  // stateStorageInitialized_ specifies that the state storage has been initialized.
  // This initialization may be postponed if the linear problem was generated without 
  // the right-hand side or solution vectors.
  bool stateStorageInitialized_;

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
    stateStorageInitialized_(false),
    iter_(0),
    assertPositiveDefiniteness_( params.get("Assert Positive Definiteness", true) ),
    numEntriesForCondEst_(params.get("Max Size For Condest",0) ),
    doCondEst_(false)
  {
    foldConvergenceDetectionIntoAllreduce_ = params.get<bool>("Fold Convergence Detection Into Allreduce",false);
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Setup the state storage.
  template <class ScalarType, class MV, class OP>
  void CGIter<ScalarType,MV,OP>::setStateSize ()
  {
    if (!stateStorageInitialized_) {

      // Check if there is any multivector to clone from.
      Teuchos::RCP<const MV> lhsMV = lp_->getLHS();
      Teuchos::RCP<const MV> rhsMV = lp_->getRHS();
      if (lhsMV == Teuchos::null && rhsMV == Teuchos::null) {
	stateStorageInitialized_ = false;
	return;
      }
      else {
	
	// Initialize the state storage
	// If the subspace has not be initialized before, generate it using the LHS or RHS from lp_.
	if (R_ == Teuchos::null) {
	  // Get the multivector that is not null.
	  Teuchos::RCP<const MV> tmp = ( (rhsMV!=Teuchos::null)? rhsMV: lhsMV );
	  TEUCHOS_TEST_FOR_EXCEPTION(tmp == Teuchos::null,std::invalid_argument,
			     "Belos::CGIter::setStateSize(): linear problem does not specify multivectors to clone from.");
          S_ = MVT::Clone( *tmp, 2 );
          std::vector<int> index(1,0);
          index[0] = 0;
          R_ = MVT::CloneViewNonConst( *S_, index );
          index[0] = 1;
          Z_ = MVT::CloneViewNonConst( *S_, index );
	  P_ = MVT::Clone( *tmp, 1 );
	  AP_ = MVT::Clone( *tmp, 1 );

        }

        // Tracking information for condition number estimation
        if(numEntriesForCondEst_ > 0) {
          diag_.resize(numEntriesForCondEst_);
          offdiag_.resize(numEntriesForCondEst_-1);
        }
        	
	// State storage has now been initialized.
	stateStorageInitialized_ = true;
      }
    }
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Initialize this iteration object
  template <class ScalarType, class MV, class OP>
  void CGIter<ScalarType,MV,OP>::initializeCG(CGIterationState<ScalarType,MV>& newstate)
  {
    // Initialize the state storage if it isn't already.
    if (!stateStorageInitialized_) 
      setStateSize();

    TEUCHOS_TEST_FOR_EXCEPTION(!stateStorageInitialized_,std::invalid_argument,
		       "Belos::CGIter::initialize(): Cannot initialize state storage!");
    
    // NOTE:  In CGIter R_, the initial residual, is required!!!  
    //
    std::string errstr("Belos::CGIter::initialize(): Specified multivectors must have a consistent length and width.");

    if (newstate.R != Teuchos::null) {

      TEUCHOS_TEST_FOR_EXCEPTION( MVT::GetGlobalLength(*newstate.R) != MVT::GetGlobalLength(*R_),
                          std::invalid_argument, errstr );
      TEUCHOS_TEST_FOR_EXCEPTION( MVT::GetNumberVecs(*newstate.R) != 1,
                          std::invalid_argument, errstr );

      // Copy basis vectors from newstate into V
      if (newstate.R != R_) {
        // copy over the initial residual (unpreconditioned).
	MVT::Assign( *newstate.R, *R_ );
      }

      // Compute initial direction vectors
      // Initially, they are set to the preconditioned residuals
      //
      if ( lp_->getLeftPrec() != Teuchos::null ) {
        lp_->applyLeftPrec( *R_, *Z_ );
        if ( lp_->getRightPrec() != Teuchos::null ) {
          Teuchos::RCP<MV> tmp = MVT::CloneCopy( *Z_ );
          lp_->applyRightPrec( *tmp, *Z_ );
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
    else {

      TEUCHOS_TEST_FOR_EXCEPTION(newstate.R == Teuchos::null,std::invalid_argument,
                         "Belos::CGIter::initialize(): CGIterationState does not have initial residual.");
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
    if (initialized_ == false) {
      initialize();
    }

    // Allocate memory for scalars.
    std::vector<ScalarType> alpha(1);
    std::vector<ScalarType> beta(1);
    std::vector<ScalarType> rHz(1), rHz_old(1), pAp(1);
    Teuchos::SerialDenseMatrix<int,ScalarType> rHs( 1, 2 );

    // Create convenience variables for zero and one.
    const ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();
    const MagnitudeType zero = Teuchos::ScalarTraits<MagnitudeType>::zero();

    // Scalars for condition estimation (if needed) - These will always use entry zero, for convenience
    ScalarType pAp_old = one, beta_old = one, rHz_old2 = one;
             
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
