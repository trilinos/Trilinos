// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef BELOS_PSEUDO_BLOCK_CG_ITER_MP_VECTOR_HPP
#define BELOS_PSEUDO_BLOCK_CG_ITER_MP_VECTOR_HPP

/*!
  \file BelosPseudoBlockCGIter.hpp
  \brief Belos concrete class for performing the pseudo-block CG iteration.

  Specialized for Sacado::MP::Vector scalar type to handle break-down when
  ensemble components are over-solved.
*/

#include "Stokhos_Sacado_Kokkos_MP_Vector.hpp"
#include "BelosPseudoBlockCGIter.hpp"

/*!
  \class Belos::PseudoBlockCGIter

  \brief This class implements the pseudo-block CG iteration, where the basic CG
  algorithm is performed on all of the linear systems simultaneously.

  \ingroup belos_solver_framework

  \author Heidi Thornquist

   This implementation is the same as the default, but is specialized for
   the Sacado::MP::Vector scalar type to handle break-down when
   ensemble components are over-solved.
*/

namespace Belos {

  template<class Storage, class MV, class OP>
  class PseudoBlockCGIter<Sacado::MP::Vector<Storage>, MV, OP> :
    virtual public CGIteration<Sacado::MP::Vector<Storage>, MV, OP> {

  public:

    //
    // Convenience typedefs
    //
    typedef Sacado::MP::Vector<Storage> ScalarType;
    typedef MultiVecTraits<ScalarType,MV> MVT;
    typedef OperatorTraits<ScalarType,MV,OP> OPT;
    typedef Teuchos::ScalarTraits<ScalarType> SCT;
    typedef typename SCT::magnitudeType MagnitudeType;
    typedef Teuchos::ScalarTraits<typename Storage::value_type> SVT;

    //! @name Constructors/Destructor
    //@{

    /*! \brief %PseudoBlockCGIter constructor with linear problem, solver utilities, and parameter list of solver options.
     *
     * This constructor takes pointers required by the linear solver, in addition
     * to a parameter list of options for the linear solver.
     */
    PseudoBlockCGIter( const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem,
                          const Teuchos::RCP<OutputManager<ScalarType> > &printer,
                          const Teuchos::RCP<StatusTest<ScalarType,MV,OP> > &tester,
                          Teuchos::ParameterList &params );

    //! Destructor.
    virtual ~PseudoBlockCGIter() {};
    //@}


    //! @name Solver methods
    //@{

    /*! \brief This method performs CG iterations on each linear system until the status
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
     *
     */
    void iterate();

    /*! \brief Initialize the solver to an iterate, providing a complete state.
     *
     * The %PseudoBlockCGIter contains a certain amount of state, consisting of the current
     * direction vectors and residuals.
     *
     * initialize() gives the user the opportunity to manually set these,
     * although this must be done with caution, abiding by the rules given
     * below.
     *
     * \post
     * <li>isInitialized() == \c true (see post-conditions of isInitialize())
     *
     * The user has the option of specifying any component of the state using
     * initialize(). However, these arguments are assumed to match the
     * post-conditions specified under isInitialized(). Any necessary component of the
     * state not given to initialize() will be generated.
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
     * \returns A CGIterationState object containing const pointers to the current
     * solver state.
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
    Teuchos::RCP<const MV> getNativeResiduals( std::vector<MagnitudeType> *norms ) const { return R_; }

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

    //! \brief Set the blocksize.
    void setBlockSize(int blockSize) {
      TEUCHOS_TEST_FOR_EXCEPTION(blockSize!=1,std::invalid_argument,
                         "Belos::PseudoBlockCGIter::setBlockSize(): Cannot use a block size that is not one.");
    }

    //! States whether the solver has been initialized or not.
    bool isInitialized() { return initialized_; }

    //@}

    //! Sets whether or not to store the diagonal for condition estimation
    void setDoCondEst(bool val){doCondEst_=val;}

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

  private:

    //
    // Classes inputed through constructor that define the linear problem to be solved.
    //
    const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> >    lp_;
    const Teuchos::RCP<OutputManager<ScalarType> >          om_;
    const Teuchos::RCP<StatusTest<ScalarType,MV,OP> >       stest_;

    //
    // Algorithmic parameters
    //
    // numRHS_ is the current number of linear systems being solved.
    int numRHS_;

    //
    // Current solver state
    //
    // initialized_ specifies that the basis vectors have been initialized and the iterate() routine
    // is capable of running; _initialize is controlled  by the initialize() member method
    // For the implications of the state of initialized_, please see documentation for initialize()
    bool initialized_;

    // Current number of iterations performed.
    int iter_;

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

    // Tolerance for ensemble breakdown
    typename SVT::magnitudeType breakDownTol_;

  };

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor.
  template<class StorageType, class MV, class OP>
  PseudoBlockCGIter<Sacado::MP::Vector<StorageType>,MV,OP>::
  PseudoBlockCGIter(const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem,
                                                               const Teuchos::RCP<OutputManager<ScalarType> > &printer,
                                                               const Teuchos::RCP<StatusTest<ScalarType,MV,OP> > &tester,
                                                               Teuchos::ParameterList &params ):
    lp_(problem),
    om_(printer),
    stest_(tester),
    numRHS_(0),
    initialized_(false),
    iter_(0),
    assertPositiveDefiniteness_( params.get("Assert Positive Definiteness", true) ),
    numEntriesForCondEst_(params.get("Max Size For Condest",0) ),
    breakDownTol_(params.get("Ensemble Breakdown Tolerance", 0.0))
  {
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Initialize this iteration object
  template <class StorageType, class MV, class OP>
  void PseudoBlockCGIter<Sacado::MP::Vector<StorageType>,MV,OP>::
  initializeCG(CGIterationState<ScalarType,MV>& newstate)
  {
    // Check if there is any mltivector to clone from.
    Teuchos::RCP<const MV> lhsMV = lp_->getCurrLHSVec();
    Teuchos::RCP<const MV> rhsMV = lp_->getCurrRHSVec();
    TEUCHOS_TEST_FOR_EXCEPTION((lhsMV==Teuchos::null && rhsMV==Teuchos::null),std::invalid_argument,
                       "Belos::PseudoBlockCGIter::initialize(): Cannot initialize state storage!");

    // Get the multivector that is not null.
    Teuchos::RCP<const MV> tmp = ( (rhsMV!=Teuchos::null)? rhsMV: lhsMV );

    // Get the number of right-hand sides we're solving for now.
    int numRHS = MVT::GetNumberVecs(*tmp);
    numRHS_ = numRHS;

    // Initialize the state storage
    // If the subspace has not be initialized before or has changed sizes, generate it using the LHS or RHS from lp_.
    if (Teuchos::is_null(R_) || MVT::GetNumberVecs(*R_)!=numRHS_) {
      R_ = MVT::Clone( *tmp, numRHS_ );
      Z_ = MVT::Clone( *tmp, numRHS_ );
      P_ = MVT::Clone( *tmp, numRHS_ );
      AP_ = MVT::Clone( *tmp, numRHS_ );
    }

    // Tracking information for condition number estimation
    if(numEntriesForCondEst_ > 0) {
      diag_.resize(numEntriesForCondEst_);
      offdiag_.resize(numEntriesForCondEst_-1);
    }

    // NOTE:  In CGIter R_, the initial residual, is required!!!
    //
    std::string errstr("Belos::BlockPseudoCGIter::initialize(): Specified multivectors must have a consistent length and width.");

    // Create convenience variables for zero and one.
    const ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();
    const MagnitudeType zero = Teuchos::ScalarTraits<MagnitudeType>::zero();

    if (!Teuchos::is_null(newstate.R)) {

      TEUCHOS_TEST_FOR_EXCEPTION( MVT::GetGlobalLength(*newstate.R) != MVT::GetGlobalLength(*R_),
                          std::invalid_argument, errstr );
      TEUCHOS_TEST_FOR_EXCEPTION( MVT::GetNumberVecs(*newstate.R) != numRHS_,
                          std::invalid_argument, errstr );

      // Copy basis vectors from newstate into V
      if (newstate.R != R_) {
        // copy over the initial residual (unpreconditioned).
        MVT::MvAddMv( one, *newstate.R, zero, *newstate.R, *R_ );
      }

      // Compute initial direction vectors
      // Initially, they are set to the preconditioned residuals
      //
      if ( lp_->getLeftPrec() != Teuchos::null ) {
        lp_->applyLeftPrec( *R_, *Z_ );
        if ( lp_->getRightPrec() != Teuchos::null ) {
          Teuchos::RCP<MV> tmp1 = MVT::Clone( *Z_, numRHS_ );
          lp_->applyRightPrec( *Z_, *tmp1 );
          Z_ = tmp1;
        }
      }
      else if ( lp_->getRightPrec() != Teuchos::null ) {
        lp_->applyRightPrec( *R_, *Z_ );
      }
      else {
        Z_ = R_;
      }
      MVT::MvAddMv( one, *Z_, zero, *Z_, *P_ );
    }
    else {

      TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::is_null(newstate.R),std::invalid_argument,
                         "Belos::CGIter::initialize(): CGStateIterState does not have initial residual.");
    }

    // The solver is initialized
    initialized_ = true;
  }


 //////////////////////////////////////////////////////////////////////////////////////////////////
  // Iterate until the status test informs us we should stop.
  template <class StorageType, class MV, class OP>
  void PseudoBlockCGIter<Sacado::MP::Vector<StorageType>,MV,OP>::
  iterate()
  {
    //
    // Allocate/initialize data structures
    //
    if (initialized_ == false) {
      initialize();
    }

    // Allocate memory for scalars.
    int i=0;
    std::vector<int> index(1);
    std::vector<ScalarType> rHz( numRHS_ ), rHz_old( numRHS_ ), pAp( numRHS_ );
    Teuchos::SerialDenseMatrix<int, ScalarType> alpha( numRHS_,numRHS_ ), beta( numRHS_,numRHS_ );

    // Create convenience variables for zero and one.
    const ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();
    const MagnitudeType zero = Teuchos::ScalarTraits<MagnitudeType>::zero();

    // Scalars for condition estimation (if needed) - These will always use entry zero, for convenience
    ScalarType pAp_old=one, beta_old=one ,rHz_old2=one;

    // Get the current solution std::vector.
    Teuchos::RCP<MV> cur_soln_vec = lp_->getCurrLHSVec();

    // Compute first <r,z> a.k.a. rHz
    MVT::MvDot( *R_, *Z_, rHz );

    if ( assertPositiveDefiniteness_ )
        for (i=0; i<numRHS_; ++i)
            TEUCHOS_TEST_FOR_EXCEPTION( SCT::real(rHz[i]) < zero,
                                CGIterateFailure,
                                "Belos::PseudoBlockCGIter::iterate(): negative value for r^H*M*r encountered!" );

    ////////////////////////////////////////////////////////////////
    // Iterate until the status test tells us to stop.
    //
    while (stest_->checkStatus(this) != Passed) {

      // Increment the iteration
      iter_++;

      // Multiply the current direction std::vector by A and store in AP_
      lp_->applyOp( *P_, *AP_ );

      // Compute alpha := <R_,Z_> / <P_,AP_>
      MVT::MvDot( *P_, *AP_, pAp );

      for (i=0; i<numRHS_; ++i) {
        if ( assertPositiveDefiniteness_ )
            // Check that pAp[i] is a positive number!
            TEUCHOS_TEST_FOR_EXCEPTION( SCT::real(pAp[i]) <= zero,
                                CGIterateFailure,
                                "Belos::PseudoBlockCGIter::iterate(): non-positive value for p^H*A*p encountered!" );

        const int ensemble_size = pAp[i].size();
        for (int k=0; k<ensemble_size; ++k) {
          if (SVT::magnitude(pAp[i].coeff(k)) <= breakDownTol_)
            alpha(i,i).fastAccessCoeff(k) = 0.0;
          else
            alpha(i,i).fastAccessCoeff(k) = rHz[i].coeff(k) / pAp[i].coeff(k);
        }
      }

      //
      // Update the solution std::vector x := x + alpha * P_
      //
      MVT::MvTimesMatAddMv( one, *P_, alpha, one, *cur_soln_vec );
      lp_->updateSolution();
      //
      // Save the denominator of beta before residual is updated [ old <R_, Z_> ]
      //
      for (i=0; i<numRHS_; ++i) {
        rHz_old[i] = rHz[i];
      }
      //
      // Compute the new residual R_ := R_ - alpha * AP_
      //
      MVT::MvTimesMatAddMv( -one, *AP_, alpha, one, *R_ );
      //
      // Compute beta := [ new <R_, Z_> ] / [ old <R_, Z_> ],
      // and the new direction std::vector p.
      //
      if ( lp_->getLeftPrec() != Teuchos::null ) {
        lp_->applyLeftPrec( *R_, *Z_ );
        if ( lp_->getRightPrec() != Teuchos::null ) {
          Teuchos::RCP<MV> tmp = MVT::Clone( *Z_, numRHS_ );
          lp_->applyRightPrec( *Z_, *tmp );
          Z_ = tmp;
        }
      }
      else if ( lp_->getRightPrec() != Teuchos::null ) {
        lp_->applyRightPrec( *R_, *Z_ );
      }
      else {
        Z_ = R_;
      }
      //
      MVT::MvDot( *R_, *Z_, rHz );
      if ( assertPositiveDefiniteness_ )
          for (i=0; i<numRHS_; ++i)
              TEUCHOS_TEST_FOR_EXCEPTION( SCT::real(rHz[i]) < zero,
                                  CGIterateFailure,
                                  "Belos::PseudoBlockCGIter::iterate(): negative value for r^H*M*r encountered!" );
      //
      // Update the search directions.
      for (i=0; i<numRHS_; ++i) {
        const int ensemble_size = rHz[i].size();
        for (int k=0; k<ensemble_size; ++k) {
          if (SVT::magnitude(rHz[i].coeff(k)) <= breakDownTol_)
            beta(i,i).fastAccessCoeff(k) = 0.0;
          else
            beta(i,i).fastAccessCoeff(k) = rHz[i].coeff(k) / rHz_old[i].coeff(k);
        }
        index[0] = i;
        Teuchos::RCP<const MV> Z_i = MVT::CloneView( *Z_, index );
        Teuchos::RCP<MV> P_i = MVT::CloneViewNonConst( *P_, index );
        MVT::MvAddMv( one, *Z_i, beta(i,i), *P_i, *P_i );
      }

      // Condition estimate (if needed)
      if(doCondEst_ > 0) {
        if(iter_ > 1 ) {
          diag_[iter_-1]    = Teuchos::ScalarTraits<ScalarType>::real((beta_old * beta_old * pAp_old + pAp[0]) / rHz_old[0]);
          offdiag_[iter_-2] = -Teuchos::ScalarTraits<ScalarType>::real(beta_old * pAp_old / (sqrt( rHz_old[0] * rHz_old2)));
        }
        else {
          diag_[iter_-1]    = Teuchos::ScalarTraits<ScalarType>::real(pAp[0] / rHz_old[0]);
        }
        rHz_old2 = rHz_old[0];
        beta_old = beta(0,0);
        pAp_old = pAp[0];
      }


      //
    } // end while (sTest_->checkStatus(this) != Passed)
  }

} // end Belos namespace

#endif /* BELOS_PSEUDO_BLOCK_CG_ITER_HPP */
