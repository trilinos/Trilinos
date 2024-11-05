// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef BELOS_BICGSTAB_ITER_HPP
#define BELOS_BICGSTAB_ITER_HPP

/*! \file BelosBiCGStabIter.hpp
    \brief Belos concrete class for performing the pseudo-block BiCGStab iteration.
*/

#include "BelosConfigDefs.hpp"
#include "BelosTypes.hpp"
#include "BelosCGIteration.hpp"

#include "BelosLinearProblem.hpp"
#include "BelosMatOrthoManager.hpp"
#include "BelosOutputManager.hpp"
#include "BelosStatusTest.hpp"
#include "BelosOperatorTraits.hpp"
#include "BelosMultiVecTraits.hpp"

#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TimeMonitor.hpp"

/*!
  \class Belos::BiCGStabIter

  \brief This class implements the pseudo-block BiCGStab iteration, where the basic BiCGStab
  algorithm is performed on all of the linear systems simultaneously.

  \ingroup belos_solver_framework

  \author Alicia Klinvex
*/

namespace Belos {

  //! @name BiCGStabIteration Structures
  //@{

  /** \brief Structure to contain pointers to BiCGStabIteration state variables.
   *
   * This struct is utilized by BiCGStabIteration::initialize() and BiCGStabIteration::getState().
   */
  template <class ScalarType, class MV>
  struct BiCGStabIterationState {

    /*! \brief The current residual. */
    Teuchos::RCP<const MV> R;

    /*! \brief The initial residual. */
    Teuchos::RCP<const MV> Rhat;

    /*! \brief The first decent direction vector */
    Teuchos::RCP<const MV> P;

    /*! \brief A * M * the first decent direction vector */
    Teuchos::RCP<const MV> V;

    std::vector<ScalarType> rho_old, alpha, omega;

    BiCGStabIterationState() : R(Teuchos::null), Rhat(Teuchos::null),
                    P(Teuchos::null), V(Teuchos::null)
    {
      rho_old.clear();
      alpha.clear();
      omega.clear();
    }
  };

  template<class ScalarType, class MV, class OP>
  class BiCGStabIter : virtual public Iteration<ScalarType,MV,OP> {

  public:

    //
    // Convenience typedefs
    //
    typedef MultiVecTraits<ScalarType,MV> MVT;
    typedef OperatorTraits<ScalarType,MV,OP> OPT;
    typedef Teuchos::ScalarTraits<ScalarType> SCT;
    typedef typename SCT::magnitudeType MagnitudeType;
    typedef Teuchos::ScalarTraits<MagnitudeType> MT;

    //! @name Constructors/Destructor
    //@{

    /*! \brief %BiCGStabIter constructor with linear problem, solver utilities, and parameter list of solver options.
     *
     * This constructor takes pointers required by the linear solver, in addition
     * to a parameter list of options for the linear solver.
     */
    BiCGStabIter( const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem,
                          const Teuchos::RCP<OutputManager<ScalarType> > &printer,
                          const Teuchos::RCP<StatusTest<ScalarType,MV,OP> > &tester,
                          Teuchos::ParameterList &params );

    //! Destructor.
    virtual ~BiCGStabIter() {};
    //@}


    //! @name Solver methods
    //@{

    /*! \brief This method performs BiCGStab iterations on each linear system until the status
     * test indicates the need to stop or an error occurs (in which case, an
     * std::exception is thrown).
     *
     * iterate() will first determine whether the solver is initialized; if
     * not, it will call initialize() using default arguments. After
     * initialization, the solver performs BiCGStab iterations until the
     * status test evaluates as ::Passed, at which point the method returns to
     * the caller.
     *
     * The status test is queried at the beginning of the iteration.
     *
     */
    void iterate();

    /*! \brief Initialize the solver to an iterate, providing a complete state.
     *
     * The %BiCGStabIter contains a certain amount of state, consisting of the current
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
    void initializeBiCGStab(BiCGStabIterationState<ScalarType,MV>& newstate);

    /*! \brief Initialize the solver with the initial vectors from the linear problem
     *  or random data.
     */
    void initialize()
    {
      BiCGStabIterationState<ScalarType,MV> empty;
      initializeBiCGStab(empty);
    }

    /*! \brief Get the current state of the linear solver.
     *
     * The data is only valid if isInitialized() == \c true.
     *
     * \returns A BiCGStabIterationState object containing const pointers to the current
     * solver state.
     */
    BiCGStabIterationState<ScalarType,MV> getState() const {
      BiCGStabIterationState<ScalarType,MV> state;
      state.R = R_;
      state.Rhat = Rhat_;
      state.P = P_;
      state.V = V_;
      state.rho_old = rho_old_;
      state.alpha = alpha_;
      state.omega = omega_;
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
    // amk TODO: are the residuals actually being set?  What is a native residual?
    Teuchos::RCP<const MV> getNativeResiduals( std::vector<MagnitudeType> * /* norms */ ) const { return R_; }

    //! Get the current update to the linear system.
    /*! \note This method returns a null pointer because the linear problem is current.
    */
    // amk TODO: what is this supposed to be doing?
    Teuchos::RCP<MV> getCurrentUpdate() const { return Teuchos::null; }

    //! Has breakdown been detected in any linear system.
    bool breakdownDetected() { return breakdown_; }

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
                         "Belos::BiCGStabIter::setBlockSize(): Cannot use a block size that is not one.");
    }

    //! States whether the solver has been initialized or not.
    bool isInitialized() { return initialized_; }

    //@}

  private:

    void axpy(const ScalarType alpha, const MV & A,
              const std::vector<ScalarType> beta, const MV& B, MV& mv, bool minus=false);

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

    // Breakdown has been observed for at least one of the linear systems
    bool breakdown_;

    // Current number of iterations performed.
    int iter_;

    //
    // State Storage
    //
    // Initial residual
    Teuchos::RCP<MV> Rhat_;
    //
    // Residual
    Teuchos::RCP<MV> R_;
    //
    // Direction vector 1
    Teuchos::RCP<MV> P_;
    //
    // Operator applied to preconditioned direction vector 1
    Teuchos::RCP<MV> V_;
    //
    std::vector<ScalarType> rho_old_, alpha_, omega_;
  };

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor.
  template<class ScalarType, class MV, class OP>
  BiCGStabIter<ScalarType,MV,OP>::BiCGStabIter(const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem,
                                                               const Teuchos::RCP<OutputManager<ScalarType> > &printer,
                                                               const Teuchos::RCP<StatusTest<ScalarType,MV,OP> > &tester,
                                                               Teuchos::ParameterList &/* params */ ):
    lp_(problem),
    om_(printer),
    stest_(tester),
    numRHS_(0),
    initialized_(false),
    breakdown_(false),
    iter_(0)
  {
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Initialize this iteration object
  template <class ScalarType, class MV, class OP>
  void BiCGStabIter<ScalarType,MV,OP>::initializeBiCGStab(BiCGStabIterationState<ScalarType,MV>& newstate)
  {
    // Check if there is any multivector to clone from.
    Teuchos::RCP<const MV> lhsMV = lp_->getCurrLHSVec();
    Teuchos::RCP<const MV> rhsMV = lp_->getCurrRHSVec();
    TEUCHOS_TEST_FOR_EXCEPTION((lhsMV==Teuchos::null && rhsMV==Teuchos::null),std::invalid_argument,
                       "Belos::BiCGStabIter::initialize(): Cannot initialize state storage!");

    // Get the multivector that is not null.
    Teuchos::RCP<const MV> tmp = ( (rhsMV!=Teuchos::null)? rhsMV: lhsMV );

    // Get the number of right-hand sides we're solving for now.
    int numRHS = MVT::GetNumberVecs(*tmp);
    numRHS_ = numRHS;

    // Initialize the state storage
    // If the subspace has not be initialized before or has changed sizes, generate it using the LHS or RHS from lp_.
    if (Teuchos::is_null(R_) || MVT::GetNumberVecs(*R_)!=numRHS_) {
      R_ = MVT::Clone( *tmp, numRHS_ );
      Rhat_ = MVT::Clone( *tmp, numRHS_ );
      P_ = MVT::Clone( *tmp, numRHS_ );
      V_ = MVT::Clone( *tmp, numRHS_ );

      rho_old_.resize(numRHS_);
      alpha_.resize(numRHS_);
      omega_.resize(numRHS_);
    }

    // Reset breakdown to false before initializing iteration
    breakdown_ = false;

    // NOTE:  In BiCGStabIter R_, the initial residual, is required!!!
    //
    std::string errstr("Belos::BlockPseudoCGIter::initialize(): Specified multivectors must have a consistent length and width.");

    // Create convenience variable for one.
    const ScalarType one = SCT::one();

    if (!Teuchos::is_null(newstate.R)) {

      TEUCHOS_TEST_FOR_EXCEPTION( MVT::GetGlobalLength(*newstate.R) != MVT::GetGlobalLength(*R_),
                          std::invalid_argument, errstr );
      TEUCHOS_TEST_FOR_EXCEPTION( MVT::GetNumberVecs(*newstate.R) != numRHS_,
                          std::invalid_argument, errstr );

      // Copy residual vectors from newstate into R
      if (newstate.R != R_) {
        // Assigned by the new state
        MVT::Assign(*newstate.R, *R_);
      }
      else {
        // Computed
        lp_->computeCurrResVec(R_.get());
      }

      // Set Rhat
      if (!Teuchos::is_null(newstate.Rhat) && newstate.Rhat != Rhat_) {
        // Assigned by the new state
        MVT::Assign(*newstate.Rhat, *Rhat_);
      }
      else {
        // Set to be the initial residual
        MVT::Assign(*R_, *Rhat_);
      }

      // Set V
      if (!Teuchos::is_null(newstate.V) && newstate.V != V_) {
        // Assigned by the new state
        MVT::Assign(*newstate.V, *V_);
      }
      else {
        // Initial V = 0
        MVT::MvInit(*V_);
      }

      // Set P
      if (!Teuchos::is_null(newstate.P) && newstate.P != P_) {
        // Assigned by the new state
        MVT::Assign(*newstate.P, *P_);
      }
      else {
        // Initial P = 0
        MVT::MvInit(*P_);
      }

      // Set rho_old
      if (newstate.rho_old.size () == static_cast<size_t> (numRHS_)) {
        // Assigned by the new state
        rho_old_ = newstate.rho_old;
      }
      else {
        // Initial rho = 1
        rho_old_.assign(numRHS_,one);
      }

      // Set alpha
      if (newstate.alpha.size() == static_cast<size_t> (numRHS_)) {
        // Assigned by the new state
        alpha_ = newstate.alpha;
      }
      else {
        // Initial rho = 1
        alpha_.assign(numRHS_,one);
      }

      // Set omega
      if (newstate.omega.size() == static_cast<size_t> (numRHS_)) {
        // Assigned by the new state
        omega_ = newstate.omega;
      }
      else {
        // Initial rho = 1
        omega_.assign(numRHS_,one);
      }

    }
    else {

      TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::is_null(newstate.R),std::invalid_argument,
                         "Belos::BiCGStabIter::initialize(): BiCGStabStateIterState does not have initial residual.");
    }

    // The solver is initialized
    initialized_ = true;
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Iterate until the status test informs us we should stop.
  template <class ScalarType, class MV, class OP>
  void BiCGStabIter<ScalarType,MV,OP>::iterate()
  {
    using Teuchos::RCP;

    //
    // Allocate/initialize data structures
    //
    if (initialized_ == false) {
      initialize();
    }

    // Allocate memory for scalars.
    int i=0;
    std::vector<ScalarType> rho_new( numRHS_ ), beta( numRHS_ );
    std::vector<ScalarType> rhatV( numRHS_ ), tT( numRHS_ ), tS( numRHS_ );

    // Create convenience variable for one.
    const ScalarType one = SCT::one();

    // TODO: We may currently be using more space than is required
    RCP<MV> leftPrecVec, leftPrecVec2;

    RCP<MV> Y, Z, S, T;
    S = MVT::Clone( *R_, numRHS_ );
    T = MVT::Clone( *R_, numRHS_ );
    if (lp_->isLeftPrec() || lp_->isRightPrec()) {
      Y = MVT::Clone( *R_, numRHS_ );
      Z = MVT::Clone( *R_, numRHS_ );
    }
    else {
      Y = P_;
      Z = S;
    }

    // Get the current solution std::vector.
    Teuchos::RCP<MV> X = lp_->getCurrLHSVec();

    ////////////////////////////////////////////////////////////////
    // Iterate until the status test tells us to stop.
    //
    while (stest_->checkStatus(this) != Passed && !breakdown_) {

      // Increment the iteration
      iter_++;

      // rho_new = <R_, Rhat_>
      MVT::MvDot(*R_,*Rhat_,rho_new);

      // beta = ( rho_new / rho_old ) (alpha / omega )
      // TODO: None of these loops are currently threaded
      for(i=0; i<numRHS_; i++) {
        // Catch breakdown in rho_old here, since
        // it is just rho_new from the previous iteration.
        if (SCT::magnitude(rho_new[i]) < MT::sfmin())
          breakdown_ = true;

        beta[i] = (rho_new[i] / rho_old_[i]) * (alpha_[i] / omega_[i]);
      }

      // p = r + beta (p - omega v)
      // TODO: Is it safe to call MvAddMv with A or B = mv?
      // TODO: Not all of these things have to be part of the state
      axpy(one, *P_, omega_, *V_, *P_, true); // p = p - omega v
      axpy(one, *R_, beta, *P_, *P_); // p = r + beta (p - omega v)

      // y = K\p, unless K does not exist
      // TODO: There may be a more efficient way to apply the preconditioners
      if(lp_->isLeftPrec()) {
        if(lp_->isRightPrec()) {
          if(leftPrecVec == Teuchos::null) {
            leftPrecVec = MVT::Clone( *R_, numRHS_ );
          }
          lp_->applyLeftPrec(*P_,*leftPrecVec);
          lp_->applyRightPrec(*leftPrecVec,*Y);
        }
        else {
          lp_->applyLeftPrec(*P_,*Y);
        }
      }
      else if(lp_->isRightPrec()) {
        lp_->applyRightPrec(*P_,*Y);
      }

      // v = Ay
      lp_->applyOp(*Y,*V_);

      // alpha = rho_new / <Rhat, V>
      MVT::MvDot(*V_,*Rhat_,rhatV);
      for(i=0; i<numRHS_; i++) {
        if (SCT::magnitude(rhatV[i]) < MT::sfmin())
        {
          breakdown_ = true;
          return;
        }
        else 
          alpha_[i] = rho_new[i] / rhatV[i];
      }

      // s = r - alpha v
      axpy(one, *R_, alpha_, *V_, *S, true);

      // z = K\s, unless K does not exist
      if(lp_->isLeftPrec()) {
        if(lp_->isRightPrec()) {
          if(leftPrecVec == Teuchos::null) {
            leftPrecVec = MVT::Clone( *R_, numRHS_ );
          }
          lp_->applyLeftPrec(*S,*leftPrecVec);
          lp_->applyRightPrec(*leftPrecVec,*Z);
        }
        else {
          lp_->applyLeftPrec(*S,*Z);
        }
      }
      else if(lp_->isRightPrec()) {
        lp_->applyRightPrec(*S,*Z);
      }

      // t = Az
      lp_->applyOp(*Z,*T);

      // omega = <K1\t,K1\s> / <K1\t,K1\t>
      if(lp_->isLeftPrec()) {
        if(leftPrecVec == Teuchos::null) {
          leftPrecVec = MVT::Clone( *R_, numRHS_ );
        }
        if(leftPrecVec2 == Teuchos::null) {
          leftPrecVec2 = MVT::Clone( *R_, numRHS_ );
        }
        lp_->applyLeftPrec(*T,*leftPrecVec2);
        MVT::MvDot(*leftPrecVec2,*leftPrecVec2,tT);
        MVT::MvDot(*leftPrecVec,*leftPrecVec2,tS);
      }
      else {
        MVT::MvDot(*T,*T,tT);
        MVT::MvDot(*S,*T,tS);
      }
      for(i=0; i<numRHS_; i++) {
        if (SCT::magnitude(tT[i]) < MT::sfmin())
        {
          omega_[i] = SCT::zero();
          breakdown_ = true;
        }
        else
          omega_[i] = tS[i] / tT[i];
      }

      // x = x + alpha y + omega z
      axpy(one, *X, alpha_, *Y, *X); // x = x + alpha y
      axpy(one, *X, omega_, *Z, *X); // x = x + alpha y + omega z

      // r = s - omega t
      axpy(one, *S, omega_, *T, *R_, true);

      // Update rho_old
      rho_old_ = rho_new;
    } // end while (sTest_->checkStatus(this) != Passed)
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Iterate until the status test informs us we should stop.
  template <class ScalarType, class MV, class OP>
  void BiCGStabIter<ScalarType,MV,OP>::axpy(const ScalarType alpha, const MV & A,
                                            const std::vector<ScalarType> beta, const MV& B, MV& mv, bool minus)
  {
    Teuchos::RCP<const MV> A1, B1;
    Teuchos::RCP<MV> mv1;
    std::vector<int> index(1);

    for(int i=0; i<numRHS_; i++) {
      index[0] = i;
      A1 = MVT::CloneView(A,index);
      B1 = MVT::CloneView(B,index);
      mv1 = MVT::CloneViewNonConst(mv,index);
      if(minus) {
        MVT::MvAddMv(alpha,*A1,-beta[i],*B1,*mv1);
      }
      else {
        MVT::MvAddMv(alpha,*A1,beta[i],*B1,*mv1);
      }
    }
  }

} // end Belos namespace

#endif /* BELOS_BICGSTAB_ITER_HPP */
