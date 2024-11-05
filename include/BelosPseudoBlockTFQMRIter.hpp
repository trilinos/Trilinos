// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
// This file contains an implementation of the TFQMR iteration
// for solving non-Hermitian linear systems of equations Ax = b, 
// where b is a single-vector and x is the corresponding solution.
//
// The implementation is a slight modification on the TFQMR iteration
// found in Saad's "Iterative Methods for Sparse Linear Systems".
//

#ifndef BELOS_PSEUDO_BLOCK_TFQMR_ITER_HPP
#define BELOS_PSEUDO_BLOCK_TFQMR_ITER_HPP

/*!
  \file BelosPseudoBlockTFQMRIter.hpp

  \brief Belos concrete class for generating iterations with the
  preconditioned tranpose-free QMR (TFQMR) method.
*/

#include "BelosConfigDefs.hpp"
#include "BelosIteration.hpp"
#include "BelosTypes.hpp"	

#include "BelosLinearProblem.hpp"
#include "BelosOutputManager.hpp"
#include "BelosStatusTest.hpp"
#include "BelosOperatorTraits.hpp"
#include "BelosMultiVecTraits.hpp"

#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TimeMonitor.hpp"

/*!	\class Belos::PseudoBlockTFQMRIter

	\brief This class implements the preconditioned transpose-free QMR algorithm for
	solving non-Hermitian linear systems of equations Ax = b, where b is the right-hand 
	side vector and x is the corresponding solution.

        \ingroup belos_solver_framework

	\author Heidi Thornquist
*/

namespace Belos {

  /** \brief Structure to contain pointers to PseudoBlockTFQMRIter state variables.
   *
   * This struct is utilized by PseudoBlockTFQMRIter::initialize() and PseudoBlockTFQMRIter::getState().
   */
  template <class ScalarType, class MV>
  struct PseudoBlockTFQMRIterState {
 
    typedef Teuchos::ScalarTraits<ScalarType> SCT;
    typedef typename SCT::magnitudeType MagnitudeType;

    /*! \brief The current residual basis. */
    Teuchos::RCP<const MV> W;
    Teuchos::RCP<const MV> U;
    Teuchos::RCP<const MV> AU;
    Teuchos::RCP<const MV> Rtilde;
    Teuchos::RCP<const MV> D;
    Teuchos::RCP<const MV> V;
    std::vector<ScalarType> alpha, eta, rho;
    std::vector<MagnitudeType> tau, theta;


    PseudoBlockTFQMRIterState() : W(Teuchos::null), U(Teuchos::null), AU(Teuchos::null),
                                  Rtilde(Teuchos::null), D(Teuchos::null), V(Teuchos::null)
    {}
  };

  template <class ScalarType, class MV, class OP>
  class PseudoBlockTFQMRIter : public Iteration<ScalarType,MV,OP> { 
  public:
    //
    // Convenience typedefs
    //
    typedef MultiVecTraits<ScalarType,MV> MVT;
    typedef OperatorTraits<ScalarType,MV,OP> OPT;
    typedef Teuchos::ScalarTraits<ScalarType> SCT;
    typedef typename SCT::magnitudeType MagnitudeType;
    
    //! @name Constructor/Destructor.
    //@{ 

    //! %Belos::PseudoBlockTFQMRIter constructor.
    PseudoBlockTFQMRIter( const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem, 
	       const Teuchos::RCP<OutputManager<ScalarType> > &printer,
	       const Teuchos::RCP<StatusTest<ScalarType,MV,OP> > &tester,
	       Teuchos::ParameterList &params );
    
    //! %Belos::PseudoBlockTFQMRIter destructor.
    virtual ~PseudoBlockTFQMRIter() {};
    //@}
  

    //! @name Solver methods
    //@{ 
    
    /*! \brief This method performs pseudo-block TFQMR iterations until the status
     * test indicates the need to stop or an error occurs (in which case, an
     * std::exception is thrown).
     *
     * iterate() will first determine whether the solver is inintialized; if
     * not, it will call initialize() using default arguments. After
     * initialization, the solver performs pseudo-block TFQMR iterations until the
     * status test evaluates as ::Passed, at which point the method returns to
     * the caller. 
     */
    void iterate();

    /*! \brief Initialize the solver to an iterate, providing a complete state.
     *
     * The %PseudoBlockTFQMRIter contains a certain amount of state, consisting of the current 
     * Krylov basis and the associated Hessenberg matrix.
     *
     * initialize() gives the user the opportunity to manually set these,
     * although this must be done with caution, abiding by the rules given
     * below. All notions of orthogonality and orthonormality are derived from
     * the inner product specified by the orthogonalization manager.
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
    void initializeTFQMR(const PseudoBlockTFQMRIterState<ScalarType,MV> & newstate);
    
    /*! \brief Initialize the solver with the initial vectors from the linear problem
     *  or random data.
     */
    void initialize()
    {
      PseudoBlockTFQMRIterState<ScalarType,MV> empty;
      initializeTFQMR(empty);
    }
    
    /*! \brief Get the current state of the linear solver.
     *
     * The data is only valid if isInitialized() == \c true.
     *
     * \returns A PseudoBlockTFQMRIterState object containing const pointers to the current
     * solver state.
     */
    PseudoBlockTFQMRIterState<ScalarType,MV> getState() const {
      PseudoBlockTFQMRIterState<ScalarType,MV> state;
 
      // Copy over the vectors.
      state.W = W_;
      state.U = U_;
      state.AU = AU_;
      state.Rtilde = Rtilde_;
      state.D = D_;
      state.V = V_;

      // Copy over the scalars.
      state.alpha = alpha_;
      state.eta = eta_;
      state.rho = rho_;
      state.tau = tau_;
      state.theta = theta_;

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
    Teuchos::RCP<const MV> getNativeResiduals( std::vector<MagnitudeType> *norms ) const;
    
    //! Get the current update to the linear system.
    /*! \note This method returns the accumulated update to the solution instead of updating
              the linear problem, since it may incur an additional preconditioner application each iteration.
     */
    Teuchos::RCP<MV> getCurrentUpdate() const { return solnUpdate_; }

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
                       "Belos::PseudoBlockTFQMRIter::setBlockSize(): Cannot use a block size that is not one.");
    }
    
    //! States whether the solver has been initialized or not.
    bool isInitialized() { return initialized_; }
    
    //@}
    
  
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

    // Storage for QR factorization of the least squares system.
    std::vector<ScalarType> alpha_, eta_, rho_, rho_old_;
    std::vector<MagnitudeType> tau_, theta_;
    
    // 
    // Current solver state
    //
    // initialized_ specifies that the basis vectors have been initialized and the iterate() routine
    // is capable of running; _initialize is controlled  by the initialize() member method
    // For the implications of the state of initialized_, please see documentation for initialize()
    bool initialized_;
    
     // Current subspace dimension, and number of iterations performed.
    int iter_;
    
    // 
    // State Storage
    //
    Teuchos::RCP<MV> W_;
    Teuchos::RCP<MV> U_, AU_;
    Teuchos::RCP<MV> Rtilde_;
    Teuchos::RCP<MV> D_;
    Teuchos::RCP<MV> V_;
    Teuchos::RCP<MV> solnUpdate_;

  };
  
  
  //
  // Implementation
  //
  
  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor.
  template <class ScalarType, class MV, class OP>
  PseudoBlockTFQMRIter<ScalarType,MV,OP>::PseudoBlockTFQMRIter(const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem,
					 const Teuchos::RCP<OutputManager<ScalarType> > &printer,
					 const Teuchos::RCP<StatusTest<ScalarType,MV,OP> > &tester,
					 Teuchos::ParameterList &/* params */ 
					 ) : 
    lp_(problem), 
    om_(printer),
    stest_(tester),
    numRHS_(0),
    initialized_(false),
    iter_(0)
  { 
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Compute native residual from TFQMR recurrence.
  template <class ScalarType, class MV, class OP>
  Teuchos::RCP<const MV> 
  PseudoBlockTFQMRIter<ScalarType,MV,OP>::getNativeResiduals( std::vector<MagnitudeType> *normvec ) const 
  {
    MagnitudeType one = Teuchos::ScalarTraits<MagnitudeType>::one();
    if (normvec) {
      // Resize the vector passed in, if it is too small.
      if ((int) normvec->size() < numRHS_)
        normvec->resize( numRHS_ );

      // Compute the native residuals.
      for (int i=0; i<numRHS_; i++) {
        (*normvec)[i] = Teuchos::ScalarTraits<MagnitudeType>::squareroot( 2*iter_ + one )*tau_[i];
      }
    }

    return Teuchos::null;
  }
	
  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Initialize this iteration object
  template <class ScalarType, class MV, class OP>
  void PseudoBlockTFQMRIter<ScalarType,MV,OP>::initializeTFQMR(const PseudoBlockTFQMRIterState<ScalarType,MV> & newstate)
  {
    // Create convenience variables for zero and one.
    const ScalarType STone = Teuchos::ScalarTraits<ScalarType>::one();
    const ScalarType STzero = Teuchos::ScalarTraits<ScalarType>::zero();
    const MagnitudeType MTzero = Teuchos::ScalarTraits<MagnitudeType>::zero();

    // NOTE:  In PseudoBlockTFQMRIter Rtilde_, the initial residual, is required!!!
    TEUCHOS_TEST_FOR_EXCEPTION(newstate.Rtilde == Teuchos::null,std::invalid_argument,
                       "Belos::PseudoBlockTFQMRIter::initialize(): PseudoBlockTFQMRIterState does not have initial residual.");

    // Get the number of right-hand sides we're solving for now.
    int numRHS = MVT::GetNumberVecs(*newstate.Rtilde);
    numRHS_ = numRHS;

    // Initialize the state storage
    // If the subspace has not be initialized before or we are reusing this solver object, generate it using Rtilde.
    if ( Teuchos::is_null(Rtilde_) || (MVT::GetNumberVecs(*Rtilde_) == numRHS_) )
    {
      // Create and/or initialize D_.
      if ( Teuchos::is_null(D_) )
        D_ = MVT::Clone( *newstate.Rtilde, numRHS_ );
      MVT::MvInit( *D_, STzero );

      // Create and/or initialize solnUpdate_;
      if ( Teuchos::is_null(solnUpdate_) )
        solnUpdate_ = MVT::Clone( *newstate.Rtilde, numRHS_ );
      MVT::MvInit( *solnUpdate_, STzero );
     
      // Create Rtilde_. 
      if (newstate.Rtilde != Rtilde_) 
        Rtilde_ = MVT::CloneCopy( *newstate.Rtilde );
      W_ = MVT::CloneCopy( *Rtilde_ );
      U_ = MVT::CloneCopy( *Rtilde_ );
      V_ = MVT::Clone( *Rtilde_, numRHS_ );

      // Multiply the current residual by Op and store in V_
      //       V_ = Op * R_ 
      lp_->apply( *U_, *V_ );
      AU_ = MVT::CloneCopy( *V_ ); 

      // Resize work vectors.
      alpha_.resize( numRHS_, STone );
      eta_.resize( numRHS_, STzero );
      rho_.resize( numRHS_ );
      rho_old_.resize( numRHS_ );
      tau_.resize( numRHS_ );
      theta_.resize( numRHS_, MTzero );

      MVT::MvNorm( *Rtilde_, tau_ );                     // tau = ||r_0||
      MVT::MvDot( *Rtilde_, *Rtilde_, rho_ );            // rho = (r_tilde, r0)
    }   
    else 
    {
      // If the subspace has changed sizes, clone it from the incoming state.
      Rtilde_ = MVT::CloneCopy( *newstate.Rtilde );
      W_ = MVT::CloneCopy( *newstate.W );
      U_ = MVT::CloneCopy( *newstate.U );
      AU_ = MVT::CloneCopy( *newstate.AU );
      D_ = MVT::CloneCopy( *newstate.D );
      V_ = MVT::CloneCopy( *newstate.V );

      // The solution update is just set to zero, since the current update has already
      // been added to the solution during deflation.
      solnUpdate_ = MVT::Clone( *Rtilde_, numRHS_ );
      MVT::MvInit( *solnUpdate_, STzero );

      // Copy work vectors.
      alpha_ = newstate.alpha;
      eta_ = newstate.eta;
      rho_ = newstate.rho;
      tau_ = newstate.tau;
      theta_ = newstate.theta;
    }

    // The solver is initialized
    initialized_ = true;
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Iterate until the status test informs us we should stop.
  template <class ScalarType, class MV, class OP>
  void PseudoBlockTFQMRIter<ScalarType,MV,OP>::iterate() 
  {
    //
    // Allocate/initialize data structures
    //
    if (initialized_ == false) {
      initialize();
    }
 
    // Create convenience variables for zero and one.
    const ScalarType STone = Teuchos::ScalarTraits<ScalarType>::one();
    const ScalarType STzero = Teuchos::ScalarTraits<ScalarType>::zero();
    const MagnitudeType MTone = Teuchos::ScalarTraits<MagnitudeType>::one();
    const MagnitudeType MTzero = Teuchos::ScalarTraits<MagnitudeType>::zero();
    std::vector< ScalarType > beta(numRHS_,STzero);
    std::vector<int> index(1);
    //
    //  Start executable statements. 
    //

    ////////////////////////////////////////////////////////////////
    // Iterate until the status test tells us to stop.
    //
    while (stest_->checkStatus(this) != Passed) {

      for (int iIter=0; iIter<2; iIter++)
      {
        //
        //--------------------------------------------------------
        // Compute the new alpha if we need to
        //--------------------------------------------------------
        //
        if (iIter == 0) {
  	  MVT::MvDot( *V_, *Rtilde_, alpha_ );      //   alpha = rho / (r_tilde, v) 
          for (int i=0; i<numRHS_; i++) {
	    rho_old_[i] = rho_[i];                   // rho_old = rho
	    alpha_[i] = rho_old_[i]/alpha_[i];
          }
        }
        //
        //--------------------------------------------------------
        // Loop over all RHS and compute updates.
        //--------------------------------------------------------
        //
        for (int i=0; i<numRHS_; ++i) {
          index[0] = i;

          //
          //--------------------------------------------------------
          // Update w.
          //   w = w - alpha*Au
          //--------------------------------------------------------
          //
          Teuchos::RCP<const MV> AU_i = MVT::CloneView( *AU_, index );
          Teuchos::RCP<MV> W_i = MVT::CloneViewNonConst( *W_, index );
          MVT::MvAddMv( STone, *W_i, -alpha_[i], *AU_i, *W_i );
          //
          //--------------------------------------------------------
          // Update d.
          //   d = u + (theta^2/alpha)eta*d
          //--------------------------------------------------------
          //
          Teuchos::RCP<const MV> U_i = MVT::CloneView( *U_, index );
          Teuchos::RCP<MV> D_i = MVT::CloneViewNonConst( *D_, index );
          MVT::MvAddMv( STone, *U_i, (theta_[i]*theta_[i]/alpha_[i])*eta_[i], *D_i, *D_i );
          //
          //--------------------------------------------------------
          // Update u if we need to.
          //   u = u - alpha*v
          //   
          // Note: This is usually computed with alpha (above), but we're trying be memory efficient.
          //--------------------------------------------------------
          //
          if (iIter == 0) {
            // Compute new U.
            Teuchos::RCP<const MV> V_i = MVT::CloneView( *V_, index );
            Teuchos::RCP<MV> U2_i = MVT::CloneViewNonConst( *U_, index );
  	    MVT::MvAddMv( STone, *U2_i, -alpha_[i], *V_i, *U2_i );
          }
        }
        //
        //--------------------------------------------------------
        // Update Au for the next iteration.
        //--------------------------------------------------------
        //
        if (iIter == 0) {
  	  lp_->apply( *U_, *AU_ );                       
        }
        //
        //--------------------------------------------------------
        // Compute the new theta, c, eta, tau; i.e. the update to the least squares solution.
        //--------------------------------------------------------
        //
        MVT::MvNorm( *W_, theta_ );     // theta = ||w|| / tau

        bool breakdown=false;
        for (int i=0; i<numRHS_; ++i) {
          theta_[i] /= tau_[i];
          // cs = 1.0 / sqrt(1.0 + theta^2)
          MagnitudeType cs = MTone / Teuchos::ScalarTraits<MagnitudeType>::squareroot(MTone + theta_[i]*theta_[i]);  
          tau_[i] *= theta_[i]*cs;     // tau = tau * theta * cs
          if ( tau_[i] == MTzero ) {
            breakdown = true;
          }
          eta_[i] = cs*cs*alpha_[i];     // eta = cs^2 * alpha
        }
        //
        //--------------------------------------------------------
        // Accumulate the update for the solution x := x + eta*D_
        //--------------------------------------------------------
        //
        for (int i=0; i<numRHS_; ++i) {
          index[0]=i;
          Teuchos::RCP<const MV> D_i = MVT::CloneView( *D_, index );
          Teuchos::RCP<MV> update_i = MVT::CloneViewNonConst( *solnUpdate_, index );
	  MVT::MvAddMv( STone, *update_i, eta_[i], *D_i, *update_i );
        } 
        //
        //--------------------------------------------------------
        // Breakdown was detected above, return to status test to
        // remove converged solutions.
        //--------------------------------------------------------
        if ( breakdown ) {
          break;
        }
        //
        if (iIter == 1) {
  	  //
	  //--------------------------------------------------------
	  // Compute the new rho, beta if we need to.
	  //--------------------------------------------------------
	  //
	  MVT::MvDot( *W_, *Rtilde_, rho_ );         // rho = (r_tilde, w)
      
          for (int i=0; i<numRHS_; ++i) {
	    beta[i] = rho_[i]/rho_old_[i];           // beta = rho / rho_old

  	    //
	    //--------------------------------------------------------
	    // Update u, v, and Au if we need to.
	    // Note: We are updating v in two stages to be memory efficient
	    //--------------------------------------------------------
	    //
            index[0]=i;
            Teuchos::RCP<const MV> W_i = MVT::CloneView( *W_, index );
            Teuchos::RCP<MV> U_i = MVT::CloneViewNonConst( *U_, index );
	    MVT::MvAddMv( STone, *W_i, beta[i], *U_i, *U_i );       // u = w + beta*u
	
	    // First stage of v update.
            Teuchos::RCP<const MV> AU_i = MVT::CloneView( *AU_, index );
            Teuchos::RCP<MV> V_i = MVT::CloneViewNonConst( *V_, index );
	    MVT::MvAddMv( STone, *AU_i, beta[i], *V_i, *V_i );      // v = Au + beta*v 
	  }

	  // Update Au.
	  lp_->apply( *U_, *AU_ );                          // Au = A*u
	
	  // Second stage of v update.
          for (int i=0; i<numRHS_; ++i) {
            index[0]=i;
            Teuchos::RCP<const MV> AU_i = MVT::CloneView( *AU_, index );
            Teuchos::RCP<MV> V_i = MVT::CloneViewNonConst( *V_, index );
	    MVT::MvAddMv( STone, *AU_i, beta[i], *V_i, *V_i );      // v = Au + beta*v 
          }
        }
      }

      // Increment the iteration
      iter_++;
      
    } // end while (sTest_->checkStatus(this) != Passed)
  } 

} // namespace Belos
//
#endif // BELOS_PSEUDO_BLOCK_TFQMR_ITER_HPP
//
// End of file BelosPseudoBlockTFQMRIter.hpp


