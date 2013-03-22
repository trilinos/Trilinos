//@HEADER
// ************************************************************************
//
//                 Belos: Block Linear Solvers Package
//                  Copyright 2004 Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER
//
// This file contains an implementation of the TFQMR iteration
// for solving non-Hermitian linear systems of equations Ax = b, 
// where b is a single-vector and x is the corresponding solution.
//
// The implementation is a slight modification on the TFQMR iteration
// found in Saad's "Iterative Methods for Sparse Linear Systems".
//

#ifndef BELOS_TFQMR_ITER_HPP
#define BELOS_TFQMR_ITER_HPP

/*!
  \file BelosTFQMRIter.hpp

  \brief Belos concrete class for generating iterations with the
  preconditioned tranpose-free QMR (TFQMR) method.
*/

#include "BelosConfigDefs.hpp"
#include "BelosIteration.hpp"
#include "BelosTypes.hpp"	

#include "BelosLinearProblem.hpp"
#include "BelosMatOrthoManager.hpp"
#include "BelosOutputManager.hpp"
#include "BelosStatusTest.hpp"
#include "BelosOperatorTraits.hpp"
#include "BelosMultiVecTraits.hpp"

#include "Teuchos_BLAS.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TimeMonitor.hpp"

/*!	\class Belos::TFQMRIter

	\brief This class implements the preconditioned transpose-free QMR algorithm for
	solving non-Hermitian linear systems of equations Ax = b, where b is the right-hand 
	side vector and x is the corresponding solution.

        \ingroup belos_solver_framework

	\author Heidi Thornquist
*/

namespace Belos {

  /** \brief Structure to contain pointers to TFQMRIter state variables.
   *
   * This struct is utilized by TFQMRIter::initialize() and TRQMRIter::getState().
   */
  template <class ScalarType, class MV>
  struct TFQMRIterState {

    /*! \brief The current residual basis. */
    Teuchos::RCP<const MV> R;
    Teuchos::RCP<const MV> W;
    Teuchos::RCP<const MV> U;
    Teuchos::RCP<const MV> Rtilde;
    Teuchos::RCP<const MV> D;
    Teuchos::RCP<const MV> V;

    TFQMRIterState() : R(Teuchos::null), W(Teuchos::null), U(Teuchos::null),
                       Rtilde(Teuchos::null), D(Teuchos::null), V(Teuchos::null)
    {}
  };

  
  //! @name TFQMRIter Exceptions
  //@{
  
  /** \brief TFQMRIterInitFailure is thrown when the TFQMRIter object is unable to
   * generate an initial iterate in the TFQMRIter::initialize() routine.
   *
   * This std::exception is thrown from the TFQMRIter::initialize() method, which is
   * called by the user or from the TFQMRIter::iterate() method if isInitialized()
   * == \c false.
   *
   * In the case that this std::exception is thrown,
   * TFQMRIter::isInitialized() will be \c false and the user will need to provide
   * a new initial iterate to the iteration.
   */
  class TFQMRIterInitFailure : public BelosError {public:
    TFQMRIterInitFailure(const std::string& what_arg) : BelosError(what_arg)
    {}};
  
  /** \brief TFQMRIterateFailure is thrown when the TFQMRIter object is unable to
   * compute the next iterate in the TFQMRIter::iterate() routine.
   *
   * This std::exception is thrown from the TFQMRIter::iterate() method.
   *
   */
  class TFQMRIterateFailure : public BelosError {public:
    TFQMRIterateFailure(const std::string& what_arg) : BelosError(what_arg)
    {}};
  
  //@}


  template <class ScalarType, class MV, class OP>
  class TFQMRIter : public Iteration<ScalarType,MV,OP> { 
  public:
    //
    // Convenience typedefs
    //
    typedef MultiVecTraits<ScalarType,MV> MVT;
    typedef MultiVecTraitsExt<ScalarType,MV> MVText;
    typedef OperatorTraits<ScalarType,MV,OP> OPT;
    typedef Teuchos::ScalarTraits<ScalarType> SCT;
    typedef typename SCT::magnitudeType MagnitudeType;
    
    //! @name Constructor/Destructor.
    //@{ 

    //! %Belos::TFQMRIter constructor.
    TFQMRIter( const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem, 
	       const Teuchos::RCP<OutputManager<ScalarType> > &printer,
	       const Teuchos::RCP<StatusTest<ScalarType,MV,OP> > &tester,
	       Teuchos::ParameterList &params );
    
    //! %Belos::TFQMRIter destructor.
    virtual ~TFQMRIter() {};
    //@}
  

    //! @name Solver methods
    //@{ 
    
    /*! \brief This method performs block TFQMR iterations until the status
     * test indicates the need to stop or an error occurs (in which case, an
     * std::exception is thrown).
     *
     * iterate() will first determine whether the solver is inintialized; if
     * not, it will call initialize() using default arguments. After
     * initialization, the solver performs block TFQMR iterations until the
     * status test evaluates as ::Passed, at which point the method returns to
     * the caller. 
     *
     * The block TFQMR iteration proceeds as follows:
     * -# The operator problem->applyOp() is applied to the newest \c blockSize vectors in the Krylov basis.
     * -# The resulting vectors are orthogonalized against the previous basis vectors, and made orthonormal.
     * -# The Hessenberg matrix is updated.
     * -# The least squares system is updated.
     *
     * The status test is queried at the beginning of the iteration.
     *
     * Possible exceptions thrown include the TFQMRIterOrthoFailure.
     *
     */
    void iterate();

    /*! \brief Initialize the solver to an iterate, providing a complete state.
     *
     * The %BlockTFQMRIter contains a certain amount of state, consisting of the current 
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
    void initializeTFQMR(TFQMRIterState<ScalarType,MV> newstate);
    
    /*! \brief Initialize the solver with the initial vectors from the linear problem
     *  or random data.
     */
    void initialize()
    {
      TFQMRIterState<ScalarType,MV> empty;
      initializeTFQMR(empty);
    }
    
    /*! \brief Get the current state of the linear solver.
     *
     * The data is only valid if isInitialized() == \c true.
     *
     * \returns A TFQMRIterState object containing const pointers to the current
     * solver state.
     */
    TFQMRIterState<ScalarType,MV> getState() const {
      TFQMRIterState<ScalarType,MV> state;
      state.R = R_;
      state.W = W_;
      state.U = U_;
      state.Rtilde = Rtilde_;
      state.D = D_;
      state.V = V_;
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
    /*! \note This method returns a null pointer because the linear problem is updated every iteration.
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
                       "Belos::TFQMRIter::setBlockSize(): Cannot use a block size that is not one.");
    }
    
    //! States whether the solver has been initialized or not.
    bool isInitialized() { return initialized_; }
    
    //@}
    
  
  private:

    //
    // Internal methods
    //
    //! Method for initalizing the state storage needed by TFQMR.
    void setStateSize();
    
    //
    // Classes inputed through constructor that define the linear problem to be solved.
    //
    const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> >    lp_;
    const Teuchos::RCP<OutputManager<ScalarType> >          om_;
    const Teuchos::RCP<StatusTest<ScalarType,MV,OP> >       stest_;
    
    //
    // Algorithmic parameters
    //      

    // Storage for QR factorization of the least squares system.
    Teuchos::SerialDenseMatrix<int,ScalarType> alpha_, rho_, rho_old_;
    std::vector<MagnitudeType> tau_, cs_, theta_;
    
    // 
    // Current solver state
    //
    // initialized_ specifies that the basis vectors have been initialized and the iterate() routine
    // is capable of running; _initialize is controlled  by the initialize() member method
    // For the implications of the state of initialized_, please see documentation for initialize()
    bool initialized_;
    
    // stateStorageInitialized_ specifies that the state storage has be initialized to the current
    // blockSize_ and numBlocks_.  This initialization may be postponed if the linear problem was
    // generated without the right-hand side or solution vectors.
    bool stateStorageInitialized_;
    
     // Current subspace dimension, and number of iterations performed.
    int iter_;
    
    // 
    // State Storage
    //
    Teuchos::RCP<MV> R_;
    Teuchos::RCP<MV> W_;
    Teuchos::RCP<MV> U_, AU_;
    Teuchos::RCP<MV> Rtilde_;
    Teuchos::RCP<MV> D_;
    Teuchos::RCP<MV> V_;

  };
  
  
  //
  // Implementation
  //
  
  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor.
  template <class ScalarType, class MV, class OP>
  TFQMRIter<ScalarType,MV,OP>::TFQMRIter(const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem,
					 const Teuchos::RCP<OutputManager<ScalarType> > &printer,
					 const Teuchos::RCP<StatusTest<ScalarType,MV,OP> > &tester,
					 Teuchos::ParameterList &params 
					 ) : 
    lp_(problem), 
    om_(printer),
    stest_(tester),
    alpha_(1,1),
    rho_(1,1),
    rho_old_(1,1),
    tau_(1),
    cs_(1),
    theta_(1),
    initialized_(false),
    stateStorageInitialized_(false),
    iter_(0)
  { 
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Compute native residual from TFQMR recurrence.
  template <class ScalarType, class MV, class OP>
  Teuchos::RCP<const MV> 
  TFQMRIter<ScalarType,MV,OP>::getNativeResiduals( std::vector<MagnitudeType> *normvec ) const 
  {
    MagnitudeType one = Teuchos::ScalarTraits<MagnitudeType>::one();
    if (normvec)
      (*normvec)[0] = Teuchos::ScalarTraits<MagnitudeType>::squareroot( iter_ + one )*tau_[0];

    return Teuchos::null;
  }
	

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Setup the state storage.
  template <class ScalarType, class MV, class OP>
  void TFQMRIter<ScalarType,MV,OP>::setStateSize ()
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
                             "Belos::TFQMRIter::setStateSize(): linear problem does not specify multivectors to clone from.");
          R_ = MVT::Clone( *tmp, 1 );
	  AU_ = MVT::Clone( *tmp, 1 );
	  D_ = MVT::Clone( *tmp, 1 );
          V_ = MVT::Clone( *tmp, 1 );
        }

        // State storage has now been initialized.
        stateStorageInitialized_ = true;
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Initialize this iteration object
  template <class ScalarType, class MV, class OP>
  void TFQMRIter<ScalarType,MV,OP>::initializeTFQMR(TFQMRIterState<ScalarType,MV> newstate)
  {
    // Initialize the state storage if it isn't already.
    if (!stateStorageInitialized_)
      setStateSize();

    TEUCHOS_TEST_FOR_EXCEPTION(!stateStorageInitialized_,std::invalid_argument,
                       "Belos::TFQMRIter::initialize(): Cannot initialize state storage!");

    // NOTE:  In TFQMRIter R_, the initial residual, is required!!!
    //
    std::string errstr("Belos::TFQMRIter::initialize(): Specified multivectors must have a consistent length and width.");

    // Create convenience variables for zero and one.
    const ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();
    const ScalarType STzero = Teuchos::ScalarTraits<ScalarType>::zero();
    const MagnitudeType MTzero = Teuchos::ScalarTraits<MagnitudeType>::zero();

    if (newstate.R != Teuchos::null) {

      TEUCHOS_TEST_FOR_EXCEPTION( MVText::GetGlobalLength(*newstate.R) != MVText::GetGlobalLength(*R_),
                          std::invalid_argument, errstr );
      TEUCHOS_TEST_FOR_EXCEPTION( MVT::GetNumberVecs(*newstate.R) != 1,
                          std::invalid_argument, errstr );

      // Copy basis vectors from newstate into V
      if (newstate.R != R_) {
        // copy over the initial residual (unpreconditioned).
        MVT::MvAddMv( one, *newstate.R, STzero, *newstate.R, *R_ );
      }

      // Compute initial vectors
      // Initially, they are set to the preconditioned residuals
      //
      W_ = MVT::CloneCopy( *R_ );
      U_ = MVT::CloneCopy( *R_ );
      Rtilde_ = MVT::CloneCopy( *R_ );
      MVT::MvInit( *D_ );
      // Multiply the current residual by Op and store in V_
      //       V_ = Op * R_ 
      //
      lp_->apply( *U_, *V_ );
      AU_ = MVT::CloneCopy( *V_ ); 
      //
      // Compute initial scalars: theta, eta, tau, rho_old
      //
      theta_[0] = MTzero;
      MVT::MvNorm( *R_, tau_ );                         // tau = ||r_0||
      MVT::MvTransMv( one, *Rtilde_, *R_, rho_old_ );   // rho = (r_tilde, r0)
    }
    else {

      TEUCHOS_TEST_FOR_EXCEPTION(newstate.R == Teuchos::null,std::invalid_argument,
                         "Belos::TFQMRIter::initialize(): TFQMRIterState does not have initial residual.");
    }

    // The solver is initialized
    initialized_ = true;
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Iterate until the status test informs us we should stop.
  template <class ScalarType, class MV, class OP>
  void TFQMRIter<ScalarType,MV,OP>::iterate() 
  {
    //
    // Allocate/initialize data structures
    //
    if (initialized_ == false) {
      initialize();
    }
 
    // Create convenience variables for zero and one.
    const ScalarType STone = Teuchos::ScalarTraits<ScalarType>::one();
    const MagnitudeType MTone = Teuchos::ScalarTraits<MagnitudeType>::one();
    const ScalarType STzero = Teuchos::ScalarTraits<ScalarType>::zero();
    ScalarType eta = STzero, beta = STzero;
    //
    //  Start executable statements. 
    //
    // Get the current solution vector.
    Teuchos::RCP<MV> cur_soln_vec = lp_->getCurrLHSVec();

    // Check that the current solution vector only has one column.
    TEUCHOS_TEST_FOR_EXCEPTION( MVT::GetNumberVecs(*cur_soln_vec) != 1, TFQMRIterateFailure,
                        "Belos::TFQMRIter::iterate(): current linear system has more than one vector!" );


    ////////////////////////////////////////////////////////////////
    // Iterate until the status test tells us to stop.
    //
    while (stest_->checkStatus(this) != Passed) {

      //
      //--------------------------------------------------------
      // Compute the new alpha if we need to
      //--------------------------------------------------------
      //
      if (iter_%2 == 0) {
	MVT::MvTransMv( STone, *Rtilde_, *V_, alpha_ );      //   alpha = rho / (r_tilde, v) 
	alpha_(0,0) = rho_old_(0,0)/alpha_(0,0);
      }
      //
      //--------------------------------------------------------
      // Update w.
      //   w = w - alpha*Au
      //--------------------------------------------------------
      //
      MVT::MvAddMv( STone, *W_, -alpha_(0,0), *AU_, *W_ );
      //
      //--------------------------------------------------------
      // Update d.
      //   d = u + (theta^2/alpha)eta*d
      //--------------------------------------------------------
      //
      MVT::MvAddMv( STone, *U_, (theta_[0]*theta_[0]/alpha_(0,0))*eta, *D_, *D_ );
      //
      //--------------------------------------------------------
      // Update u if we need to.
      //   u = u - alpha*v
      //   
      // Note: This is usually computed with alpha (above), but we're trying be memory efficient.
      //--------------------------------------------------------
      //
      if (iter_%2 == 0) {
        // Compute new U.
	MVT::MvAddMv( STone, *U_, -alpha_(0,0), *V_, *U_ );

	// Update Au for the next iteration.
	lp_->apply( *U_, *AU_ );                       
      }
      //
      //--------------------------------------------------------
      // Compute the new theta, c, eta, tau; i.e. the update to the least squares solution.
      //--------------------------------------------------------
      //
      MVT::MvNorm( *W_, theta_ );     // theta = ||w|| / tau
      theta_[0] /= tau_[0];
      // cs = 1.0 / sqrt(1.0 + theta^2)
      cs_[0] = MTone / Teuchos::ScalarTraits<MagnitudeType>::squareroot(MTone + theta_[0]*theta_[0]);  
      tau_[0] *= theta_[0]*cs_[0];     // tau = tau * theta * cs
      eta = cs_[0]*cs_[0]*alpha_(0,0);     // eta = cs^2 * alpha
      
      //
      //--------------------------------------------------------
      // Update the solution.
      //--------------------------------------------------------
      //
      lp_->updateSolution( D_, true, eta );
      //
      if (iter_%2) {
	//
	//--------------------------------------------------------
	// Compute the new rho, beta if we need to.
	//--------------------------------------------------------
	//
	MVT::MvTransMv( STone, *Rtilde_, *W_, rho_ );     // rho = (r_tilde, w)
	beta = rho_(0,0)/rho_old_(0,0);                   // beta = rho / rho_old
	rho_old_(0,0) = rho_(0,0);                        // rho_old = rho
	//
	//--------------------------------------------------------
	// Update u, v, and Au if we need to.
	// Note: We are updating v in two stages to be memory efficient
	//--------------------------------------------------------
	//
	MVT::MvAddMv( STone, *W_, beta, *U_, *U_ );       // u = w + beta*u
	
	// First stage of v update.
	MVT::MvAddMv( STone, *AU_, beta, *V_, *V_ );      // v = Au + beta*v 
	
	// Update Au.
	lp_->apply( *U_, *AU_ );                          // Au = A*u
	
	// Second stage of v update.
	MVT::MvAddMv( STone, *AU_, beta, *V_, *V_ );      // v = Au + beta*v
      }

      // Increment the iteration
      iter_++;
      
    } // end while (sTest_->checkStatus(this) != Passed)
  } 

} // namespace Belos
//
#endif // BELOS_TFQMR_ITER_HPP
//
// End of file BelosTFQMRIter.hpp


