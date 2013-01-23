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

#ifndef BELOS_MINRES_ITER_HPP
#define BELOS_MINRES_ITER_HPP

/// \file BelosMinresIter.hpp
/// \brief MINRES iteration implementation 
///
/// The Minimal Residual Method (MINRES) is a Krylov subspace method
/// for solving symmetric (in real arithmetic, or Hermitian in complex
/// arithmetic), nonsingular, but possibly indefinite linear systems
/// \fn$Ax=b\fn$.  It works on one right-hand side \fn$b\fn$ at a
/// time.
///
/// References:
///
/// C. Paige and M. Saunders.  "Solution of sparse indefinite systems
/// of linear equations."  SIAM J. Numer. Anal., vol. 12, pp. 617-629,
/// 1975.
///
/// http://www.stanford.edu/group/SOL/software/minres/matlab/minres.m

#include "BelosConfigDefs.hpp"
#include "BelosTypes.hpp"
#include "BelosMinresIteration.hpp"

#include "BelosLinearProblem.hpp"
#include "BelosOutputManager.hpp"
#include "BelosStatusTest.hpp"
#include "BelosOperatorTraits.hpp"
#include "BelosMultiVecTraits.hpp"

#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_BLAS.hpp"

namespace Belos {

///
/// \class MinresIter
/// \brief MINRES implementation
/// \author Nico Schl\"omer
///
/// Implementation of the preconditioned Minimal Residual Method
/// (MINRES) iteration.  This a bilinear form implementation, that
/// uses inner products of the form <x,My> to solve the preconditioned
/// linear system M^{-1}*A x = b.  Thus, it is necessary that the
/// left preconditioner M is positive definite.
///
/// \ingroup belos_solver_framework
///
template<class ScalarType, class MV, class OP>
class MinresIter : virtual public MinresIteration<ScalarType,MV,OP> {

  public:

  //
  // Convenience typedefs
  //
  typedef MultiVecTraits< ScalarType, MV > MVT;
  typedef MultiVecTraitsExt< ScalarType, MV > MVText;
  typedef OperatorTraits< ScalarType, MV, OP > OPT;
  typedef Teuchos::ScalarTraits< ScalarType > SCT;
  typedef typename SCT::magnitudeType MagnitudeType;
  typedef Teuchos::ScalarTraits< MagnitudeType > SMT;

  //! @name Constructors/Destructor
  //@{ 

  /// \brief Constructor
  ///
  /// \params problem The linear problem to solve
  /// \params printer Output manager, for intermediate solver output
  /// \params tester Status test for determining when the current
  ///   approximate solution has converged 
  /// \params params Parameter list of solver options 
  ///
  MinresIter (const Teuchos::RCP< LinearProblem< ScalarType, MV, OP > >& problem,
	      const Teuchos::RCP< OutputManager< ScalarType > > &        printer,
	      const Teuchos::RCP< StatusTest< ScalarType, MV, OP > >&    tester,
	      const Teuchos::ParameterList& params);

  //! Destructor
  virtual ~MinresIter() {};
  //@}


  //! @name Solver methods
  //@{ 

  /// \brief Perform MINRES iterations until convergence or error
  ///
  /// Perform MINRES iterations until the status test indicates the
  /// need to stop, or until an error occurs.  In the latter case, a
  /// (subclass of) std::exception is thrown.
  ///
  /// iterate() will first determine whether the solver is
  /// initialized; if not, it will call initialize() using default
  /// arguments. After initialization, the solver performs MINRES
  /// iterations until the status test evaluates as ::Passed, at which
  /// point the method returns to the caller.
  ///
  /// The status test is queried at the beginning of the iteration.
  ///
  void iterate();

  /*! \brief Initialize the solver to an iterate, providing a complete state.
   *
   * The %MinresIter contains a certain amount of state, consisting of the current
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
  void initializeMinres (MinresIterationState<ScalarType,MV> newstate);

  /// \brief Initialize the solver
  ///
  /// Initialize the solver.  If a starting guess is provided in the
  /// linear problem, use that.  Otherwise, choose a random starting
  /// guess.
  void initialize()
  {
    MinresIterationState<ScalarType,MV> empty;
    initializeMinres(empty);
  }

  /// Get the current state of the linear solver.
  ///
  /// The returned current state is only valid if isInitialized() == \c true.
  ///
  /// \return A MinresIterationState object containing const pointers
  ///         to the current solver state.
  MinresIterationState<ScalarType,MV> getState() const {
    if (! isInitialized())
      throw std::logic_error("getState() cannot be called unless "
			     "the state has been initialized");
    MinresIterationState<ScalarType,MV> state;
    state.Y = Y_;
    state.R1 = R1_;
    state.R2 = R2_;
    state.W = W_;
    state.W1 = W1_;
    state.W2 = W2_;
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
  Teuchos::RCP<const MV> 
  getNativeResiduals( std::vector<MagnitudeType> *norms ) const 
  { 
    if (norms != NULL)
      {
	std::vector<MagnitudeType>& theNorms = *norms;
	if (theNorms.size() < 1)
	  theNorms.resize(1);
	theNorms[0] = phibar_;
      }
    return Teuchos::null; 
  }

  //! Get the current update to the linear system.
  /*! \note This method returns a null pointer because the linear problem is current.
  */
  Teuchos::RCP<MV> getCurrentUpdate() const { return Teuchos::null; }

  //!
  void symOrtho( ScalarType a, ScalarType b, ScalarType *c, ScalarType *s, ScalarType *r );

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
                       "Belos::MinresIter::setBlockSize(): Cannot use a block size that is not one.");
  }

  //! States whether the solver has been initialized or not.
  bool isInitialized() const { return initialized_; }
  bool isInitialized() { return initialized_; }

  //@}

  private:

  //
  // Internal methods
  //
  //! Method for initalizing the state storage needed by MINRES.
  void setStateSize();

  //
  // Classes inputed through constructor that define the linear problem to be solved.
  //
  const Teuchos::RCP< LinearProblem< ScalarType, MV, OP > > lp_;
  const Teuchos::RCP< OutputManager< ScalarType > >         om_;
  const Teuchos::RCP< StatusTest< ScalarType, MV, OP > >    stest_;


  /// \brief Whether the solver has been initialized
  ///
  /// If initialized_ == true, then the basis vectors have been
  /// initialized and the iterate() routine is capable of running.
  /// _initialize is set by the initialize() member method.  For the
  /// implications of the state of initialized_, please see
  /// documentation for initialize().
  bool initialized_;

  /// \brief Whether the state storage has been initialized
  ///
  /// If stateStorageInitialized_ == true, then the state storage has
  /// been initialized.  This initialization may be postponed if the
  /// linear problem was generated without the right-hand side or
  /// solution vectors.
  bool stateStorageInitialized_;

  //! Current number of iterations performed.
  int iter_;

  /// \brief Current "native" residual
  ///
  /// Current "native" residual (not the "exact" residual \fn$\|b -
  /// Ax\|_2\fn$).
  MagnitudeType phibar_;

  // 
  // State Storage
  //

  //! Preconditioned residual
  Teuchos::RCP< MV > Y_;
  //! Previous residual
  Teuchos::RCP< MV > R1_;
  //! Previous residual
  Teuchos::RCP< MV > R2_;
  //! Direction vector
  Teuchos::RCP< MV > W_;
  //! Previous direction vector
  Teuchos::RCP< MV > W1_;
  //! Previous direction vector
  Teuchos::RCP< MV > W2_;

  /// Coefficient in the MINRES iteration
  ///
  /// \note If we could be sure that the preconditioner is Hermitian
  ///   in complex arithmetic (which must be true anyway, in order for
  ///   MINRES to work), we could make beta1_ a MagnitudeType.  This
  ///   would certainly be cleaner, considering it will be copied into
  ///   beta (which is of MagnitudeType).
  Teuchos::SerialDenseMatrix<int,ScalarType> beta1_;

};

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor.
  template<class ScalarType, class MV, class OP>
  MinresIter<ScalarType,MV,OP>::MinresIter(const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem,
                                                   const Teuchos::RCP<OutputManager<ScalarType> > &printer,
                                                   const Teuchos::RCP<StatusTest<ScalarType,MV,OP> > &tester,
                                                   const Teuchos::ParameterList &params ):
    lp_(problem),
    om_(printer),
    stest_(tester),
    initialized_(false),
    stateStorageInitialized_(false),
    iter_(0)
  {
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Setup the state storage.
  template <class ScalarType, class MV, class OP>
  void MinresIter<ScalarType,MV,OP>::setStateSize ()
  {
    if (!stateStorageInitialized_) {

      // Check if there is any multivector to clone from.
      Teuchos::RCP< const MV > lhsMV = lp_->getLHS();
      Teuchos::RCP< const MV > rhsMV = lp_->getRHS();
      if (lhsMV == Teuchos::null && rhsMV == Teuchos::null) {
        stateStorageInitialized_ = false;
        return;
      }
      else {

        // Initialize the state storage
        // If the subspace has not be initialized before, generate it using the LHS or RHS from lp_.
        if (Y_ == Teuchos::null) {
          // Get the multivector that is not null.
          Teuchos::RCP< const MV > tmp = ( (rhsMV!=Teuchos::null)? rhsMV: lhsMV );
          TEUCHOS_TEST_FOR_EXCEPTION( tmp == Teuchos::null,
                              std::invalid_argument,
                              "Belos::MinresIter::setStateSize(): linear problem does not specify multivectors to clone from.");
          Y_  = MVT::Clone( *tmp, 1 );
          R1_ = MVT::Clone( *tmp, 1 );
          R2_ = MVT::Clone( *tmp, 1 );
          W_  = MVT::Clone( *tmp, 1 );
          W1_ = MVT::Clone( *tmp, 1 );
          W2_ = MVT::Clone( *tmp, 1 );
        }
        // State storage has now been initialized.
        stateStorageInitialized_ = true;
      }
    }
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Initialize this iteration object
  template <class ScalarType, class MV, class OP>
  void MinresIter<ScalarType,MV,OP>::initializeMinres(MinresIterationState<ScalarType,MV> newstate)
  {
    // Initialize the state storage if it isn't already.
    if (!stateStorageInitialized_) 
      setStateSize();

    TEUCHOS_TEST_FOR_EXCEPTION( !stateStorageInitialized_,
                        std::invalid_argument,
                        "Belos::MinresIter::initialize(): Cannot initialize state storage!" );

    TEUCHOS_TEST_FOR_EXCEPTION( newstate.Y == Teuchos::null,
                        std::invalid_argument,
                        "Belos::MinresIter::initialize(): MinresIterationState does not have initial residual.");

    std::string errstr("Belos::MinresIter::initialize(): Specified multivectors must have a consistent length and width.");
    TEUCHOS_TEST_FOR_EXCEPTION( MVText::GetGlobalLength(*newstate.Y) != MVText::GetGlobalLength(*Y_),
                        std::invalid_argument,
                        errstr );
    TEUCHOS_TEST_FOR_EXCEPTION( MVT::GetNumberVecs(*newstate.Y) != 1,
                        std::invalid_argument,
                        errstr );

    // Create convenience variables for zero, one.
    const ScalarType one = SCT::one();
    const MagnitudeType zero = SMT::zero();

    // Set up y and v for the first Lanczos vector v_1.
    // y  =  beta1_ P' v1,  where  P = C**(-1).
    // v is really P' v1.
    MVT::MvAddMv( one, *newstate.Y, zero, *newstate.Y, *R2_ );
    MVT::MvAddMv( one, *newstate.Y, zero, *newstate.Y, *R1_ );

    // Initialize the W's to 0.
    MVT::MvInit ( *W_ );
    MVT::MvInit ( *W2_ );

    if ( lp_->getLeftPrec() != Teuchos::null ) {
      lp_->applyLeftPrec( *newstate.Y, *Y_ );
    } 
    else {
      if (newstate.Y != Y_) {
        // copy over the initial residual (unpreconditioned).
        MVT::MvAddMv( one, *newstate.Y, zero, *newstate.Y, *Y_ );
      }
    }

    // beta1_ = b'*y;
    beta1_ = Teuchos::SerialDenseMatrix<int,ScalarType>( 1, 1 );
    MVT::MvTransMv( one, *newstate.Y, *Y_, beta1_ );

    TEUCHOS_TEST_FOR_EXCEPTION( SCT::real(beta1_(0,0)) < zero,
                        std::invalid_argument,
                        "The preconditioner is not positive definite." );

    if( SCT::magnitude(beta1_(0,0)) == zero )
    {
        // X = 0
        Teuchos::RCP<MV> cur_soln_vec = lp_->getCurrLHSVec();
        MVT::MvInit( *cur_soln_vec );
    }

    beta1_(0,0) = SCT::squareroot( beta1_(0,0) );

    // The solver is initialized
    initialized_ = true;
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Iterate until the status test informs us we should stop.
  template <class ScalarType, class MV, class OP>
  void MinresIter<ScalarType,MV,OP>::iterate()
  {
    //
    // Allocate/initialize data structures
    //
    if (initialized_ == false) {
      initialize();
    }

    Teuchos::BLAS<int,ScalarType> blas;

    // Create convenience variables for zero and one.
    const ScalarType one = SCT::one();
    const MagnitudeType zero = SMT::zero();

    // Allocate memory for scalars.
    Teuchos::SerialDenseMatrix<int,ScalarType> alpha( 1, 1 );
    Teuchos::SerialDenseMatrix<int,ScalarType> beta( beta1_ );
    phibar_ = Teuchos::ScalarTraits<ScalarType>::magnitude( beta1_(0,0) );
    ScalarType shift = zero; // TODO Allow for proper shift.

    // Initialize a few variables.
    ScalarType oldBeta = zero;
    ScalarType epsln = zero;
    ScalarType cs = -one;
    ScalarType sn = zero;
    ScalarType dbar = zero;

    // Declare a few others that will be initialized in the loop.
    ScalarType oldeps;
    ScalarType delta;
    ScalarType gbar;
    ScalarType phi;
    ScalarType gamma;

    // Allocate workspace.
    Teuchos::RCP<MV> V    = MVT::Clone( *Y_, 1 );
    Teuchos::RCP<MV> tmpY, tmpW;  // Not allocated, just used to transfer ownership.

    // Get the current solution vector.
    Teuchos::RCP<MV> cur_soln_vec = lp_->getCurrLHSVec();

    // Check that the current solution vector only has one column. 
    TEUCHOS_TEST_FOR_EXCEPTION( MVT::GetNumberVecs(*cur_soln_vec) != 1,
                        MinresIterateFailure,
                        "Belos::MinresIter::iterate(): current linear system has more than one vector!" );

    ////////////////////////////////////////////////////////////////
    // Iterate until the status test tells us to stop.
    //
    while (stest_->checkStatus(this) != Passed) {

      // Increment the iteration
      iter_++;

      // Normalize previous vector.
      //   v = y / beta(0,0);
      MVT::MvAddMv (one / beta(0,0), *Y_, zero, *Y_, *V);

      // Apply operator.
      lp_->applyOp (*V, *Y_);

      // Apply shift
      if (shift != zero)
	MVT::MvAddMv (one, *Y_, -shift, *V, *Y_);

      if (iter_ > 1)
	MVT::MvAddMv (one, *Y_, -beta(0,0)/oldBeta, *R1_, *Y_);

      // alpha := dot(V, Y_)
      MVT::MvTransMv (one, *V, *Y_, alpha);

      // y := y - alpha/beta r2
      MVT::MvAddMv (one, *Y_, -alpha(0,0)/beta(0,0), *R2_, *Y_);

      // r1 = r2;
      // r2 = y;
      tmpY = R1_;
      R1_ = R2_;
      R2_ = Y_;
      Y_ = tmpY;

      // apply left preconditioner
      if ( lp_->getLeftPrec() != Teuchos::null ) {
        lp_->applyLeftPrec( *R2_, *Y_ );
      } // else "y = r2"
      else {
        MVT::MvAddMv( one, *R2_, zero, *R2_, *Y_ );
      }

      // Get new beta.
      oldBeta = beta(0,0);
      MVT::MvTransMv( one, *R2_, *Y_, beta );

      // Intercept beta <= 0.
      //
      // Note: we don't try to test for nonzero imaginary component of
      // beta, because (a) it could be small and nonzero due to
      // rounding error in computing the inner product, and (b) it's
      // hard to tell how big "not small" should be, without computing
      // some error bounds (for example, by modifying the linear
      // algebra library to compute a posteriori rounding error bounds
      // for the inner product, and then changing
      // Belos::MultiVecTraits to make this information available).
      TEUCHOS_TEST_FOR_EXCEPTION( SCT::real(beta(0,0)) <= zero,
                          MinresIterateFailure,
                          "Belos::MinresIter::iterate(): Encountered nonpositi"
			  "ve value " << beta(0,0) << " for r2^H*M*r2 at itera"
			  "tion " << iter_ << ": MINRES cannot continue." );
      beta(0,0) = SCT::squareroot( beta(0,0) );

      // Apply previous rotation Q_{k-1} to get
      //
      //    [delta_k epsln_{k+1}] = [cs  sn][dbar_k  0         ]
      //    [gbar_k  dbar_{k+1} ]   [-sn cs][alpha_k beta_{k+1}].
      //
      oldeps = epsln;
      delta  = cs*dbar + sn*alpha(0,0);
      gbar   = sn*dbar - cs*alpha(0,0);
      epsln  =           sn*beta(0,0);
      dbar   =         - cs*beta(0,0);

      // Compute the next plane rotation Q_k.
      this->symOrtho(gbar, beta(0,0), &cs, &sn, &gamma);

      phi    = cs * phibar_; // phi_k
      phibar_ = Teuchos::ScalarTraits<ScalarType>::magnitude( sn * phibar_ ); // phibar_{k+1}

      //  w1 = w2;
      //  w2 = w;
      MVT::MvAddMv( one, *W_, zero, *W_, *W1_ );
      tmpW = W1_;
      W1_ = W2_;
      W2_ = W_;
      W_ = tmpW;

      //  w = (v - oldeps*w1 - delta*w2) / gamma;
      MVT::MvAddMv( one, *V, -oldeps, *W1_, *W_ );
      MVT::MvAddMv( one, *W_, -delta,  *W2_, *W_ );
      MVT::MvScale( *W_, one / gamma );

      // Update x:
      // x = x + phi*w;
      MVT::MvAddMv( one, *cur_soln_vec, phi, *W_, *cur_soln_vec );
      lp_->updateSolution();
    } // end while (sTest_->checkStatus(this) != Passed)
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Compute the next plane rotation Qk.
  //   r = norm([a b]);
  //   c = a / r;
  //   s = b / r;
  template <class ScalarType, class MV, class OP>
  void MinresIter<ScalarType,MV,OP>::symOrtho( ScalarType a, ScalarType b,
                                               ScalarType *c, ScalarType *s, ScalarType *r
                                             )
  {
    const ScalarType one = SCT::one();
    const ScalarType zero = SCT::zero();
    const MagnitudeType m_zero = SMT::zero();
    const MagnitudeType absA = SCT::magnitude( a );
    const MagnitudeType absB = SCT::magnitude( b );
    if ( absB == m_zero ) {
        *s = zero;
        *r = absA;
        if ( absA == m_zero )
            *c = one;
        else
            *c = a / absA;
    } else if ( absA == m_zero ) {
        *c = zero;
        *s = b / absB;
        *r = absB;
    } else if ( absB >= absA ) { // && a!=0 && b!=0
        ScalarType tau = a / b;
        if ( Teuchos::ScalarTraits<ScalarType>::real(b) < m_zero )
            *s = -one / SCT::squareroot( one+tau*tau );
        else
            *s =  one / SCT::squareroot( one+tau*tau );
        *c = *s * tau;
        *r = b / *s;
    } else { // (absA > absB) && a!=0 && b!=0
        ScalarType tau = b / a;
        if ( Teuchos::ScalarTraits<ScalarType>::real(a) < m_zero )
            *c = -one / SCT::squareroot( one+tau*tau );
        else
            *c =  one / SCT::squareroot( one+tau*tau );
        *s = *c * tau;
        *r = a / *c;
    }
  }

} // end Belos namespace

#endif /* BELOS_MINRES_ITER_HPP */
