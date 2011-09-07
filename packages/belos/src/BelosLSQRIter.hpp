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

#ifndef BELOS_LSQR_ITER_HPP
#define BELOS_LSQR_ITER_HPP

/*! \file LSQRIter.hpp
    \brief Belos concrete class that iterates LSQR.
*/

#include "BelosConfigDefs.hpp"
#include "BelosTypes.hpp"
#include "BelosLSQRIteration.hpp"

#include "BelosLinearProblem.hpp"
//#include "BelosMatOrthoManager.hpp"
#include "BelosOutputManager.hpp"
#include "BelosStatusTest.hpp"
#include "BelosOperatorTraits.hpp"
#include "BelosMultiVecTraits.hpp"

//#include "Teuchos_BLAS.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TimeMonitor.hpp"

/*!	
  \class Belos::LSQRIter
  
  \brief This class implements the LSQR iteration.

  \ingroup belos_solver_framework 

  \author David Day

*/
  
namespace Belos {

template<class ScalarType, class MV, class OP>
class LSQRIter : virtual public Belos::Iteration<ScalarType,MV,OP> {

  public:
    
  //
  // Convenience typedefs
  //
  typedef Belos::MultiVecTraits<ScalarType,MV> MVT;
  typedef Belos::OperatorTraits<ScalarType,MV,OP> OPT;
  typedef Teuchos::ScalarTraits<ScalarType> SCT;
  typedef typename SCT::magnitudeType MagnitudeType;

  //! @name Constructors/Destructor
  //@{ 

  /*! \brief %LSQRIter constructor with linear problem, solver utilities, and parameter
   * list of solver options.
   *
   * This constructor takes pointers required by the linear solver iteration, in addition
   * to a parameter list of options for the linear solver.
   */
  LSQRIter( const Teuchos::RCP<Belos::LinearProblem<ScalarType,MV,OP> > &problem, 
	    const Teuchos::RCP<Belos::OutputManager<ScalarType> > &printer,
	    const Teuchos::RCP<Belos::StatusTest<ScalarType,MV,OP> > &tester,
		  Teuchos::ParameterList &params );

//                  const Teuchos::RCP<MatOrthoManager<ScalarType,MV,OP> > &ortho,


  //! Destructor.
  virtual ~LSQRIter() {};
  //@}


  //! @name Solver methods
  //@{ 
  
  /*! \brief This method performs LSQR iterations until the status
   * test indicates the need to stop or an error occurs (in which case, an
   * std::exception is thrown).
   *
   * iterate() will first determine whether the solver is initialized; if
   * not, it will call initialize() using default arguments. After
   * initialization, the solver performs LSQR iterations until the
   * status test evaluates as ::Passed, at which point the method returns to
   * the caller. 
   *
   * The status test is queried at the beginning of the iteration.
   */
  void iterate();

  /*! \brief Initialize the solver to an iterate, providing a complete state.
   *
   * The %LSQRIter contains a certain amount of state, consisting of two bidiagonalization 
   * vectors, a descent direction, a damping value, and various estimates of errors and
   * the like.
   *
   * \note For any pointer in \c newstate which directly points to the multivectors in 
   * the solver, the data is not copied.
   */
  void initializeLSQR(LSQRIterationState<ScalarType,MV> newstate);

  /*! \brief Initialize the solver.
   */
  void initialize()
  {
    LSQRIterationState<ScalarType,MV> empty;
    initializeLSQR(empty);
  }
  
  /*! \brief Get the current state of the linear solver.
   *
   * The data is only valid if isInitialized() == \c true.
   *
   * \returns A LSQRIterationState object containing const pointers to the current solver
   * state.
   */
  LSQRIterationState<ScalarType,MV> getState() const {
    LSQRIterationState<ScalarType,MV> state;
    state.U = U_;
    state.V = V_;
    state.W = W_;
    state.lambda = lambda_;
    state.resid_norm = resid_norm_;
    state.frob_mat_norm = frob_mat_norm_;
    state.mat_cond_num = mat_cond_num_;
    state.mat_resid_norm = mat_resid_norm_; 
    state.sol_norm = sol_norm_;
    state.bnorm = bnorm_;
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
  //! \return This method returns a null pointer because residuals aren't used with LSQR.
  Teuchos::RCP<const MV> getNativeResiduals( std::vector<MagnitudeType> *norms ) const { return Teuchos::null; }

  //! Get the current update to the linear system.
  /*! \note This method returns a null pointer because the linear problem is current.
  */
  Teuchos::RCP<MV> getCurrentUpdate() const { return Teuchos::null; }

  //@}
  
  //! @name Accessor methods
  //@{ 

  //! Get a constant reference to the linear problem.
  const Belos::LinearProblem<ScalarType,MV,OP>& getProblem() const { return *lp_; }

  //! Get the blocksize to be used by the iterative solver in solving this linear problem.
  int getBlockSize() const { return 1; }

  //! \brief Set the blocksize to be used by the iterative solver to solve this linear problem.
  //This is unique to single vector methods.
  void setBlockSize(int blockSize) {
    TEST_FOR_EXCEPTION(blockSize!=1,std::invalid_argument,
		       "LSQRIter::setBlockSize(): Cannot use a block size that is not one.");
  }

  //! States whether the solver has been initialized or not.
  bool isInitialized() { return initialized_; }

  //@}

  private:

  //
  // Internal methods
  //
  //! Method for initalizing the state storage needed by LSQR.
  void setStateSize();
  
  //
  // Classes inputed through constructor that define the linear problem to be solved.
  //
  const Teuchos::RCP<Belos::LinearProblem<ScalarType,MV,OP> >    lp_;
  const Teuchos::RCP<Belos::OutputManager<ScalarType> >          om_;
  const Teuchos::RCP<Belos::StatusTest<ScalarType,MV,OP> >       stest_;
//  const Teuchos::RCP<OrthoManager<ScalarType,MV> >        ortho_;

  //  
  // Current solver state
  //
  // initialized_ specifies that the basis vectors have been initialized and the iterate()
  // routine is capable of running; _initialize is controlled  by the initialize() member
  // method.  For the implications of the state of initialized_, please see documentation
  // for initialize()
  bool initialized_;

  // stateStorageInitialized_ specifies that the state storage has been initialized.
  // This initialization may be postponed if the linear problem was generated without 
  // the right-hand side or solution vectors.
  bool stateStorageInitialized_;

  // Current number of iterations performed.
  int iter_;
  
  // 
  // State Storage
  // 
  //
  // Bidiagonalization vector
  Teuchos::RCP<MV> U_;
  //
  // Bidiagonalization vector
  Teuchos::RCP<MV> V_;
  //
  // Direction vector
  Teuchos::RCP<MV> W_;
  //
  // Damping value
  ScalarType lambda_;
  //
  // Residual norm estimate
  ScalarType resid_norm_;
  //
  // Frobenius norm estimate
  ScalarType frob_mat_norm_;
  //
  // Condition number estimate
  ScalarType mat_cond_num_;
  //
  // A^T*resid norm estimate
  ScalarType mat_resid_norm_;
  //
  // Solution norm estimate
  ScalarType sol_norm_;
  //
  // RHS norm
  ScalarType bnorm_;

};

  /////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor.

//                                     const Teuchos::RCP<MatOrthoManager<ScalarType,MV,OP> > &ortho,

  template<class ScalarType, class MV, class OP>
  LSQRIter<ScalarType,MV,OP>::LSQRIter(const Teuchos::RCP<Belos::LinearProblem<ScalarType,MV,OP> > &problem, 
				       const Teuchos::RCP<Belos::OutputManager<ScalarType> > &printer,
				       const Teuchos::RCP<Belos::StatusTest<ScalarType,MV,OP> > &tester,
						   Teuchos::ParameterList &params ):
    lp_(problem),
    om_(printer),
    stest_(tester),
    initialized_(false),
    stateStorageInitialized_(false),
    iter_(0),
    lambda_(params.get("Lambda", 0.0))
  {
  }

  /////////////////////////////////////////////////////////////////////////////////////////////
  // Setup the state storage.
  template <class ScalarType, class MV, class OP>
  void LSQRIter<ScalarType,MV,OP>::setStateSize ()
  {
    if (!stateStorageInitialized_) {
      // Check if there is any multivector to clone from.
      Teuchos::RCP<const MV> lhsMV = lp_->getLHS();
      Teuchos::RCP<const MV> rhsMV = lp_->getRHS();
      if (lhsMV == Teuchos::null || rhsMV == Teuchos::null) {
	stateStorageInitialized_ = false;
	return;
      }
      else {
	// Initialize the state storage
	// If the subspace has not been initialized before, generate it
        // using the LHS and RHS from lp_.
	if (U_ == Teuchos::null) {
	  // Get the multivectors.
	  TEST_FOR_EXCEPTION(rhsMV == Teuchos::null, std::invalid_argument, "LSQRIter::setStateSize(): linear problem does not specify right hand multivector to clone from.");
	  TEST_FOR_EXCEPTION(lhsMV == Teuchos::null, std::invalid_argument, "LSQRIter::setStateSize(): linear problem does not specify left hand multivector to clone from.");

	  U_ = MVT::Clone( *rhsMV, 1 );
	  V_ = MVT::Clone( *lhsMV, 1 );
	  W_ = MVT::Clone( *lhsMV, 1 );
	}
	
	// State storage has now been initialized.
	stateStorageInitialized_ = true;
      }
    }
  }


  /////////////////////////////////////////////////////////////////////////////////////////////
  // Initialize this iteration object
  template <class ScalarType, class MV, class OP>
  void LSQRIter<ScalarType,MV,OP>::initializeLSQR(LSQRIterationState<ScalarType,MV> newstate)
  {
    // Initialize the state storage if it isn't already.
    if (!stateStorageInitialized_) 
      setStateSize();

    TEST_FOR_EXCEPTION(!stateStorageInitialized_,std::invalid_argument,
		       "LSQRIter::initialize(): Cannot initialize state storage!");
    
    std::string errstr("LSQRIter::initialize(): Specified multivectors must have a consistent length and width.");

    // Create convenience variables for zero and one.
    const ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();
    const ScalarType zero = Teuchos::ScalarTraits<ScalarType>::zero();

    // Compute initial bidiagonalization vectors and search direction
    //
    Teuchos::RCP<const MV> lhsMV = lp_->getLHS();
    Teuchos::RCP<const MV> rhsMV = lp_->getRHS();

    OPT::Apply(*(lp_->getOperator()), *lhsMV, *U_);
    MVT::MvAddMv( one, *rhsMV, -one, *U_, *U_);
    OPT::Apply(*(lp_->getOperator()), *U_, *V_, CONJTRANS);
    MVT::MvAddMv( one, *V_, zero, *V_, *W_);
    
    frob_mat_norm_ = zero;
    mat_cond_num_ = zero;
    sol_norm_ = zero;
    
    // The solver is initialized
    initialized_ = true;
  }


  /////////////////////////////////////////////////////////////////////////////////////////////
  // Iterate until the status test informs us we should stop.
  template <class ScalarType, class MV, class OP>
  void LSQRIter<ScalarType,MV,OP>::iterate()
  {
    //
    // Allocate/initialize data structures
    //
    if (initialized_ == false) {
      initialize();
    }
    
    // Create convenience variables for zero and one.
    const ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();
    const MagnitudeType zero = Teuchos::ScalarTraits<MagnitudeType>::zero();

    // Allocate memory for scalars.
    std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> alpha(1);
    std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> beta(1);
    std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> wnorm2(1);
    ScalarType rhobar, phibar, cs1, phi, rho, cs, sn, theta, xxnorm = zero, common;
    ScalarType zetabar, sn1, psi, res = zero, bbnorm = zero, ddnorm = zero, gamma, tau;
    ScalarType cs2 = -one, sn2 = zero, gammabar, zeta = zero, delta;
    
    // Allocate memory for working vectors.
    // Operator applied to bidiagonalization vector
    Teuchos::RCP<MV> AV;
    // Transpose of operator applied to bidiagonalization vector
    Teuchos::RCP<MV> AtU;
    AV = MVT::Clone( *U_, 1);
    AtU = MVT::Clone( *V_, 1);

    // Get the current solution vector.
    Teuchos::RCP<MV> cur_soln_vec = lp_->getCurrLHSVec();

    // Check that the current solution vector only has one column. 
    TEST_FOR_EXCEPTION( MVT::GetNumberVecs(*cur_soln_vec) != 1, LSQRIterateFailure,
                        "LSQRIter::iterate(): current linear system has more than one vector!" );

    // Compute alpha and beta and scale bidiagonalization vectors
    MVT::MvNorm( *U_, beta );
    MVT::MvScale( *U_, one / beta[0] );
    MVT::MvScale( *V_, one / beta[0] );
    MVT::MvNorm( *V_, alpha );
    MVT::MvScale( *V_, one / alpha[0] );
    MVT::MvScale( *W_, one / (beta[0] * alpha[0]) );

    rhobar = alpha[0];
    phibar = beta[0];
    
    resid_norm_ = beta[0];
    mat_resid_norm_ = alpha[0] * beta[0];
    bnorm_ = beta[0];

    ////////////////////////////////////////////////////////////////
    // Iterate until the status test tells us to stop.
    //
    while (stest_->checkStatus(this) != Belos::Passed) {
      // Increment the iteration
      iter_++;

      // Perform the next step of the bidiagonalization to obtain the next U_ and V_ vectors
      // and the scalars alpha and beta.  These satisfy the relations
      // beta*U_ = AV - alpha*U_
      // alpha*V_ = AtU - beta*V_

      // AV := A * V_
      OPT::Apply(*(lp_->getOperator()), *V_, *AV);
      MVT::MvAddMv( one, *AV, -alpha[0], *U_, *U_ );
      MVT::MvNorm( *U_, beta);
      // Check that beta is a positive number!
      TEST_FOR_EXCEPTION( SCT::real(beta[0]) <= zero, LSQRIterateFailure, "LSQRIter::iterate(): non-positive value for beta encountered!");
      bbnorm += alpha[0]*alpha[0] + beta[0]*beta[0] + lambda_*lambda_;
      MVT::MvScale( *U_, one / beta[0] );

      // AtU := A^T * U_
      OPT::Apply(*(lp_->getOperator()), *U_, *AtU, CONJTRANS);
      MVT::MvAddMv( one, *AtU, -beta[0], *V_, *V_ );
      MVT::MvNorm( *V_, alpha );
      // Check that alpha is a positive number!
      TEST_FOR_EXCEPTION( SCT::real(alpha[0]) <= zero, LSQRIterateFailure, "LSQRIter::iterate(): non-positive value for alpha encountered!");
      MVT::MvScale( *V_, one / alpha[0] );

      // Use a plane rotation to eliminate the damping parameter.  
      // This alters the diagonal (rhobar) of the lower-bidiagonal matrix.
      common = Teuchos::ScalarTraits< ScalarType >::squareroot(rhobar*rhobar + lambda_*lambda_);
      cs1 = rhobar / common;
      sn1 = lambda_ / common;
      psi = sn1 * phibar;
      phibar = cs1 * phibar;

      // Use a plane rotation to eliminate the subdiagonal element (beta)
      // of the lower-bidiagonal matrix, giving an upper-bidiagonal matrix.
      rho = Teuchos::ScalarTraits< ScalarType >::squareroot(rhobar*rhobar + lambda_*lambda_ + beta[0]*beta[0]);
      cs = common / rho; 
      sn = beta[0] / rho;
      theta = sn * alpha[0];
      rhobar = -cs * alpha[0];
      phi = cs * phibar;
      phibar = sn * phibar;
      tau = sn * phi;

      delta = sn2 * rho;
      gammabar = -cs2 * rho;
      zetabar = (phi - delta*zeta) / gammabar;
      sol_norm_ = Teuchos::ScalarTraits< ScalarType >::squareroot(xxnorm + zetabar*zetabar);
      gamma = Teuchos::ScalarTraits< ScalarType >::squareroot(gammabar*gammabar + theta*theta);
      cs2 = gammabar / gamma;
      sn2 = theta / gamma;
      zeta = (phi - delta*zeta) / gamma;
      xxnorm += zeta*zeta;

      // Update the solution vector and search direction vector
      MVT::MvAddMv( phi / rho, *W_, one, *cur_soln_vec, *cur_soln_vec);
      lp_->updateSolution();
      MVT::MvNorm( *W_, wnorm2 );
      ddnorm += (one / rho)*(one / rho) * wnorm2[0]*wnorm2[0];
      MVT::MvAddMv( one, *V_, -theta / rho, *W_, *W_ );

      frob_mat_norm_ = Teuchos::ScalarTraits< ScalarType >::squareroot(bbnorm);
      mat_cond_num_ = frob_mat_norm_ * Teuchos::ScalarTraits< ScalarType >::squareroot(ddnorm);
      res+= psi*psi;
      resid_norm_ = Teuchos::ScalarTraits< ScalarType >::squareroot(phibar*phibar + res);
      mat_resid_norm_ = alpha[0] * Teuchos::ScalarTraits< ScalarType >::magnitude(tau);

    } // end while (sTest_->checkStatus(this) != Passed)
  } // iterate()

} // end Belos namespace

#endif /* BELOS_LSQR_ITER_HPP */
