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

#ifndef BELOS_BLOCK_CG_ITER_HPP
#define BELOS_BLOCK_CG_ITER_HPP

/*! \file BelosBlockCGIter.hpp
    \brief Belos concrete class for performing the block conjugate-gradient (CG) iteration.
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

#include "Teuchos_BLAS.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_SerialSymDenseMatrix.hpp"
#include "Teuchos_SerialSpdDenseSolver.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TimeMonitor.hpp"

/*!	
  \class Belos::BlockCGIter
  
  \brief This class implements the block, preconditioned Conjugate Gradient (CG) iteration.

  \ingroup belos_solver_framework
 
  \author Teri Barth and Heidi Thornquist
*/

namespace Belos {
  
template<class ScalarType, class MV, class OP>
class BlockCGIter : virtual public CGIteration<ScalarType,MV,OP> {

  public:
    
  //
  // Convenience typedefs
  //
  typedef MultiVecTraits<ScalarType,MV> MVT;
  typedef MultiVecTraitsExt<ScalarType,MV> MVText;
  typedef OperatorTraits<ScalarType,MV,OP> OPT;
  typedef Teuchos::ScalarTraits<ScalarType> SCT;
  typedef typename SCT::magnitudeType MagnitudeType;

  //! @name Constructors/Destructor
  //@{ 

  /*! \brief %BlockCGIter constructor with linear problem, solver utilities, and parameter list of solver options.
   *
   * This constructor takes pointers required by the linear solver iteration, in addition
   * to a parameter list of options for the linear solver.
   */
  BlockCGIter( const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem, 
	       const Teuchos::RCP<OutputManager<ScalarType> > &printer,
	       const Teuchos::RCP<StatusTest<ScalarType,MV,OP> > &tester,
	       const Teuchos::RCP<MatOrthoManager<ScalarType,MV,OP> > &ortho,
	       Teuchos::ParameterList &params );

  //! Destructor.
  virtual ~BlockCGIter() {};
  //@}


  //! @name Solver methods
  //@{ 
  
  /*! \brief This method performs BlockCG iterations until the status
   * test indicates the need to stop or an error occurs (in which case, an
   * std::exception is thrown).
   *
   * iterate() will first determine whether the solver is initialized; if
   * not, it will call initialize() using default arguments. After
   * initialization, the solver performs BlockCG iterations until the
   * status test evaluates as ::Passed, at which point the method returns to
   * the caller. 
   *
   * The status test is queried at the beginning of the iteration.
   */
  void iterate();

  /*! \brief Initialize the solver to an iterate, providing a complete state.
   *
   * The %BlockCGIter contains a certain amount of state, consisting of the current 
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
  void initializeCG(CGIterationState<ScalarType,MV> newstate);

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
  void resetNumIters( int iter=0 ) { iter_ = iter; }

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

  //! Get the block size to be used by the iterative solver in solving this linear problem.
  int getBlockSize() const { return blockSize_; }

  //! \brief Set the block size to be used by the iterative solver in solving this linear problem.
  void setBlockSize(int blockSize);

  //! States whether the solver has been initialized or not.
  bool isInitialized() { return initialized_; }

  //@}

  private:

  //
  // Internal methods
  //
  //! Method for initalizing the state storage needed by block CG.
  void setStateSize();
  
  //
  // Classes inputed through constructor that define the linear problem to be solved.
  //
  const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> >    lp_;
  const Teuchos::RCP<OutputManager<ScalarType> >          om_;
  const Teuchos::RCP<StatusTest<ScalarType,MV,OP> >       stest_;
  const Teuchos::RCP<OrthoManager<ScalarType,MV> >        ortho_;

  //
  // Algorithmic parameters
  //
  // blockSize_ is the solver block size.
  int blockSize_;

  //  
  // Current solver state
  //
  // initialized_ specifies that the basis vectors have been initialized and the iterate() routine
  // is capable of running; _initialize is controlled  by the initialize() member method
  // For the implications of the state of initialized_, please see documentation for initialize()
  bool initialized_;

  // stateStorageInitialized_ specified that the state storage has be initialized.
  // This initialization may be postponed if the linear problem was generated without 
  // the right-hand side or solution vectors.
  bool stateStorageInitialized_;

  // Current subspace dimension, and number of iterations performed.
  int iter_;
  
  // 
  // State Storage
  // 
  // Residual
  Teuchos::RCP<MV> R_;
  //
  // Preconditioned residual
  Teuchos::RCP<MV> Z_;
  //
  // Direction std::vector
  Teuchos::RCP<MV> P_;
  //
  // Operator applied to direction std::vector
  Teuchos::RCP<MV> AP_;

};

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor.
  template<class ScalarType, class MV, class OP>
  BlockCGIter<ScalarType,MV,OP>::BlockCGIter(const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem, 
					     const Teuchos::RCP<OutputManager<ScalarType> > &printer,
					     const Teuchos::RCP<StatusTest<ScalarType,MV,OP> > &tester,
					     const Teuchos::RCP<MatOrthoManager<ScalarType,MV,OP> > &ortho,
					     Teuchos::ParameterList &params ):
    lp_(problem),
    om_(printer),
    stest_(tester),
    ortho_(ortho),
    blockSize_(0),
    initialized_(false),
    stateStorageInitialized_(false),
    iter_(0)
  {
    // Set the block size and allocate data
    int bs = params.get("Block Size", 1);
    setBlockSize( bs );
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Setup the state storage.
  template <class ScalarType, class MV, class OP>
  void BlockCGIter<ScalarType,MV,OP>::setStateSize ()
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
	if (R_ == Teuchos::null || MVT::GetNumberVecs(*R_)!=blockSize_) {
	  // Get the multivector that is not null.
	  Teuchos::RCP<const MV> tmp = ( (rhsMV!=Teuchos::null)? rhsMV: lhsMV );
	  TEUCHOS_TEST_FOR_EXCEPTION(tmp == Teuchos::null,std::invalid_argument,
			     "Belos::BlockCGIter::setStateSize(): linear problem does not specify multivectors to clone from.");
	  R_ = MVT::Clone( *tmp, blockSize_ );
	  Z_ = MVT::Clone( *tmp, blockSize_ );
	  P_ = MVT::Clone( *tmp, blockSize_ );
	  AP_ = MVT::Clone( *tmp, blockSize_ );
	}
	
	// State storage has now been initialized.
	stateStorageInitialized_ = true;
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Set the block size and make necessary adjustments.
  template <class ScalarType, class MV, class OP>
  void BlockCGIter<ScalarType,MV,OP>::setBlockSize (int blockSize)
  {
    // This routine only allocates space; it doesn't not perform any computation
    // any change in size will invalidate the state of the solver.

    TEUCHOS_TEST_FOR_EXCEPTION(blockSize <= 0, std::invalid_argument, "Belos::BlockGmresIter::setBlockSize was passed a non-positive argument.");
    if (blockSize == blockSize_) {
      // do nothing
      return;
    }

    if (blockSize!=blockSize_)
      stateStorageInitialized_ = false;

    blockSize_ = blockSize;

    initialized_ = false;

    // Use the current blockSize_ to initialize the state storage.
    setStateSize();

  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Initialize this iteration object
  template <class ScalarType, class MV, class OP>
  void BlockCGIter<ScalarType,MV,OP>::initializeCG(CGIterationState<ScalarType,MV> newstate)
  {
    // Initialize the state storage if it isn't already.
    if (!stateStorageInitialized_) 
      setStateSize();

    TEUCHOS_TEST_FOR_EXCEPTION(!stateStorageInitialized_,std::invalid_argument,
		       "Belos::BlockCGIter::initialize(): Cannot initialize state storage!");
    
    // NOTE:  In BlockCGIter R_, the initial residual, is required!!!  
    //
    std::string errstr("Belos::BlockCGIter::initialize(): Specified multivectors must have a consistent length and width.");

    // Create convenience variables for zero and one.
    const ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();
    const MagnitudeType zero = Teuchos::ScalarTraits<MagnitudeType>::zero();

    if (newstate.R != Teuchos::null) {

      TEUCHOS_TEST_FOR_EXCEPTION( MVText::GetGlobalLength(*newstate.R) != MVText::GetGlobalLength(*R_),
                          std::invalid_argument, errstr );
      TEUCHOS_TEST_FOR_EXCEPTION( MVT::GetNumberVecs(*newstate.R) != blockSize_,
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
          Teuchos::RCP<MV> tmp = MVT::Clone( *Z_, blockSize_ );
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
      MVT::MvAddMv( one, *Z_, zero, *Z_, *P_ );
    }
    else {

      TEUCHOS_TEST_FOR_EXCEPTION(newstate.R == Teuchos::null,std::invalid_argument,
                         "Belos::BlockCGIter::initialize(): BlockCGStateIterState does not have initial residual.");
    }

    // The solver is initialized
    initialized_ = true;
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Iterate until the status test informs us we should stop.
  template <class ScalarType, class MV, class OP>
  void BlockCGIter<ScalarType,MV,OP>::iterate()
  {
    //
    // Allocate/initialize data structures
    //
    if (initialized_ == false) {
      initialize();
    }
    // Allocate data needed for LAPACK work.
    int info = 0;
    char UPLO = 'U';
    bool uplo = true;
    Teuchos::LAPACK<int,ScalarType> lapack;

    // Allocate memory for scalars.
    Teuchos::SerialDenseMatrix<int,ScalarType> alpha( blockSize_, blockSize_ );
    Teuchos::SerialDenseMatrix<int,ScalarType> beta( blockSize_, blockSize_ );
    Teuchos::SerialDenseMatrix<int,ScalarType> rHz( blockSize_, blockSize_ ), 
      rHz_old( blockSize_, blockSize_ ), pAp( blockSize_, blockSize_ );
    Teuchos::SerialSymDenseMatrix<int,ScalarType> pApHerm(Teuchos::View, uplo, pAp.values(), blockSize_, blockSize_);

    // Create dense spd solver.
    Teuchos::SerialSpdDenseSolver<int,ScalarType> lltSolver;

    // Create convenience variables for zero and one.
    const ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();
    
    // Get the current solution std::vector.
    Teuchos::RCP<MV> cur_soln_vec = lp_->getCurrLHSVec();

    // Check that the current solution std::vector has blockSize_ columns. 
    TEUCHOS_TEST_FOR_EXCEPTION( MVT::GetNumberVecs(*cur_soln_vec) != blockSize_, CGIterateFailure,
                        "Belos::BlockCGIter::iterate(): current linear system does not have the right number of vectors!" );
    int rank = ortho_->normalize( *P_, Teuchos::null );
    TEUCHOS_TEST_FOR_EXCEPTION(rank != blockSize_,CGIterationOrthoFailure,
                         "Belos::BlockCGIter::iterate(): Failed to compute initial block of orthonormal direction vectors.");


    ////////////////////////////////////////////////////////////////
    // Iterate until the status test tells us to stop.
    //
    while (stest_->checkStatus(this) != Passed) {
        
      // Increment the iteration
      iter_++;
    
      // Multiply the current direction std::vector by A and store in Ap_
      lp_->applyOp( *P_, *AP_ );
      
      // Compute alpha := <P_,R_> / <P_,AP_>
      // 1) Compute P^T * A * P = pAp and P^T * R 
      // 2) Compute the Cholesky Factorization of pAp
      // 3) Back and forward solves to compute alpha
      //
      MVT::MvTransMv( one, *P_, *R_, alpha );
      MVT::MvTransMv( one, *P_, *AP_, pAp );      
     
      // Compute Cholesky factorization of pAp
      lltSolver.setMatrix( Teuchos::rcp(&pApHerm, false) );
      lltSolver.factorWithEquilibration( true );
      info = lltSolver.factor();
      TEUCHOS_TEST_FOR_EXCEPTION(info != 0,CGIterationLAPACKFailure,
                         "Belos::BlockCGIter::iterate(): Failed to compute Cholesky factorization using LAPACK routine POTRF.");

      // Compute alpha by performing a back and forward solve with the Cholesky factorization in pAp.
      lltSolver.setVectors( Teuchos::rcp( &alpha, false ), Teuchos::rcp( &alpha, false ) );
      info = lltSolver.solve();
      TEUCHOS_TEST_FOR_EXCEPTION(info != 0,CGIterationLAPACKFailure,
                         "Belos::BlockCGIter::iterate(): Failed to compute alpha using Cholesky factorization (POTRS).");
      
      //
      // Update the solution std::vector X := X + alpha * P_
      //
      MVT::MvTimesMatAddMv( one, *P_, alpha, one, *cur_soln_vec );
      lp_->updateSolution();
      //
      // Compute the new residual R_ := R_ - alpha * AP_
      //
      MVT::MvTimesMatAddMv( -one, *AP_, alpha, one, *R_ );
      //
      // Compute the new preconditioned residual, Z_.
      if ( lp_->getLeftPrec() != Teuchos::null ) {
        lp_->applyLeftPrec( *R_, *Z_ );
        if ( lp_->getRightPrec() != Teuchos::null ) {
          Teuchos::RCP<MV> tmp = MVT::Clone( *Z_, blockSize_ );
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
      // Compute beta := <AP_,Z_> / <P_,AP_> 
      // 1) Compute AP_^T * Z_ 
      // 2) Compute the Cholesky Factorization of pAp (already have)
      // 3) Back and forward solves to compute beta

      // Compute <AP_,Z>
      MVT::MvTransMv( -one, *AP_, *Z_, beta );
      //
      lltSolver.setVectors( Teuchos::rcp( &beta, false ), Teuchos::rcp( &beta, false ) );
      info = lltSolver.solve();
      TEUCHOS_TEST_FOR_EXCEPTION(info != 0,CGIterationLAPACKFailure,
                         "Belos::BlockCGIter::iterate(): Failed to compute beta using Cholesky factorization (POTRS).");
      //
      // Compute the new direction vectors P_ = Z_ + P_ * beta 
      //
      Teuchos::RCP<MV> Pnew = MVT::CloneCopy( *Z_ );
      MVT::MvTimesMatAddMv(one, *P_, beta, one, *Pnew);
      P_ = Pnew;

      // Compute orthonormal block of new direction vectors.
      rank = ortho_->normalize( *P_, Teuchos::null );
      TEUCHOS_TEST_FOR_EXCEPTION(rank != blockSize_,CGIterationOrthoFailure,
                         "Belos::BlockCGIter::iterate(): Failed to compute block of orthonormal direction vectors.");
      
    } // end while (sTest_->checkStatus(this) != Passed)
  }

} // end Belos namespace

#endif /* BELOS_BLOCK_CG_ITER_HPP */
