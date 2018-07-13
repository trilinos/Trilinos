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

#ifndef BELOS_CG_SINGLE_RED_ITER_HPP
#define BELOS_CG_SINGLE_RED_ITER_HPP

/*! \file BelosCGSingleRedIter.hpp
    \brief Belos concrete class for performing a single-reduction conjugate-gradient (CG) iteration.
*/

#include "BelosConfigDefs.hpp"
#include "BelosTypes.hpp"
#include "BelosCGIteration.hpp"

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

/*!	
  \class Belos::CGSingleRedIter
  
  \brief This class implements the preconditioned single-reduction Conjugate Gradient (CG) iteration.

  \ingroup belos_solver_framework
 
  \author Heidi Thornquist
*/

namespace Belos {
  
template<class ScalarType, class MV, class OP>
class CGSingleRedIter : virtual public CGIteration<ScalarType,MV,OP> {

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

  /*! \brief %CGSingleRedIter constructor with linear problem, solver utilities, and parameter list of solver options.
   *
   * This constructor takes pointers required by the linear solver iteration, in addition
   * to a parameter list of options for the linear solver.
   */
  CGSingleRedIter( const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem, 
		  const Teuchos::RCP<OutputManager<ScalarType> > &printer,
		  const Teuchos::RCP<StatusTest<ScalarType,MV,OP> > &tester,
		  Teuchos::ParameterList &params );

  //! Destructor.
  virtual ~CGSingleRedIter() {};
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
   * The %CGSingleRedIter contains a certain amount of state, consisting of the current 
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

  //! \brief Set the blocksize to be used by the iterative solver in solving this linear problem.
  void setBlockSize(int blockSize) {
    TEUCHOS_TEST_FOR_EXCEPTION(blockSize!=1,std::invalid_argument,
		       "Belos::CGSingleRedIter::setBlockSize(): Cannot use a block size that is not one.");
  }

  //! States whether the solver has been initialized or not.
  bool isInitialized() { return initialized_; }

  //! Sets whether or not to store the diagonal for condition estimation
  void setDoCondEst(bool val){/*ignored*/}

  //! Gets the diagonal for condition estimation (NOT_IMPLEMENTED)
  Teuchos::ArrayView<MagnitudeType> getDiag() { 
    Teuchos::ArrayView<MagnitudeType> temp;
    return temp;
  }

  //! Gets the off-diagonal for condition estimation (NOT_IMPLEMENTED)
  Teuchos::ArrayView<MagnitudeType> getOffDiag() {
    Teuchos::ArrayView<MagnitudeType> temp;
    return temp;
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
  
  // 
  // State Storage
  // 
  // Residual
  Teuchos::RCP<MV> R_;
  //
  // Preconditioned residual
  Teuchos::RCP<MV> Z_;
  //
  // Operator applied to preconditioned residual
  Teuchos::RCP<MV> AZ_;
  //
  // This is the additional storage needed for single-reduction CG (Chronopoulos/Gear)
  // R_ views the first vector and AZ_ views the second vector
  Teuchos::RCP<MV> S_;
  //
  // Direction vector
  Teuchos::RCP<MV> P_;
  //
  // Operator applied to direction vector (S)
  Teuchos::RCP<MV> AP_;

};

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor.
  template<class ScalarType, class MV, class OP>
  CGSingleRedIter<ScalarType,MV,OP>::CGSingleRedIter(const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem, 
						     const Teuchos::RCP<OutputManager<ScalarType> > &printer,
						     const Teuchos::RCP<StatusTest<ScalarType,MV,OP> > &tester,
						     Teuchos::ParameterList &params ):
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
  void CGSingleRedIter<ScalarType,MV,OP>::setStateSize ()
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
			     "Belos::CGSingleRedIter::setStateSize(): linear problem does not specify multivectors to clone from.");
	  S_ = MVT::Clone( *tmp, 2 );
	  Z_ = MVT::Clone( *tmp, 1 );
	  P_ = MVT::Clone( *tmp, 1 );
	  AP_ = MVT::Clone( *tmp, 1 );

          // R_ will view the first column of S_, AZ_ will view the second.
          std::vector<int> index(1,0);
	  R_ = MVT::CloneViewNonConst( *S_, index );
          index[0] = 1;
          AZ_ = MVT::CloneViewNonConst( *S_, index );
	}
	
	// State storage has now been initialized.
	stateStorageInitialized_ = true;
      }
    }
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Initialize this iteration object
  template <class ScalarType, class MV, class OP>
  void CGSingleRedIter<ScalarType,MV,OP>::initializeCG(CGIterationState<ScalarType,MV>& newstate)
  {
    // Initialize the state storage if it isn't already.
    if (!stateStorageInitialized_) 
      setStateSize();

    TEUCHOS_TEST_FOR_EXCEPTION(!stateStorageInitialized_,std::invalid_argument,
		       "Belos::CGSingleRedIter::initialize(): Cannot initialize state storage!");
    
    // NOTE:  In CGSingleRedIter R_, the initial residual, is required!!!  
    //
    std::string errstr("Belos::CGSingleRedIter::initialize(): Specified multivectors must have a consistent length and width.");

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
          Teuchos::RCP<MV> tmp = MVT::Clone( *Z_, 1 );
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
      MVT::Assign( *Z_, *P_ );

      // Multiply the current preconditioned residual vector by A and store in AZ_
      lp_->applyOp( *Z_, *AZ_ );

      // Logically, AP_ = AZ_
      MVT::Assign( *AZ_, *AP_ );
    }
    else {

      TEUCHOS_TEST_FOR_EXCEPTION(newstate.R == Teuchos::null,std::invalid_argument,
                         "Belos::CGSingleRedIter::initialize(): CGIterationState does not have initial residual.");
    }

    // The solver is initialized
    initialized_ = true;
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Iterate until the status test informs us we should stop.
  template <class ScalarType, class MV, class OP>
  void CGSingleRedIter<ScalarType,MV,OP>::iterate()
  {
    //
    // Allocate/initialize data structures
    //
    if (initialized_ == false) {
      initialize();
    }

    // Allocate memory for scalars.
    Teuchos::SerialDenseMatrix<int,ScalarType> sHz( 2, 1 );
    ScalarType rHz, rHz_old, alpha, beta, delta;

    // Create convenience variables for zero and one.
    const ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();
    const MagnitudeType zero = Teuchos::ScalarTraits<MagnitudeType>::zero();
    
    // Get the current solution vector.
    Teuchos::RCP<MV> cur_soln_vec = lp_->getCurrLHSVec();

    // Check that the current solution vector only has one column. 
    TEUCHOS_TEST_FOR_EXCEPTION( MVT::GetNumberVecs(*cur_soln_vec) != 1, CGIterateFailure,
                        "Belos::CGSingleRedIter::iterate(): current linear system has more than one vector!" );

    // Compute first <s,z> a.k.a. <r,z> and <Az,z> combined
    MVT::MvTransMv( one, *S_, *Z_, sHz );
    rHz = sHz(0,0);
    delta = sHz(1,0);
    alpha = rHz / delta;   

    // Check that alpha is a positive number!
    TEUCHOS_TEST_FOR_EXCEPTION( SCT::real(alpha) <= zero, CGIterateFailure,
      "Belos::CGSingleRedIter::iterate(): non-positive value for p^H*A*p encountered!" );
 
    ////////////////////////////////////////////////////////////////
    // Iterate until the status test tells us to stop.
    //
    while (1) {
      
      // Update the solution vector x := x + alpha * P_
      //
      MVT::MvAddMv( one, *cur_soln_vec, alpha, *P_, *cur_soln_vec );
      lp_->updateSolution();
      //
      // Compute the new residual R_ := R_ - alpha * AP_
      //
      MVT::MvAddMv( one, *R_, -alpha, *AP_, *R_ );
      //
      // Check the status test, now that the solution and residual have been updated
      //
      if (stest_->checkStatus(this) == Passed) {
        break; 
      }
      // Increment the iteration
      iter_++;
      //
      // Apply preconditioner to new residual to update Z_
      //
      if ( lp_->getLeftPrec() != Teuchos::null ) {
        lp_->applyLeftPrec( *R_, *Z_ );
        if ( lp_->getRightPrec() != Teuchos::null ) {
          Teuchos::RCP<MV> tmp = MVT::Clone( *Z_, 1 );
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
      // Multiply the current preconditioned residual vector by A and store in AZ_
      lp_->applyOp( *Z_, *AZ_ );
      //
      // Compute <S_,Z_> a.k.a. <R_,Z_> and <AZ_,Z_> combined
      MVT::MvTransMv( one, *S_, *Z_, sHz );
      //
      // Update scalars.
      rHz_old = rHz;
      rHz = sHz(0,0);
      delta = sHz(1,0);
      //
      beta = rHz / rHz_old;
      alpha = rHz / (delta - (beta*rHz / alpha));
      //
      // Check that alpha is a positive number!
      TEUCHOS_TEST_FOR_EXCEPTION( SCT::real(alpha) <= zero, CGIterateFailure,
        "Belos::CGSingleRedIter::iterate(): non-positive value for p^H*A*p encountered!" );
      //
      // Update the direction vector P_ := Z_ + beta * P_
      //
      MVT::MvAddMv( one, *Z_, beta, *P_, *P_ );
      // 
      // Update AP_ through recurrence relation AP_ := AZ_ + beta * AP_
      // NOTE: This increases the number of vector updates by 1, in exchange for
      //       reducing the collectives from 2 to 1.
      //
      MVT::MvAddMv( one, *AZ_, beta, *AP_, *AP_ );
      // 
    } // end while (sTest_->checkStatus(this) != Passed)
  }

} // end Belos namespace

#endif /* BELOS_CG_SINGLE_RED_ITER_HPP */
