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

#ifndef BELOS_PSEUDO_BLOCK_CG_ITER_HPP
#define BELOS_PSEUDO_BLOCK_CG_ITER_HPP

/*! \file BelosPseudoBlockCGIter.hpp
    \brief Belos concrete class for performing the pseudo-block CG iteration.
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
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TimeMonitor.hpp"

/*!	
  \class Belos::PseudoBlockCGIter
  
  \brief This class implements the pseudo-block CG iteration, where the basic CG
  algorithm is performed on all of the linear systems simultaneously.  

  \ingroup belos_solver_framework
 
  \author Heidi Thornquist
*/

namespace Belos {
  
  template<class ScalarType, class MV, class OP>
  class PseudoBlockCGIter : virtual public CGIteration<ScalarType,MV,OP> {
    
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

    // Current number of iterations performed.
    bool assertPositiveDefiniteness_;

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

  };
  
  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor.
  template<class ScalarType, class MV, class OP>
  PseudoBlockCGIter<ScalarType,MV,OP>::PseudoBlockCGIter(const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem, 
							       const Teuchos::RCP<OutputManager<ScalarType> > &printer,
							       const Teuchos::RCP<StatusTest<ScalarType,MV,OP> > &tester,
							       Teuchos::ParameterList &params ):
    lp_(problem),
    om_(printer),
    stest_(tester),
    numRHS_(0),
    initialized_(false),
    iter_(0),
    assertPositiveDefiniteness_( params.get("Assert Positive Definiteness", true) )
  {
  }
  

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Initialize this iteration object
  template <class ScalarType, class MV, class OP>
  void PseudoBlockCGIter<ScalarType,MV,OP>::initializeCG(CGIterationState<ScalarType,MV> newstate)
  {
    // Check if there is any multivector to clone from.
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
	
    // NOTE:  In CGIter R_, the initial residual, is required!!!  
    //
    std::string errstr("Belos::BlockPseudoCGIter::initialize(): Specified multivectors must have a consistent length and width.");

    // Create convenience variables for zero and one.
    const ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();
    const MagnitudeType zero = Teuchos::ScalarTraits<MagnitudeType>::zero();

    if (!Teuchos::is_null(newstate.R)) {

      TEUCHOS_TEST_FOR_EXCEPTION( MVText::GetGlobalLength(*newstate.R) != MVText::GetGlobalLength(*R_),
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
  template <class ScalarType, class MV, class OP>
  void PseudoBlockCGIter<ScalarType,MV,OP>::iterate()
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

        alpha(i,i) = rHz[i] / pAp[i];
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
        beta(i,i) = rHz[i] / rHz_old[i];
	index[0] = i;
	Teuchos::RCP<const MV> Z_i = MVT::CloneView( *Z_, index );
	Teuchos::RCP<MV> P_i = MVT::CloneViewNonConst( *P_, index );
        MVT::MvAddMv( one, *Z_i, beta(i,i), *P_i, *P_i );
      }
      //      
    } // end while (sTest_->checkStatus(this) != Passed)
  }

} // end Belos namespace

#endif /* BELOS_PSEUDO_BLOCK_CG_ITER_HPP */
