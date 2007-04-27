// @HEADER
// ***********************************************************************
//
//                 Belos: Block Linear Solvers Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef BELOS_CG_ITER_HPP
#define BELOS_CG_ITER_HPP

/*! \file BelosCGIter.hpp
    \brief Belos concrete class for performing the conjugate-gradient (CG) iteration.
*/

#include "BelosConfigDefs.hpp"
#include "BelosTypes.hpp"
#include "BelosIteration.hpp"

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
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TimeMonitor.hpp"

/*!	
  \class Belos::CGIter
  
  \brief This class implements the preconditioned Conjugate Gradient (CG) iteration.
 
  \author Teri Barth and Heidi Thornquist
*/

namespace Belos {
  
  //! @name CGIter Structures 
  //@{ 
  
  /** \brief Structure to contain pointers to CGIter state variables.
   *
   * This struct is utilized by CGIter::initialize() and CGIter::getState().
   */
  template <class ScalarType, class MV>
  struct CGIterState {

    /*! \brief The current residual. */
    Teuchos::RefCountPtr<MV> r;

    /*! \brief The current preconditioned residual. */
    Teuchos::RefCountPtr<MV> z;

    /*! \brief The current decent direction vector */
    Teuchos::RefCountPtr<MV> p;

    /*! \brief The matrix A applied to current decent direction vector */
    Teuchos::RefCountPtr<MV> Ap;
    
    CGIterState() : r(Teuchos::null), z(Teuchos::null), 
		    p(Teuchos::null), Ap(Teuchos::null)
    {}
  };
  
  //! @name CGIter Exceptions
  //@{ 
  
  /** \brief CGIterInitFailure is thrown when the CGIter object is unable to
   * generate an initial iterate in the CGIter::initialize() routine. 
   *
   * This exception is thrown from the CGIter::initialize() method, which is
   * called by the user or from the CGIter::iterate() method if isInitialized()
   * == \c false.
   *
   * In the case that this exception is thrown, 
   * CGIter::isInitialized() will be \c false and the user will need to provide
   * a new initial iterate to the iteration.
   */
  class CGIterInitFailure : public BelosError {public:
    CGIterInitFailure(const std::string& what_arg) : BelosError(what_arg)
    {}};

  /** \brief CGIterateFailure is thrown when the CGIter object is unable to
   * compute the next iterate in the CGIter::iterate() routine. 
   *
   * This exception is thrown from the CGIter::iterate() method.
   *
   */
  class CGIterateFailure : public BelosError {public:
    CGIterateFailure(const std::string& what_arg) : BelosError(what_arg)
    {}};
  
  //@}


template<class ScalarType, class MV, class OP>
class CGIter : virtual public Iteration<ScalarType,MV,OP> {

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
  CGIter( const Teuchos::RefCountPtr<LinearProblem<ScalarType,MV,OP> > &problem, 
		  const Teuchos::RefCountPtr<OutputManager<ScalarType> > &printer,
		  const Teuchos::RefCountPtr<StatusTest<ScalarType,MV,OP> > &tester,
		  Teuchos::ParameterList &params );

  //! Destructor.
  virtual ~CGIter() {};
  //@}


  //! @name Solver methods
  //@{ 
  
  /*! \brief This method performs CG iterations until the status
   * test indicates the need to stop or an error occurs (in which case, an
   * exception is thrown).
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
  void initialize(CGIterState<ScalarType,MV> newstate);

  /*! \brief Initialize the solver with the initial vectors from the linear problem
   *  or random data.
   */
  void initialize()
  {
    CGIterState<ScalarType,MV> empty;
    initialize(empty);
  }
  
  /*! \brief Get the current state of the linear solver.
   *
   * The data is only valid if isInitialized() == \c true.
   *
   * \returns A CGIterState object containing const pointers to the current solver state.
   */
  CGIterState<ScalarType,MV> getState() const {
    CGIterState<ScalarType,MV> state;
    state.r = r_;
    state.p = p_;
    state.Ap = Ap_;
    state.z = z_;
    return state;
  }

  //@}

  
  //! @name Status methods
  //@{ 

  //! \brief Get the current iteration count.
  int getNumIters() const { return iter_; }
  
  //! \brief Reset the iteration count.
  void resetNumIters() { iter_ = 0; }

  //! Get the norms of the residuals native to the solver.
  //! \return A vector of length blockSize containing the native residuals.
  Teuchos::RefCountPtr<const MV> getNativeResiduals( std::vector<MagnitudeType> *norms ) const;

  //! Get the current update to the linear system.
  /*! \note This method returns a null pointer because the linear problem is current.
  */
  Teuchos::RefCountPtr<MV> getCurrentUpdate() const { return Teuchos::null; }

  //@}
  
  //! @name Accessor methods
  //@{ 

  //! Get a constant reference to the linear problem.
  const LinearProblem<ScalarType,MV,OP>& getProblem() const { return *lp_; }

  //! States whether the solver has been initialized or not.
  bool isInitialized() { return initialized_; }

  //@}

  private:

  //
  // Internal methods
  //
  //! Method for initalizing the state storage needed by block GMRES.
  void setStateSize();
  
  //
  // Classes inputed through constructor that define the linear problem to be solved.
  //
  const Teuchos::RefCountPtr<LinearProblem<ScalarType,MV,OP> >    lp_;
  const Teuchos::RefCountPtr<OutputManager<ScalarType> >          om_;
  const Teuchos::RefCountPtr<StatusTest<ScalarType,MV,OP> >       stest_;

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
  Teuchos::RefCountPtr<MV> r_;
  //
  // Preconditioned residual
  Teuchos::RefCountPtr<MV> z_;
  //
  // Direction vector
  Teuchos::RefCountPtr<MV> p_;
  //
  // Operator applied to direction vector
  Teuchos::RefCountPtr<MV> Ap_;

};

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor.
  template<class ScalarType, class MV, class OP>
  CGIter<ScalarType,MV,OP>::CGIter(const Teuchos::RefCountPtr<LinearProblem<ScalarType,MV,OP> > &problem, 
						   const Teuchos::RefCountPtr<OutputManager<ScalarType> > &printer,
						   const Teuchos::RefCountPtr<StatusTest<ScalarType,MV,OP> > &tester,
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
  void CGIter<ScalarType,MV,OP>::setStateSize ()
  {
    if (!stateStorageInitialized_) {

      // Check if there is any multivector to clone from.
      Teuchos::RefCountPtr<const MV> lhsMV = lp_->getLHS();
      Teuchos::RefCountPtr<const MV> rhsMV = lp_->getRHS();
      if (lhsMV == Teuchos::null && rhsMV == Teuchos::null) {
	stateStorageInitialized_ = false;
	return;
      }
      else {
	
	// Initialize the state storage
	// If the subspace has not be initialized before, generate it using the LHS or RHS from lp_.
	if (r_ == Teuchos::null) {
	  // Get the multivector that is not null.
	  Teuchos::RefCountPtr<const MV> tmp = ( (rhsMV!=Teuchos::null)? rhsMV: lhsMV );
	  TEST_FOR_EXCEPTION(tmp == Teuchos::null,std::invalid_argument,
			     "Belos::CGIter::setStateSize(): linear problem does not specify multivectors to clone from.");
	  r_ = MVT::Clone( *tmp, 1 );
	  z_ = MVT::Clone( *tmp, 1 );
	  p_ = MVT::Clone( *tmp, 1 );
	  Ap_ = MVT::Clone( *tmp, 1 );
	}
	
	// State storage has now been initialized.
	stateStorageInitialized_ = true;
      }
    }
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Get the native residuals stored in this iteration.  
  // Note:  No residual vector will be returned by Gmres.
  template <class ScalarType, class MV, class OP>
  Teuchos::RefCountPtr<const MV> CGIter<ScalarType,MV,OP>::getNativeResiduals( std::vector<MagnitudeType> *norms ) const 
  {
    std::vector<int> index( 1 );
    index[0] = 0;
    RefCountPtr<MV> curR = MVT::CloneView( *r_, index );
    return curR;
  }
  

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Initialize this iteration object
  template <class ScalarType, class MV, class OP>
  void CGIter<ScalarType,MV,OP>::initialize(CGIterState<ScalarType,MV> newstate)
  {
    // Initialize the state storage if it isn't already.
    if (!stateStorageInitialized_) 
      setStateSize();

    TEST_FOR_EXCEPTION(!stateStorageInitialized_,std::invalid_argument,
		       "Belos::CGIter::initialize(): Cannot initialize state storage!");
    
    // NOTE:  In CGIter r_, the initial residual, is required!!!  
    //
    std::string errstr("Belos::CGIter::initialize(): Specified multivectors must have a consistent length and width.");

    // Create convenience variables for zero and one.
    const ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();
    const MagnitudeType zero = Teuchos::ScalarTraits<MagnitudeType>::zero();

    if (newstate.r != Teuchos::null) {

      TEST_FOR_EXCEPTION( MVT::GetVecLength(*newstate.r) != MVT::GetVecLength(*r_),
                          std::invalid_argument, errstr );
      TEST_FOR_EXCEPTION( MVT::GetNumberVecs(*newstate.r) != 1,
                          std::invalid_argument, errstr );

      // Copy basis vectors from newstate into V
      if (newstate.r != r_) {
        // copy over the initial residual (unpreconditioned).
	MVT::MvAddMv( one, *newstate.r, zero, *newstate.r *r_ );
      }

      // Compute initial direction vectors
      // Initially, they are set to the preconditioned residuals
      //
      if ( lp_->getLeftPrec() != Teuchos::null ) {
	lp_->applyLeftPrec( *r_, *z_ ); 
	MVT::MvAddMv( one, *z_, zero, *z_, *p_ );
      } else {
	z_ = r_;
	MVT::MvAddMv( one, *r_, zero, *r_, *p_ );
      }
    }
    else {

      TEST_FOR_EXCEPTION(newstate.r == Teuchos::null,std::invalid_argument,
                         "Belos::CGIter::initialize(): CGStateIterState does not have initial residual.");
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
    Teuchos::SerialDenseMatrix<int,ScalarType> alpha( 1, 1 );
    Teuchos::SerialDenseMatrix<int,ScalarType> beta( 1, 1 );
    Teuchos::SerialDenseMatrix<int,ScalarType> rHz( 1, 1 ), rHz_old( 1, 1 ), pAp( 1, 1 );

    // Create convenience variables for zero and one.
    const ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();
    const MagnitudeType zero = Teuchos::ScalarTraits<MagnitudeType>::zero();
    
    // Get the current solution vector.
    Teuchos::RefCountPtr<MV> cur_soln_vec = lp_->getCurrLHSVec();

    // Compute first <r,z> a.k.a. rHz
    MVT::MvTransMv( one, *r_, *z_, rHz );
    
    ////////////////////////////////////////////////////////////////
    // Iterate until the status test tells us to stop.
    //
    while (stest_->checkStatus(this) != Passed) {
      
      // Increment the iteration
      iter_++;

      // Multiply the current direction vector by A and store in Ap_
      lp_->applyOp( *_p, *Ap_ );
      
      // Compute alpha := <r_,z_> / <p_,Ap_>
      MVT::MvTransMv( one, *p_, *Ap_, pAp );
      alpha(0,0) = rHz(0,0) / pAp(0,0);
      
      // Check that alpha is a positive number!
      TEST_FOR_EXCEPTION( SCT::real(alpha(0,0)) <= zero, CGIterateFailure,
			  "Belos::CGIter::iterate(): non-positive value for p^H*A*p encountered!" );
      //
      // Update the solution vector x := x + alpha * p_
      //
      MVT::MvAddMv( one, *cur_soln_vec, alpha(0,0), *p_, *cur_soln_vec );
      lp_->updateSolution();
      //
      // Save the denominator of beta before residual is updated [ old <r_, z_> ]
      //
      rHz_old(0,0) = rHz(0,0);
      //
      // Compute the new residual r_ := r_ - alpha * Ap_
      //
      MVT::MvAddMv( one, *r_, -alpha(0,0), *Ap_, *r_ );
      //
      // Compute beta := [ new <r_, z_> ] / [ old <r_, z_> ], 
      // and the new direction vector p.
      //
      if ( lp_->getLeftPrec() != Teuchos::null ) {
	lp_->applyLeftPrec( *r_, *z_ );
      } else {
	z_ = r_;
      }
      //
      MVT::MvTransMv( one, *r_, *z_, rHz );
      //
      beta(0,0) = rHz(0,0) / rHz_old(0,0);
      //
      MVT::MvAddMv( one, *z_, beta(0,0), *p_, *p_ );
      
    } // end while (sTest_->checkStatus(this) != Passed)
  }

} // end Belos namespace

#endif /* BELOS_CG_ITER_HPP */
