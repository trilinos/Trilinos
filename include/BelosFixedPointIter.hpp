// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef BELOS_FIXEDPOINT_ITER_HPP
#define BELOS_FIXEDPOINT_ITER_HPP

/*! \file BelosFixedPointIter.hpp
    \brief Belos concrete class for performing fixed point iteration iteration.
*/

#include "BelosConfigDefs.hpp"
#include "BelosTypes.hpp"
#include "BelosFixedPointIteration.hpp"

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
  \class Belos::FixedPointIter

  \brief This class implements the preconditioned fixed point iteration.

  \ingroup belos_solver_framework

*/

namespace Belos {

template<class ScalarType, class MV, class OP>
class FixedPointIter : virtual public FixedPointIteration<ScalarType,MV,OP> {

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

  /*! \brief %FixedPointIter constructor with linear problem, solver utilities, and parameter list of solver options.
   *
   * This constructor takes pointers required by the linear solver iteration, in addition
   * to a parameter list of options for the linear solver.
   */
  FixedPointIter( const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem,
                  const Teuchos::RCP<OutputManager<ScalarType> > &printer,
                  const Teuchos::RCP<StatusTest<ScalarType,MV,OP> > &tester,
                  Teuchos::ParameterList &params );

  //! Destructor.
  virtual ~FixedPointIter() {};
  //@}


  //! @name Solver methods
  //@{

  /*! \brief This method performs Fixed Point iterations until the status
   * test indicates the need to stop or an error occurs (in which case, an
   * std::exception is thrown).
   *
   * iterate() will first determine whether the solver is initialized; if
   * not, it will call initialize() using default arguments. After
   * initialization, the solver performs FixedPoint iterations until the
   * status test evaluates as ::Passed, at which point the method returns to
   * the caller.
   *
   * The status test is queried at the beginning of the iteration.
   */
  void iterate();

  /*! \brief Initialize the solver to an iterate, providing a complete state.
   *
   * The %FixedPointIter contains a certain amount of state, consisting of the current
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
  void initializeFixedPoint(FixedPointIterationState<ScalarType,MV>& newstate);

  /*! \brief Initialize the solver with the initial vectors from the linear problem
   *  or random data.
   */
  void initialize()
  {
    FixedPointIterationState<ScalarType,MV> empty;
    initializeFixedPoint(empty);
  }

  /*! \brief Get the current state of the linear solver.
   *
   * The data is only valid if isInitialized() == \c true.
   *
   * \returns A FixedPointIterationState object containing const pointers to the current solver state.
   */
  FixedPointIterationState<ScalarType,MV> getState() const {
    FixedPointIterationState<ScalarType,MV> state;
    state.R = R_;
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
  Teuchos::RCP<const MV> getNativeResiduals( std::vector<MagnitudeType> * /* norms */ ) const { return R_; }

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
  int getBlockSize() const { return numRHS_; }

  //! \brief Set the blocksize to be used by the iterative solver in solving this linear problem.
  void setBlockSize(int blockSize);

  //! States whether the solver has been initialized or not.
  bool isInitialized() { return initialized_; }

  //@}

  private:

  //
  // Internal methods
  //
  //! Method for initalizing the state storage needed by FixedPoint.
  void setStateSize();

  //
  // Classes inputed through constructor that define the linear problem to be solved.
  //
  const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> >    lp_;
  const Teuchos::RCP<OutputManager<ScalarType> >          om_;
  const Teuchos::RCP<StatusTest<ScalarType,MV,OP> >       stest_;

  // Algorithmic parameters
  //
  // blockSize_ is the solver block size.
  int numRHS_;

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

};

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor.
  template<class ScalarType, class MV, class OP>
  FixedPointIter<ScalarType,MV,OP>::FixedPointIter(const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem,
                                                   const Teuchos::RCP<OutputManager<ScalarType> > &printer,
                                                   const Teuchos::RCP<StatusTest<ScalarType,MV,OP> > &tester,
                                                   Teuchos::ParameterList &params ):
    lp_(problem),
    om_(printer),
    stest_(tester),
    numRHS_(0),
    initialized_(false),
    stateStorageInitialized_(false),
    iter_(0)
  {
    setBlockSize(params.get("Block Size",MVT::GetNumberVecs(*problem->getCurrRHSVec())));
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Setup the state storage.
  template <class ScalarType, class MV, class OP>
  void FixedPointIter<ScalarType,MV,OP>::setStateSize ()
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
                             "Belos::FixedPointIter::setStateSize(): linear problem does not specify multivectors to clone from.");
          R_ = MVT::Clone( *tmp, numRHS_ );
          Z_ = MVT::Clone( *tmp, numRHS_ );
        }

        // State storage has now been initialized.
        stateStorageInitialized_ = true;
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Set the block size and make necessary adjustments.
  template <class ScalarType, class MV, class OP>
  void FixedPointIter<ScalarType,MV,OP>::setBlockSize(int blockSize)
  {
    // This routine only allocates space; it doesn't not perform any computation
    // any change in size will invalidate the state of the solver.

    TEUCHOS_TEST_FOR_EXCEPTION(blockSize != MVT::GetNumberVecs(*lp_->getCurrRHSVec()), std::invalid_argument, "Belos::FixedPointIter::setBlockSize size must match # RHS.");

    TEUCHOS_TEST_FOR_EXCEPTION(blockSize <= 0, std::invalid_argument, "Belos::FixedPointIter::setBlockSize was passed a non-positive argument.");
    if (blockSize == numRHS_) {
      // do nothing
      return;
    }

    if (blockSize!=numRHS_)
      stateStorageInitialized_ = false;

    numRHS_ = blockSize;

    initialized_ = false;

    // Use the current blockSize_ to initialize the state storage.
    setStateSize();

  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Initialize this iteration object
  template <class ScalarType, class MV, class OP>
  void FixedPointIter<ScalarType,MV,OP>::initializeFixedPoint(FixedPointIterationState<ScalarType,MV>& newstate)
  {
    // Initialize the state storage if it isn't already.
    if (!stateStorageInitialized_)
      setStateSize();

    TEUCHOS_TEST_FOR_EXCEPTION(!stateStorageInitialized_,std::invalid_argument,
                       "Belos::FixedPointIter::initialize(): Cannot initialize state storage!");

    // NOTE:  In FixedPointIter R_, the initial residual, is required!!!
    //
    std::string errstr("Belos::FixedPointIter::initialize(): Specified multivectors must have a consistent length and width.");

    // Create convenience variables for zero and one.
    const ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();
    const ScalarType zero = Teuchos::ScalarTraits<ScalarType>::zero();

    if (newstate.R != Teuchos::null) {
      TEUCHOS_TEST_FOR_EXCEPTION( MVT::GetNumberVecs(*R_) != MVT::GetNumberVecs(*newstate.R),
                                  std::invalid_argument, errstr );

      TEUCHOS_TEST_FOR_EXCEPTION( MVT::GetGlobalLength(*newstate.R) != MVT::GetGlobalLength(*R_),
                          std::invalid_argument, errstr );
      TEUCHOS_TEST_FOR_EXCEPTION( MVT::GetNumberVecs(*newstate.R) != numRHS_,
                          std::invalid_argument, errstr );

      // Copy basis vectors from newstate into V
      if (newstate.R != R_) {
        // copy over the initial residual (unpreconditioned).
        MVT::MvAddMv( one, *newstate.R, zero, *newstate.R, *R_ );
      }

    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(newstate.R == Teuchos::null,std::invalid_argument,
                         "Belos::FixedPointIter::initialize(): FixedPointIterationState does not have initial residual.");
    }

    // The solver is initialized
    initialized_ = true;
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Iterate until the status test informs us we should stop.
  template <class ScalarType, class MV, class OP>
  void FixedPointIter<ScalarType,MV,OP>::iterate()
  {
    //
    // Allocate/initialize data structures
    //
    if (initialized_ == false) {
      initialize();
    }

    // Create convenience variables for zero and one.
    const ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();
    const MagnitudeType zero = Teuchos::ScalarTraits<MagnitudeType>::zero(); // unused

    // Get the current solution vector.
    Teuchos::RCP<MV> cur_soln_vec = lp_->getCurrLHSVec();

    // Temp vector
    Teuchos::RCP<MV> tmp = MVT::Clone( *R_, numRHS_ );

    if (lp_->getRightPrec() != Teuchos::null) {
      // Set rhs to initial residual
      Teuchos::RCP<MV> rhs = MVT::CloneCopy( *R_ );

      // Zero initial guess
      MVT::MvInit( *Z_, zero );

      ////////////////////////////////////////////////////////////////
      // Iterate until the status test tells us to stop.
      //
      while (stest_->checkStatus(this) != Passed) {

        // Increment the iteration
        iter_++;

        // Apply preconditioner
        lp_->applyRightPrec( *R_, *tmp );

        // Update solution vector
        MVT::MvAddMv( one, *cur_soln_vec, one, *tmp, *cur_soln_vec );
        lp_->updateSolution();

        // Update solution vector
        MVT::MvAddMv( one, *Z_, one, *tmp, *Z_ );

        // Compute new residual
        lp_->applyOp (*Z_, *tmp );
        MVT::MvAddMv( one, *rhs, -one, *tmp, *R_ );

      } // end while (sTest_->checkStatus(this) != Passed)

    } else {
      Teuchos::RCP<const MV> rhs = lp_->getCurrRHSVec();

      ////////////////////////////////////////////////////////////////
      // Iterate until the status test tells us to stop.
      //
      while (stest_->checkStatus(this) != Passed) {

        // Increment the iteration
        iter_++;

        // Compute initial preconditioned residual
        if ( lp_->getLeftPrec() != Teuchos::null ) {
          lp_->applyLeftPrec( *R_, *Z_ );
        }
        else {
          Z_ = R_;
        }

        // Update solution vector
        MVT::MvAddMv(one,*cur_soln_vec,one,*Z_,*cur_soln_vec);
        lp_->updateSolution();

        // Compute new residual
        lp_->applyOp(*cur_soln_vec,*tmp);
        MVT::MvAddMv(one,*rhs,-one,*tmp,*R_);

      } // end while (sTest_->checkStatus(this) != Passed)
    }
  }

} // end Belos namespace

#endif /* BELOS_FIXEDPOINT_ITER_HPP */
