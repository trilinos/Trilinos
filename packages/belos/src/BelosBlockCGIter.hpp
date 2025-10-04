// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
#include "BelosDenseMatTraits.hpp"

#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TimeMonitor.hpp"

namespace Belos {

//! @name BlockCGIteration Structures
  //@{

  /** \brief Structure to contain pointers to BlockCGIteration state variables.
   *
   * This struct is utilized by BlockCGIteration::initialize() and BlockCGIteration::getState().
   */
  template <class ScalarType, class MV>
  class BlockCGIterationState : public CGIterationStateBase<ScalarType, MV> {

  public:
    BlockCGIterationState() = default;

    BlockCGIterationState(Teuchos::RCP<const MV> tmp) {
      initialize(tmp);
    }

    virtual ~BlockCGIterationState() = default;

    void initialize(Teuchos::RCP<const MV> tmp, int _numVectors) {
      using MVT = MultiVecTraits<ScalarType, MV>;
      this->R = MVT::Clone( *tmp, _numVectors );
      this->Z = MVT::Clone( *tmp, _numVectors );
      this->P = MVT::Clone( *tmp, _numVectors );
      this->AP = MVT::Clone(*tmp, _numVectors );

      CGIterationStateBase<ScalarType, MV>::initialize(tmp, _numVectors);
    }

    bool matches(Teuchos::RCP<const MV> tmp, int _numVectors=1) const {
      return CGIterationStateBase<ScalarType, MV>::matches(tmp, _numVectors);
    }

};

/// \class BlockCGIter
/// \brief Implementation of the block preconditioned Conjugate
///   Gradient (CG) iteration.
/// \ingroup belos_solver_framework
/// \author Teri Barth and Heidi Thornquist

/// \brief Stub implementation of BlockCGIter, for ScalarType types
///   for which Teuchos::LAPACK does NOT have a valid implementation.
template<class ScalarType, class MV, class OP, class DM,
         const bool lapackSupportsScalarType =
         Belos::Details::LapackSupportsScalar<ScalarType>::value>
class BlockCGIter : virtual public CGIteration<ScalarType, MV, OP, DM> {
public:
  typedef MultiVecTraits<ScalarType,MV, DM> MVT;
  typedef OperatorTraits<ScalarType,MV,OP> OPT;
  typedef Teuchos::ScalarTraits<ScalarType> SCT;
  typedef typename SCT::magnitudeType MagnitudeType;

  BlockCGIter( const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > & /* problem */,
               const Teuchos::RCP<OutputManager<ScalarType> > & /* printer */,
               const Teuchos::RCP<StatusTest<ScalarType,MV,OP,DM> > & /* tester */,
               const Teuchos::RCP<MatOrthoManager<ScalarType,MV,OP,DM> > & /* ortho */,
               Teuchos::ParameterList & /* params */ )
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Stub");
  }

  virtual ~BlockCGIter() {}

  void iterate () {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Stub");
  }

  void initializeCG (Teuchos::RCP<BlockCGIterationState<ScalarType,MV> > /* newstate */, Teuchos::RCP<MV> /* R_0 */) {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Stub");
  }

  void initialize () {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Stub");
  }

  Teuchos::RCP<CGIterationStateBase<ScalarType,MV> > getState () const {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Stub");
  }

  void setState(Teuchos::RCP<CGIterationStateBase<ScalarType,MV> > state) {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Stub");
  }

  int getNumIters() const {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Stub");
  }

  void resetNumIters( int iter=0 ) {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Stub");
  }

  Teuchos::RCP<const MV>
  getNativeResiduals (std::vector<MagnitudeType>* /* norms */) const {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Stub");
  }

  Teuchos::RCP<MV> getCurrentUpdate() const {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Stub");
  }

  const LinearProblem<ScalarType,MV,OP>& getProblem() const {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Stub");
  }

  int getBlockSize() const {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Stub");
  }

  void setBlockSize(int blockSize) {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Stub");
  }

  bool isInitialized() {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Stub");
  }

  void setDoCondEst(bool val){
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Stub");
  }

};

/// \brief Partial specialization for ScalarType types for which
///   Teuchos::LAPACK has a valid implementation.
///
/// This is the (non-stub) actual implementation of BlockCGIter.
template<class ScalarType, class MV, class OP, class DM>
class BlockCGIter<ScalarType, MV, OP, DM, true> :
    virtual public CGIteration<ScalarType,MV,OP,DM>
{
public:
  //
  // Convenience typedefs
  //
  using MVT = MultiVecTraits<ScalarType, MV, DM>;
  using OPT = OperatorTraits<ScalarType, MV, OP>;
  using SCT = Teuchos::ScalarTraits<ScalarType>;
  using MagnitudeType = typename SCT::magnitudeType;

  //! @name Constructors/Destructor
  //@{

  /*! \brief %BlockCGIter constructor with linear problem, solver utilities, and parameter list of solver options.
   *
   * This constructor takes pointers required by the linear solver iteration, in addition
   * to a parameter list of options for the linear solver.
   */
  BlockCGIter( const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem,
               const Teuchos::RCP<OutputManager<ScalarType> > &printer,
               const Teuchos::RCP<StatusTest<ScalarType,MV,OP,DM> > &tester,
               const Teuchos::RCP<MatOrthoManager<ScalarType,MV,OP,DM> > &ortho,
               Teuchos::ParameterList &params );

  //! Destructor.
  virtual ~BlockCGIter() = default;
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
  void initializeCG(Teuchos::RCP<CGIterationStateBase<ScalarType,MV> > newstate, Teuchos::RCP<MV> R_0);

  /*! \brief Initialize the solver with the initial vectors from the linear problem
   *  or random data.
   */
  void initialize()
  {
    initializeCG(Teuchos::null, Teuchos::null);
  }

  /*! \brief Get the current state of the linear solver.
   *
   * The data is only valid if isInitialized() == \c true.
   *
   * \returns A BlockCGIterationState object containing const pointers to the current solver state.
   */
  Teuchos::RCP<CGIterationStateBase<ScalarType,MV> > getState() const {
    auto state = Teuchos::rcp(new BlockCGIterationState<ScalarType,MV>());
    state->R = R_;
    state->P = P_;
    state->AP = AP_;
    state->Z = Z_;
    return state;
  }

  void setState(Teuchos::RCP<CGIterationStateBase<ScalarType,MV> > state) {
    auto s = Teuchos::rcp_dynamic_cast<BlockCGIterationState<ScalarType,MV> >(state, true);
    R_ = s->R;
    Z_ = s->Z;
    P_ = s->P;
    AP_ = s->AP;
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

  //! Get the block size to be used by the iterative solver in solving this linear problem.
  int getBlockSize() const { return blockSize_; }

  //! \brief Set the block size to be used by the iterative solver in solving this linear problem.
  void setBlockSize(int blockSize);

  //! States whether the solver has been initialized or not.
  bool isInitialized() { return initialized_; }

  //! Sets whether or not to store the diagonal for condition estimation
  void setDoCondEst(bool /* val */){/*ignored*/}

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
  // Classes inputed through constructor that define the linear problem to be solved.
  //
  const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> >    lp_;
  const Teuchos::RCP<OutputManager<ScalarType> >          om_;
  const Teuchos::RCP<StatusTest<ScalarType,MV,OP,DM> >       stest_;
  const Teuchos::RCP<OrthoManager<ScalarType,MV,DM> >        ortho_;

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

  template<class ScalarType, class MV, class OP, class DM>
  BlockCGIter<ScalarType,MV,OP,DM,true>::
  BlockCGIter (const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> >& problem,
               const Teuchos::RCP<OutputManager<ScalarType> >& printer,
               const Teuchos::RCP<StatusTest<ScalarType,MV,OP,DM> >& tester,
               const Teuchos::RCP<MatOrthoManager<ScalarType,MV,OP,DM> >& ortho,
               Teuchos::ParameterList& params) :
    lp_(problem),
    om_(printer),
    stest_(tester),
    ortho_(ortho),
    blockSize_(0),
    initialized_(false),
    iter_(0)
  {
    // Set the block size and allocate data
    int bs = params.get("Block Size", 1);
    setBlockSize( bs );
  }

  template<class ScalarType, class MV, class OP, class DM>
  void BlockCGIter<ScalarType,MV,OP,DM,true>::setBlockSize (int blockSize)
  {
    // This routine only allocates space; it doesn't not perform any computation
    // any change in size will invalidate the state of the solver.
    TEUCHOS_TEST_FOR_EXCEPTION
      (blockSize <= 0, std::invalid_argument, "Belos::BlockGmresIter::"
       "setBlockSize: blockSize = " << blockSize << " <= 0.");
    if (blockSize == blockSize_) {
      return; // do nothing
    }
    blockSize_ = blockSize;
    initialized_ = false;
  }

  template <class ScalarType, class MV, class OP, class DM>
  void BlockCGIter<ScalarType,MV,OP,DM,true>::
  initializeCG (Teuchos::RCP<CGIterationStateBase<ScalarType,MV> > newstate, Teuchos::RCP<MV> R_0)
  {
    const char prefix[] = "Belos::BlockCGIter::initialize: ";

    // Initialize the state storage if it isn't already.
    Teuchos::RCP<const MV> lhsMV = lp_->getLHS();
    Teuchos::RCP<const MV> rhsMV = lp_->getRHS();
    Teuchos::RCP<const MV> tmp = ( (rhsMV!=Teuchos::null)? rhsMV: lhsMV );
    TEUCHOS_ASSERT(!newstate.is_null());
    if (!Teuchos::rcp_dynamic_cast<BlockCGIterationState<ScalarType,MV> >(newstate, true)->matches(tmp, blockSize_))
      newstate->initialize(tmp, blockSize_);
    setState(newstate);

    // NOTE:  In BlockCGIter R_, the initial residual, is required!!!
    const char errstr[] = "Specified multivectors must have a consistent "
      "length and width.";

    {

      TEUCHOS_TEST_FOR_EXCEPTION
        (MVT::GetGlobalLength(*R_0) != MVT::GetGlobalLength(*R_),
         std::invalid_argument, prefix << errstr );
      TEUCHOS_TEST_FOR_EXCEPTION
        (MVT::GetNumberVecs(*R_0) != blockSize_,
         std::invalid_argument, prefix << errstr );

      // Copy basis vectors from newstate into V
      if (R_0 != R_) {
        // copy over the initial residual (unpreconditioned).
        MVT::Assign( *R_0, *R_ );
      }
      // Compute initial direction vectors
      // Initially, they are set to the preconditioned residuals
      //
      if ( lp_->getLeftPrec() != Teuchos::null ) {
        lp_->applyLeftPrec( *R_, *Z_ );
        if ( lp_->getRightPrec() != Teuchos::null ) {
          Teuchos::RCP<MV> tmp2 = MVT::Clone( *Z_, blockSize_ );
          lp_->applyRightPrec( *Z_, *tmp2 );
          Z_ = tmp2;
        }
      }
      else if ( lp_->getRightPrec() != Teuchos::null ) {
        lp_->applyRightPrec( *R_, *Z_ );
      }
      else {
        Z_ = R_;
      }
      MVT::Assign( *Z_, *P_ );
    }

    // The solver is initialized
    initialized_ = true;
  }

  template<class ScalarType, class MV, class OP, class DM>
  void BlockCGIter<ScalarType,MV,OP,DM,true>::iterate()
  {
    const char prefix[] = "Belos::BlockCGIter::iterate: ";
    typedef DenseMatTraits<ScalarType,DM>    DMT;

    //
    // Allocate/initialize data structures
    //
    if (initialized_ == false) {
      initialize();
    }

    // Allocate memory for scalars.
    Teuchos::RCP<DM> alpha = DMT::Create( blockSize_, blockSize_ );
    Teuchos::RCP<DM> beta = DMT::Create( blockSize_, blockSize_ );
    Teuchos::RCP<DM> pAp = DMT::Create( blockSize_, blockSize_ );

    // Create dense spd solver.
    Teuchos::RCP<DenseSolver<ScalarType,DM>> lltSolver = DMT::createDenseSolver();
    lltSolver->setSPD( true );

    // Create convenience variable for one.
    const ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();

    // Get the current solution std::vector.
    Teuchos::RCP<MV> cur_soln_vec = lp_->getCurrLHSVec();

    // Check that the current solution std::vector has blockSize_ columns.
    TEUCHOS_TEST_FOR_EXCEPTION
      (MVT::GetNumberVecs(*cur_soln_vec) != blockSize_, CGIterateFailure,
       prefix << "Current linear system does not have the right number of vectors!" );
    int rank = ortho_->normalize( *P_, Teuchos::null );
    TEUCHOS_TEST_FOR_EXCEPTION
      (rank != blockSize_, CGIterationOrthoFailure,
       prefix << "Failed to compute initial block of orthonormal direction vectors.");

    //
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
      MVT::MvTransMv( one, *P_, *R_, *alpha );
      MVT::MvTransMv( one, *P_, *AP_, *pAp );

      // Compute Cholesky factorization of pAp
      lltSolver->setMatrix( pAp );
      lltSolver->factorWithEquilibration( true );
      int info = lltSolver->factor();
      TEUCHOS_TEST_FOR_EXCEPTION
        (info != 0, CGIterationLAPACKFailure,
         prefix << "Failed to compute Cholesky factorization using LAPACK routine POTRF.");

      // Compute alpha by performing a back and forward solve with the
      // Cholesky factorization in pAp.
      lltSolver->setVectors (alpha, alpha);
      info = lltSolver->solve();
      TEUCHOS_TEST_FOR_EXCEPTION
        (info != 0, CGIterationLAPACKFailure,
         prefix << "Failed to compute alpha using Cholesky factorization (POTRS).");

      // Update the solution std::vector X := X + alpha * P_
      MVT::MvTimesMatAddMv( one, *P_, *alpha, one, *cur_soln_vec );
      lp_->updateSolution();

      // Compute the new residual R_ := R_ - alpha * AP_
      MVT::MvTimesMatAddMv( -one, *AP_, *alpha, one, *R_ );

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

      // Compute beta := <AP_,Z_> / <P_,AP_>
      // 1) Compute AP_^T * Z_
      // 2) Compute the Cholesky Factorization of pAp (already have)
      // 3) Back and forward solves to compute beta

      // Compute <AP_,Z>
      MVT::MvTransMv( -one, *AP_, *Z_, *beta );

      lltSolver->setVectors( beta, beta );
      info = lltSolver->solve();
      TEUCHOS_TEST_FOR_EXCEPTION
        (info != 0, CGIterationLAPACKFailure,
         prefix << "Failed to compute beta using Cholesky factorization (POTRS).");

      // Compute the new direction vectors P_ = Z_ + P_ * beta
      Teuchos::RCP<MV> Pnew = MVT::CloneCopy( *Z_ );
      MVT::MvTimesMatAddMv(one, *P_, *beta, one, *Pnew);
      P_ = Pnew;

      // Compute orthonormal block of new direction vectors.
      rank = ortho_->normalize( *P_, Teuchos::null );
      TEUCHOS_TEST_FOR_EXCEPTION
        (rank != blockSize_, CGIterationOrthoFailure,
         prefix << "Failed to compute block of orthonormal direction vectors.");

    } // end while (sTest_->checkStatus(this) != Passed)
  }

} // namespace Belos

#endif /* BELOS_BLOCK_CG_ITER_HPP */
