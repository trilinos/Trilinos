// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef BELOS_RCG_ITER_HPP
#define BELOS_RCG_ITER_HPP

/*! \file BelosRCGIter.hpp
    \brief Belos concrete class for performing the RCG iteration.
*/

#include "BelosConfigDefs.hpp"
#include "BelosTypes.hpp"

#include "BelosLinearProblem.hpp"
#include "BelosMatOrthoManager.hpp"
#include "BelosOutputManager.hpp"
#include "BelosStatusTest.hpp"
#include "BelosOperatorTraits.hpp"
#include "BelosMultiVecTraits.hpp"
#include "BelosCGIteration.hpp"

#include "Teuchos_LAPACK.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TimeMonitor.hpp"

// MLP Remove after debugging
#include <fstream>
#include <iomanip>

/*!
  \class Belos::RCGIter

  \brief This class implements the RCG iteration, where a
  single-std::vector Krylov subspace is constructed.

  \ingroup belos_solver_framework

  \author Michael Parks and Heidi Thornquist
*/

namespace Belos {

  //! @name RCGIter Structures
  //@{

  /** \brief Structure to contain pointers to RCGIter state variables.
   *
   * This struct is utilized by RCGIter::initialize()
   */
  template <class ScalarType, class MV>
  struct RCGIterState {
    /*! \brief The current dimension of the reduction.
     *
     * This should always be equal to BlockGmresIter::getCurSubspaceDim()
     */
    int curDim;

    /*! \brief The current Krylov basis. */
    Teuchos::RCP<MV> P;

    /*! \brief A times current search vector */
    Teuchos::RCP<MV> Ap;

    /*! \brief The current residual. */
    Teuchos::RCP<MV> r;

    /*! \brief The current preconditioned residual. */
    Teuchos::RCP<MV> z;

    /*! \brief Flag to indicate the recycle space should be used */
    bool existU;

    /*! \brief The recycled subspace and its image. */
    Teuchos::RCP<MV> U, AU;

    /*! \brief Coefficients arising in RCG iteration
     */
    Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > Alpha;
    Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > Beta;
    Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > D;
    Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > rTz_old;

    /*! \brief Solutions to local least-squares problems */
    Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > Delta;

    /*! \brief The LU factorization of the matrix U^T A U  */
    Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > LUUTAU;
    /*! \brief Data from LU factorization of U^T A U */
    Teuchos::RCP<std::vector<int> > ipiv;


    RCGIterState() : curDim(0), P(Teuchos::null), Ap(Teuchos::null), r(Teuchos::null),
                     z(Teuchos::null),
                     existU(false),
		     U(Teuchos::null), AU(Teuchos::null),
		     Alpha(Teuchos::null), Beta(Teuchos::null), D(Teuchos::null), rTz_old(Teuchos::null),
		     Delta(Teuchos::null), LUUTAU(Teuchos::null), ipiv(Teuchos::null)
    {}
  };

  //@}

  template<class ScalarType, class MV, class OP>
  class RCGIter : virtual public Iteration<ScalarType,MV,OP> {

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

    /*! \brief %RCGIter constructor with linear problem, solver utilities, and parameter list of solver options.
     *
     * This constructor takes pointers required by the linear solver, in addition
     * to a parameter list of options for the linear solver. These options include the following:
     *   - "Num Blocks" - an \c int specifying the maximum number of blocks allocated for the solver basis. Default: 25
     *   - "Restart Timers" = a \c bool specifying whether the timers should be restarted each time iterate() is called. Default: false
     */
    RCGIter( const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem,
		const Teuchos::RCP<OutputManager<ScalarType> > &printer,
		const Teuchos::RCP<StatusTest<ScalarType,MV,OP> > &tester,
		Teuchos::ParameterList &params );

    //! Destructor.
    virtual ~RCGIter() {};
    //@}


    //! @name Solver methods
    //@{

   /*! \brief This method performs RCG iterations until the status
    * test indicates the need to stop or an error occurs (in which case, an
    * std::exception is thrown).
    *
    * iterate() will first determine whether the solver is initialized; if
    * not, it will call initialize() using default arguments. After
    * initialization, the solver performs RCG iterations until the
    * status test evaluates as ::Passed, at which point the method returns to
    * the caller.
    *
    * The status test is queried at the beginning of the iteration.
    */
    void iterate();

   /*! \brief Initialize the solver to an iterate, providing a complete state.
    *
    * The %RCGIter contains a certain amount of state, consisting of the current
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
    void initialize(RCGIterState<ScalarType,MV> &newstate);

   /*! \brief Initialize the solver with the initial vectors from the linear problem
    *  or random data.
    */
    void initialize()
    {
      RCGIterState<ScalarType,MV> empty;
      initialize(empty);
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
    Teuchos::RCP<const MV> getNativeResiduals( std::vector<MagnitudeType> * /* norms */ ) const { return r_; }

    //! Get the current update to the linear system.
    /*! \note Some solvers, like GMRES, do not compute updates to the solution every iteration.
      This method forces its computation.  Other solvers, like CG and RCG update the solution
      each iteration, so this method will return a zero std::vector indicating that the linear
      problem contains the current solution.
    */
    Teuchos::RCP<MV> getCurrentUpdate() const { return Teuchos::null; }

    //! Get the dimension of the search subspace used to generate the current solution to the linear problem.
    int getCurSubspaceDim() const {
      if (!initialized_) return 0;
      return curDim_;
    };

    //! Get the maximum dimension allocated for the search subspace.
    int getMaxSubspaceDim() const { return numBlocks_+1; }

    //@}


    //! @name Accessor methods
    //@{

    //! Get a constant reference to the linear problem.
    const LinearProblem<ScalarType,MV,OP>& getProblem() const { return *lp_; }

    //! Get the maximum number of blocks used by the iterative solver in solving this linear problem.
    int getNumBlocks() const { return numBlocks_; }

    //! \brief Set the maximum number of blocks used by the iterative solver.
    void setNumBlocks(int numBlocks) { setSize( recycleBlocks_, numBlocks ); };

    //! Get the maximum number of recycled blocks used by the iterative solver in solving this linear problem.
    int getRecycledBlocks() const { return recycleBlocks_; }

    //! \brief Set the maximum number of recycled blocks used by the iterative solver.
    void setRecycledBlocks(int recycleBlocks) { setSize( recycleBlocks, numBlocks_ ); };

    //! Get the blocksize to be used by the iterative solver in solving this linear problem.
    int getBlockSize() const { return 1; }

    //! \brief Set the blocksize.
    void setBlockSize(int blockSize) {
      TEUCHOS_TEST_FOR_EXCEPTION(blockSize!=1,std::invalid_argument,
			 "Belos::RCGIter::setBlockSize(): Cannot use a block size that is not one.");
    }

    //! \brief Set the maximum number of blocks used by the iterative solver and the number of recycled vectors.
    void setSize( int recycleBlocks, int numBlocks );

    //! States whether the solver has been initialized or not.
    bool isInitialized() { return initialized_; }

    //@}

  private:

    //
    // Internal methods
    //

    //
    // Classes input through constructor that define the linear problem to be solved.
    //
    const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> >    lp_;
    const Teuchos::RCP<OutputManager<ScalarType> >          om_;
    const Teuchos::RCP<StatusTest<ScalarType,MV,OP> >       stest_;

    //
    // Algorithmic parameters
    //
    // numBlocks_ is the size of the allocated space for the Krylov basis, in blocks.
    int numBlocks_;

    // recycleBlocks_ is the size of the allocated space for the recycled subspace, in blocks.
    int recycleBlocks_;

    //
    // Current solver state
    //
    // initialized_ specifies that the basis vectors have been initialized and the iterate() routine
    // is capable of running; _initialize is controlled  by the initialize() member method
    // For the implications of the state of initialized_, please see documentation for initialize()
    bool initialized_;

    // Current subspace dimension, and number of iterations performed.
    int curDim_, iter_;

    //
    // State Storage
    //
    // Search vectors
    Teuchos::RCP<MV> P_;
    //
    // A times current search vector
    Teuchos::RCP<MV> Ap_;
    //
    // Residual vector
    Teuchos::RCP<MV> r_;
    //
    // Preconditioned residual
    Teuchos::RCP<MV> z_;
    //
    // Flag to indicate that the recycle space should be used
    bool existU_;
    // Recycled subspace and its image
    Teuchos::RCP<MV> U_, AU_;
    //
    // Coefficients arising in RCG iteration
    Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > Alpha_,Beta_,D_;
    //
    // Solutions to local least-squares problems
    Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > Delta_;
    //
    // The LU factorization of the matrix U^T A U
    Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > LUUTAU_;
    //
    // Data from LU factorization of UTAU
    Teuchos::RCP<std::vector<int> > ipiv_;
    //
    // The scalar r'*z
    Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > rTz_old_;
  };

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor.
  template<class ScalarType, class MV, class OP>
  RCGIter<ScalarType,MV,OP>::RCGIter(const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem,
				     const Teuchos::RCP<OutputManager<ScalarType> > &printer,
				     const Teuchos::RCP<StatusTest<ScalarType,MV,OP> > &tester,
					   Teuchos::ParameterList &params ):
    lp_(problem),
    om_(printer),
    stest_(tester),
    numBlocks_(0),
    recycleBlocks_(0),
    initialized_(false),
    curDim_(0),
    iter_(0),
    existU_(false)
  {
    // Get the maximum number of blocks allowed for this Krylov subspace
    TEUCHOS_TEST_FOR_EXCEPTION(!params.isParameter("Num Blocks"), std::invalid_argument,
                       "Belos::RCGIter::constructor: mandatory parameter \"Num Blocks\" is not specified.");
    int nb = Teuchos::getParameter<int>(params, "Num Blocks");

    TEUCHOS_TEST_FOR_EXCEPTION(!params.isParameter("Recycled Blocks"), std::invalid_argument,
                       "Belos::RCGIter::constructor: mandatory parameter \"Recycled Blocks\" is not specified.");
    int rb = Teuchos::getParameter<int>(params, "Recycled Blocks");

    // Set the number of blocks and allocate data
    setSize( rb, nb );
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Set the block size and make necessary adjustments.
  template <class ScalarType, class MV, class OP>
  void RCGIter<ScalarType,MV,OP>::setSize( int recycleBlocks, int numBlocks )
  {

    TEUCHOS_TEST_FOR_EXCEPTION(numBlocks <= 0, std::invalid_argument, "Belos::RCGIter::setSize() was passed a non-positive argument for \"Num Blocks\".");
    TEUCHOS_TEST_FOR_EXCEPTION(recycleBlocks <= 0, std::invalid_argument, "Belos::RCGIter::setSize() was passed a non-positive argument for \"Recycled Blocks\".");
    TEUCHOS_TEST_FOR_EXCEPTION(recycleBlocks >= numBlocks, std::invalid_argument, "Belos::RCGIter::setSize() the number of recycled blocks is larger than the allowable subspace.");

    numBlocks_ = numBlocks;
    recycleBlocks_ = recycleBlocks;

  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Initialize this iteration object
  template <class ScalarType, class MV, class OP>
  void RCGIter<ScalarType,MV,OP>::initialize(RCGIterState<ScalarType,MV> &newstate)
  {

    if (newstate.P != Teuchos::null &&
        newstate.Ap != Teuchos::null &&
        newstate.r != Teuchos::null &&
        newstate.z != Teuchos::null &&
        newstate.U != Teuchos::null &&
        newstate.AU != Teuchos::null &&
        newstate.Alpha != Teuchos::null &&
        newstate.Beta != Teuchos::null &&
        newstate.D != Teuchos::null &&
        newstate.Delta != Teuchos::null &&
        newstate.LUUTAU != Teuchos::null &&
        newstate.ipiv != Teuchos::null &&
        newstate.rTz_old != Teuchos::null) {

      curDim_ = newstate.curDim;
      P_ = newstate.P;
      Ap_ = newstate.Ap;
      r_ = newstate.r;
      z_ = newstate.z;
      existU_ = newstate.existU;
      U_ = newstate.U;
      AU_ = newstate.AU;
      Alpha_ = newstate.Alpha;
      Beta_ = newstate.Beta;
      D_ = newstate.D;
      Delta_ = newstate.Delta;
      LUUTAU_ = newstate.LUUTAU;
      ipiv_ = newstate.ipiv;
      rTz_old_ = newstate.rTz_old;
    }
    else {

      TEUCHOS_TEST_FOR_EXCEPTION(newstate.P == Teuchos::null,std::invalid_argument,
                         "Belos::RCGIter::initialize(): RCGIterState does not have P initialized.");

      TEUCHOS_TEST_FOR_EXCEPTION(newstate.Ap == Teuchos::null,std::invalid_argument,
                         "Belos::RCGIter::initialize(): RCGIterState does not have Ap initialized.");

      TEUCHOS_TEST_FOR_EXCEPTION(newstate.r == Teuchos::null,std::invalid_argument,
                         "Belos::RCGIter::initialize(): RCGIterState does not have r initialized.");

      TEUCHOS_TEST_FOR_EXCEPTION(newstate.z == Teuchos::null,std::invalid_argument,
                         "Belos::RCGIter::initialize(): RCGIterState does not have z initialized.");

      TEUCHOS_TEST_FOR_EXCEPTION(newstate.U == Teuchos::null,std::invalid_argument,
                         "Belos::RCGIter::initialize(): RCGIterState does not have U initialized.");

      TEUCHOS_TEST_FOR_EXCEPTION(newstate.AU == Teuchos::null,std::invalid_argument,
                         "Belos::RCGIter::initialize(): RCGIterState does not have AU initialized.");

      TEUCHOS_TEST_FOR_EXCEPTION(newstate.Alpha == Teuchos::null,std::invalid_argument,
                         "Belos::RCGIter::initialize(): RCGIterState does not have Alpha initialized.");

      TEUCHOS_TEST_FOR_EXCEPTION(newstate.Beta == Teuchos::null,std::invalid_argument,
                         "Belos::RCGIter::initialize(): RCGIterState does not have Beta initialized.");

      TEUCHOS_TEST_FOR_EXCEPTION(newstate.D == Teuchos::null,std::invalid_argument,
                         "Belos::RCGIter::initialize(): RCGIterState does not have D initialized.");

      TEUCHOS_TEST_FOR_EXCEPTION(newstate.Delta == Teuchos::null,std::invalid_argument,
                         "Belos::RCGIter::initialize(): RCGIterState does not have Delta initialized.");

      TEUCHOS_TEST_FOR_EXCEPTION(newstate.LUUTAU == Teuchos::null,std::invalid_argument,
                         "Belos::RCGIter::initialize(): RCGIterState does not have LUUTAU initialized.");

      TEUCHOS_TEST_FOR_EXCEPTION(newstate.ipiv == Teuchos::null,std::invalid_argument,
                         "Belos::RCGIter::initialize(): RCGIterState does not have ipiv initialized.");

      TEUCHOS_TEST_FOR_EXCEPTION(newstate.rTz_old == Teuchos::null,std::invalid_argument,
                         "Belos::RCGIter::initialize(): RCGIterState does not have rTz_old initialized.");

    }

    // the solver is initialized
    initialized_ = true;

  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Iterate until the status test informs us we should stop.
  template <class ScalarType, class MV, class OP>
  void RCGIter<ScalarType,MV,OP>::iterate()
  {
    TEUCHOS_TEST_FOR_EXCEPTION( initialized_ == false, CGIterateFailure,
                        "Belos::RCGIter::iterate(): RCGIter class not initialized." );

    // We'll need LAPACK
    Teuchos::LAPACK<int,ScalarType> lapack;

    // Create convenience variables for zero and one.
    ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();
    ScalarType zero = Teuchos::ScalarTraits<ScalarType>::zero();

    // Allocate memory for scalars
    std::vector<int> index(1);
    Teuchos::SerialDenseMatrix<int,ScalarType> pAp(1,1), rTz(1,1);

    // Get the current solution std::vector.
    Teuchos::RCP<MV> cur_soln_vec = lp_->getCurrLHSVec();

    // Check that the current solution std::vector only has one column.
    TEUCHOS_TEST_FOR_EXCEPTION( MVT::GetNumberVecs(*cur_soln_vec) != 1, CGIterateFailure,
                        "Belos::RCGIter::iterate(): current linear system has more than one std::vector!" );

    // Compute the current search dimension.
    int searchDim = numBlocks_+1;

    // index of iteration within current cycle
    int i_ = 0;

    ////////////////////////////////////////////////////////////////
    // iterate until the status test tells us to stop.
    //
    // also break if our basis is full
    //
    Teuchos::RCP<const MV> p_ = Teuchos::null;
    Teuchos::RCP<MV> pnext_ = Teuchos::null;
    while (stest_->checkStatus(this) != Passed && curDim_+1 <= searchDim) {

      // Ap = A*p;
      index.resize( 1 );
      index[0] = i_;
      p_  = MVT::CloneView( *P_,  index );
      lp_->applyOp( *p_, *Ap_ );

      // d = p'*Ap;
      MVT::MvTransMv( one, *p_, *Ap_, pAp );
      (*D_)(i_,0) = pAp(0,0);

      // alpha = rTz_old / pAp
      (*Alpha_)(i_,0) = (*rTz_old_)(0,0) / pAp(0,0);

      // Check that alpha is a positive number
      TEUCHOS_TEST_FOR_EXCEPTION( SCT::real(pAp(0,0)) <= zero, CGPositiveDefiniteFailure,
                                  "Belos::RCGIter::iterate(): non-positive value for p^H*A*p encountered!" );

      // x = x + (alpha * p);
      MVT::MvAddMv( one, *cur_soln_vec, (*Alpha_)(i_,0), *p_, *cur_soln_vec );
      lp_->updateSolution();

      // r = r - (alpha * Ap);
      MVT::MvAddMv( one, *r_, -(*Alpha_)(i_,0), *Ap_, *r_ );

      std::vector<MagnitudeType> norm(1);
      MVT::MvNorm( *r_, norm );
//printf("i = %i\tnorm(r) = %e\n",i_,norm[0]);

      // z = M\r
      if ( lp_->getLeftPrec() != Teuchos::null ) {
        lp_->applyLeftPrec( *r_, *z_ );
      }
      else if ( lp_->getRightPrec() != Teuchos::null ) {
        lp_->applyRightPrec( *r_, *z_ );
      }
      else {
        z_ = r_;
      }

      // rTz_new = r'*z;
      MVT::MvTransMv( one, *r_, *z_, rTz );

      // beta = rTz_new/rTz_old;
      (*Beta_)(i_,0) = rTz(0,0) / (*rTz_old_)(0,0);

      // rTz_old = rTz_new;
      (*rTz_old_)(0,0) = rTz(0,0);

      // get pointer for next p
      index.resize( 1 );
      index[0] = i_+1;
      pnext_ = MVT::CloneViewNonConst( *P_,  index );

      if (existU_) {
        // mu = UTAU \ (AU'*z);
        Teuchos::SerialDenseMatrix<int,ScalarType> mu( Teuchos::View, *Delta_, recycleBlocks_, 1, 0, i_ );
        MVT::MvTransMv( one, *AU_, *z_, mu );
        char TRANS = 'N';
        int info;
        lapack.GETRS( TRANS, recycleBlocks_, 1, LUUTAU_->values(), LUUTAU_->stride(), &(*ipiv_)[0], mu.values(), mu.stride(), &info );
        TEUCHOS_TEST_FOR_EXCEPTION(info != 0, CGIterationLAPACKFailure,
                           "Belos::RCGIter::solve(): LAPACK GETRS failed to compute a solution.");
        // p = -(U*mu) + (beta*p) + z (in two steps)
        // p = (beta*p) + z;
        MVT::MvAddMv( (*Beta_)(i_,0), *p_, one, *z_, *pnext_ );
        // pnext = -(U*mu) + (one)*pnext;
        MVT::MvTimesMatAddMv( -one, *U_, mu, one, *pnext_ );
      }
      else {
        // p = (beta*p) + z;
        MVT::MvAddMv( (*Beta_)(i_,0), *p_, one, *z_, *pnext_ );
      }

      // Done with this view; release pointer
      p_ = Teuchos::null;
      pnext_ = Teuchos::null;

      // increment iteration count and dimension index
      i_++;
      iter_++;
      curDim_++;

    } // end while (statusTest == false)

   }

} // end Belos namespace

#endif /* BELOS_RCG_ITER_HPP */
