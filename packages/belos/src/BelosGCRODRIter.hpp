// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef BELOS_GCRODR_ITER_HPP
#define BELOS_GCRODR_ITER_HPP

/*! \file BelosGCRODRIter.hpp
    \brief Belos concrete class for performing the GCRO-DR iteration.
*/

#include "BelosConfigDefs.hpp"
#include "BelosTypes.hpp"

#include "BelosLinearProblem.hpp"
#include "BelosMatOrthoManager.hpp"
#include "BelosOutputManager.hpp"
#include "BelosStatusTest.hpp"
#include "BelosOperatorTraits.hpp"
#include "BelosMultiVecTraits.hpp"
#include "BelosDenseMatTraits.hpp"

#include "Teuchos_BLAS.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TimeMonitor.hpp"

/*!	
  \class Belos::GCRODRIter
  
  \brief This class implements the GCRODR iteration, where a
  single-vector Krylov subspace is constructed.  The QR decomposition of 
  block, upper Hessenberg matrix is performed each iteration to update
  the least squares system and give the current linear system residuals.
  \ingroup belos_solver_framework
 
  \author Michael Parks and Heidi Thornquist
*/

namespace Belos {
  
  //! @name GCRODRIter Structures
  //@{
  
  /** \brief Structure to contain pointers to GCRODRIter state variables.
   *
   * This struct is utilized by GCRODRIter::initialize() and GCRODRIter::getState().
   */
  template <class ScalarType, class MV, class DM>
  struct GCRODRIterState {
    /*! \brief The current dimension of the reduction.
     *
     * This should always be equal to BlockGmresIter::getCurSubspaceDim()
     */
    int curDim;
    
    /*! \brief The current Krylov basis. */
    Teuchos::RCP<MV> V;
   
    /*! \brief The recycled subspace and its projection. */
    Teuchos::RCP<MV> U, C;

    /*! \brief The global projection matrix including Krylov subpace and recycled subspace
     */
    Teuchos::RCP<DM> H2;

    /*! \brief The current Hessenberg matrix.
     *
     * The \c curDim by \c curDim leading submatrix of H is the
     * projection of problem->getOperator() by the first \c curDim vectors in V.
     */
    Teuchos::RCP<DM> H;

    /*! \brief The projection of the Krylov subspace against the recycled subspace
     */
    Teuchos::RCP<DM> B;

    GCRODRIterState() : curDim(0), V(Teuchos::null), 
			U(Teuchos::null), C(Teuchos::null),
			H2(Teuchos::null), H(Teuchos::null), 
			B(Teuchos::null)
    {}
  };
  
  //@}
  
  //! @name GCRODRIter Exceptions
  //@{
  
  /** \brief GCRODRIterInitFailure is thrown when the GCRODRIter object is unable to
   * generate an initial iterate in the GCRODRIter::initialize() routine.
   *
   * This std::exception is thrown from the GCRODRIter::initialize() method, which is
   * called by the user or from the GCRODRIter::iterate() method if isInitialized()
   * == \c false.
   *
   * In the case that this std::exception is thrown,
   * GCRODRIter::isInitialized() will be \c false and the user will need to provide
   * a new initial iterate to the iteration.
   */
  class GCRODRIterInitFailure : public BelosError {
    public:
      GCRODRIterInitFailure(const std::string& what_arg) : BelosError(what_arg) {}
  };
  
  /** \brief GCRODRIterOrthoFailure is thrown when the GCRODRIter object is unable to
   * compute independent direction vectors in the GCRODRIter::iterate() routine.
   *
   * This std::exception is thrown from the GCRODRIter::iterate() method.
   *
   */
  class GCRODRIterOrthoFailure : public BelosError {
    public:
      GCRODRIterOrthoFailure(const std::string& what_arg) : BelosError(what_arg) {}
  };
  
  //@}
  
 
  template<class ScalarType, class MV, class OP, class DM>
  class GCRODRIter : virtual public Iteration<ScalarType,MV,OP,DM> {
    
  public:
    
    //
    // Convenience typedefs
    //
    typedef MultiVecTraits<ScalarType,MV,DM> MVT;
    typedef OperatorTraits<ScalarType,MV,OP> OPT;
    typedef DenseMatTraits<ScalarType,DM>    DMT;
    typedef Teuchos::ScalarTraits<ScalarType> SCT;
    typedef typename SCT::magnitudeType MagnitudeType;
    
    //! @name Constructors/Destructor
    //@{ 
    
    /*! \brief %GCRODRIter constructor with linear problem, solver utilities, and parameter list of solver options.
     *
     * This constructor takes pointers required by the linear solver, in addition
     * to a parameter list of options for the linear solver. These options include the following:
     *   - "Num Blocks" - an \c int specifying the maximum number of blocks allocated for the solver basis. Default: 25
     *   - "Restart Timers" = a \c bool specifying whether the timers should be restarted each time iterate() is called. Default: false
     *   - "Keep Hessenberg" = a \c bool specifying whether the upper Hessenberg should be stored separately from the least squares system. Default: false
     */
    GCRODRIter( const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem, 
		const Teuchos::RCP<OutputManager<ScalarType> > &printer,
		const Teuchos::RCP<StatusTest<ScalarType,MV,OP,DM> > &tester,
		const Teuchos::RCP<MatOrthoManager<ScalarType,MV,OP,DM> > &ortho,
		Teuchos::ParameterList &params );
    
    //! Destructor.
    virtual ~GCRODRIter() {};
    //@}
    
    
    //! @name Solver methods
    //@{ 
    
    /*! \brief This method performs block Gmres iterations until the status
     * test indicates the need to stop or an error occurs (in which case, an
     * std::exception is thrown).
     *
     * iterate() will first determine whether the solver is inintialized; if
     * not, it will call initialize() using default arguments. After
     * initialization, the solver performs block Gmres iterations until the
     * status test evaluates as ::Passed, at which point the method returns to
     * the caller. 
     *
     * The block Gmres iteration proceeds as follows:
     * -# The operator problem->applyOp() is applied to the newest \c blockSize vectors in the Krylov basis.
     * -# The resulting vectors are orthogonalized against the previous basis vectors, and made orthonormal.
     * -# The Hessenberg matrix is updated.
     * -# The least squares system is updated.
     *
     * The status test is queried at the beginning of the iteration.
     *
     * Possible exceptions thrown include the GCRODRIterOrthoFailure.
     *
     */
    void iterate();
    
    /*! \brief Initialize the solver to an iterate, providing a complete state.
     *
     * The %GCRODRIter contains a certain amount of state, consisting of the current 
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
    void initialize(GCRODRIterState<ScalarType,MV,DM>& newstate);
    
    /*! \brief Initialize the solver with empty data. Calling this method will result in error,
     *  as GCRODRIter must be initialized with a valid state.
     */
    void initialize() {
      GCRODRIterState<ScalarType,MV,DM> empty;
      initialize(empty);
    }
    
    /*! \brief Get the current state of the linear solver.
     *
     * The data is only valid if isInitialized() == \c true.
     *
     * \returns A GCRODRIterState object containing const pointers to the current
     * solver state.
     */
    GCRODRIterState<ScalarType,MV,DM> getState() const {
      GCRODRIterState<ScalarType,MV,DM> state;
      state.curDim = curDim_;
      state.V = V_;
      state.U = U_;
      state.C = C_;
      state.H2 = H2_;
      state.H = H_;
      state.B = B_;

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
    /*! \note Some solvers, like GMRES, do not compute updates to the solution every iteration.
      This method forces its computation.  Other solvers, like CG, update the solution
      each iteration, so this method will return a zero vector indicating that the linear
      problem contains the current solution.
    */
    Teuchos::RCP<MV> getCurrentUpdate() const;
    
    //! Method for updating QR factorization of upper Hessenberg matrix
    /*! \note If \c dim >= \c getCurSubspaceDim() and \c dim < \c getMaxSubspaceDim(), then 
      the \c dim-th equations of the least squares problem will be updated.
    */
    void updateLSQR( int dim = -1 );
    
    //! Get the dimension of the search subspace used to generate the current solution to the linear problem.
    int getCurSubspaceDim() const { 
      if (!initialized_) return 0;
      return curDim_;
    };
    
    //! Get the maximum dimension allocated for the search subspace.
    int getMaxSubspaceDim() const { return numBlocks_; }
    
    //@}
    
    
    //! @name Accessor methods
    //@{ 
    
    //! Get a constant reference to the linear problem.
    const LinearProblem<ScalarType,MV,OP>& getProblem() const { return *lp_; }
    
    //! Get the maximum number of blocks used by the iterative solver in solving this linear problem.
    int getNumBlocks() const { return numBlocks_; }
    
    //! \brief Set the maximum number of blocks used by the iterative solver.
    void setNumBlocks(int numBlocks) { setSize( recycledBlocks_, numBlocks ); };
    
    //! Get the blocksize to be used by the iterative solver in solving this linear problem.
    int getBlockSize() const { return 1; }
    
    //! \brief Set the blocksize.
    void setBlockSize(int blockSize) {
      TEUCHOS_TEST_FOR_EXCEPTION(blockSize!=1,std::invalid_argument,"Belos::GCRODRIter::setBlockSize(): Cannot use a block size that is not one.");
    }

    //! \brief Set the maximum number of blocks used by the iterative solver and the number of recycled vectors.
    void setSize( int recycledBlocks, int numBlocks ) {
      // only call resize if size changed
      if ( recycledBlocks_ != recycledBlocks )
        recycledBlocks_ = recycledBlocks;
      if ( numBlocks_ != numBlocks ) {
        numBlocks_ = numBlocks;
        cs_.resize( numBlocks_+1 );
        sn_.resize( numBlocks_+1 );
        z_ = DMT::Create( numBlocks_+1, 1, false );
        R_ = DMT::Create( numBlocks_+1, numBlocks, false );
      }
    }

    //! States whether the solver has been initialized or not.
    bool isInitialized() { return initialized_; }
    
    //@}
    
  private:
    
    //
    // Internal methods
    //
    
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
    // numBlocks_ is the size of the allocated space for the Krylov basis, in blocks.
    int numBlocks_; 

    // recycledBlocks_ is the size of the allocated space for the recycled subspace, in blocks.
    int recycledBlocks_; 
   
    // Storage for QR factorization of the least squares system.
    std::vector<ScalarType> sn_;
    std::vector<MagnitudeType> cs_;

    // 
    // Current solver state
    //
    // initialized_ specifies that the basis vectors have been initialized and the iterate() routine
    // is capable of running; _initialize is controlled  by the initialize() member method
    // For the implications of the state of initialized_, please see documentation for initialize()
    bool initialized_;
    
    // Current subspace dimension, and number of iterations performed.
    int curDim_, iter_;
   
    // Pointer to the (0,0) position of the upper Hessenberg matrix.
    int ptrH00_; 
    // 
    // State Storage
    //
    // Krylov vectors.
    Teuchos::RCP<MV> V_;
    //
    // Recycled subspace vectors.
    Teuchos::RCP<MV> U_, C_;
    //
    // Global projected matrix, including H_ and B_
    Teuchos::RCP<DM> H2_;
    //
    // Projected matrices
    // H_ : Projected matrix from the Krylov factorization AV = VH + FE^T
    Teuchos::RCP<DM> H_;
    //
    // B_ : Projected matrix from the recycled subspace B = C^H*A*V
    Teuchos::RCP<DM> B_;
    //
    // QR decomposition of Projected matrices for solving the least squares system HY = B.
    // R_: Upper triangular reduction of H
    // z_: Q applied to right-hand side of the least squares system
    Teuchos::RCP<DM> R_;
    Teuchos::RCP<DM> z_;  
  };
  
  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor.
  template<class ScalarType, class MV, class OP, class DM>
  GCRODRIter<ScalarType,MV,OP,DM>::GCRODRIter(const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem, 
					   const Teuchos::RCP<OutputManager<ScalarType> > &printer,
					   const Teuchos::RCP<StatusTest<ScalarType,MV,OP,DM> > &tester,
					   const Teuchos::RCP<MatOrthoManager<ScalarType,MV,OP,DM> > &ortho,
					   Teuchos::ParameterList &params ):
    lp_(problem), om_(printer), stest_(tester), ortho_(ortho) {
    numBlocks_      = 0;
    recycledBlocks_ = 0;
    initialized_    = false;
    curDim_         = 0;
    iter_           = 0;
    ptrH00_         = 0;
    V_              = Teuchos::null;
    U_              = Teuchos::null;
    C_              = Teuchos::null;
    H2_             = Teuchos::null;
    H_              = Teuchos::null;
    B_              = Teuchos::null;

    // Get the maximum number of blocks allowed for this Krylov subspace
    TEUCHOS_TEST_FOR_EXCEPTION(!params.isParameter("Num Blocks"), std::invalid_argument, "Belos::GCRODRIter::constructor: mandatory parameter \"Num Blocks\" is not specified.");
    int nb = Teuchos::getParameter<int>(params, "Num Blocks");

    TEUCHOS_TEST_FOR_EXCEPTION(!params.isParameter("Recycled Blocks"), std::invalid_argument,"Belos::GCRODRIter::constructor: mandatory parameter \"Recycled Blocks\" is not specified.");
    int rb = Teuchos::getParameter<int>(params, "Recycled Blocks");

    TEUCHOS_TEST_FOR_EXCEPTION(nb <= 0, std::invalid_argument, "Belos::GCRODRIter() was passed a non-positive argument for \"Num Blocks\".");
    TEUCHOS_TEST_FOR_EXCEPTION(rb >= nb, std::invalid_argument, "Belos::GCRODRIter() the number of recycled blocks is larger than the allowable subspace.");

    numBlocks_ = nb;
    recycledBlocks_ = rb;
    cs_.resize( numBlocks_+1 );
    sn_.resize( numBlocks_+1 );
    z_ = DMT::Create( numBlocks_+1, 1, false );
    R_ = DMT::Create( numBlocks_+1, numBlocks_, false );

  }
  
  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Get the current update from this subspace.
  template<class ScalarType, class MV, class OP, class DM>
  Teuchos::RCP<MV> GCRODRIter<ScalarType,MV,OP,DM>::getCurrentUpdate() const {
    //
    // If this is the first iteration of the Arnoldi factorization, 
    // there is no update, so return Teuchos::null. 
    //
    Teuchos::RCP<MV> currentUpdate = Teuchos::null;
    if (curDim_==0) { 
      return currentUpdate; 
    } else {
      const ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();
      const ScalarType zero = Teuchos::ScalarTraits<ScalarType>::zero();
      Teuchos::BLAS<int,ScalarType> blas;
      currentUpdate = MVT::Clone( *V_, 1 );
      //
      //  Make a view and then copy the RHS of the least squares problem.  DON'T OVERWRITE IT!
      //
      Teuchos::RCP<DM> y = DMT::SubviewCopy(*z_, curDim_, 1);
      DMT::SyncDeviceToHost( *y );
      DMT::SyncDeviceToHost( *R_ );
      //
      //  Solve the least squares problem.
      //
      blas.TRSM( Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
		 Teuchos::NON_UNIT_DIAG, curDim_, 1, one,  
		 DMT::GetConstRawHostPtr(*R_), DMT::GetStride(*R_), 
		 DMT::GetRawHostPtr(*y), DMT::GetStride(*y) );
      //
      //  Compute the current update from the Krylov basis; V(:,1:curDim_)*y.
      //
      std::vector<int> index(curDim_);
      for ( int i=0; i<curDim_; i++ ) index[i] = i;
      Teuchos::RCP<const MV> Vjp1 = MVT::CloneView( *V_, index );
      MVT::MvTimesMatAddMv( one, *Vjp1, *y, zero, *currentUpdate );
      //
      //  Add in portion of update from recycled subspace U; U(:,1:recycledBlocks_)*B*y.
      //
      if (U_ != Teuchos::null) {
	Teuchos::RCP<DM> z = DMT::Create( recycledBlocks_, 1 );
	DMT::SyncDeviceToHost( *H2_ );
        blas.GEMM( Teuchos::NO_TRANS, Teuchos::NO_TRANS, recycledBlocks_, 1, curDim_, one, 
		   DMT::GetConstRawHostPtr(*B_), DMT::GetStride(*B_),
		   DMT::GetConstRawHostPtr(*y), DMT::GetStride(*y),
		   zero, DMT::GetRawHostPtr(*z), DMT::GetStride(*z));
        DMT::SyncHostToDevice( *z );
        MVT::MvTimesMatAddMv( -one, *U_, *z, one, *currentUpdate );
      }
    }
    std::vector<MagnitudeType> normUpdate( 1 );
    MVT::MvNorm( *currentUpdate, normUpdate );

    return currentUpdate;
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Get the native residuals stored in this iteration.  
  // Note: This method does not return a MultiVector of the residual vectors, only return Teuchos::null
  template<class ScalarType, class MV, class OP, class DM>
  Teuchos::RCP<const MV> GCRODRIter<ScalarType,MV,OP,DM>::getNativeResiduals( std::vector<MagnitudeType> *norms ) const {
    //
    // NOTE: Make sure the incoming std::vector is the correct size!
    //
    if ( norms && (int)norms->size()==0 )                         
      norms->resize( 1 );                                          
    
    if (norms) {
      DMT::SyncDeviceToHost( *z_ );
      const ScalarType curNativeResid = DMT::ValueConst(*z_,curDim_,0);
      (*norms)[0] = SCT::magnitude (curNativeResid);
    }
    return Teuchos::null;
  }
  
  

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Initialize this iteration object
  template<class ScalarType, class MV, class OP, class DM>
  void GCRODRIter<ScalarType,MV,OP,DM>::initialize(GCRODRIterState<ScalarType,MV,DM>& newstate) {
   
    if (newstate.V != Teuchos::null && newstate.H2 != Teuchos::null) {
      curDim_ = newstate.curDim;
      V_      = newstate.V;
      U_      = newstate.U;
      C_      = newstate.C;
      H2_     = newstate.H2;
      // There is no recycled space, so this iteration is creating the first one.
      if (newstate.U == Teuchos::null) {
        ptrH00_ = recycledBlocks_+1;
        H_ = DMT::Subview( *H2_, numBlocks_+1, numBlocks_, ptrH00_, ptrH00_ );
        B_ = Teuchos::null;
      }
      else {	      
        ptrH00_ = recycledBlocks_;
        H_ = DMT::Subview( *H2_, numBlocks_+1, numBlocks_, ptrH00_, ptrH00_ );
        B_ = DMT::Subview( *H2_, recycledBlocks_, numBlocks_, 0, ptrH00_ );
      }
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(newstate.V == Teuchos::null,std::invalid_argument,"Belos::GCRODRIter::initialize(): GCRODRIterState does not have V initialized.");
      TEUCHOS_TEST_FOR_EXCEPTION(newstate.H2 == Teuchos::null,std::invalid_argument,"Belos::GCRODRIter::initialize(): GCRODRIterState does not have H2 initialized.");
    }

    // the solver is initialized
    initialized_ = true;

  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Iterate until the status test informs us we should stop.
  template<class ScalarType, class MV, class OP, class DM>
  void GCRODRIter<ScalarType,MV,OP,DM>::iterate() {

    TEUCHOS_TEST_FOR_EXCEPTION( initialized_ == false, GCRODRIterInitFailure,"Belos::GCRODRIter::iterate(): GCRODRIter class not initialized." );

    // Force call to setsize to ensure internal storage is correct dimension
    setSize( recycledBlocks_, numBlocks_ );
    
    Teuchos::RCP<MV> Vnext;
    Teuchos::RCP<const MV> Vprev;
    std::vector<int> curind(1);

    // z_ must be zeroed out in order to compute Givens rotations correctly
    DMT::PutScalar(*z_);

    // Orthonormalize the new V_0
    curind[0] = 0;
    Vnext = MVT::CloneViewNonConst(*V_,curind);
    // Orthonormalize first column
    // First, get a view of the first element of z_ to hold the orthonormalization coefficients
    Teuchos::RCP<DM> z0 = DMT::Subview( *z_, 1, 1 );
    int rank = ortho_->normalize( *Vnext, z0 );
    TEUCHOS_TEST_FOR_EXCEPTION(rank != 1,GCRODRIterOrthoFailure, "Belos::GCRODRIter::iterate(): couldn't generate basis of full rank.");

    std::vector<int> prevind(numBlocks_+1);

    ////////////////////////////////////////////////////////////////
    // iterate until the status test tells us to stop.
    //
    // also break if our basis is full
    //
    if (U_ == Teuchos::null) { // iterate without recycle space (because none is available yet)
      while (stest_->checkStatus(this) != Passed && curDim_+1 <= numBlocks_) {
        iter_++;
        int lclDim = curDim_ + 1;

        // Get next index in basis
        curind[0] = lclDim;
        Vnext = MVT::CloneViewNonConst(*V_,curind);

        // Get previous index in basis
        curind[0] = curDim_;
        Vprev = MVT::CloneView(*V_,curind);

        // Compute the next vector in the Krylov basis:  Vnext = Op*Vprev
        lp_->apply(*Vprev,*Vnext);

        // Remove all previous Krylov basis vectors from Vnext and put coefficients in H and R.

        // Get a view of all the previous vectors
        prevind.resize(lclDim);
        for (int i=0; i<lclDim; i++) { prevind[i] = i; }
        Vprev = MVT::CloneView(*V_,prevind);
        Teuchos::Array<Teuchos::RCP<const MV> > AVprev(1, Vprev);

        // Get a view of the part of the Hessenberg matrix needed to hold the ortho coeffs.
	Teuchos::RCP<DM> subH = DMT::Subview(*H2_, lclDim, 1, ptrH00_, ptrH00_+curDim_);
        Teuchos::Array<Teuchos::RCP<DM> > AsubH( 1, subH );

        // Get a view of the part of the Hessenberg matrix needed to hold the norm coeffs.
        Teuchos::RCP<DM> subR = DMT::Subview(*H2_, 1, 1, ptrH00_+lclDim, ptrH00_+curDim_);

        // Project out the previous Krylov vectors and normalize the next vector.
	rank = ortho_->projectAndNormalize(*Vnext, AsubH, subR, AVprev);

	// Copy over the coefficients to R just in case we run into an error.
	Teuchos::RCP<DM> subR2 = DMT::Subview(*R_, lclDim+1, 1, 0, curDim_);
        Teuchos::RCP<const DM> subH2 = DMT::SubviewConst(*H2_, lclDim+1, 1, ptrH00_, ptrH00_+curDim_);
	DMT::Assign(*subR2, *subH2);
	subR2 = Teuchos::null;

        TEUCHOS_TEST_FOR_EXCEPTION(rank != 1,GCRODRIterOrthoFailure, "Belos::GCRODRIter::iterate(): couldn't generate basis of full rank.");

        // Update the QR factorization of the upper Hessenberg matrix
        updateLSQR();

        // Update basis dim
        curDim_++;
      } // end while (statusTest == false)
    }
    else { // iterate with recycle space
      while (stest_->checkStatus(this) != Passed && curDim_+1 <= numBlocks_) {
        iter_++;
        int lclDim = curDim_ + 1; 

        // Get the current part of the basis.
        curind[0] = lclDim;
        Vnext = MVT::CloneViewNonConst(*V_,curind);

        // Get a view of the previous vectors.
        // This is used for orthogonalization and for computing V^H K H
        curind[0] = curDim_;
        Vprev = MVT::CloneView(*V_,curind);

        // Compute the next std::vector in the Krylov basis:  Vnext = Op*Vprev
        lp_->apply(*Vprev,*Vnext);
        Vprev = Teuchos::null;

	DMT::SyncHostToDevice(*H2_);

        // First, remove the recycled subspace (C) from Vnext and put coefficients in B.
        Teuchos::Array<Teuchos::RCP<const MV> > C(1, C_);
        Teuchos::RCP<DM> subB = DMT::Subview(*H2_, recycledBlocks_, 1, 0, ptrH00_+curDim_);
        Teuchos::RCP<DM> tmpB = DMT::Create(recycledBlocks_, 1);
	Teuchos::Array<Teuchos::RCP<DM> > AsubB( 1, subB );

        // Project out the recycled subspace.
        ortho_->project( *Vnext, AsubB, C );

        // Now, remove all previous Krylov basis vectors from Vnext and put coefficients in H and R.          
        // Get a view of all the previous vectors
        prevind.resize(lclDim);
        for (int i=0; i<lclDim; i++) { prevind[i] = i; }
        Vprev = MVT::CloneView(*V_,prevind);
        Teuchos::Array<Teuchos::RCP<const MV> > AVprev(1, Vprev);
	
        // Get a view of the part of the Hessenberg matrix needed to hold the ortho coeffs.
        Teuchos::RCP<DM> subH = DMT::Subview(*H2_, lclDim, 1, ptrH00_, ptrH00_+curDim_);
        Teuchos::Array<Teuchos::RCP<DM> > AsubH(1, subH);
      
        // Get a view of the part of the Hessenberg matrix needed to hold the norm coeffs.
        Teuchos::RCP<DM> subR = DMT::Subview(*H2_, 1, 1, ptrH00_+lclDim, ptrH00_+curDim_);

        // Project out the previous Krylov vectors and normalize the next vector.
        rank = ortho_->projectAndNormalize(*Vnext, AsubH, subR, AVprev);
 
        // Copy over the coefficients to R just in case we run into an error.
	Teuchos::RCP<DM> subR2 = DMT::Subview(*R_, lclDim+1, 1, 0, curDim_);
        Teuchos::RCP<const DM> subH2 = DMT::SubviewConst(*H2_, lclDim+1, 1, ptrH00_, ptrH00_+curDim_);
	DMT::Assign(*subR2, *subH2);
	subR2 = Teuchos::null;
      
        TEUCHOS_TEST_FOR_EXCEPTION(rank != 1,GCRODRIterOrthoFailure, "Belos::GCRODRIter::iterate(): couldn't generate basis of full rank.");

        // Update the QR factorization of the upper Hessenberg matrix
        updateLSQR();

        // Update basis dim
        curDim_++;
      
      } // end while (statusTest == false)
    } // end if (U_ == Teuchos::null)
   
  }

  
  template<class ScalarType, class MV, class OP, class DM>
  void GCRODRIter<ScalarType,MV,OP,DM>::updateLSQR( int dim ) {

    int i;
    const ScalarType zero = Teuchos::ScalarTraits<ScalarType>::zero();
   
    // Get correct dimension based on input "dim"
    // Remember that ortho failures result in an exit before updateLSQR() is called.
    // Therefore, it is possible that dim == curDim_.
    int curDim = curDim_;
    if ( (dim >= curDim_) && (dim < getMaxSubspaceDim()) )
      curDim = dim;
    
    Teuchos::BLAS<int, ScalarType> blas;
    //
    // Apply previous transformations and compute new transformation to reduce upper-Hessenberg
    // system to upper-triangular form.
    //
    // QR factorization of Least-Squares system with Givens rotations
    //
    DMT::SyncDeviceToHost(*R_);
    DMT::SyncDeviceToHost(*z_);
    //
    for (i=0; i<curDim; i++) {
      //
      // Apply previous Givens rotations to new column of Hessenberg matrix
      //
      blas.ROT( 1, &(DMT::Value(*R_,i,curDim)), 1, &(DMT::Value(*R_,i+1, curDim)), 1, &cs_[i], &sn_[i] );
 
    }
    //
    // Calculate new Givens rotation
    //
    blas.ROTG( &(DMT::Value(*R_,curDim,curDim)), &(DMT::Value(*R_,curDim+1,curDim)), &cs_[curDim], &sn_[curDim] );
    DMT::Value(*R_,curDim+1,curDim) = zero;
    //
    // Update RHS w/ new transformation
    //
    blas.ROT( 1, &(DMT::Value(*z_,curDim,0)), 1, &(DMT::Value(*z_,curDim+1,0)), 1, &cs_[curDim], &sn_[curDim] );
    //
    DMT::SyncHostToDevice(*R_);
    DMT::SyncHostToDevice(*z_);
    //
  } // end updateLSQR()

} // end Belos namespace

#endif /* BELOS_GCRODR_ITER_HPP */
