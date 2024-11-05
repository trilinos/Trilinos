// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef BELOS_PCPG_ITER_HPP
#define BELOS_PCPG_ITER_HPP

/*! \file BelosPCPGIter.hpp
    \brief Belos concrete class to iterate Preconditioned Conjugate Projected Gradients
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

#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TimeMonitor.hpp"

/*!	
  \class Belos::PCPGIter
  
  \brief This class implements the PCPG iteration, where a
  single-std::vector Krylov subspace is constructed.  The documentation
  refers to blocks, but note that at this point, all blocks have unit
  dimension.
 
  \author David Day
*/

namespace Belos {
  
  //! @name PCPGIter Structures
  //@{
  
  /** \brief Structure to contain pointers to PCPGIter state variables.
   *
   * The structure is utilized by initialize() and getState().
   */
  template <class ScalarType, class MV>
  struct PCPGIterState {
    /*! \brief The current dimension of the reduction.
     *
     * This ought always to equal PCPGIter::getCurSubspaceDim()
     */
    /*! \brief Number of block columns in matrices C and U */
    int curDim;
    /*! \brief Number of block columns in matrices C and U before current iteration */
    int prevUdim;

    /*! \brief The current residual. */
    Teuchos::RCP<MV> R;

    /*! \brief The current preconditioned residual. */
    Teuchos::RCP<MV> Z;

    /*! \brief The current decent direction std::vector */
    Teuchos::RCP<MV> P;

    /*! \brief The matrix A applied to current decent direction std::vector */
    Teuchos::RCP<MV> AP;

    /*! \brief The recycled subspace */
    Teuchos::RCP<MV> U;

    /*! \brief C = AU, U spans recycled subspace */
    Teuchos::RCP<MV> C;

    /*! \brief The current Hessenberg matrix.
     *
     * The \c curDim by \c curDim D = diag(P'*AP) = U' * C
     */
    Teuchos::RCP<const Teuchos::SerialDenseMatrix<int,ScalarType> > D;

    PCPGIterState() : curDim(0), 
                      prevUdim(0), 
                      R(Teuchos::null), Z(Teuchos::null), 
                      P(Teuchos::null), AP(Teuchos::null),
                      U(Teuchos::null), C(Teuchos::null),
		      D(Teuchos::null)
    {}
  };
  
  //@}
  
  template<class ScalarType, class MV, class OP>
  class PCPGIter : virtual public Iteration<ScalarType,MV,OP> {
    
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
    
    /*! \brief %PCPGIter constructor with linear problem, solver utilities, and parameter list of solver options.
     *
     * This constructor takes pointers required by the linear solver, in addition
     * to a parameter list of options for the linear solver. These options include the following:
     *   - "Restart Timers" = a \c bool specifying whether the timers should be restarted each time iterate() is called. Default: false
     *   - "Keep Diagonal" = a \c bool specifying whether the upper Hessenberg should be stored separately from the least squares system. Default: false
     */
    PCPGIter( const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem, 
		const Teuchos::RCP<OutputManager<ScalarType> > &printer,
		const Teuchos::RCP<StatusTest<ScalarType,MV,OP> > &tester,
		const Teuchos::RCP<MatOrthoManager<ScalarType,MV,OP> > &ortho,
		Teuchos::ParameterList &params );
    
    //! Destructor.
    virtual ~PCPGIter() {};
    //@}
    
    
    //! @name Solver methods
    //@{ 
    
    /*! \brief PCPGIter iterates CG until the status test either requests a stop or detects an error.  
     * In the latter case, std::exception is thrown.
     *
     * iterate() will first determine whether or not the solver is inintialized; if not, iterate()
     * will call initialize() using default arguments. After initialization, the solver performs CG
     * iterations until the status test evaluates as ::Passed, at which point the method returns to
     * the caller. 
     *
     * The PCPG iteration proceeds as follows:
     * -# the operator problem->applyOp() is applied to the newest vector in the Krylov basis,
     * -# the result is (approximately) A-orthogonalized to the previous basis vectors,
     * -# the coupled two-term recurrence is iterated,
     * -# the search direction P is projected into a complement of the seed space U.
     *
     * The status test is queried at the beginning of the iteration.
     * Potential CG exceptions include IterationInit, Iterate
     *
     */
    void iterate();
    
    /*! \brief Initialize the solver to an iterate, providing a complete state.
     *
     * The %PCPGIter state consists of the %CGIter state, the stored search directions
     * and the seed space.  The constructor calls setSize() which calls setStateSize().
     *
     * initialize(IterState) also calls setStateSize(), passing the current seed space to Iterate.
     * initialize(IterState) sets the state to the specified IterState and then 
     * initialized_ := true.  Fundamental state changes cause initialized_ := false.
     *
     * \post 
     * <li>isInitialized() == \c true
     *
     * Optionally, the user may specify any component of the state using initialize().  
     * Any component of the state not given to initialize() will be generated.
     *
     * \note For any pointer in \c newstate which directly points to the multivectors in 
     * the solver, the data is not (supposed to be) copied.
     */
    void initialize(PCPGIterState<ScalarType,MV>& newstate);
    
    /*! \brief Initialize the solver with the initial vectors from the linear problem.
     *  An exception is thrown if initialzed is called and newstate.R is null.
     */
    void initialize()
    {
      PCPGIterState<ScalarType,MV> empty;
      initialize(empty);
    }
    
    /*! \brief Get the current state of the linear solver.
     *
     * The data is only valid if isInitialized() == \c true.
     *
     * \returns A PCPGIterState object containing const pointers to the current
     * solver state.
     */
    PCPGIterState<ScalarType,MV> getState() const {
      PCPGIterState<ScalarType,MV> state;
      state.Z = Z_;         // CG state
      state.P = P_;
      state.AP = AP_;
      state.R = R_;
      state.U = U_;         // seed state 
      state.C = C_; 
      state.D = D_;
      state.curDim = curDim_;
      state.prevUdim = prevUdim_;
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

    //! Get the current update to the linear system solution?.
    /*! \note getCurrentUpdate returns a null pointer indicating that the linear problem
        contains the current solution.
    */
    Teuchos::RCP<MV> getCurrentUpdate() const { return Teuchos::null; }

    //! Get the current dimension of the whole seed subspace.
    int getCurSubspaceDim() const { 
      if (!initialized_) return 0;
      return curDim_;
    };

    //! Get the dimension of the search subspace used to solve the current solution to the linear problem.
    int getPrevSubspaceDim() const { 
      if (!initialized_) return 0;
      return prevUdim_;
    };
    
    //@}
    
    
    //! @name Accessor methods
    //@{ 
    
    //! Get a constant reference to the linear problem.
    const LinearProblem<ScalarType,MV,OP>& getProblem() const { return *lp_; }
    
    //! Get the maximum number of blocks used by the iterative solver in solving this linear problem.
    int getBlockSize() const { return 1; }
    
    //! Get the maximum number of recycled blocks used by the iterative solver in solving this linear problem.
    int getNumRecycledBlocks() const { return savedBlocks_; }

    //! Get the blocksize to be used by the iterative solver in solving this linear problem.
    
    //! \brief Set the blocksize.
    void setBlockSize(int blockSize) {
      TEUCHOS_TEST_FOR_EXCEPTION(blockSize!=1,std::invalid_argument,
			 "Belos::PCPGIter::setBlockSize(): Cannot use a block size that is not one.");
    }

    //! \brief Set the maximum number of saved or recycled blocks used by the iterative solver
    void setSize( int savedBlocks );

    //! States whether the solver has been initialized or not.
    bool isInitialized() { return initialized_; }

    //! tell the Iterator to "reset" itself;  delete and rebuild the seed space.
    void resetState(); 
    
    //@}
    
  private:
    
    //
    // Internal methods
    //
    //! Method for initalizing the state storage needed by PCPG
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
    // savedBlocks_ is the number of blocks allocated for the reused subspace
    int savedBlocks_; 
    //
    // 
    // Current solver state
    //
    // initialized_ specifies that the basis vectors have been initialized and the iterate() routine
    // is capable of running; _initialize is controlled  by the initialize() member method
    // For the implications of the state of initialized_, please see documentation for initialize()
    bool initialized_;
    
    // stateStorageInitialized_ indicates that the state storage has be initialized to the current
    // savedBlocks_.  State storage initialization may be postponed if the linear problem was
    // generated without either the right-hand side or solution vectors.
    bool stateStorageInitialized_;

    // keepDiagonal_ specifies that the iteration must keep the diagonal matrix of pivots 
    bool keepDiagonal_;

    // initDiagonal_ specifies that the iteration will reinitialize the diagonal matrix by zeroing
    // out all entries before an iteration is started.
    bool initDiagonal_;
    
    // Current subspace dimension
    int curDim_;

    // Dimension of seed space used to solve current linear system
    int prevUdim_;
    
    // Number of iterations performed
    int iter_;
    // 
    // State Storage    ... of course this part is different for CG
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
    //
    // Recycled subspace vectors.
    Teuchos::RCP<MV> U_;
    // 
    // C = A * U,  linear system is Ax=b 
    Teuchos::RCP<MV> C_;
    //
    // Projected matrices
    // D_ : Diagonal matrix of pivots D = P'AP 
    Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > D_;
  };
  
  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor.
  template<class ScalarType, class MV, class OP>
  PCPGIter<ScalarType,MV,OP>::PCPGIter(const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem, 
					   const Teuchos::RCP<OutputManager<ScalarType> > &printer,
					   const Teuchos::RCP<StatusTest<ScalarType,MV,OP> > &tester,
					   const Teuchos::RCP<MatOrthoManager<ScalarType,MV,OP> > &ortho,
					   Teuchos::ParameterList &params ):
    lp_(problem),
    om_(printer),
    stest_(tester),
    ortho_(ortho),
    savedBlocks_(0),
    initialized_(false),
    stateStorageInitialized_(false),
    keepDiagonal_(false), 
    initDiagonal_(false),
    curDim_(0),
    prevUdim_(0),
    iter_(0)
  {
    // Get the maximum number of blocks allowed for this Krylov subspace

    TEUCHOS_TEST_FOR_EXCEPTION(!params.isParameter("Saved Blocks"), std::invalid_argument,
                       "Belos::PCPGIter::constructor: mandatory parameter \"Saved Blocks\" is not specified.");
    int rb = Teuchos::getParameter<int>(params, "Saved Blocks");

    // Find out whether we are saving the Diagonal matrix.
    keepDiagonal_ = params.get("Keep Diagonal", false);

    // Find out whether we are initializing the Diagonal matrix.
    initDiagonal_ = params.get("Initialize Diagonal", false);

    // Set the number of blocks and allocate data
    setSize( rb );
  }
  
  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Set the block size and adjust as necessary
  template <class ScalarType, class MV, class OP>
  void PCPGIter<ScalarType,MV,OP>::setSize( int savedBlocks )
  {
    // allocate space only; perform no computation
    // Any change in size invalidates the state of the solver as implemented here.

    TEUCHOS_TEST_FOR_EXCEPTION(savedBlocks <= 0, std::invalid_argument, "Belos::PCPGIter::setSize() was passed a non-positive argument for \"Num Saved Blocks\".");

    if ( savedBlocks_ != savedBlocks) {
      stateStorageInitialized_ = false;
      savedBlocks_ = savedBlocks;
      initialized_ = false;
      curDim_ = 0;
      prevUdim_ = 0;
      setStateSize(); // Use the current savedBlocks_ to initialize the state storage.    
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Enable the reuse of a single solver object for completely different linear systems
  template <class ScalarType, class MV, class OP>
  void PCPGIter<ScalarType,MV,OP>::resetState()
  {
      stateStorageInitialized_ = false;
      initialized_ = false;
      curDim_ = 0;
      prevUdim_ = 0;
      setStateSize();
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Setup the state storage.  Called by either initialize or, if savedBlocks_ changes, setSize.
  template <class ScalarType, class MV, class OP>
  void PCPGIter<ScalarType,MV,OP>::setStateSize ()
  {
    if (!stateStorageInitialized_) {

      // Check if there is any multivector to clone from.
      Teuchos::RCP<const MV> lhsMV = lp_->getLHS();
      Teuchos::RCP<const MV> rhsMV = lp_->getRHS();
      if (lhsMV == Teuchos::null && rhsMV == Teuchos::null) {
	return;  // postpone exception 
      }
      else {
	
	//////////////////////////////////
	// blockSize*recycledBlocks dependent
        int newsd = savedBlocks_ ; //int newsd = blockSize_* savedBlocks_ ;
	//
	// Initialize the CG state storage
        // If the subspace is not initialized, generate it using the LHS or RHS from lp_.
	// Generate CG state only if it does not exist, otherwise resize it.
        if (Z_ == Teuchos::null) {
          Teuchos::RCP<const MV> tmp = ( (rhsMV!=Teuchos::null)? rhsMV: lhsMV );
          Z_ = MVT::Clone( *tmp, 1 );
        }
        if (P_ == Teuchos::null) {
          Teuchos::RCP<const MV> tmp = ( (rhsMV!=Teuchos::null)? rhsMV: lhsMV );
          P_ = MVT::Clone( *tmp, 1 );
        }
        if (AP_ == Teuchos::null) {
          Teuchos::RCP<const MV> tmp = ( (rhsMV!=Teuchos::null)? rhsMV: lhsMV );
          AP_ = MVT::Clone( *tmp, 1 );
        }

	if (C_ == Teuchos::null) {        

	  // Get the multivector that is not null. 
	  Teuchos::RCP<const MV> tmp = ( (rhsMV!=Teuchos::null)? rhsMV: lhsMV );
	  TEUCHOS_TEST_FOR_EXCEPTION(tmp == Teuchos::null,std::invalid_argument,
			     "Belos::PCPGIter::setStateSize(): linear problem does not specify multivectors to clone from.");
	  TEUCHOS_TEST_FOR_EXCEPTION( 0 != prevUdim_,std::invalid_argument,
			     "Belos::PCPGIter::setStateSize(): prevUdim not zero and C is null.");
	  C_ = MVT::Clone( *tmp, savedBlocks_ );
	}
	else {
	  // Generate C_ by cloning itself ONLY if more space is needed.
	  if (MVT::GetNumberVecs(*C_) < savedBlocks_ ) {
	    Teuchos::RCP<const MV> tmp = C_;
	    C_ = MVT::Clone( *tmp, savedBlocks_ );
	  }
	}
	if (U_ == Teuchos::null) {        
	  Teuchos::RCP<const MV> tmp = ( (rhsMV!=Teuchos::null)? rhsMV: lhsMV );
	  TEUCHOS_TEST_FOR_EXCEPTION( 0 != prevUdim_,std::invalid_argument,
			     "Belos::PCPGIter::setStateSize(): prevUdim not zero and U is null.");
	  U_ = MVT::Clone( *tmp, savedBlocks_ );
	}
	else {
	  // Generate U_ by cloning itself ONLY if more space is needed.
	  if (MVT::GetNumberVecs(*U_) < savedBlocks_ ) {
	    Teuchos::RCP<const MV> tmp = U_;
	    U_ = MVT::Clone( *tmp, savedBlocks_ );
	  }
	}
        if (keepDiagonal_) {
          if (D_ == Teuchos::null) {
            D_ = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>() );
          }
          if (initDiagonal_) {
            D_->shape( newsd, newsd );
          }
          else {
            if (D_->numRows() < newsd || D_->numCols() < newsd) {
              D_->shapeUninitialized( newsd, newsd );
            }
          }
        }
	// State storage has now been initialized.
	stateStorageInitialized_ = true;
      } // if there is a vector to clone from
    } // if !stateStorageInitialized_ 
  } // end of setStateSize

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Initialize the iteration object
  template <class ScalarType, class MV, class OP>
  void PCPGIter<ScalarType,MV,OP>::initialize(PCPGIterState<ScalarType,MV>& newstate)
  {

    TEUCHOS_TEST_FOR_EXCEPTION(!stateStorageInitialized_,std::invalid_argument,
		       "Belos::PCPGIter::initialize(): Cannot initialize state storage!");
    
    // Requirements: R_ and consistent multivectors widths and lengths
    //
    std::string errstr("Belos::PCPGIter::initialize(): Specified multivectors must have a consistent length and width.");

    if (newstate.R != Teuchos::null){ 

      R_ = newstate.R; // SolverManager::R_ == newstate.R == Iterator::R_
      if (newstate.U == Teuchos::null){ 
        prevUdim_ = 0;
        newstate.U = U_;
        newstate.C = C_;
      }
      else {
        prevUdim_ =  newstate.curDim;
        if (newstate.C == Teuchos::null){  // Stub for new feature
          std::vector<int> index(prevUdim_);
          for (int i=0; i< prevUdim_; ++i)  
            index[i] = i; 
          Teuchos::RCP<const MV> Ukeff = MVT::CloneView( *newstate.U, index );
          newstate.C = MVT::Clone( *newstate.U, prevUdim_ ); 
          Teuchos::RCP<MV> Ckeff = MVT::CloneViewNonConst( *newstate.C, index );    
          lp_->apply( *Ukeff, *Ckeff );
        }
        curDim_ = prevUdim_ ;
      }

      // Initialize the state storage if not already allocated in the constructor
      if (!stateStorageInitialized_) 
        setStateSize();

      //TEUCHOS_TEST_FOR_EXCEPTION( MVT::GetGlobalLength(*newstate.V) != MVT::GetGlobalLength(*V_), std::invalid_argument, errstr );
      //TEUCHOS_TEST_FOR_EXCEPTION( MVT::GetNumberVecs(*newstate.V) < 1, std::invalid_argument, errstr );

      newstate.prevUdim =  prevUdim_; // big change in functionality from GCRODR 
      newstate.curDim =  curDim_; 

      //TEUCHOS_TEST_FOR_EXCEPTION(newstate.z->numRows() < curDim_ || newstate.z->numCols() < 1, std::invalid_argument, errstr);

      std::vector<int> zero_index(1);
      zero_index[0] = 0;
      if ( lp_->getLeftPrec() != Teuchos::null ) { // Compute the initial search direction 
        lp_->applyLeftPrec( *R_, *Z_ );
        MVT::SetBlock( *Z_,  zero_index , *P_ );  // P(:,zero_index) := Z
      } else {                                      
        Z_ = R_;
        MVT::SetBlock( *R_, zero_index, *P_ );
      }

      std::vector<int> nextind(1);
      nextind[0] = curDim_;

      MVT::SetBlock( *P_,  nextind, *newstate.U ); // U(:,curDim_ ) := P_

      ++curDim_;
      newstate.curDim = curDim_; 

      TEUCHOS_TEST_FOR_EXCEPTION( MVT::GetNumberVecs(*newstate.U) != savedBlocks_ ,
                          std::invalid_argument, errstr );
      if (newstate.U != U_) { // Why this is needed?
	U_ = newstate.U;
      }

      TEUCHOS_TEST_FOR_EXCEPTION( MVT::GetNumberVecs(*newstate.C) != savedBlocks_ ,
                          std::invalid_argument, errstr );
      if (newstate.C != C_) {
	C_ = newstate.C;
      }
    }
    else {

      TEUCHOS_TEST_FOR_EXCEPTION(newstate.R == Teuchos::null,std::invalid_argument,
                         "Belos::PCPGIter::initialize(): PCPGStateIterState does not have initial kernel R_0.");
    }

    // the solver is initialized
    initialized_ = true;
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Iterate until the status test informs us we should stop.
  template <class ScalarType, class MV, class OP>
  void PCPGIter<ScalarType,MV,OP>::iterate()
  {
    //
    // Allocate/initialize data structures
    //
    if (initialized_ == false) {
      initialize();
    }
    const bool debug = false;

    // Allocate memory for scalars.
    Teuchos::SerialDenseMatrix<int,ScalarType> alpha( 1, 1 );
    Teuchos::SerialDenseMatrix<int,ScalarType> pAp( 1, 1 );
    Teuchos::SerialDenseMatrix<int,ScalarType> beta( 1, 1 );
    Teuchos::SerialDenseMatrix<int,ScalarType> rHz( 1, 1 ), rHz_old( 1, 1 );

    if( iter_ != 0 )
      std::cout << " Iterate Warning: begin from nonzero iter_ ?" << std::endl;  //DMD

    // GenOrtho Project Stubs
    std::vector<int> prevInd;
    Teuchos::RCP<const MV> Uprev;
    Teuchos::RCP<const MV> Cprev;
    Teuchos::SerialDenseMatrix<int,ScalarType> CZ;

    if( prevUdim_ ){
      prevInd.resize( prevUdim_ );
      for( int i=0; i<prevUdim_ ; i++) prevInd[i] = i;
      CZ.reshape( prevUdim_ , 1 );
      Uprev = MVT::CloneView(*U_, prevInd);
      Cprev = MVT::CloneView(*C_, prevInd);
    }

    // Get the current solution std::vector.
    Teuchos::RCP<MV> cur_soln_vec = lp_->getCurrLHSVec();

    // Check that the current solution std::vector only has one column.
    TEUCHOS_TEST_FOR_EXCEPTION( MVT::GetNumberVecs(*cur_soln_vec) != 1, CGIterationInitFailure,
                        "Belos::PCPGIter::iterate(): current linear system has more than one std::vector!" );

    //Check that the input is correctly set up 
    TEUCHOS_TEST_FOR_EXCEPTION( curDim_  != prevUdim_ + 1, CGIterationInitFailure,
                        "Belos::PCPGIter::iterate(): mistake in initialization !" );


    const ScalarType zero = Teuchos::ScalarTraits<ScalarType>::zero();
    const ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();


    std::vector<int> curind(1);
    std::vector<ScalarType> rnorm(MVT::GetNumberVecs(*cur_soln_vec));
    if (prevUdim_ > 0){                 // A-orthonalize P=Z to Uprev
      Teuchos::RCP<MV> P; 
      curind[0] = curDim_ - 1;          // column = dimension - 1 
      P = MVT::CloneViewNonConst(*U_,curind); 
      MVT::MvTransMv( one, *Cprev, *P, CZ );
      MVT::MvTimesMatAddMv( -one, *Uprev, CZ, one, *P );       // P -= U*(C'Z)

      if( debug ){
        MVT::MvTransMv( one, *Cprev, *P, CZ );
        std::cout << " Input CZ post ortho " << std::endl;
        CZ.print( std::cout );
      }
      if( curDim_ == savedBlocks_ ){
        std::vector<int> zero_index(1);
        zero_index[0] = 0;
        MVT::SetBlock( *P, zero_index, *P_ );
      }
      P = Teuchos::null;
    }

    // Compute first <r,z> a.k.a. rHz
    MVT::MvTransMv( one, *R_, *Z_, rHz );

    ////////////////////////////////////////////////////////////////
    // iterate until the status test is satisfied
    //
    while (stest_->checkStatus(this) != Passed ) {
      Teuchos::RCP<const MV> P; 
      Teuchos::RCP<MV> AP;
      iter_++;                          // The next iteration begins.
      //std::vector<int> curind(1);
      curind[0] = curDim_ - 1;          // column = dimension - 1 
      if( debug ){
        MVT::MvNorm(*R_, rnorm);
        std::cout << iter_ << "  " << curDim_ <<  "   " << rnorm[0] << std::endl;
      }
      if( prevUdim_ + iter_ < savedBlocks_ ){
        P = MVT::CloneView(*U_,curind); 
        AP = MVT::CloneViewNonConst(*C_,curind); 
        lp_->applyOp( *P, *AP );
        MVT::MvTransMv( one, *P, *AP, pAp );
      }else{
        if( prevUdim_ + iter_ == savedBlocks_ ){
          AP = MVT::CloneViewNonConst(*C_,curind); 
          lp_->applyOp( *P_, *AP );
          MVT::MvTransMv( one, *P_, *AP, pAp );
        }else{
          lp_->applyOp( *P_, *AP_ );
          MVT::MvTransMv( one, *P_, *AP_, pAp );
        }
      }

      if( keepDiagonal_  && prevUdim_ + iter_ <= savedBlocks_ )
        (*D_)(iter_ -1 ,iter_ -1 ) = pAp(0,0);

      // positive pAp required 
      TEUCHOS_TEST_FOR_EXCEPTION( pAp(0,0) <= zero, CGPositiveDefiniteFailure,
                          "Belos::PCPGIter::iterate(): non-positive value for p^H*A*p encountered!" );

      // alpha := <R_,Z_> / <P,AP>
      alpha(0,0) = rHz(0,0) / pAp(0,0);

      // positive alpha required 
      TEUCHOS_TEST_FOR_EXCEPTION( alpha(0,0) <= zero, CGPositiveDefiniteFailure,
                          "Belos::PCPGIter::iterate(): non-positive value for alpha encountered!" );

      // solution update  x += alpha * P
      if( curDim_ < savedBlocks_ ){
         MVT::MvAddMv( one, *cur_soln_vec, alpha(0,0), *P, *cur_soln_vec );
      }else{
         MVT::MvAddMv( one, *cur_soln_vec, alpha(0,0), *P_, *cur_soln_vec );
      }
      //lp_->updateSolution(); ... does nothing.
      //
      // The denominator of beta is saved before residual is updated [ old <R_, Z_> ].
      //
      rHz_old(0,0) = rHz(0,0);
      //
      // residual update R_ := R_ - alpha * AP
      //
      if( prevUdim_ + iter_ <= savedBlocks_ ){
         MVT::MvAddMv( one, *R_, -alpha(0,0), *AP, *R_ );
         AP = Teuchos::null;
      }else{
         MVT::MvAddMv( one, *R_, -alpha(0,0), *AP_, *R_ );
      }
      //
      // update beta := [ new <R_, Z_> ] / [ old <R_, Z_> ] and the search direction p.
      //
      if ( lp_->getLeftPrec() != Teuchos::null ) {
        lp_->applyLeftPrec( *R_, *Z_ );
      } else {
        Z_ = R_;
      }
      //
      MVT::MvTransMv( one, *R_, *Z_, rHz );
      //
      beta(0,0) = rHz(0,0) / rHz_old(0,0);
      //
      if( curDim_ < savedBlocks_ ){
         curDim_++;                                                         // update basis dim
         curind[0] = curDim_ - 1;
         Teuchos::RCP<MV> Pnext = MVT::CloneViewNonConst(*U_,curind);
         MVT::MvAddMv( one, *Z_, beta(0,0), *P, *Pnext );
         if( prevUdim_ ){ // Deflate seed space 
             MVT::MvTransMv( one, *Cprev, *Z_, CZ );
             MVT::MvTimesMatAddMv( -one, *Uprev, CZ, one, *Pnext ); // Pnext -= U*(C'Z)
             if( debug ){
               std::cout << " Check CZ " << std::endl;
               MVT::MvTransMv( one, *Cprev, *Pnext, CZ );
               CZ.print( std::cout );
             }
         }
         P = Teuchos::null;
         if( curDim_ == savedBlocks_ ){
           std::vector<int> zero_index(1);
           zero_index[0] = 0;
           MVT::SetBlock( *Pnext, zero_index, *P_ );
         }
         Pnext = Teuchos::null;
      }else{
         MVT::MvAddMv( one, *Z_, beta(0,0), *P_, *P_ );
         if( prevUdim_ ){ // Deflate seed space
             MVT::MvTransMv( one, *Cprev, *Z_, CZ );
             MVT::MvTimesMatAddMv( -one, *Uprev, CZ, one, *P_ );       // P_ -= U*(C'Z)

             if( debug ){
               std::cout << " Check CZ " << std::endl;
               MVT::MvTransMv( one, *Cprev, *P_, CZ );
               CZ.print( std::cout );
             }
         }
      }
      // CGB: 5/26/2010
      // this RCP<const MV> P was previously a variable outside the loop. however, it didn't appear to be see any use between
      // loop iterations. therefore, I moved it inside to avoid scoping errors with previously used variables named P.
      // to ensure that this wasn't a bug, I verify below that we have set P == null, i.e., that we are not going to use it again
      // same for AP
      TEUCHOS_TEST_FOR_EXCEPTION( AP != Teuchos::null || P != Teuchos::null, std::logic_error, "Loop recurrence violated. Please contact Belos team.");
    } // end coupled two-term recursion
    if( prevUdim_ + iter_ < savedBlocks_ ) --curDim_; // discard negligible search direction
  }

} // end Belos namespace

#endif /* BELOS_PCPG_ITER_HPP */
