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

#ifndef BELOS_PSEUDO_BLOCK_GMRES_ITER_HPP
#define BELOS_PSEUDO_BLOCK_GMRES_ITER_HPP

/*! \file BelosBlockGmresIter.hpp
    \brief Belos concrete class for performing the pseudo-block GMRES iteration.
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
  \class Belos::PseudoBlockGmresIter
  
  \brief This class implements the pseudo-block GMRES iteration, where a
  block Krylov subspace is constructed for all of the linear systems simultaneously.  
  The QR decomposition of each block, upper Hessenberg matrix is performed each iteration 
  to update the least squares system and give the current linear system residuals.
 
  \author Heidi Thornquist
*/

namespace Belos {
  
  //! @name PseudoBlockGmresIter Structures 
  //@{ 
  
  /** \brief Structure to contain pointers to PseudoBlockGmresIter state variables.
   *
   * This struct is utilized by PseudoBlockGmresIter::initialize() and PseudoBlockGmresIter::getState().
   */
  template <class ScalarType, class MV>
  struct PseudoBlockGmresIterState {
    /*! \brief The current dimension of the reduction.
     *
     * This should always be equal to PseudoBlockGmresIter::getCurSubspaceDim()
     */
    int curDim;
    /*! \brief The current Krylov basis. */
    std::vector<Teuchos::RefCountPtr<const MV> > V;
    /*! \brief The current Hessenberg matrix. 
     *
     * The \c curDim by \c curDim leading submatrix of H is the 
     * projection of problem->getOperator() by the first \c curDim vectors in V. 
     */
    std::vector<Teuchos::RefCountPtr<const Teuchos::SerialDenseMatrix<int,ScalarType> > > H;
    /*! \brief The current upper-triangular matrix from the QR reduction of H. */
    std::vector<Teuchos::RefCountPtr<const Teuchos::SerialDenseMatrix<int,ScalarType> > > R;
    /*! \brief The current right-hand side of the least squares system RY = Z. */
    std::vector<Teuchos::RefCountPtr<const Teuchos::SerialDenseVector<int,ScalarType> > > Z;
    
    PseudoBlockGmresIterState() : curDim(0), V(0),
				  H(0), R(0), Z(0)
    {}
  };
  
  //! @name PseudoBlockGmresIter Exceptions
  //@{ 
  
  /** \brief PseudoBlockGmresIterInitFailure is thrown when the PseudoBlockGmresIter object is unable to
   * generate an initial iterate in the PseudoBlockGmresIter::initialize() routine. 
   *
   * This exception is thrown from the PseudoBlockGmresIter::initialize() method, which is
   * called by the user or from the PseudoBlockGmresIter::iterate() method if isInitialized()
   * == \c false.
   *
   * In the case that this exception is thrown, 
   * PseudoBlockGmresIter::isInitialized() will be \c false and the user will need to provide
   * a new initial iterate to the iteration.
   *
   */
  class PseudoBlockGmresIterInitFailure : public BelosError {public:
    PseudoBlockGmresIterInitFailure(const std::string& what_arg) : BelosError(what_arg)
    {}};
  
  /** \brief PseudoBlockGmresIterOrthoFailure is thrown when the orthogonalization manager is
   * unable to generate orthonormal columns from the new basis vectors.
   *
   * This exception is thrown from the PseudoBlockGmresIter::iterate() method.
   *
   */
  class PseudoBlockGmresIterOrthoFailure : public BelosError {public:
    PseudoBlockGmresIterOrthoFailure(const std::string& what_arg) : BelosError(what_arg)
    {}};
  
  //@}
  
  
  template<class ScalarType, class MV, class OP>
  class PseudoBlockGmresIter : virtual public Iteration<ScalarType,MV,OP> {
    
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
    
    /*! \brief %PseudoBlockGmresIter constructor with linear problem, solver utilities, and parameter list of solver options.
     *
     * This constructor takes pointers required by the linear solver, in addition
     * to a parameter list of options for the linear solver. These options include the following:
     *   - "Block Size" - an \c int specifying the block size used by the algorithm. This can also be specified using the setBlockSize() method. Default: 1
     *   - "Num Blocks" - an \c int specifying the maximum number of blocks allocated for the solver basis. Default: 25
     *   - "Restart Timers" = a \c bool specifying whether the timers should be restarted each time iterate() is called. Default: false
     */
    PseudoBlockGmresIter( const Teuchos::RefCountPtr<LinearProblem<ScalarType,MV,OP> > &problem, 
			  const Teuchos::RefCountPtr<OutputManager<ScalarType> > &printer,
			  const Teuchos::RefCountPtr<StatusTest<ScalarType,MV,OP> > &tester,
			  const Teuchos::RefCountPtr<MatOrthoManager<ScalarType,MV,OP> > &ortho,
			  Teuchos::ParameterList &params );
    
    //! Destructor.
    virtual ~PseudoBlockGmresIter() {};
    //@}
    
    
    //! @name Solver methods
    //@{ 
    
    /*! \brief This method performs block Gmres iterations until the status
     * test indicates the need to stop or an error occurs (in which case, an
     * exception is thrown).
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
     * Possible exceptions thrown include the PseudoBlockGmresIterOrthoFailure.
     *
     */
    void iterate();
    
    /*! \brief Initialize the solver to an iterate, providing a complete state.
     *
     * The %PseudoBlockGmresIter contains a certain amount of state, consisting of the current 
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
    void initialize(PseudoBlockGmresIterState<ScalarType,MV> newstate);
    
    /*! \brief Initialize the solver with the initial vectors from the linear problem
     *  or random data.
     */
    void initialize()
    {
      PseudoBlockGmresIterState<ScalarType,MV> empty;
      initialize(empty);
    }
    
    /*! \brief Get the current state of the linear solver.
     *
     * The data is only valid if isInitialized() == \c true.
     *
     * \returns A PseudoBlockGmresIterState object containing const pointers to the current
     * solver state.
     */
    PseudoBlockGmresIterState<ScalarType,MV> getState() const {
      PseudoBlockGmresIterState<ScalarType,MV> state;
      state.curDim = curDim_;
      state.V.resize(numRHS_);
      state.H.resize(numRHS_);
      state.Z.resize(numRHS_);
      for (int i=0; i<numRHS_; ++i) {
	state.V[i] = V_[i];
	state.H[i] = H_[i];
	state.Z[i] = Z_[i];
      }
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
    /*! \note Some solvers, like GMRES, do not compute updates to the solution every iteration.
      This method forces its computation.  Other solvers, like CG, update the solution
      each iteration, so this method will return a zero vector indicating that the linear
      problem contains the current solution.
    */
    Teuchos::RefCountPtr<MV> getCurrentUpdate() const;
    
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
    
    //! Get the blocksize to be used by the iterative solver in solving this linear problem.
    int getBlockSize() const { return 1; }
    
    //! \brief Set the blocksize.
    void setBlockSize(int blockSize) { 
      TEST_FOR_EXCEPTION(blockSize!=1,std::invalid_argument,
			 "Belos::CGIter::setBlockSize(): Cannot use a block size that is not one.");
    }
    
    //! Get the maximum number of blocks used by the iterative solver in solving this linear problem.
    int getNumBlocks() const { return numBlocks_; }
    
    //! \brief Set the maximum number of blocks used by the iterative solver.
    void setNumBlocks(int numBlocks);
    
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
    const Teuchos::RefCountPtr<OrthoManager<ScalarType,MV> >        ortho_;
    
    //
    // Algorithmic parameters
    //  
    // numRHS_ is the current number of linear systems being solved.
    int numRHS_;
    // numBlocks_ is the size of the allocated space for the Krylov basis, in blocks.
    int numBlocks_; 
    
    // Storage for QR factorization of the least squares system.
    Teuchos::Array<Teuchos::SerialDenseVector<int,ScalarType> > sn;
    Teuchos::Array<Teuchos::SerialDenseVector<int,MagnitudeType> > cs;
    
    // Pointers to a work vector used to improve aggregate performance.
    RefCountPtr<MV> U_vec_, AU_vec_;    

    // Pointers to the current right-hand side and solution multivecs being solved for.
    RefCountPtr<MV> cur_block_rhs_, cur_block_sol_;

    // 
    // Current solver state
    //
    // initialized_ specifies that the basis vectors have been initialized and the iterate() routine
    // is capable of running; _initialize is controlled  by the initialize() member method
    // For the implications of the state of initialized_, please see documentation for initialize()
    bool initialized_;
    
    // stateStorageInitialized_ specified that the state storage has be initialized to the current
    // blockSize_ and numBlocks_.  This initialization may be postponed if the linear problem was
    // generated without the right-hand side or solution vectors.
    bool stateStorageInitialized_;
    
    // Current subspace dimension, and number of iterations performed.
    int curDim_, iter_;

    // 
    // State Storage
    //
    std::vector<Teuchos::RefCountPtr<MV> > V_;
    //
    // Projected matrices
    // H_ : Projected matrix from the Krylov factorization AV = VH + FE^T
    //
    std::vector<Teuchos::RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > > H_;
    // 
    // QR decomposition of Projected matrices for solving the least squares system HY = B.
    // R_: Upper triangular reduction of H
    // Z_: Q applied to right-hand side of the least squares system
    std::vector<Teuchos::RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > > R_;
    std::vector<Teuchos::RefCountPtr<Teuchos::SerialDenseVector<int,ScalarType> > > Z_;  
  };
  
  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor.
  template<class ScalarType, class MV, class OP>
  PseudoBlockGmresIter<ScalarType,MV,OP>::PseudoBlockGmresIter(const Teuchos::RefCountPtr<LinearProblem<ScalarType,MV,OP> > &problem, 
							       const Teuchos::RefCountPtr<OutputManager<ScalarType> > &printer,
							       const Teuchos::RefCountPtr<StatusTest<ScalarType,MV,OP> > &tester,
							       const Teuchos::RefCountPtr<MatOrthoManager<ScalarType,MV,OP> > &ortho,
							       Teuchos::ParameterList &params ):
    lp_(problem),
    om_(printer),
    stest_(tester),
    ortho_(ortho),
    numRHS_(0),
    numBlocks_(0),
    initialized_(false),
    stateStorageInitialized_(false),
    curDim_(0),
    iter_(0)
  {
    // Get the maximum number of blocks allowed for each Krylov subspace
    TEST_FOR_EXCEPTION(!params.isParameter("Num Blocks"), std::invalid_argument,
                       "Belos::PseudoBlockGmresIter::constructor: mandatory parameter 'Num Blocks' is not specified.");
    int nb = Teuchos::getParameter<int>(params, "Num Blocks");

    setNumBlocks( nb );
  }
  
  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Set the block size and make necessary adjustments.
  template <class ScalarType, class MV, class OP>
  void PseudoBlockGmresIter<ScalarType,MV,OP>::setNumBlocks (int numBlocks)
  {
    // This routine only allocates space; it doesn't not perform any computation
    // any change in size will invalidate the state of the solver.
    
    TEST_FOR_EXCEPTION(numBlocks <= 0, std::invalid_argument, "Belos::PseudoBlockGmresIter::setNumBlocks was passed a non-positive argument.");

    numBlocks_ = numBlocks;
    curDim_ = 0;

    initialized_ = false;
    stateStorageInitialized_ = false;
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Setup the state storage.
  template <class ScalarType, class MV, class OP>
  void PseudoBlockGmresIter<ScalarType,MV,OP>::setStateSize ()
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
	
	//////////////////////////////////
	// Reinitialize storage for least squares solve
	//
	cs.resize(numRHS_);
	sn.resize(numRHS_);
	for (int i=0; i<numRHS_; ++i) {
	  cs[i].resize(numBlocks_+1);
	  sn[i].resize(numBlocks_+1);
	}
	
	// Get the multivector that is not null.
	Teuchos::RefCountPtr<const MV> tmp = ( (rhsMV!=Teuchos::null)? rhsMV: lhsMV );
	TEST_FOR_EXCEPTION(tmp == Teuchos::null,std::invalid_argument,
			   "Belos::PseudoBlockGmresIter::setStateSize(): linear problem does not specify multivectors to clone from.");

	// Initialize the state storage
	// If the subspace has not be initialized before, generate it using the LHS or RHS from lp_.
	if (V_.size() == 0) {
          V_.resize(numRHS_);
	  for (int i=0; i<numRHS_; ++i) {
	    V_[i] = MVT::Clone(*tmp,numBlocks_+1);
	  }
	}
	else {
	  if ( ((int)V_.size() < numRHS_) || (MVT::GetNumberVecs(*V_[0]) < numBlocks_+1) ) {
	    for (int i=0; i<numRHS_; ++i) {
	      V_[i] = MVT::Clone(*tmp,numBlocks_+1);
	    }
	  }
	}
	//
	// Create the rectangular Hessenberg matrix and right-hand side of least squares problem.
	//
	if ((int)H_.size() != numRHS_) {
	  H_.resize(numRHS_);
	  Z_.resize(numRHS_);
	  for (int i=0; i<numRHS_; ++i) {
	    if (H_[i] == Teuchos::null)
	      H_[i] = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>() );
	    if (Z_[i] == Teuchos::null)
	      Z_[i] = Teuchos::rcp( new Teuchos::SerialDenseVector<int,ScalarType>() );
	  }
	}
	if ( (int)Z_[0]->length() != numBlocks_+1 ) {
	  for (int i=0; i<numRHS_; ++i) {
	    H_[i]->shapeUninitialized(numBlocks_+1, numBlocks_);
	    Z_[i]->shapeUninitialized(numBlocks_+1, 1); 
	  }
	}

	// State storage has now been initialized.
	stateStorageInitialized_ = true;
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Get the current update from this subspace.
  template <class ScalarType, class MV, class OP>
  Teuchos::RefCountPtr<MV> PseudoBlockGmresIter<ScalarType,MV,OP>::getCurrentUpdate() const
  {
    //
    // If this is the first iteration of the Arnoldi factorization, 
    // there is no update, so return Teuchos::null. 
    //
    RefCountPtr<MV> currentUpdate = Teuchos::null;
    if (curDim_==0) { 
      return currentUpdate; 
    } else {
      currentUpdate = MVT::CloneCopy(*(lp_->getCurrLHSVec()));
      std::vector<int> index(1), index2(curDim_);
      for (int i=0; i<curDim_; ++i) {
        index2[i] = i;
      }
      const ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();
      Teuchos::BLAS<int,ScalarType> blas;
      
      for (int i=0; i<numRHS_; ++i) {
        index[0] = i;
        RefCountPtr<MV> cur_block_copy_vec = MVT::CloneView( *currentUpdate, index );
        //
        //  Make a view and then copy the RHS of the least squares problem.  DON'T OVERWRITE IT!
        //
        Teuchos::SerialDenseVector<int,ScalarType> y( Teuchos::Copy, Z_[i]->values(), curDim_ );
        //
        //  Solve the least squares problem and compute current solutions.
        //
        blas.TRSM( Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
	           Teuchos::NON_UNIT_DIAG, curDim_, 1, one,  
		   H_[i]->values(), H_[i]->stride(), y.values(), y.stride() );
	
	RefCountPtr<const MV> Vjp1 = MVT::CloneView( *V_[i], index2 );
	MVT::MvTimesMatAddMv( one, *Vjp1, y, one, *cur_block_copy_vec );
      }
    }
    return currentUpdate;
  }
  

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Get the native residuals stored in this iteration.  
  // Note:  No residual vector will be returned by Gmres.
  template <class ScalarType, class MV, class OP>
  Teuchos::RefCountPtr<const MV> PseudoBlockGmresIter<ScalarType,MV,OP>::getNativeResiduals( std::vector<MagnitudeType> *norms ) const 
  {
    //
    // NOTE: Make sure the incoming vector is the correct size!
    //
    if ( norms && (int)norms->size() < numRHS_ )                         
      norms->resize( numRHS_ );                                          
    
    if (norms) {
      Teuchos::BLAS<int,ScalarType> blas;
      for (int j=0; j<numRHS_; j++) {
        (*norms)[j] = blas.NRM2( numRHS_, &(*Z_[j])(curDim_), 1);
      }
    }
    return Teuchos::null;
  }

  

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Initialize this iteration object
  template <class ScalarType, class MV, class OP>
  void PseudoBlockGmresIter<ScalarType,MV,OP>::initialize(PseudoBlockGmresIterState<ScalarType,MV> newstate)
  {
    // Get the number of right-hand sides we're solving for now.
    int numRHS = MVT::GetNumberVecs(*(lp_->getCurrLHSVec()));
    if (numRHS != numRHS_)
      stateStorageInitialized_ = false;
    numRHS_ = numRHS;

    // Initialize the state storage if it isn't already.
    if (!stateStorageInitialized_) 
      setStateSize();

      TEST_FOR_EXCEPTION(!stateStorageInitialized_,std::invalid_argument,
			 "Belos::PseudoBlockGmresIter::initialize(): Cannot initialize state storage!");      
          
    // NOTE:  In PseudoBlockGmresIter, V and Z are required!!!  
    // inconsitent multivectors widths and lengths will not be tolerated, and
    // will be treated with exceptions.
    //
    std::string errstr("Belos::PseudoBlockGmresIter::initialize(): Specified multivectors must have a consistent length and width.");

    if ((int)newstate.V.size() != 0 && (int)newstate.Z.size() != 0) {

      // initialize V_,Z_, and curDim_
      
      TEST_FOR_EXCEPTION( MVT::GetVecLength(*newstate.V[0]) != MVT::GetVecLength(*V_[0]),
                          std::invalid_argument, errstr );
      TEST_FOR_EXCEPTION( MVT::GetNumberVecs(*newstate.V[0]) < numRHS_,
                          std::invalid_argument, errstr );
      TEST_FOR_EXCEPTION( newstate.curDim > numBlocks_+1,
                          std::invalid_argument, errstr );

      curDim_ = newstate.curDim;

      // check size of Z
      for (int i=0; i<numRHS_; ++i)
	TEST_FOR_EXCEPTION(newstate.Z[i]->numRows() < curDim_, std::invalid_argument, errstr);
      

      // copy basis vectors from newstate into V
      for (int i=0; i<numRHS_; ++i) {
	int lclDim = MVT::GetNumberVecs(*newstate.V[i]);
	
	if (newstate.V[i] != V_[i]) {
	  // only copy over the first block and print a warning.
	  if (curDim_ == 0 && lclDim > 1) {
	    om_->stream(Warnings) << "Belos::PseudoBlockGmresIter::initialize(): the solver was initialized with a kernel of " 
				  << lclDim << endl << "The block size however is only " << 1 << endl
				  << "The last " << lclDim - 1 << " vectors will be discarded." << endl;
	  }
	  std::vector<int> nevind(curDim_+1);
	  for (int j=0; j<curDim_+1; j++) nevind[j] = j;
	  Teuchos::RefCountPtr<const MV> newV = MVT::CloneView( *newstate.V[i], nevind );
	  Teuchos::RefCountPtr<MV> lclV = MVT::CloneView( *V_[i], nevind );
	  MVT::MvAddMv( 1.0, *newV, 0.0, *newV, *lclV );
	  
	  // done with local pointers
	  lclV = Teuchos::null;
	}
	
	// put data into Z_, make sure old information is not still hanging around.
	if (newstate.Z[i] != Z_[i]) {
	  Z_[i]->putScalar();
	  Teuchos::SerialDenseVector<int,ScalarType> newZ(Teuchos::View,newstate.Z[i]->values(),curDim_+1);
	  Teuchos::RefCountPtr<Teuchos::SerialDenseVector<int,ScalarType> > lclZ;
	  lclZ = Teuchos::rcp( new Teuchos::SerialDenseVector<int,ScalarType>(Teuchos::View,Z_[i]->values(),curDim_+1) );
	  lclZ->assign(newZ);
	  
	  // done with local pointers
	  lclZ = Teuchos::null;
	}
      }
    }
    else {
      
      TEST_FOR_EXCEPTION(newstate.V.size() == 0, std::invalid_argument,
                         "Belos::PseudoBlockGmresIter::initialize(): BlockGmresStateIterState does not have initial kernel V_0.");

      TEST_FOR_EXCEPTION(newstate.Z.size() == 0, std::invalid_argument,
                         "Belos::PseudoBlockGmresIter::initialize(): BlockGmresStateIterState does not have initial norms Z_0.");
    }
    
    // the solver is initialized
    initialized_ = true;

    /*
    if (om_->isVerbosity( Debug ) ) {
      // Check almost everything here
      CheckList chk;
      chk.checkV = true;
      chk.checkArn = true;
      chk.checkAux = true;
      om_->print( Debug, accuracyCheck(chk, ": after initialize()") );
    }
    */

  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Iterate until the status test informs us we should stop.
  template <class ScalarType, class MV, class OP>
  void PseudoBlockGmresIter<ScalarType,MV,OP>::iterate()
  {
    //
    // Allocate/initialize data structures
    //
    if (initialized_ == false) {
      initialize();
    }

    const ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();
    const ScalarType zero = Teuchos::ScalarTraits<ScalarType>::zero();
    
    // Compute the current search dimension. 
    int searchDim = numBlocks_;
    //
    // Associate each initial block of V_[i] with U_vec[i]
    // Reset the index vector (this might have been changed if there was a restart)
    //
    std::vector<int> index(1);
    std::vector<int> index2(1);
    index[0] = 0;
    Teuchos::RefCountPtr<MV> tmp_vec;
    Teuchos::RefCountPtr<MV> U_vec = MVT::Clone( *V_[0], numRHS_ );

    // Create AU_vec to hold A*U_vec.
    Teuchos::RefCountPtr<MV> AU_vec = MVT::Clone( *V_[0], numRHS_ );

    for (int i=0; i<numRHS_; ++i) {
      index2[0] = i;
      tmp_vec = MVT::CloneView( *V_[i], index );
      MVT::MvAddMv( one, *tmp_vec, zero, *tmp_vec, *MVT::CloneView( *U_vec, index2 ) );
    }
    
    ////////////////////////////////////////////////////////////////
    // iterate until the status test tells us to stop.
    //
    // also break if our basis is full
    //
    while (stest_->checkStatus(this) != Passed && curDim_+1 <= searchDim) {

      iter_++;

      //
      // Apply the operator to _work_vector
      //
      lp_->apply( *U_vec, *AU_vec );
      //
      //
      // Resize index.
      //
      int num_prev = curDim_+1;
      index.resize( num_prev );
      for (int i=0; i<num_prev; ++i) { 
	index[i] = i; 
      }
      //
      // Orthogonalize next Krylov vector for each right-hand side.
      //
      for (int i=0; i<numRHS_; ++i) {
	//
	// Get previous Krylov vectors.
	//
	RefCountPtr<MV> V_prev = MVT::CloneView( *V_[i], index );
	Teuchos::Array< RefCountPtr<const MV> > V_array( 1, V_prev );
	//
	// Get a view of the new candidate vector.
	//
	index2[0] = i;
	RefCountPtr<MV> V_new = MVT::CloneView( *AU_vec, index2 );
	//
	// Get a view of the current part of the upper-hessenberg matrix.
	//
	RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > h_new 
	  = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>
			  ( Teuchos::View, *H_[i], num_prev, 1, 0, curDim_ ) );
	Teuchos::Array< RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > > h_array( 1, h_new );
	
	RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > r_new
	  = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>
			  ( Teuchos::View, *H_[i], 1, 1, num_prev, curDim_ ) );
	//
	// Orthonormalize the new block of the Krylov expansion
	// NOTE:  Rank deficiencies are not checked because this is a single-vector Krylov method.
	//
	ortho_->projectAndNormalize( *V_new, h_array, r_new, V_array );
	//
	// NOTE:  V_new is a copy of the iter+1 vector in V_[i], so the normalized vector has to be
	// be copied back in when V_new is changed.  
	//
	index2[0] = curDim_+1;
	tmp_vec = MVT::CloneView( *V_[i], index2 );
	MVT::MvAddMv( one, *V_new, zero, *V_new, *tmp_vec );
      }
      // 
      // Now _AU_vec is the new _U_vec, so swap these two vectors.
      // NOTE: This alleviates the need for allocating a vector for AU_vec each iteration.
      // 
      RefCountPtr<MV> tmp_AU_vec = U_vec;
      U_vec = AU_vec;
      AU_vec = tmp_AU_vec;
      //
      // V has been extended, and H has been extended. 
      //
      // Update the QR factorization of the upper Hessenberg matrix
      //
      updateLSQR();
      //
      // Update basis dim and release all pointers.
      //
      curDim_ += 1;
      //        
      /*      
      // When required, monitor some orthogonalities
      if (om_->isVerbosity( Debug ) ) {
      // Check almost everything here
      CheckList chk;
      chk.checkV = true;
      chk.checkArn = true;
      om_->print( Debug, accuracyCheck(chk, ": after local update") );
      }
      else if (om_->isVerbosity( OrthoDetails ) ) {
        CheckList chk;
        chk.checkV = true;
        om_->print( OrthoDetails, accuracyCheck(chk, ": after local update") );
      }
      */ 
      
    } // end while (statusTest == false)
   
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Update the least squares solution for each right-hand side.
  template<class ScalarType, class MV, class OP>
  void PseudoBlockGmresIter<ScalarType,MV,OP>::updateLSQR( int dim )
  {
    // Get correct dimension based on input "dim"
    // Remember that ortho failures result in an exit before updateLSQR() is called.
    // Therefore, it is possible that dim == curDim_.
    int curDim = curDim_;
    if (dim >= curDim_ && dim < getMaxSubspaceDim()) {
      curDim = dim;
    }
    
    int i, j;
    const ScalarType zero = Teuchos::ScalarTraits<ScalarType>::zero();

    Teuchos::LAPACK<int, ScalarType> lapack;
    Teuchos::BLAS<int, ScalarType> blas;
    
    for (i=0; i<numRHS_; ++i) {
      //
      // Update the least-squares QR for each linear system.
      //
      // QR factorization of Least-Squares system with Givens rotations
      //
      for (j=0; j<curDim_; j++) {
	//
	// Apply previous Givens rotations to new column of Hessenberg matrix
	//
	blas.ROT( 1, &(*H_[i])(j,curDim_), 1, &(*H_[i])(j+1, curDim_), 1, &(cs[i])[j], &(sn[i])[j] );
      }
      //
      // Calculate new Givens rotation
      //
      blas.ROTG( &(*H_[i])(curDim_,curDim_), &(*H_[i])(curDim_+1,curDim_), &(cs[i])[curDim_], &(sn[i])[curDim_] );
      (*H_[i])(curDim_+1,curDim_) = zero;
      //
      // Update RHS w/ new transformation
      //
      blas.ROT( 1, &(*Z_[i])(curDim_), 1, &(*Z_[i])(curDim_+1), 1, &cs[i][curDim_], &sn[i][curDim_] );
    }

  } // end updateLSQR()

} // end Belos namespace

#endif /* BELOS_PSEUDO_BLOCK_GMRES_ITER_HPP */
