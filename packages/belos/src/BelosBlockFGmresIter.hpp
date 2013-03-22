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

#ifndef BELOS_BLOCK_FGMRES_ITER_HPP
#define BELOS_BLOCK_FGMRES_ITER_HPP

/*! \file BelosBlockFGmresIter.hpp
    \brief Belos concrete class for performing the block, flexible GMRES iteration.
*/

#include "BelosConfigDefs.hpp"
#include "BelosTypes.hpp"
#include "BelosGmresIteration.hpp"

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
  \class Belos::BlockFGmresIter
  
  \brief This class implements the block flexible GMRES iteration, where a
  block Krylov subspace is constructed.  The QR decomposition of 
  block, upper Hessenberg matrix is performed each iteration to update
  the least squares system and give the current linear system residuals.

  \ingroup belos_solver_framework
 
  \author Teri Barth and Heidi Thornquist
*/

namespace Belos {
  
template<class ScalarType, class MV, class OP>
class BlockFGmresIter : virtual public GmresIteration<ScalarType,MV,OP> {

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

  /*! \brief %BlockFGmresIter constructor with linear problem, solver utilities, and parameter list of solver options.
   *
   * This constructor takes pointers required by the linear solver, in addition
   * to a parameter list of options for the linear solver. These options include the following:
   *   - "Block Size" - an \c int specifying the block size used by the algorithm. This can also be specified using the setBlockSize() method. Default: 1
   *   - "Num Blocks" - an \c int specifying the maximum number of blocks allocated for the solver basis. Default: 25
   *   - "Restart Timers" = a \c bool specifying whether the timers should be restarted each time iterate() is called. Default: false
   *   - "Keep Hessenberg" = a \c bool specifying whether the upper Hessenberg should be stored separately from the least squares system. Default: false
   */
  BlockFGmresIter( const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem, 
		   const Teuchos::RCP<OutputManager<ScalarType> > &printer,
		   const Teuchos::RCP<StatusTest<ScalarType,MV,OP> > &tester,
		   const Teuchos::RCP<MatOrthoManager<ScalarType,MV,OP> > &ortho,
		   Teuchos::ParameterList &params );
  
  //! Destructor.
  virtual ~BlockFGmresIter() {};
  //@}
  
  
  //! @name Solver methods
  //@{ 
  
  /*! \brief This method performs block FGmres iterations until the status
   * test indicates the need to stop or an error occurs (in which case, an
   * std::exception is thrown).
   *
   * iterate() will first determine whether the solver is inintialized; if
   * not, it will call initialize() using default arguments. After
   * initialization, the solver performs block FGmres iterations until the
   * status test evaluates as ::Passed, at which point the method returns to
   * the caller. 
   *
   * The block FGmres iteration proceeds as follows:
   * -# The operator problem->applyOp() is applied to the newest \c blockSize vectors in the Krylov basis.
   * -# The resulting vectors are orthogonalized against the previous basis vectors, and made orthonormal.
   * -# The Hessenberg matrix is updated.
   * -# The least squares system is updated.
   *
   * The status test is queried at the beginning of the iteration.
   *
   * Possible exceptions thrown include the GmresIterationOrthoFailure.
   *
   */
  void iterate();

  /*! \brief Initialize the solver to an iterate, providing a complete state.
   *
   * The %BlockFGmresIter contains a certain amount of state, consisting of the current 
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
  void initializeGmres(GmresIterationState<ScalarType,MV> newstate);

  /*! \brief Initialize the solver with the initial vectors from the linear problem
   *  or random data.
   */
  void initialize()
  {
    GmresIterationState<ScalarType,MV> empty;
    initializeGmres(empty);
  }
  
  /*! \brief Get the current state of the linear solver.
   *
   * The data is only valid if isInitialized() == \c true.
   *
   * \returns A GmresIterationState object containing const pointers to the current
   * solver state.
   */
  GmresIterationState<ScalarType,MV> getState() const {
    GmresIterationState<ScalarType,MV> state;
    state.curDim = curDim_;
    state.V = V_;
    state.Z = Z_;
    state.H = H_;
    state.R = R_;
    state.z = z_;
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
  /*! \note Some solvers, like flexible GMRES, do not compute updates to the solution every iteration.
            This method forces its computation.  Other solvers, like CG, update the solution
            each iteration, so this method will return a zero std::vector indicating that the linear
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
  int getMaxSubspaceDim() const { return blockSize_*numBlocks_; }

  //@}

  
  //! @name Accessor methods
  //@{ 

  //! Get a constant reference to the linear problem.
  const LinearProblem<ScalarType,MV,OP>& getProblem() const { return *lp_; }

  //! Get the blocksize to be used by the iterative solver in solving this linear problem.
  int getBlockSize() const { return blockSize_; }
  
  //! \brief Set the blocksize.
  void setBlockSize(int blockSize) { setSize( blockSize, numBlocks_ ); }
  
  //! Get the maximum number of blocks used by the iterative solver in solving this linear problem.
  int getNumBlocks() const { return numBlocks_; }
  
  //! \brief Set the maximum number of blocks used by the iterative solver.
  void setNumBlocks(int numBlocks) { setSize( blockSize_, numBlocks ); }
  
  /*! \brief Set the blocksize and number of blocks to be used by the
   * iterative solver in solving this linear problem.
   *
   *  Changing either the block size or the number of blocks will reset the
   *  solver to an uninitialized state.
   */
  void setSize(int blockSize, int numBlocks);

  //! States whether the solver has been initialized or not.
  bool isInitialized() { return initialized_; }

  //@}

  private:

  //
  // Internal methods
  //
  //! Method for initalizing the state storage needed by block flexible GMRES.
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
  // It controls the number of vectors added to the basis on each iteration.
  int blockSize_;
  // numBlocks_ is the size of the allocated space for the Krylov basis, in blocks.
  int numBlocks_; 
  
  // Storage for QR factorization of the least squares system.
  Teuchos::SerialDenseVector<int,ScalarType> beta, sn;
  Teuchos::SerialDenseVector<int,MagnitudeType> cs;
  
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
  Teuchos::RCP<MV> V_;
  Teuchos::RCP<MV> Z_;
  //
  // Projected matrices
  // H_ : Projected matrix from the Krylov factorization AV = VH + FE^T
  //
  Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > H_;
  // 
  // QR decomposition of Projected matrices for solving the least squares system HY = B.
  // R_: Upper triangular reduction of H
  // z_: Q applied to right-hand side of the least squares system
  Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > R_;
  Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > z_;  
};

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor.
  template<class ScalarType, class MV, class OP>
  BlockFGmresIter<ScalarType,MV,OP>::BlockFGmresIter(const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem, 
						   const Teuchos::RCP<OutputManager<ScalarType> > &printer,
						   const Teuchos::RCP<StatusTest<ScalarType,MV,OP> > &tester,
						   const Teuchos::RCP<MatOrthoManager<ScalarType,MV,OP> > &ortho,
						   Teuchos::ParameterList &params ):
    lp_(problem),
    om_(printer),
    stest_(tester),
    ortho_(ortho),
    blockSize_(0),
    numBlocks_(0),
    initialized_(false),
    stateStorageInitialized_(false),
    curDim_(0),
    iter_(0)
  {
    // Get the maximum number of blocks allowed for this Krylov subspace
    TEUCHOS_TEST_FOR_EXCEPTION(!params.isParameter("Num Blocks"), std::invalid_argument,
                       "Belos::BlockFGmresIter::constructor: mandatory parameter 'Num Blocks' is not specified.");
    int nb = Teuchos::getParameter<int>(params, "Num Blocks");

    // Set the block size and allocate data
    int bs = params.get("Block Size", 1);
    setSize( bs, nb );
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Set the block size and make necessary adjustments.
  template <class ScalarType, class MV, class OP>
  void BlockFGmresIter<ScalarType,MV,OP>::setSize (int blockSize, int numBlocks)
  {
    // This routine only allocates space; it doesn't not perform any computation
    // any change in size will invalidate the state of the solver.

    TEUCHOS_TEST_FOR_EXCEPTION(numBlocks <= 0 || blockSize <= 0, std::invalid_argument, "Belos::BlockFGmresIter::setSize was passed a non-positive argument.");
    if (blockSize == blockSize_ && numBlocks == numBlocks_) {
      // do nothing
      return;
    }

    if (blockSize!=blockSize_ || numBlocks!=numBlocks_)
      stateStorageInitialized_ = false;

    blockSize_ = blockSize;
    numBlocks_ = numBlocks;

    initialized_ = false;
    curDim_ = 0;

    // Use the current blockSize_ and numBlocks_ to initialize the state storage.    
    setStateSize();
    
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Setup the state storage.
  template <class ScalarType, class MV, class OP>
  void BlockFGmresIter<ScalarType,MV,OP>::setStateSize ()
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
	
	//////////////////////////////////
	// blockSize*numBlocks dependent
	//
	int newsd = blockSize_*(numBlocks_+1);
	
	if (blockSize_==1) {
	  cs.resize( newsd );
	  sn.resize( newsd );
	}
	else {
	  beta.resize( newsd );
	}
	
	// Initialize the state storage
        TEUCHOS_TEST_FOR_EXCEPTION( blockSize_*static_cast<ptrdiff_t>(numBlocks_) > MVText::GetGlobalLength(*rhsMV),std::invalid_argument,
                           "Belos::BlockFGmresIter::setStateSize(): Cannot generate a Krylov basis with dimension larger the operator!");

	// If the subspace has not be initialized before, generate it using the LHS or RHS from lp_.
	if (V_ == Teuchos::null) {
	  // Get the multivector that is not null.
	  Teuchos::RCP<const MV> tmp = ( (rhsMV!=Teuchos::null)? rhsMV: lhsMV );
	  TEUCHOS_TEST_FOR_EXCEPTION(tmp == Teuchos::null,std::invalid_argument,
			     "Belos::BlockFGmresIter::setStateSize(): linear problem does not specify multivectors to clone from.");
	  V_ = MVT::Clone( *tmp, newsd );
	}
	else {
	  // Generate V_ by cloning itself ONLY if more space is needed.
	  if (MVT::GetNumberVecs(*V_) < newsd) {
	    Teuchos::RCP<const MV> tmp = V_;
	    V_ = MVT::Clone( *tmp, newsd );
	  }
	}

	if (Z_ == Teuchos::null) {
	  // Get the multivector that is not null.
	  Teuchos::RCP<const MV> tmp = ( (rhsMV!=Teuchos::null)? rhsMV: lhsMV );
	  TEUCHOS_TEST_FOR_EXCEPTION(tmp == Teuchos::null,std::invalid_argument,
			     "Belos::BlockFGmresIter::setStateSize(): linear problem does not specify multivectors to clone from.");
	  Z_ = MVT::Clone( *tmp, newsd );
	}
	else {
	  // Generate Z_ by cloning itself ONLY if more space is needed.
	  if (MVT::GetNumberVecs(*Z_) < newsd) {
	    Teuchos::RCP<const MV> tmp = Z_;
	    Z_ = MVT::Clone( *tmp, newsd );
	  }
	} 

	// Generate H_ only if it doesn't exist, otherwise resize it.
	if (H_ == Teuchos::null)
	  H_ = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>( newsd, newsd-blockSize_ ) );	
	else
	  H_->shapeUninitialized( newsd, newsd-blockSize_ );
	
	// TODO:  Insert logic so that Hessenberg matrix can be saved and reduced matrix is stored in R_
	//R_ = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>( newsd, newsd-blockSize_ ) );
	// Generate z_ only if it doesn't exist, otherwise resize it.
	if (z_ == Teuchos::null)
	  z_ = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(newsd,blockSize_) );
	else
	  z_->shapeUninitialized( newsd, blockSize_ );
	
	// State storage has now been initialized.
	stateStorageInitialized_ = true;
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Get the current update from this subspace.
  template <class ScalarType, class MV, class OP>
  Teuchos::RCP<MV> BlockFGmresIter<ScalarType,MV,OP>::getCurrentUpdate() const
  {
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
      currentUpdate = MVT::Clone( *Z_, blockSize_ );
      //
      //  Make a view and then copy the RHS of the least squares problem.  DON'T OVERWRITE IT!
      //
      Teuchos::SerialDenseMatrix<int,ScalarType> y( Teuchos::Copy, *z_, curDim_, blockSize_ );
      //
      //  Solve the least squares problem.
      //
      blas.TRSM( Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
		 Teuchos::NON_UNIT_DIAG, curDim_, blockSize_, one,  
		 H_->values(), H_->stride(), y.values(), y.stride() );
      //
      //  Compute the current update.
      //
      std::vector<int> index(curDim_);
      for ( int i=0; i<curDim_; i++ ) {   
        index[i] = i;
      }
      Teuchos::RCP<const MV> Zjp1 = MVT::CloneView( *Z_, index );
      MVT::MvTimesMatAddMv( one, *Zjp1, y, zero, *currentUpdate );
    }
    return currentUpdate;
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Get the native residuals stored in this iteration.  
  // Note:  No residual std::vector will be returned by FGmres.
  template <class ScalarType, class MV, class OP>
  Teuchos::RCP<const MV> BlockFGmresIter<ScalarType,MV,OP>::getNativeResiduals( std::vector<MagnitudeType> *norms ) const 
  {
    //
    // NOTE: Make sure the incoming std::vector is the correct size!
    //
    if ( norms && (int)norms->size() < blockSize_ )                         
      norms->resize( blockSize_ );                                          
    
    if (norms) {
      Teuchos::BLAS<int,ScalarType> blas;
      for (int j=0; j<blockSize_; j++) {
        (*norms)[j] = blas.NRM2( blockSize_, &(*z_)(curDim_, j), 1);
      }
    }
    return Teuchos::null;
  }

  

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Initialize this iteration object
  template <class ScalarType, class MV, class OP>
  void BlockFGmresIter<ScalarType,MV,OP>::initializeGmres(GmresIterationState<ScalarType,MV> newstate)
  {
    // Initialize the state storage if it isn't already.
    if (!stateStorageInitialized_) 
      setStateSize();

    TEUCHOS_TEST_FOR_EXCEPTION(!stateStorageInitialized_,std::invalid_argument,
		       "Belos::BlockFGmresIter::initialize(): Cannot initialize state storage!");
    
    // NOTE:  In BlockFGmresIter, V and Z are required!!!  
    // inconsitent multivectors widths and lengths will not be tolerated, and
    // will be treated with exceptions.
    //
    std::string errstr("Belos::BlockFGmresIter::initialize(): Specified multivectors must have a consistent length and width.");

    if (newstate.V != Teuchos::null && newstate.z != Teuchos::null) {

      // initialize V_,z_, and curDim_

      TEUCHOS_TEST_FOR_EXCEPTION( MVText::GetGlobalLength(*newstate.V) != MVText::GetGlobalLength(*V_),
                          std::invalid_argument, errstr );
      TEUCHOS_TEST_FOR_EXCEPTION( MVT::GetNumberVecs(*newstate.V) < blockSize_,
                          std::invalid_argument, errstr );
      TEUCHOS_TEST_FOR_EXCEPTION( newstate.curDim > blockSize_*(numBlocks_+1),
                          std::invalid_argument, errstr );

      curDim_ = newstate.curDim;
      int lclDim = MVT::GetNumberVecs(*newstate.V);

      // check size of Z
      TEUCHOS_TEST_FOR_EXCEPTION(newstate.z->numRows() < curDim_ || newstate.z->numCols() < blockSize_, std::invalid_argument, errstr);
      

      // copy basis vectors from newstate into V
      if (newstate.V != V_) {
        // only copy over the first block and print a warning.
        if (curDim_ == 0 && lclDim > blockSize_) {
	  om_->stream(Warnings) << "Belos::BlockFGmresIter::initialize(): the solver was initialized with a kernel of " << lclDim << std::endl
				<< "The block size however is only " << blockSize_ << std::endl
				<< "The last " << lclDim - blockSize_ << " vectors will be discarded." << std::endl;
	}
        std::vector<int> nevind(curDim_+blockSize_);
        for (int i=0; i<curDim_+blockSize_; i++) nevind[i] = i;
	Teuchos::RCP<const MV> newV = MVT::CloneView( *newstate.V, nevind );
	Teuchos::RCP<MV> lclV = MVT::CloneViewNonConst( *V_, nevind );
        MVT::MvAddMv( 1.0, *newV, 0.0, *newV, *lclV );

        // done with local pointers
        lclV = Teuchos::null;
      }

      // put data into z_, make sure old information is not still hanging around.
      if (newstate.z != z_) {
        z_->putScalar();
        Teuchos::SerialDenseMatrix<int,ScalarType> newZ(Teuchos::View,*newstate.z,curDim_+blockSize_,blockSize_);
        Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > lclz;
        lclz = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(Teuchos::View,*z_,curDim_+blockSize_,blockSize_) );
        lclz->assign(newZ);

        // done with local pointers
        lclz = Teuchos::null;
      }

    }
    else {

      TEUCHOS_TEST_FOR_EXCEPTION(newstate.V == Teuchos::null,std::invalid_argument,
                         "Belos::BlockFGmresIter::initialize(): BlockFGmresStateIterState does not have initial kernel V_0.");

      TEUCHOS_TEST_FOR_EXCEPTION(newstate.z == Teuchos::null,std::invalid_argument,
                         "Belos::BlockFGmresIter::initialize(): BlockFGmresStateIterState does not have initial norms z_0.");
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
  void BlockFGmresIter<ScalarType,MV,OP>::iterate()
  {
    //
    // Allocate/initialize data structures
    //
    if (initialized_ == false) {
      initialize();
    }
    
    // Compute the current search dimension. 
    int searchDim = blockSize_*numBlocks_;

    ////////////////////////////////////////////////////////////////
    // iterate until the status test tells us to stop.
    //
    // also break if our basis is full
    //
    while (stest_->checkStatus(this) != Passed && curDim_+blockSize_ <= searchDim) {

      iter_++;

      // F can be found at the curDim_ block, but the next block is at curDim_ + blockSize_.
      int lclDim = curDim_ + blockSize_; 

      // Get the current part of the basis.
      std::vector<int> curind(blockSize_);
      for (int i=0; i<blockSize_; i++) { curind[i] = lclDim + i; }
      Teuchos::RCP<MV> Vnext = MVT::CloneViewNonConst(*V_,curind);

      // Get a view of the previous vectors.
      // This is used for orthogonalization and for computing V^H K H
      for (int i=0; i<blockSize_; i++) { curind[i] = curDim_ + i; }
      Teuchos::RCP<const MV> Vprev = MVT::CloneView(*V_,curind);
      Teuchos::RCP<MV> Znext = MVT::CloneViewNonConst(*Z_,curind);

      // Compute the next std::vector in the Krylov basis:  Znext = M*Vprev
      lp_->applyRightPrec(*Vprev,*Znext);
      Vprev = Teuchos::null;

      // Compute the next std::vector in the Krylov basis:  Vnext = A*Znext
      lp_->applyOp(*Znext,*Vnext);
      Znext = Teuchos::null;
      
      // Remove all previous Krylov basis vectors from Vnext      
      // Get a view of all the previous vectors
      std::vector<int> prevind(lclDim);
      for (int i=0; i<lclDim; i++) { prevind[i] = i; }
      Vprev = MVT::CloneView(*V_,prevind);
      Teuchos::Array<Teuchos::RCP<const MV> > AVprev(1, Vprev);
      
      // Get a view of the part of the Hessenberg matrix needed to hold the ortho coeffs.
      Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> >
	subH = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>
			     ( Teuchos::View,*H_,lclDim,blockSize_,0,curDim_ ) );
      Teuchos::Array<Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > > AsubH;
      AsubH.append( subH );
      
      // Get a view of the part of the Hessenberg matrix needed to hold the norm coeffs.
      Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> >
	subR = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>
			     ( Teuchos::View,*H_,blockSize_,blockSize_,lclDim,curDim_ ) );
      int rank = ortho_->projectAndNormalize(*Vnext,AsubH,subR,AVprev);
      TEUCHOS_TEST_FOR_EXCEPTION(rank != blockSize_,GmresIterationOrthoFailure,
			 "Belos::BlockFGmresIter::iterate(): couldn't generate basis of full rank.");
      //
      // V has been extended, and H has been extended. 
      //
      // Update the QR factorization of the upper Hessenberg matrix
      //
      updateLSQR();
      //
      // Update basis dim and release all pointers.
      //
      Vnext = Teuchos::null;
      curDim_ += blockSize_;
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

  
  template<class ScalarType, class MV, class OP>
  void BlockFGmresIter<ScalarType,MV,OP>::updateLSQR( int dim )
  {
    int i, j, maxidx;
    ScalarType sigma, mu, vscale, maxelem;
    const ScalarType zero = Teuchos::ScalarTraits<ScalarType>::zero();
    
    // Get correct dimension based on input "dim"
    // Remember that ortho failures result in an exit before updateLSQR() is called.
    // Therefore, it is possible that dim == curDim_.
    int curDim = curDim_;
    if (dim >= curDim_ && dim < getMaxSubspaceDim()) {
      curDim = dim;
    }
    
    Teuchos::BLAS<int, ScalarType> blas;
    //
    // Apply previous transformations and compute new transformation to reduce upper-Hessenberg
    // system to upper-triangular form.
    //
    if (blockSize_ == 1) {
      //
      // QR factorization of Least-Squares system with Givens rotations
      //
      for (i=0; i<curDim; i++) {
	//
	// Apply previous Givens rotations to new column of Hessenberg matrix
	//
	blas.ROT( 1, &(*H_)(i,curDim), 1, &(*H_)(i+1, curDim), 1, &cs[i], &sn[i] );
      }
      //
      // Calculate new Givens rotation
      //
      blas.ROTG( &(*H_)(curDim,curDim), &(*H_)(curDim+1,curDim), &cs[curDim], &sn[curDim] );
      (*H_)(curDim+1,curDim) = zero;
      //
      // Update RHS w/ new transformation
      //
      blas.ROT( 1, &(*z_)(curDim,0), 1, &(*z_)(curDim+1,0), 1, &cs[curDim], &sn[curDim] );
    }
    else {
      //
      // QR factorization of Least-Squares system with Householder reflectors
      //
      for (j=0; j<blockSize_; j++) {
	//
	// Apply previous Householder reflectors to new block of Hessenberg matrix
	//
	for (i=0; i<curDim+j; i++) {
	  sigma = blas.DOT( blockSize_, &(*H_)(i+1,i), 1, &(*H_)(i+1,curDim+j), 1);
	  sigma += (*H_)(i,curDim+j);
	  sigma *= beta[i];
	  blas.AXPY(blockSize_, -sigma, &(*H_)(i+1,i), 1, &(*H_)(i+1,curDim+j), 1);
	  (*H_)(i,curDim+j) -= sigma;
	}
	//
	// Compute new Householder reflector
	//
	maxidx = blas.IAMAX( blockSize_+1, &(*H_)(curDim+j,curDim+j), 1 );
	maxelem = (*H_)(curDim+j+maxidx-1,curDim+j);
	for (i=0; i<blockSize_+1; i++)
	  (*H_)(curDim+j+i,curDim+j) /= maxelem;
	sigma = blas.DOT( blockSize_, &(*H_)(curDim+j+1,curDim+j), 1,
			  &(*H_)(curDim+j+1,curDim+j), 1 );
	if (sigma == zero) {
	  beta[curDim + j] = zero;
	} else {
	  mu = Teuchos::ScalarTraits<ScalarType>::squareroot((*H_)(curDim+j,curDim+j)*(*H_)(curDim+j,curDim+j)+sigma);
	  if ( Teuchos::ScalarTraits<ScalarType>::real((*H_)(curDim+j,curDim+j)) 
	       < Teuchos::ScalarTraits<MagnitudeType>::zero() ) {
	    vscale = (*H_)(curDim+j,curDim+j) - mu;
	  } else {
	    vscale = -sigma / ((*H_)(curDim+j,curDim+j) + mu);
	  }
	  beta[curDim+j] = 2.0*vscale*vscale/(sigma + vscale*vscale);
	  (*H_)(curDim+j,curDim+j) = maxelem*mu;
	  for (i=0; i<blockSize_; i++)
	    (*H_)(curDim+j+1+i,curDim+j) /= vscale;
	}
	//
	// Apply new Householder reflector to rhs
	//
	for (i=0; i<blockSize_; i++) {
	  sigma = blas.DOT( blockSize_, &(*H_)(curDim+j+1,curDim+j),
			    1, &(*z_)(curDim+j+1,i), 1);
	  sigma += (*z_)(curDim+j,i);
	  sigma *= beta[curDim+j];
	  blas.AXPY(blockSize_, -sigma, &(*H_)(curDim+j+1,curDim+j),
		    1, &(*z_)(curDim+j+1,i), 1);
	  (*z_)(curDim+j,i) -= sigma;
	}
      }
    } // end if (blockSize_ == 1)

    // If the least-squares problem is updated wrt "dim" then update the curDim_.
    if (dim >= curDim_ && dim < getMaxSubspaceDim()) {
      curDim_ = dim + blockSize_;
    }

  } // end updateLSQR()

} // end Belos namespace

#endif /* BELOS_BLOCK_FGMRES_ITER_HPP */
