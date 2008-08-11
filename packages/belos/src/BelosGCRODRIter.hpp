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

#include "Teuchos_BLAS.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TimeMonitor.hpp"

/*!	
  \class Belos::GCRODRIter
  
  \brief This class implements the GCRODR iteration, where a
  single-std::vector Krylov subspace is constructed.  The QR decomposition of 
  block, upper Hessenberg matrix is performed each iteration to update
  the least squares system and give the current linear system residuals.
 
  \author Michael Parks and Heidi Thornquist
*/

namespace Belos {
  
  //! @name GCRODRIter Structures
  //@{
  
  /** \brief Structure to contain pointers to GCRODRIter state variables.
   *
   * This struct is utilized by GCRODRIter::initialize() and GCRODRIter::getState().
   */
  template <class ScalarType, class MV>
  struct GCRODRIterState {
    /*! \brief The current dimension of the reduction.
     *
     * This should always be equal to BlockGmresIter::getCurSubspaceDim()
     */
    int curDim;
    
    /*! \brief The current Krylov basis. */
    Teuchos::RCP<const MV> V;
   
    /*! \brief The recycled subspace and its projection. */
    Teuchos::RCP<const MV> U, C;

    /*! \brief The current Hessenberg matrix.
     *
     * The \c curDim by \c curDim leading submatrix of H is the
     * projection of problem->getOperator() by the first \c curDim vectors in V.
     */
    Teuchos::RCP<const Teuchos::SerialDenseMatrix<int,ScalarType> > H;

    /*! \brief The projection of the Krylov subspace against the recycled subspace
     */
    Teuchos::RCP<const Teuchos::SerialDenseMatrix<int,ScalarType> > B;

    /*! \brief The current upper-triangular matrix from the QR reduction of H. */
    Teuchos::RCP<const Teuchos::SerialDenseMatrix<int,ScalarType> > R;

    /*! \brief The current right-hand side of the least squares system RY = Z. */
    Teuchos::RCP<const Teuchos::SerialDenseVector<int,ScalarType> > z;

    GCRODRIterState() : curDim(0), V(Teuchos::null), 
			U(Teuchos::null), C(Teuchos::null),
			H(Teuchos::null), R(Teuchos::null),
			z(Teuchos::null)
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
  class GCRODRIterInitFailure : public BelosError {public:
    GCRODRIterInitFailure(const std::string& what_arg) : BelosError(what_arg)
    {}};
  
  /** \brief GCRODRIterOrthoFailure is thrown when the GCRODRIter object is unable to
   * compute independent direction vectors in the GCRODRIter::iterate() routine.
   *
   * This std::exception is thrown from the GCRODRIter::iterate() method.
   *
   */
  class GCRODRIterOrthoFailure : public BelosError {public:
    GCRODRIterOrthoFailure(const std::string& what_arg) : BelosError(what_arg)
    {}};
  
  /** \brief GCRODRIterLAPACKFailure is thrown when a nonzero return value is passed back
   * from an LAPACK routine.
   *
   * This std::exception is thrown from the GCRODRIter::iterate() method.
   *
   */
  class GCRODRIterLAPACKFailure : public BelosError {public:
    GCRODRIterLAPACKFailure(const std::string& what_arg) : BelosError(what_arg)
    {}};
  
  //@}
  
  
  template<class ScalarType, class MV, class OP>
  class GCRODRIter : virtual public Iteration<ScalarType,MV,OP> {
    
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
		const Teuchos::RCP<StatusTest<ScalarType,MV,OP> > &tester,
		const Teuchos::RCP<MatOrthoManager<ScalarType,MV,OP> > &ortho,
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
    void initialize(GCRODRIterState<ScalarType,MV> newstate);
    
    /*! \brief Initialize the solver with the initial vectors from the linear problem
     *  or random data.
     */
    void initialize()
    {
      GCRODRIterState<ScalarType,MV> empty;
      initialize(empty);
    }
    
    /*! \brief Get the current state of the linear solver.
     *
     * The data is only valid if isInitialized() == \c true.
     *
     * \returns A GCRODRIterState object containing const pointers to the current
     * solver state.
     */
    GCRODRIterState<ScalarType,MV> getState() const {
      GCRODRIterState<ScalarType,MV> state;
      state.curDim = curDim_;
      state.V = V_;
      state.U = U_;
      state.C = C_;
      state.H = H_;
      state.B = B_;
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
    /*! \note Some solvers, like GMRES, do not compute updates to the solution every iteration.
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
    int getMaxSubspaceDim() const { return numBlocks_-recycledBlocks_; }
    
    //@}
    
    
    //! @name Accessor methods
    //@{ 
    
    //! Get a constant reference to the linear problem.
    const LinearProblem<ScalarType,MV,OP>& getProblem() const { return *lp_; }
    
    //! Get the maximum number of blocks used by the iterative solver in solving this linear problem.
    int getNumBlocks() const { return numBlocks_; }
    
    //! \brief Set the maximum number of blocks used by the iterative solver.
    void setNumBlocks(int numBlocks) { setSize( recycledBlocks_, numBlocks ); };
    
    //! Get the maximum number of recycled blocks used by the iterative solver in solving this linear problem.
    int getRecycledBlocks() const { return recycledBlocks_; }

    //! \brief Set the maximum number of recycled blocks used by the iterative solver.
    void setRecycledBlocks(int recycledBlocks) { setSize( recycledBlocks, numBlocks_ ); };
    
    //! Get the blocksize to be used by the iterative solver in solving this linear problem.
    int getBlockSize() const { return 1; }
    
    //! \brief Set the blocksize.
    void setBlockSize(int blockSize) {
      TEST_FOR_EXCEPTION(blockSize!=1,std::invalid_argument,
			 "Belos::GCRODRIter::setBlockSize(): Cannot use a block size that is not one.");
    }

    //! \brief Set the maximum number of blocks used by the iterative solver and the number of recycled vectors.
    void setSize( int recycledBlocks, int numBlocks );

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
    const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> >    lp_;
    const Teuchos::RCP<OutputManager<ScalarType> >          om_;
    const Teuchos::RCP<StatusTest<ScalarType,MV,OP> >       stest_;
    const Teuchos::RCP<OrthoManager<ScalarType,MV> >        ortho_;
    
    //
    // Algorithmic parameters
    //  
    // numBlocks_ is the size of the allocated space for the Krylov basis, in blocks.
    int numBlocks_; 

    // recycledBlocks_ is the size of the allocated space for the recycled subspace, in blocks.
    int recycledBlocks_; 
    
    // Storage for QR factorization of the least squares system.
    Teuchos::SerialDenseVector<int,ScalarType> sn;
    Teuchos::SerialDenseVector<int,MagnitudeType> cs;
    
    // 
    // Current solver state
    //
    // initialized_ specifies that the basis vectors have been initialized and the iterate() routine
    // is capable of running; _initialize is controlled  by the initialize() member method
    // For the implications of the state of initialized_, please see documentation for initialize()
    bool initialized_;
    
    // stateStorageInitialized_ specified that the state storage has be initialized to the current
    // numBlocks_ and numRecycledBlocks_.  This initialization may be postponed if the linear problem was
    // generated without the right-hand side or solution vectors.
    bool stateStorageInitialized_;
    
    // Current subspace dimension, and number of iterations performed.
    int curDim_, iter_;
    
    // 
    // State Storage
    //
    // Krylov vectors.
    Teuchos::RCP<MV> V_;
    //
    // Recycled subspace vectors.
    Teuchos::RCP<const MV> U_, C_;
    //
    // Projected matrices
    // H_ : Projected matrix from the Krylov factorization AV = VH + FE^T
    Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > H_;
    //
    // B_ : Projected matrix from the recycled subspace B = C^H*A*V
    Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > B_;
    //
    // QR decomposition of Projected matrices for solving the least squares system HY = B.
    // R_: Upper triangular reduction of H
    // z_: Q applied to right-hand side of the least squares system
    Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > R_;
    Teuchos::RCP<Teuchos::SerialDenseVector<int,ScalarType> > z_;  
  };
  
  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor.
  template<class ScalarType, class MV, class OP>
  GCRODRIter<ScalarType,MV,OP>::GCRODRIter(const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem, 
					   const Teuchos::RCP<OutputManager<ScalarType> > &printer,
					   const Teuchos::RCP<StatusTest<ScalarType,MV,OP> > &tester,
					   const Teuchos::RCP<MatOrthoManager<ScalarType,MV,OP> > &ortho,
					   Teuchos::ParameterList &params ):
    lp_(problem),
    om_(printer),
    stest_(tester),
    ortho_(ortho),
    numBlocks_(0),
    recycledBlocks_(0),
    initialized_(false),
    stateStorageInitialized_(false),
    curDim_(0),
    iter_(0)
  {
    // Get the maximum number of blocks allowed for this Krylov subspace
    TEST_FOR_EXCEPTION(!params.isParameter("Num Blocks"), std::invalid_argument,
                       "Belos::GCRODRIter::constructor: mandatory parameter \"Num Blocks\" is not specified.");
    int nb = Teuchos::getParameter<int>(params, "Num Blocks");

    TEST_FOR_EXCEPTION(!params.isParameter("Recycled Blocks"), std::invalid_argument,
                       "Belos::GCRODRIter::constructor: mandatory parameter \"Recycled Blocks\" is not specified.");
    int rb = Teuchos::getParameter<int>(params, "Recycled Blocks");

    // Set the number of blocks and allocate data
    setSize( rb, nb );
  }
  
  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Set the block size and make necessary adjustments.
  template <class ScalarType, class MV, class OP>
  void GCRODRIter<ScalarType,MV,OP>::setSize( int recycledBlocks, int numBlocks )
  {
    // This routine only allocates space; it doesn't not perform any computation
    // any change in size will invalidate the state of the solver.

    TEST_FOR_EXCEPTION(numBlocks <= 0, std::invalid_argument, "Belos::GCRODRIter::setSize() was passed a non-positive argument for \"Num Blocks\".");
    TEST_FOR_EXCEPTION(recycledBlocks <= 0, std::invalid_argument, "Belos::GCRODRIter::setSize() was passed a non-positive argument for \"Recycled Blocks\".");
    TEST_FOR_EXCEPTION(recycledBlocks >= numBlocks, std::invalid_argument, "Belos::GCRODRIter::setSize() the number of recycled blocks is larger than the allowable subspace.");

    if (numBlocks == numBlocks_ && recycledBlocks == recycledBlocks_) {
      // do nothing
      return;
    }
    else {
      stateStorageInitialized_ = false;
    }

    numBlocks_ = numBlocks;
    recycledBlocks_ = recycledBlocks;

    initialized_ = false;
    curDim_ = 0;

    // Use the current numBlocks_ to initialize the state storage.    
    setStateSize();
    
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Setup the state storage.
  template <class ScalarType, class MV, class OP>
  void GCRODRIter<ScalarType,MV,OP>::setStateSize ()
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
	int newsd = numBlocks_ - recycledBlocks_ + 1;
	
        cs.resize( newsd );
	sn.resize( newsd );
	
	// Initialize the state storage
        TEST_FOR_EXCEPTION(numBlocks_ > MVT::GetVecLength(*rhsMV),std::invalid_argument,
                           "Belos::GCRODRIter::setStateSize(): Cannot generate a Krylov basis with dimension larger the operator!");

	// If the subspace has not be initialized before, generate it using the LHS or RHS from lp_.
	if (V_ == Teuchos::null) {
	  // Get the multivector that is not null.
	  Teuchos::RCP<const MV> tmp = ( (rhsMV!=Teuchos::null)? rhsMV: lhsMV );
	  TEST_FOR_EXCEPTION(tmp == Teuchos::null,std::invalid_argument,
			     "Belos::GCRODRIter::setStateSize(): linear problem does not specify multivectors to clone from.");
	  V_ = MVT::Clone( *tmp, newsd );
	}
	else {
	  // Generate V_ by cloning itself ONLY if more space is needed.
	  if (MVT::GetNumberVecs(*V_) < newsd) {
	    Teuchos::RCP<const MV> tmp = V_;
	    V_ = MVT::Clone( *tmp, newsd );
	  }
	}

	// Generate B_ only if it doesn't exist, otherwise resize it.
	if (B_ == Teuchos::null)
	  B_ = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>( recycledBlocks_, newsd-1 ) );
	else
	  B_->shapeUninitialized( recycledBlocks_, newsd-1 );
	
	// Generate H_ only if it doesn't exist, otherwise resize it.
	if (H_ == Teuchos::null)
	  H_ = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>( newsd, newsd-1 ) );		
        else
	  H_->shape( newsd, newsd-1 );
	
	// Generate R_ only if it doesn't exist, otherwise resize it.
	if (R_ == Teuchos::null)
	  R_ = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>( newsd, newsd-1 ) );
	else
	  R_->shape( newsd, newsd-1 );

	if (z_ == Teuchos::null)
	  z_ = Teuchos::rcp( new Teuchos::SerialDenseVector<int,ScalarType>(newsd) );
	else
	  z_->shapeUninitialized( newsd, 1 );
	
	// State storage has now been initialized.
	stateStorageInitialized_ = true;
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Get the current update from this subspace.
  template <class ScalarType, class MV, class OP>
  Teuchos::RCP<MV> GCRODRIter<ScalarType,MV,OP>::getCurrentUpdate() const
  {
    //
    // If this is the first iteration of the Arnoldi factorization, 
    // there is no update, so return Teuchos::null. 
    //
    RCP<MV> currentUpdate = Teuchos::null;
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
      Teuchos::SerialDenseMatrix<int,ScalarType> y( Teuchos::Copy, *z_, curDim_, 1 );
      //
      //  Solve the least squares problem.
      //
      blas.TRSM( Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
		 Teuchos::NON_UNIT_DIAG, curDim_, 1, one,  
		 R_->values(), R_->stride(), y.values(), y.stride() );
      //
      //  Compute the current update from the Krylov basis; V(:,1:curDim_)*y.
      //
      std::vector<int> index(curDim_);
      for ( int i=0; i<curDim_; i++ ) {   
        index[i] = i;
      }
      RCP<const MV> Vjp1 = MVT::CloneView( *V_, index );
      MVT::MvTimesMatAddMv( one, *Vjp1, y, zero, *currentUpdate );
      //
      //  Add in portion of update from recycled subspace U; U(:,1:recycledBlocks_)*B*y.
      //
      Teuchos::SerialDenseMatrix<int,ScalarType> z(recycledBlocks_,1);
      Teuchos::SerialDenseMatrix<int,ScalarType> subB( Teuchos::View, *B_, recycledBlocks_, curDim_ );
      z.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS, one, subB, y, zero );
      MVT::MvTimesMatAddMv( -one, *U_, z, one, *currentUpdate );
    }
    return currentUpdate;
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Get the native residuals stored in this iteration.  
  // Note:  No residual std::vector will be returned by Gmres.
  template <class ScalarType, class MV, class OP>
  Teuchos::RCP<const MV> GCRODRIter<ScalarType,MV,OP>::getNativeResiduals( std::vector<MagnitudeType> *norms ) const 
  {
    //
    // NOTE: Make sure the incoming std::vector is the correct size!
    //
    if ( norms && (int)norms->size()==0 )                         
      norms->resize( 1 );                                          
    
    if (norms) {
      Teuchos::BLAS<int,ScalarType> blas;
      (*norms)[0] = blas.NRM2( 1, &(*z_)(curDim_), 1);
    }
    return Teuchos::null;
  }
  
  

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Initialize this iteration object
  template <class ScalarType, class MV, class OP>
  void GCRODRIter<ScalarType,MV,OP>::initialize(GCRODRIterState<ScalarType,MV> newstate)
  {
    // Initialize the state storage if it isn't already.
    if (!stateStorageInitialized_) 
      setStateSize();

    TEST_FOR_EXCEPTION(!stateStorageInitialized_,std::invalid_argument,
		       "Belos::GCRODRIter::initialize(): Cannot initialize state storage!");
    
    // NOTE:  In GCRODRIter, V and z are required!!!  
    // inconsitent multivectors widths and lengths will not be tolerated, and
    // will be treated with exceptions.
    //
    std::string errstr("Belos::GCRODRIter::initialize(): Specified multivectors must have a consistent length and width.");

    if (newstate.V != Teuchos::null && newstate.U != Teuchos::null && newstate.C != Teuchos::null && newstate.z != Teuchos::null) {

      // initialize V_,z_, and curDim_

      TEST_FOR_EXCEPTION( MVT::GetVecLength(*newstate.V) != MVT::GetVecLength(*V_),
                          std::invalid_argument, errstr );
      TEST_FOR_EXCEPTION( MVT::GetNumberVecs(*newstate.V) < 1,
                          std::invalid_argument, errstr );
      TEST_FOR_EXCEPTION( newstate.curDim > 1*(numBlocks_+1),
                          std::invalid_argument, errstr );

      curDim_ = newstate.curDim;
      int lclDim = MVT::GetNumberVecs(*newstate.V);

      // check size of Z
      TEST_FOR_EXCEPTION(newstate.z->numRows() < curDim_ || newstate.z->numCols() < 1, std::invalid_argument, errstr);
      

      // copy basis vectors from newstate into V
      if (newstate.V != V_) {
        // only copy over the first block and print a warning.
        if (curDim_ == 0 && lclDim > 1) {
	  om_->stream(Warnings) << "Belos::GCRODRIter::initialize(): the solver was initialized with a kernel of " << lclDim << std::endl
				<< "The last " << lclDim - 1 << " vectors will be discarded." << std::endl;
	}
        std::vector<int> nevind(curDim_+1);
        for (int i=0; i<curDim_+1; i++) nevind[i] = i;
	Teuchos::RCP<const MV> newV = MVT::CloneView( *newstate.V, nevind );
	Teuchos::RCP<MV> lclV = MVT::CloneView( *V_, nevind );
        MVT::MvAddMv( 1.0, *newV, 0.0, *newV, *lclV );

        // done with local pointers
        lclV = Teuchos::null;
      }

      // put data into z_, make sure old information is not still hanging around.
      if (newstate.z != z_) {
        z_->putScalar();
        Teuchos::SerialDenseMatrix<int,ScalarType> newZ(Teuchos::View,*newstate.z,curDim_+1,1);
        Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > lclZ;
        lclZ = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(Teuchos::View,*z_,curDim_+1,1) );
        lclZ->assign(newZ);

        // done with local pointers
        lclZ = Teuchos::null;
      }

      TEST_FOR_EXCEPTION( MVT::GetNumberVecs(*newstate.U) != recycledBlocks_,
                          std::invalid_argument, errstr );
      if (newstate.U != U_) {
	U_ = newstate.U;
      }

      TEST_FOR_EXCEPTION( MVT::GetNumberVecs(*newstate.C) != recycledBlocks_,
                          std::invalid_argument, errstr );
      if (newstate.C != C_) {
	C_ = newstate.C;
      }
    }
    else {

      TEST_FOR_EXCEPTION(newstate.V == Teuchos::null,std::invalid_argument,
                         "Belos::GCRODRIter::initialize(): GCRODRStateIterState does not have initial kernel V_0.");

      TEST_FOR_EXCEPTION(newstate.U == Teuchos::null,std::invalid_argument,
                         "Belos::GCRODRIter::initialize(): GCRODRStateIterState does not have recycled basis U.");

      TEST_FOR_EXCEPTION(newstate.C == Teuchos::null,std::invalid_argument,
                         "Belos::GCRODRIter::initialize(): GCRODRStateIterState does not have recycled basis C = A*U.");

      TEST_FOR_EXCEPTION(newstate.z == Teuchos::null,std::invalid_argument,
                         "Belos::GCRODRIter::initialize(): GCRODRStateIterState does not have initial norms z_0.");
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
  void GCRODRIter<ScalarType,MV,OP>::iterate()
  {
    //
    // Allocate/initialize data structures
    //
    if (initialized_ == false) {
      initialize();
    }
    
    // Compute the current search dimension. 
    int searchDim = numBlocks_ - recycledBlocks_;

    ////////////////////////////////////////////////////////////////
    // iterate until the status test tells us to stop.
    //
    // also break if our basis is full
    //
    while (stest_->checkStatus(this) != Passed && curDim_+1 <= searchDim) {

      iter_++;

      // F can be found at the curDim_ block, but the next block is at curDim_ + 1.
      int lclDim = curDim_ + 1; 

      // Get the current part of the basis.
      std::vector<int> curind(1);
      curind[0] = lclDim;
      Teuchos::RCP<MV> Vnext = MVT::CloneView(*V_,curind);

      // Get a view of the previous vectors.
      // This is used for orthogonalization and for computing V^H K H
      curind[0] = curDim_;
      Teuchos::RCP<MV> Vprev = MVT::CloneView(*V_,curind);

      // Compute the next std::vector in the Krylov basis:  Vnext = Op*Vprev
      lp_->apply(*Vprev,*Vnext);
      Vprev = Teuchos::null;

      // First, remove the recycled subspace (C) from Vnext and put coefficients in B.
      Teuchos::Array<Teuchos::RCP<const MV> > C(1, C_);
      Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> >
	subB = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>
			     ( Teuchos::View,*B_,recycledBlocks_,1,0,curDim_ ) );
      Teuchos::Array<Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > > AsubB;
      AsubB.append( subB );

      // Project out the recycled subspace.
      ortho_->project( *Vnext, AsubB, C );

      // Now, remove all previous Krylov basis vectors from Vnext and put coefficients in H and R.          
      // Get a view of all the previous vectors
      std::vector<int> prevind(lclDim);
      for (int i=0; i<lclDim; i++) { prevind[i] = i; }
      Vprev = MVT::CloneView(*V_,prevind);
      Teuchos::Array<Teuchos::RCP<const MV> > AVprev(1, Vprev);
	
      // Get a view of the part of the Hessenberg matrix needed to hold the ortho coeffs.
      Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> >
	subH = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>
			     ( Teuchos::View,*H_,lclDim,1,0,curDim_ ) );
      Teuchos::Array<Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > > AsubH;
      AsubH.append( subH );
      
      // Get a view of the part of the Hessenberg matrix needed to hold the norm coeffs.
      Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> >
	subR = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>
			     ( Teuchos::View,*H_,1,1,lclDim,curDim_ ) );

      // Project out the previous Krylov vectors and normalize the next vector.
      int rank = ortho_->projectAndNormalize(*Vnext,AsubH,subR,AVprev);

      // Copy over the coefficients to R just in case we run into an error.
      Teuchos::SerialDenseMatrix<int,ScalarType> subR2( Teuchos::View,*R_,lclDim+1,1,0,curDim_ );
      Teuchos::SerialDenseMatrix<int,ScalarType> subH2( Teuchos::View,*H_,lclDim+1,1,0,curDim_ );
      subR2.assign(subH2);
      
      TEST_FOR_EXCEPTION(rank != 1,GCRODRIterOrthoFailure,
			 "Belos::GCRODRIter::iterate(): couldn't generate basis of full rank.");
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
      curDim_++;
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
  void GCRODRIter<ScalarType,MV,OP>::updateLSQR( int dim )
  {
    int i;
    const ScalarType zero = Teuchos::ScalarTraits<ScalarType>::zero();
    
    // Get correct dimension based on input "dim"
    // Remember that ortho failures result in an exit before updateLSQR() is called.
    // Therefore, it is possible that dim == curDim_.
    int curDim = curDim_;
    if (dim >= curDim_ && dim < getMaxSubspaceDim()) {
      curDim = dim;
    }
    
    Teuchos::LAPACK<int, ScalarType> lapack;
    Teuchos::BLAS<int, ScalarType> blas;
    //
    // Apply previous transformations and compute new transformation to reduce upper-Hessenberg
    // system to upper-triangular form.
    //
    // QR factorization of Least-Squares system with Givens rotations
    //
    for (i=0; i<curDim; i++) {
      //
      // Apply previous Givens rotations to new column of Hessenberg matrix
      //
      blas.ROT( 1, &(*R_)(i,curDim), 1, &(*R_)(i+1, curDim), 1, &cs[i], &sn[i] );
    }
    //
    // Calculate new Givens rotation
    //
    blas.ROTG( &(*R_)(curDim,curDim), &(*R_)(curDim+1,curDim), &cs[curDim], &sn[curDim] );
    (*R_)(curDim+1,curDim) = zero;
    //
    // Update RHS w/ new transformation
    //
    blas.ROT( 1, &(*z_)(curDim), 1, &(*z_)(curDim+1), 1, &cs[curDim], &sn[curDim] );

  } // end updateLSQR()

} // end Belos namespace

#endif /* BELOS_GCRODR_ITER_HPP */
