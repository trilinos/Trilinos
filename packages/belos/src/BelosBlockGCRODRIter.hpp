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

#ifndef BELOS_BLOCK_GCRODR_ITER_HPP
#define BELOS_BLOCK_GCRODR_ITER_HPP


/*! \file BelosBlockGCRODRIter.hpp
 *     \brief Belos concrete class for performing the block GCRO-DR (block GMRES with recycling) iteration.
 *     */

#include "BelosConfigDefs.hpp"
#include "BelosTypes.hpp"

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

// MLP
#include <unistd.h>

/*!
  \class Belos::BlockGCRODRIter
  \brief Implementation of the Block GCRO-DR (Block Recycling GMRES) iteration.

  This class implements the Block GCRODR (block GMRES with recycling)
  iteration, wherein a block-vector Krylov subspace is constructed.
  The QR decomposition of a block upper Hessenberg matrix is performed
  each iteration to update the least squares system and give the
  current linear system residuals.  

  \ingroup belos_solver_framework
  \author Kirk M. Soodhalter and Michael Parks
*/

namespace Belos{

  //! @name BlockGCRODRIter Structures
  //@{

  /** \brief Structure to contain pointers to BlockGCRODRIter state variables.
   *
   * This struct is utilized by BlockGCRODRIter::initialize() and BlockGCRODRIter::getState().
   */
 
	template <class ScalarType, class MV>
	struct BlockGCRODRIterState {
		/*! \brief The current dimension of the reduction.
		*
		* This should always be equal to BlockGCRODRIter::getCurSubspaceDim()
		*/

		int curDim;

		/*! \brief The current Krylov basis. */
		Teuchos::RCP<MV> V; 

		/*! \brief The recycled subspace and its projection. */
		Teuchos::RCP<MV> U, C;	

		/*! \brief The current Hessenberg matrix.
		 *
		 * The \c curDim by \c curDim leading submatrix of H is the
		 *  projection of problem->getOperator() by the first \c curDim vectors in V.
		 */
		 Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > H;
		    
		/*! \brief The projection of the Krylov subspace against the recycled subspace *      
		 */
		 Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > B;

		 BlockGCRODRIterState() : curDim(0), V(Teuchos::null),
			U(Teuchos::null), C(Teuchos::null),
			H(Teuchos::null), B(Teuchos::null)
		 {}

	};

	//@}



	//! @name BlockGCRODRIter Exceptions
	//@{

	/** \brief BlockGCRODRIterInitFailure is thrown when the BlockGCRODRIter object is unable to
	* generate an initial iterate in the BlockGCRODRIter::initialize() routine.
	*
	* This std::exception is thrown from the GCRODRIter::initialize() method, which is
	* called by the user or from the BlockGCRODRIter::iterate() method if isInitialized()
	* == \c false.
	*
	* In the case that this std::exception is thrown,
	* BlockGCRODRIter::isInitialized() will be \c false and the user will need to provide
	* a new initial iterate to the iteration.
	*/

	class BlockGCRODRIterInitFailure : public BelosError {
		public:
		BlockGCRODRIterInitFailure(const std::string& what_arg) : BelosError(what_arg) {}
	};

	/** \brief BlockGCRODRIterOrthoFailure is thrown when the BlockGCRODRIter object is unable to
	* compute independent direction vectors in the BlockGCRODRIter::iterate() routine.
	*
	* This std::exception is thrown from the BlockGCRODRIter::iterate() method.
	*
	*/

	class BlockGCRODRIterOrthoFailure : public BelosError {
		public:
		BlockGCRODRIterOrthoFailure(const std::string& what_arg) : BelosError(what_arg) {}
	};

  	//@}


    template<class ScalarType, class MV, class OP>
    class BlockGCRODRIter : virtual public Iteration<ScalarType,MV,OP> {
    	public:

  	//
  	//Convenience typedefs
  	//
  	typedef MultiVecTraits<ScalarType,MV> MVT;
  	typedef OperatorTraits<ScalarType,MV,OP> OPT;
  	typedef Teuchos::ScalarTraits<ScalarType> SCT;
  	typedef typename SCT::magnitudeType MagnitudeType;
 	typedef Teuchos::SerialDenseMatrix<int,ScalarType> SDM;
    	typedef Teuchos::SerialDenseVector<int,ScalarType> SDV;

	//! @name Constructors/Destructor
	//@{

	/*! \brief %BlockGCRODRIter constructor with linear problem, solver utilities, and parameter list of solver options.
	*
	* This constructor takes pointers required by the linear solver, in addition
	* to a parameter list of options for the linear solver. These options include the following:
	*   - "Num Blocks" - an \c int specifying the maximum number of blocks allocated for the solver basis. Default: 25
	*   - "Restart Timers" = a \c bool specifying whether the timers should be restarted each time iterate() is called. Default: false
	*/

       BlockGCRODRIter( const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem,
                        const Teuchos::RCP<OutputManager<ScalarType> > &printer,
                        const Teuchos::RCP<StatusTest<ScalarType,MV,OP> > &tester,
                        const Teuchos::RCP<MatOrthoManager<ScalarType,MV,OP> > &ortho,
                              Teuchos::ParameterList &params );

       //! Destructor.
       virtual ~BlockGCRODRIter() {};
       //@}
 
	//! @name Solver methods
	//@{

	/*! \brief This method performs block GCRODR iterations until the status
	* test indicates the need to stop or an error occurs (in which case, an
	* std::exception is thrown).
	*
	* iterate() will first determine whether the solver is inintialized; if
	* not, it will call initialize() using default arguments. After
	* initialization, the solver performs block GCRODR iterations until the
	* status test evaluates as ::Passed, at which point the method returns to
	* the caller.
	*
	* The block GCRODR iteration proceeds as follows:
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
	* The %BlockGCRODRIter contains a certain amount of state, consisting of the current
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

       void initialize() {
       	  BlockGCRODRIterState<ScalarType,MV> empty;
     	  initialize(empty);
       }  

	/*! \brief Initialize the solver with empty data. Calling this method will result in error,
	*  as GCRODRIter must be initialized with a valid state.
	*/
       void initialize(BlockGCRODRIterState<ScalarType,MV> newstate);
       
	/*! \brief Get the current state of the linear solver.
	*
	* The data is only valid if isInitialized() == \c true.
	*
	* \returns A BlockGCRODRIterState object containing const pointers to the current
	* solver state.
	*/
       BlockGCRODRIterState<ScalarType,MV> getState() const{
		BlockGCRODRIterState<ScalarType,MV> state;
		state.curDim = curDim_;
	        state.V = V_;
      		state.U = U_;
      		state.C = C_;
      		state.H = H_;
      		state.B = B_;
      		return state;
       }
       //@}

	//! @name Status methods
	//@{

	//! States whether the solver has been initialized or not.
       bool isInitialized(){ return initialized_;};

       //! \brief Get the current iteration count.
       int getNumIters() const { return iter_; };

       //! \brief Reset the iteration count.
       void resetNumIters( int iter = 0 ) { iter_ = iter; };

       //! Get the norms of the residuals native to the solver.
       //! \return A vector of length blockSize containing the native residuals.
       Teuchos::RCP<const MV> getNativeResiduals( std::vector<MagnitudeType> *norms ) const;

	//! Get the current update to the linear system.
	/*! \note Some solvers, like GMRES, do not compute updates to the solution every iteration.
		This method forces its computation.  Other solvers, like CG, update the solution
		each iteration, so this method will return a zero vector indicating that the linear
		problem contains the current solution.
	*/
       Teuchos::RCP<MV> getCurrentUpdate() const;


	//@}


	//! @name Accessor methods
	//@{


       //! Get a constant reference to the linear problem.
       const LinearProblem<ScalarType,MV,OP>& getProblem() const { return *lp_; };

       //! Get the maximum number of blocks used by the iterative solver in solving this linear problem.
       int getNumBlocks() const { return numBlocks_; }

       //! Get the blocksize to be used by the iterative solver in solving this linear problem.
       int getBlockSize() const { return blockSize_; };

	//! Get the dimension of the search subspace used to generate the current solution to the linear problem.
       int getCurSubspaceDim() const {
		if (!initialized_) return 0;
      		return curDim_;
       };

	//! Get the maximum dimension allocated for the search subspace.
       int getMaxSubspaceDim() const { return numBlocks_*blockSize_; };

	//! \brief Set the maximum number of recycled blocks used by the iterative solver.
       int getRecycledBlocks() const { return recycledBlocks_; };

	//@}


	//! @name Set methods
	//@{

       void updateLSQR( int dim = -1);

       //! \brief Set the blocksize.
       void setBlockSize(int blockSize){ blockSize_ = blockSize; }

       //! \brief Set the maximum number of recycled blocks used by the iterative solver.
       void setRecycledBlocks(int recycledBlocks) { setSize( recycledBlocks, numBlocks_ ); };

       //! \brief Set the maximum number of blocks used by the iterative solver.
       void setNumBlocks(int numBlocks) { setSize( recycledBlocks_, numBlocks ); };

       //! \brief Set the maximum number of blocks used by the iterative solver and the number of recycled vectors.
       void setSize( int recycledBlocks, int numBlocks ) {
       // only call resize if size changed
       if ( (recycledBlocks_ != recycledBlocks) || (numBlocks_ != numBlocks) ) {
      	recycledBlocks_ = recycledBlocks;
       	numBlocks_ = numBlocks;
        cs_.sizeUninitialized( numBlocks_ );
        sn_.sizeUninitialized( numBlocks_ );
        Z_.shapeUninitialized( numBlocks_*blockSize_, blockSize_ );
      }
    }

	//@}

      private:

	//
	// Internal methods
	//



      //Classes inputed through constructor that define the linear problem to be solved
      //
      const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> >    lp_;
      const Teuchos::RCP<OutputManager<ScalarType> >          om_;
      const Teuchos::RCP<StatusTest<ScalarType,MV,OP> >       stest_;
      const Teuchos::RCP<OrthoManager<ScalarType,MV> >        ortho_;

      //
      //Algorithmic Parameters
      //

      //numBlocks_ is the size of the allocated space for the Krylov basis, in blocks.
      //blockSize_ is the number of columns in each block Krylov vector.
      int numBlocks_, blockSize_;

      //boolean vector indicating which right-hand sides we care about
      //when we are testing residual norms.  THIS IS NOT IMPLEMENTED.  RIGHT NOW JUST
      //SELECTS ALL RIGHT HANDS SIDES FOR NORM COMPUTATION.
      std::vector<bool> trueRHSIndices_;

      // recycledBlocks_ is the size of the allocated space for the recycled subspace, in vectors.
      int recycledBlocks_;    

      //Storage for QR factorization of the least squares system if using plane rotations.
      SDV sn_;
      Teuchos::SerialDenseVector<int,MagnitudeType> cs_;

      //Storage for QR factorization of the least squares system if using Householder reflections
      //Per block Krylov vector, we actually construct a 2*blockSize_ by 2*blockSize_ matrix which
      //is the product of all Householder transformations for that block.  This has been shown to yield
      //speed ups without losing accuracy because we can apply previous Householder transformations
      //with BLAS3 operations.
      std::vector< SDM >House_;
      SDV beta_;

      //
      //Current Solver State
      //
      //initialized_ specifies that the basis vectors have been initialized and the iterate() routine
      //is capable of running; _initialize is controlled  by the initialize() member method
      //For the implications of the state of initialized_, please see documentation for initialize()
      bool initialized_;    
    
      // Current subspace dimension, number of iterations performed, and number of iterations performed in this cycle.
      int curDim_, iter_, lclIter_;
  
      //
      // Recycled Krylov Method Storage
      //


      //! The Krylov basis vectors.
      Teuchos::RCP<MV> V_;

      //! Recycled subspace vectors.
      Teuchos::RCP<MV> U_, C_;


      //! @name Projected operators on the augmented Krylov subspace
      //@{

      /// \brief Projected matrix from the Krylov factorization.
      ///
      /// The matrix H satisfies \f$AV = VH + C*B\f$, wherein \f$B = C^H*A*V\f$.
      Teuchos::RCP<SDM > H_;

      /// \brief Projected matrix from the recycled subspace.
      ///
      /// The matrix B satisfies \f$B = C^H*A*V\f$.
      Teuchos::RCP<SDM > B_;

      /// \brief Upper triangular reduction of H_ (see above).
      ///
      /// R_ and Z_ together form the QR decomposition of the
      /// projected matrices for solving the least-squares system
      /// HY = RHS.
      ///
      Teuchos::RCP<SDM> R_;

      //! Q applied to right-hand side of the least squares system.
      SDM Z_;

      //@}

      // File stream variables to use Mike Parks' Matlab output codes.
      std::ofstream ofs;
      char filename[30];

   };//End BlockGCRODRIter Class Definition

   //////////////////////////////////////////////////////////////////////////////////////////////////
   //Constructor.
   template<class ScalarType, class MV, class OP>
   BlockGCRODRIter<ScalarType,MV,OP>::BlockGCRODRIter(const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem,
                                            const Teuchos::RCP<OutputManager<ScalarType> > &printer,
                                            const Teuchos::RCP<StatusTest<ScalarType,MV,OP> > &tester,
                                            const Teuchos::RCP<MatOrthoManager<ScalarType,MV,OP> > &ortho,
                                                  Teuchos::ParameterList &params ):lp_(problem), 
                                                  om_(printer), stest_(tester), ortho_(ortho) {
    	numBlocks_      = 0;
	blockSize_ 	= 0;
     	recycledBlocks_ = 0;
     	initialized_    = false;
    	curDim_         = 0;
    	iter_           = 0;
	lclIter_        = 0;
    	V_              = Teuchos::null;
    	U_              = Teuchos::null;
   	C_              = Teuchos::null;
    	H_              = Teuchos::null;
    	B_              = Teuchos::null;
	R_              = Teuchos::null;
	// Get the maximum number of blocks allowed for this Krylov subspace
	TEST_FOR_EXCEPTION(!params.isParameter("Num Blocks"), std::invalid_argument, "Belos::BlockGCRODRIter::constructor: mandatory parameter \"Num Blocks\" is not specified.");
	int nb = Teuchos::getParameter<int>(params, "Num Blocks");

	TEST_FOR_EXCEPTION(!params.isParameter("Recycled Blocks"), std::invalid_argument,"Belos::BlockGCRODRIter::constructor: mandatory parameter \"Recycled Blocks\" is not specified.");
	int rb = Teuchos::getParameter<int>(params, "Recycled Blocks");

	TEST_FOR_EXCEPTION(nb <= 0, std::invalid_argument, "Belos::BlockGCRODRIter() was passed a non-positive argument for \"Num Blocks\".");
	TEST_FOR_EXCEPTION(rb >= nb, std::invalid_argument, "Belos::BlockGCRODRIter() the number of recycled blocks is larger than the allowable subspace.");


	int bs = Teuchos::getParameter<int>(params, "Block Size");

	TEST_FOR_EXCEPTION(bs <= 0, std::invalid_argument, "Belos::BlockGCRODRIter() the block size was passed a non-postitive argument.");


	numBlocks_ = nb;
	recycledBlocks_ = rb;
	blockSize_ = bs;

	//INITIALIZE ENTRIES OF trueRHSIndices_ TO CHECK EVERY NORM FOR NOW.  LATER, THE USER
	//SHOULD SPECIFY WHICH RIGHT HAND SIDES ARE IMPORTANT FOR CONVERGENCE TESTING
	trueRHSIndices_.resize(blockSize_);
	int i;
	for(i=0; i<blockSize_; i++){
		trueRHSIndices_[i] = true;
	}

        //THIS MAKES SPACE FOR GIVENS ROTATIONS BUT IN REALITY WE NEED TO DO TESTING ON BLOCK SIZE
        //AND CHOOSE BETWEEN GIVENS ROTATIONS AND HOUSEHOLDER TRANSFORMATIONS.        
        cs_.sizeUninitialized( numBlocks_+1 );
    	sn_.sizeUninitialized( numBlocks_+1 );
    	Z_.shapeUninitialized( (numBlocks_+1)*blockSize_,blockSize_ );

	House_.resize(numBlocks_);

	for(i=0; i<numBlocks_;i++){
		House_[i].shapeUninitialized(2*blockSize_, 2*blockSize_);
	}
   }//End Constructor Definition

   //////////////////////////////////////////////////////////////////////////////////////////////////
   // Iterate until the status test informs us we should stop.
   template <class ScalarType, class MV, class OP>
   void BlockGCRODRIter<ScalarType,MV,OP>::iterate() {
	TEST_FOR_EXCEPTION( initialized_ == false, BlockGCRODRIterInitFailure,"Belos::BlockGCRODRIter::iterate(): GCRODRIter class not initialized." );

// MLP
sleep(1);
std::cout << "Calling setSize" << std::endl;
	// Force call to setsize to ensure internal storage is correct dimension
	setSize( recycledBlocks_, numBlocks_ );

	Teuchos::RCP<MV> Vnext;
	Teuchos::RCP<const MV> Vprev;
	std::vector<int> curind(blockSize_);

	// z_ must be zeroed out in order to compute Givens rotations correctly
	Z_.putScalar(0.0);

	// Orthonormalize the new V_0
	for(int i = 0; i<blockSize_; i++){curind[i] = i;};
	
// MLP
sleep(1);
std::cout << "Calling normalize" << std::endl;
	Vnext = MVT::CloneViewNonConst(*V_,curind);
	//Orthonormalize Initial Columns
	//Store orthogonalization coefficients in Z0
	Teuchos::RCP<SDM > Z0 =
               Teuchos::rcp( new SDM(blockSize_,blockSize_) );
	int rank = ortho_->normalize(*Vnext,Z0);

// MLP
sleep(1);
std::cout << "Assigning Z" << std::endl;
	TEST_FOR_EXCEPTION(rank != blockSize_,BlockGCRODRIterOrthoFailure, "Belos::BlockGCRODRIter::iterate(): couldn't generate basis of full rank at the initial step.");
	// Copy Z0 into the leading blockSize_ by blockSize_ block of Z_
	Teuchos::RCP<SDM > Z_block = Teuchos::rcp( new SDM(Teuchos::View, Z_, blockSize_,blockSize_) );
	Z_block->assign(*Z0);

	std::vector<int> prevind(blockSize_*(numBlocks_ + 1));

	////////////////////////////////////////////////////////////////
	// iterate until the status test tells us to stop.
	//
	// also break if the basis is full
	//
	while( (stest_->checkStatus(this) != Passed) && (curDim_+blockSize_-1) < (numBlocks_*blockSize_)) {
		lclIter_++;
		iter_++;
//KMS
//std::cout << "Iter=" << iter_ << std::endl << "lclIter=" << lclIter_ <<  std::endl;

		int HFirstCol = curDim_-blockSize_;//First column of H we need view of
		int HLastCol = HFirstCol + blockSize_-1 ;//last column of H we need a view of
		int HLastOrthRow = HLastCol;//Last row of H we will put orthog coefficients in
		int HFirstNormRow = HLastOrthRow + 1;//First row of H where normalization matrix goes
//KMS
//std::cout << "curDim_ = " << curDim_ << ", HFirstCol = " << HFirstCol << ", HLastCol =  " << HLastCol <<", HLastOrthRow =  " << HLastOrthRow << ", HFirstNormRow =  " << HFirstNormRow << std::endl;		
		// Get next basis indices
		for(int i = 0; i< blockSize_; i++){
			curind[i] = curDim_ + i;
		}
		Vnext = MVT::CloneViewNonConst(*V_,curind);

		//Get a view of the previous block Krylov vector.
		//This is used for orthogonalization and for computing V^H K H
		// Get next basis indices
		for(int i = 0; i< blockSize_; i++){
                        curind[blockSize_ - 1 - i] = curDim_ -  i - 1;
                }
		Vprev = MVT::CloneView(*V_,curind);
		// Compute the next vector in the Krylov basis:  Vnext = Op*Vprev
		lp_->apply(*Vprev,*Vnext);
		Vprev = Teuchos::null;

		//First, remove the recycled subspace (C) from Vnext and put coefficients in B.

		//Get a view of the matrix B and put the pointer into an array
		//Put a pointer to C in another array
		Teuchos::Array<Teuchos::RCP<const MV> > C(1, C_);

		Teuchos::RCP<SDM >
          	subB = Teuchos::rcp( new SDM ( Teuchos::View,*B_,recycledBlocks_,blockSize_,0, HFirstCol ) );

		Teuchos::Array<Teuchos::RCP<SDM > > AsubB;
		AsubB.append( subB );
		// Project out the recycled subspace.
                ortho_->project( *Vnext, AsubB, C );
		//Now, remove block Krylov Subspace from Vnext and store coefficients in H_ and R_
		
		// Get a view of all the previous vectors
		prevind.resize(curDim_);
		for (int i=0; i<curDim_; i++) { prevind[i] = i; }
		Vprev = MVT::CloneView(*V_,prevind);
		Teuchos::Array<Teuchos::RCP<const MV> > AVprev(1, Vprev);

		// Get a view of the part of the Hessenberg matrix needed to hold the ortho coeffs.
		Teuchos::RCP<SDM> subH = Teuchos::rcp( new SDM  ( Teuchos::View,*H_,curDim_,blockSize_,0,HFirstCol ) );
		Teuchos::Array<Teuchos::RCP<SDM > > AsubH;
		AsubH.append( subH );
		// Get a view of the part of the Hessenberg matrix needed to hold the norm coeffs.
		Teuchos::RCP<SDM >  subR = Teuchos::rcp( new SDM ( Teuchos::View,*H_,blockSize_,blockSize_,HFirstNormRow,HFirstCol ) );
		// Project out the previous Krylov vectors and normalize the next vector.
		int rank = ortho_->projectAndNormalize(*Vnext,AsubH,subR,AVprev);

		// Copy over the coefficients to R just in case we run into an error.
		SDM subR2( Teuchos::View,*R_,(lclIter_+1)*blockSize_,blockSize_,0,HFirstCol);
		SDM subH2( Teuchos::View,*H_,(lclIter_+1)*blockSize_,blockSize_,0,HFirstCol);
		subR2.assign(subH2);

		TEST_FOR_EXCEPTION(rank != blockSize_,BlockGCRODRIterOrthoFailure, "Belos::BlockGCRODRIter::iterate(): couldn't generate basis of full rank.");

		// Update the QR factorization of the upper Hessenberg matrix
		updateLSQR();

        	// Update basis dim
                curDim_ = curDim_ + blockSize_;



	}//end while(stest_->checkStatus(this) ~= Passed && curDim_+1 <= numBlocks_*blockSize_)
	
   }//end iterate() defintition

   //////////////////////////////////////////////////////////////////////////////////////////////////
   //Initialize this iteration object.
   template <class ScalarType, class MV, class OP>
   void BlockGCRODRIter<ScalarType,MV,OP>::initialize(BlockGCRODRIterState<ScalarType,MV> newstate) {
	if (newstate.V != Teuchos::null &&  newstate.H != Teuchos::null) {
      		curDim_ = newstate.curDim;
      		V_      = newstate.V;
      		U_      = newstate.U;
      		C_      = newstate.C;
      		H_      = newstate.H;
      		B_      = newstate.B;
		lclIter_ = 0;//resets the count of local iterations for the new cycle
		R_      = Teuchos::rcp(new SDM(H_->numRows(), H_->numCols() )); //R_ should look like H but point to separate memory

		//All Householder product matrices start out as identity matrices.
		//We construct an identity from which to copy.
		SDM Identity(2*blockSize_, 2*blockSize_);
		for(int i=0;i<2*blockSize_; i++){
			Identity[i][i] = 1;
		}	
		for(int i=0; i<numBlocks_;i++){
			House_[i].assign(Identity);
		}
    	}
    	else {
      		TEST_FOR_EXCEPTION(newstate.V == Teuchos::null,std::invalid_argument,"Belos::GCRODRIter::initialize(): BlockGCRODRIterState does not have V initialized.");
      		TEST_FOR_EXCEPTION(newstate.H == Teuchos::null,std::invalid_argument,"Belos::GCRODRIter::initialize(): BlockGCRODRIterState does not have H initialized.");
    	}
    	// the solver is initialized
        initialized_ = true;
   }//end initialize() defintition

   //////////////////////////////////////////////////////////////////////////////////////////////////
   //Get the native residuals stored in this iteration.
   //This function will only compute the native residuals for 
   //right-hand sides we are interested in, as dictated by 
   //std::vector<int> trueRHSIndices_ (THIS IS NOT YET IMPLEMENTED.  JUST GETS ALL RESIDUALS)
   //A norm of -1 is entered for all residuals about which we do not care.
   template <class ScalarType, class MV, class OP>
   Teuchos::RCP<const MV> 
   BlockGCRODRIter<ScalarType,MV,OP>::getNativeResiduals( std::vector<MagnitudeType> *norms ) const
   {
	//
	// NOTE: Make sure the incoming std::vector is the correct size!
	//
	if (norms != NULL) {
		if (static_cast<int> (norms->size()) < blockSize_) {
     			norms->resize( blockSize_ );
		}
      		Teuchos::BLAS<int,ScalarType> blas;
      		for (int j=0; j<blockSize_; j++) {
			if(trueRHSIndices_[j]){
        			(*norms)[j] = blas.NRM2( blockSize_, &Z_(curDim_-blockSize_+j, j), 1);
			}
			else{
				(*norms)[j] = -1;
			}
      		}
   		return Teuchos::null;
    	} else { // norms is NULL
		// FIXME If norms is NULL, return residual vectors.
   		return Teuchos::null;
	}
   }//end getNativeResiduals() definition

   //////////////////////////////////////////////////////////////////////////////////////////////////
   //Get the current update from this subspace.
   template <class ScalarType, class MV, class OP>
   Teuchos::RCP<MV> BlockGCRODRIter<ScalarType,MV,OP>::getCurrentUpdate() const {
	//
	// If this is the first iteration of the Arnoldi factorization,
	// there is no update, so return Teuchos::null.
	//
	Teuchos::RCP<MV> currentUpdate = Teuchos::null;
//KMS	if(curDim_==0) {
	if(curDim_<=blockSize_) {
		return currentUpdate;
	}
	else{
		const ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();
		const ScalarType zero = Teuchos::ScalarTraits<ScalarType>::zero();
		Teuchos::BLAS<int,ScalarType> blas;
		currentUpdate = MVT::Clone( *V_, blockSize_ );
		//
		// Make a view and then copy the RHS of the least squares problem.  DON'T OVERWRITE IT!
		//
		SDM Y( Teuchos::Copy, Z_, curDim_-blockSize_, blockSize_ );
		Teuchos::RCP<SDM> Rtmp = Teuchos::rcp(new SDM(Teuchos::View, *R_, curDim_, curDim_-blockSize_));
//KMS
sleep(1);
std::cout << "Before TRSM" << std::endl;
sleep(1);
std::cout << "The size of Rtmp is " << Rtmp -> numRows() << " by " << Rtmp -> numCols() << std::endl;
std::cout << "The size of Y is " << Y.numRows() << " by " << Y.numCols() << std::endl;
std::cout << "blockSize_ = " << blockSize_ << std::endl;
std::cout << "curDim_ =  " << curDim_ << std::endl;
std::cout << "curDim_ - blockSize_ =  " << curDim_ - blockSize_ << std::endl;
		//
		// Solve the least squares problem.
		// Observe that in calling TRSM, we use the value
		// curDim_ -blockSize_. This is because curDim_ has
		// already been incremented but the proper size of R is still 
		// based on the previous value.
		//
		blas.TRSM( Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
                Teuchos::NON_UNIT_DIAG, curDim_-blockSize_, blockSize_, one,
                Rtmp->values(), Rtmp->stride(), Y.values(), Y.stride() );
//KMS
sleep(1);
std::cout << "After TRSM" << std::endl;
      		//
              	//  Compute the current update from the Krylov basis; V(:,1:curDim_)*y.
                //
                std::vector<int> index(curDim_-blockSize_);
	      	for ( int i=0; i<curDim_-blockSize_; i++ ) index[i] = i;
      		Teuchos::RCP<const MV> Vjp1 = MVT::CloneView( *V_, index );
      		MVT::MvTimesMatAddMv( one, *Vjp1, Y, zero, *currentUpdate );



      		//
              	//  Add in portion of update from recycled subspace U; U(:,1:recycledBlocks_)*B*y.
                //
                if (U_ != Teuchos::null) {
        		SDM z(recycledBlocks_,blockSize_);
        		SDM subB( Teuchos::View, *B_, recycledBlocks_, curDim_-blockSize_ );
        		z.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS, one, subB, Y, zero );

		        //std::cout << (*U_).MyLength() << " " << (*U_).NumVectors() << " " << subB.numRows() << " " << subB.numCols() << " " << Y.numRows() << " " << Y.numCols()<< " " << curDim_ << " " << recycledBlocks_;	
        		MVT::MvTimesMatAddMv( -one, *U_, z, one, *currentUpdate );
      		}
    	}



    	return currentUpdate;
    }//end getCurrentUpdate() definition

    template<class ScalarType, class MV, class OP>
    void BlockGCRODRIter<ScalarType,MV,OP>::updateLSQR( int dim ) {
	
	int i;
	const ScalarType zero = Teuchos::ScalarTraits<ScalarType>::zero();
	const ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();

	// Get correct dimension based on input "dim"
	// Remember that ortho failures result in an exit before updateLSQR() is called.
	// Therefore, it is possible that dim == curDim_.
	int curDim = curDim_;
    	if ( (dim >= curDim_) && (dim < getMaxSubspaceDim()) ){
      		curDim = dim;
	}

	Teuchos::BLAS<int, ScalarType> blas;

	if(blockSize_ == 1){//if only one right-hand side then use Givens rotations
		//
		// Apply previous rotations and compute new rotations to reduce upper-Hessenberg
		// system to upper-triangular form.
		//
		// QR factorization of Least-Squares system with Givens rotations
		//
		for (i=0; i<curDim-1; i++) {
      			//
        		// Apply previous Givens rotations to new column of Hessenberg matrix
        		//
        		blas.ROT( 1, &(*R_)(i,curDim-1), 1, &(*R_)(i+1, curDim-1), 1, &cs_[i], &sn_[i] );
    		}
    		//
        	// Calculate new Givens rotation
        	//
        	blas.ROTG( &(*R_)(curDim-1,curDim-1), &(*R_)(curDim,curDim-1), &cs_[curDim-1], &sn_[curDim-1] );
    		(*R_)(curDim,curDim-1) = zero;
    		//
        	// Update RHS w/ new transformation
        	//
	       	blas.ROT( 1, &Z_(curDim-1,0), 1, &Z_(curDim,0), 1, &cs_[curDim-1], &sn_[curDim-1] );
	}
	else{//if multiple right-hand sides then use Householder transormations
		//
		//apply previous reflections and compute new reflections to reduce upper-Hessenberg
		//system to upper-triagular form.

		//In Matlab, applying the reflection to a matrix M would look like
		// M_refl = M - 2*v_refl*(v_refl'*M)/(norm(v_refl)^2)

		//In order to take advantage of BLAS while applying reflections to a matrix, we
		//perform it in a two step process
		//1. workvec = M'*v_refl    {using BLAS.GEMV()}
		//2. M_refl = M_refl - 2*v_refl*workvec'/(norm(v_refl)^2)    {using BLAS.GER()}

		Teuchos::RCP< SDM > workmatrix = Teuchos::null;//matrix of column vectors onto which we apply the reflection
		Teuchos::RCP< SDV > workvec = Teuchos::null;//where we store the result of the first step of the 2-step reflection process
		Teuchos::RCP<SDV> v_refl = Teuchos::null;//the reflection vector
		int R_colStart = curDim_-blockSize_;
		Teuchos::RCP< SDM >Rblock = Teuchos::null;

		//
		//Apply previous reflections
		//
		for(i=0; i<lclIter_-1; i++){
			int R_rowStart = i*blockSize_;
			//get a view of the part of R_ effected by these reflections. 
			Teuchos::RCP< SDM > RblockCopy = rcp(new SDM (Teuchos::Copy, *R_, 2*blockSize_,blockSize_, R_rowStart, R_colStart));
			Teuchos::RCP< SDM > RblockView = rcp(new SDM (Teuchos::View, *R_, 2*blockSize_,blockSize_, R_rowStart, R_colStart));
			blas.GEMM(Teuchos::NO_TRANS,Teuchos::NO_TRANS, 2*blockSize_,blockSize_,2*blockSize_,one,House_[i].values(),House_[i].stride(), RblockCopy->values(),RblockCopy -> stride(), zero, RblockView->values(),RblockView -> stride());

		}


		//Get a view of last 2*blockSize entries of entire block to 
		//generate new reflections.
		Rblock = rcp(new SDM (Teuchos::View, *R_, 2*blockSize_,blockSize_, curDim_-blockSize_, curDim_-blockSize_));

		//Calculate and apply the new reflections
		for(i=0; i<blockSize_; i++){
			//
			//Calculating Reflection
			//
			int curcol = (lclIter_ - 1)*blockSize_ + i;//current column of R_
			int lclCurcol = i;//current column of Rblock
			ScalarType signDiag = (*R_)(curcol,curcol) / Teuchos::ScalarTraits<ScalarType>::magnitude((*R_)(curcol,curcol));

                        // Norm of the vector to be reflected.
                        // BLAS returns a ScalarType, but it really should be a magnitude.
			ScalarType nvs = blas.NRM2(blockSize_+1,&((*R_)[curcol][curcol]),1);
			ScalarType alpha = -signDiag*nvs;

			//norm of reflection vector which is just the vector being reflected
			//i.e. v = R_(curcol:curcol+blockSize_,curcol))
			//v_refl = v - alpha*e1
			//norm(v_refl) = norm(v) + alpha^2 - 2*v*alpha
			//store in v_refl

                        // Beware, nvs should have a zero imaginary part (since
                        // it is a norm of a vector), but it may not due to rounding 
                        // error.
			//nvs = nvs + alpha*alpha - 2*(*R_)(curcol,curcol)*alpha;
			//(*R_)(curcol,curcol) -= alpha;

			//Copy relevant values of the current column of R_ into the reflection vector
		 	//Modify first entry
		 	//Take norm of reflection vector
		 	//Square the norm
			v_refl = rcp(new SDV(Teuchos::Copy, &((*R_)(curcol,curcol)), blockSize_ + 1 ));
			(*v_refl)[0] -= alpha;
			nvs = blas.NRM2(blockSize_+1,v_refl -> values(),1);
			nvs *= nvs;

			//
			//Apply new reflector to:
			//1. To subsequent columns of R_ in the current block
			//2. To House[iter_] to store product of reflections for this column
			//3. To the least-squares right-hand side.
			//4. The current column 
			//  		
			//


			//
			//1.
			//
			if(i < blockSize_ - 1){//only do this when there are subsquent columns in the block to apply to
				workvec = Teuchos::rcp(new SDV(blockSize_ - i -1));
				//workvec = Teuchos::rcp(new SDV(2*blockSize_));
				workmatrix = Teuchos::rcp(new SDM (Teuchos::View, *Rblock, blockSize_+1, blockSize_ - i -1, lclCurcol, lclCurcol +1 ) );
				blas.GEMV(Teuchos::TRANS, workmatrix->numRows(), workmatrix->numCols(), one, workmatrix->values(), workmatrix->stride(), v_refl->values(), 1, zero, workvec->values(), 1);
				blas.GER(workmatrix->numRows(),workmatrix->numCols(), -2.*one/nvs, v_refl->values(),1,workvec->values(),1,workmatrix->values(),workmatrix->stride());
			}


			//
			//2.
			//
			workvec = Teuchos::rcp(new SDV(2*blockSize_));
			workmatrix = Teuchos::rcp(new SDM (Teuchos::View, House_[lclIter_ -1], blockSize_+1, 2*blockSize_, i, 0 ) );
			blas.GEMV(Teuchos::TRANS,workmatrix->numRows(),workmatrix->numCols(),one,workmatrix->values(),workmatrix->stride(), v_refl->values(), 1,zero,workvec->values(),1);
			blas.GER(workmatrix->numRows(),workmatrix->numCols(), -2*one/nvs, v_refl -> values(),1,workvec->values(),1,workmatrix->values(),(*workmatrix).stride());

                        //
                        //3.
                        //
			workvec = Teuchos::rcp(new SDV(blockSize_));
			workmatrix = Teuchos::rcp(new SDM (Teuchos::View, Z_, blockSize_+1, blockSize_, curcol, 0 ) );
			blas.GEMV(Teuchos::TRANS, workmatrix->numRows(), workmatrix->numCols(), one, workmatrix-> values(), workmatrix->stride(), v_refl -> values(), 1, zero, workvec->values(), 1);
                        blas.GER((*workmatrix).numRows(),(*workmatrix).numCols(), -2*one/nvs,v_refl -> values(), 1,&((*workvec)[0]),1,(*workmatrix)[0],(*workmatrix).stride());

			//
			//4.
			//
			(*R_)[curcol][curcol] = alpha;
			for(int ii=1; ii<= blockSize_; ii++){
				(*R_)[curcol][curcol+ii] = 0;
			}
		}

	}

    } // end updateLSQR()


}//End Belos Namespace

#endif /* BELOS_BLOCK_GCRODR_ITER_HPP */
