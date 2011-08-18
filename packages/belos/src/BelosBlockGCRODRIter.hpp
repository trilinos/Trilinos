#ifndef BELOS_GCRODR_ITER_HPP
#define BELOS_GCRODR_ITER_HPP
#endif

/*! \file BelosBlockGCRODRIter.hpp
 *     \brief Belos concrete class for performing the block GCRO-DR iteration.
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

namespace Belos{


   //Structure to contain pointers to GCRODRIter state variables.  
   template <class ScalarType, class MV>
   struct BlockGCRODRIterState {
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
   };

   //NEED TO CREATE EXCEPTION CLASSES
   class BlockGCRODRIterInitFailure : public BelosError {
   	public:
   	BlockGCRODRIterInitFailure(const std::string& what_arg) : BelosError(what_arg) {}
   };

   class BlockGCRODRIterOrthoFailure : public BelosError {
  	public:
  	BlockGCRODRIterOrthoFailure(const std::string& what_arg) : BelosError(what_arg) {}
    };


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

       //CONSTRUCTORS AND DESTRUCTORS

       //ONLY ONE CONSTRUCTOR       
       //NEED COMMENTS HERE EXPLAINING WHAT PARAMETERS ARE VALID FOR PARAMETER LIST
       //RIGHT NOW THIS CONSTRUCTOR WILL IGNORE WHATEVER PARAMETER LIST IS PASSED 
       BlockGCRODRIter( const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > &problem,
                        const Teuchos::RCP<OutputManager<ScalarType> > &printer,
                        const Teuchos::RCP<StatusTest<ScalarType,MV,OP> > &tester,
                        const Teuchos::RCP<MatOrthoManager<ScalarType,MV,OP> > &ortho,
                              Teuchos::ParameterList &params );

       //DESTRUCTOR
       virtual ~BlockGCRODRIter() {};
 
       //SOLVER METHODS

       //NEED COMMENTS HERE EXPLAINING HOW THE ITERATION PROCEEDS
       void iterate();

       //Initialize the solver with empty data. THIS MUST BE IMPLEMENTED ACCORDING TO 
       //INTERFACE OF ABSTRACT CLASS iterator.  GCRODR CODE SAYS THIS FUNCTION WILL
       //CAUSE ERROR BECAUSE WE SHOULD INITIALIZE TO VALID STATE
       void initialize() {
       	  BlockGCRODRIterState<ScalarType,MV> empty;
     	  initialize(empty);
       }  

       //NEED COMMENTS EXPLAINING HOW INITIALIZATION WORKS
       void initialize(BlockGCRODRIterState<ScalarType,MV> newstate);
       
       //STATUS METHODS

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

       bool isInitialized(){ return initialized_;};

       //! \brief Get the current iteration count.
       int getNumIters() const { return iter_; };

       //! \brief Reset the iteration count.
       void resetNumIters( int iter = 0 ) { iter_ = iter; };

       //! Get the norms of the residuals native to the solver.
       //! \return A std::vector of length blockSize containing the native residuals.
       Teuchos::RCP<const MV> getNativeResiduals( std::vector<MagnitudeType> *norms ) const;

       //NEED COMMENTS EXPLAINING GET CURRENT UPDATE
       Teuchos::RCP<MV> getCurrentUpdate() const;

       //ACCESSOR METHODS

       //! Get a constant reference to the linear problem.
       const LinearProblem<ScalarType,MV,OP>& getProblem() const { return *lp_; };

       //! Get the maximum number of blocks used by the iterative solver in solving this linear problem.
       int getNumBlocks() const { return numBlocks_; }

       //! Get the blocksize to be used by the iterative solver in solving this linear problem.
       int getBlockSize() const { return blockSize_; };

       int getCurSubspaceDim() const {
		if (!initialized_) return 0;
      		return curDim_;
       };

       int getMaxSubspaceDim() const { return numBlocks_*blockSize_; };

       int getRecycledBlocks() const { return recycledBlocks_; };

       //SET METHODS

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



      private:

      //Classes inputed through constructor that define the linear problem to be solved
      //
      const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> >    lp_;
      const Teuchos::RCP<OutputManager<ScalarType> >          om_;
      const Teuchos::RCP<StatusTest<ScalarType,MV,OP> >       stest_;
      const Teuchos::RCP<OrthoManager<ScalarType,MV> >        ortho_;

      //
      //Algorithmic Paramters
      //

      //numBlocks_ is the size of the allocated space for the Krylov basis, in blocks.
      //blockSize_ is the number of columns in each block Krylov vector.
      int numBlocks_, blockSize_;

      //boolean std::vector indicating which right-hand sides we care about
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

      //
      //Krylov Vectors 
      Teuchos::RCP<MV> V_;

      //
      //Recycled subspace vectors.
      Teuchos::RCP<MV> U_, C_;

      //
      //Projected Operators on the Augmented Krylov Subspace

      //
      //H_: The Projected Matrix from Krylov Factorization AV = VH + C*B, B = C^H*A*V
      Teuchos::RCP<SDM > H_;

      //
      //B_ : Projected matrix from the recycled subspace B = C^H*A*V
      Teuchos::RCP<SDM > B_;

      // QR decomposition of Projected matrices for solving the least squares system HY = RHS
      // R_: Upper triangular reduction of H
      // Z_: Q applied to right-hand side of the least squares system
      Teuchos::RCP<SDM> R_;
      SDM Z_;

      // File stream variables to use Mike Parks matlab output codes
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

        //PARAMETER LIST EXCEPTIONS WOULD GO HERE BUT WE ARE GOING TO USE THE NEW
        //METHOD OF PARAMETER LIST VERIFICATION
        int nb = Teuchos::getParameter<int>(params, "Num Blocks");
	int rb = Teuchos::getParameter<int>(params, "Recycled Blocks");
	int bs = Teuchos::getParameter<int>(params, "Block Size");
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

	// Force call to setsize to ensure internal storage is correct dimension
	setSize( recycledBlocks_, numBlocks_ );

	Teuchos::RCP<MV> Vnext;
	Teuchos::RCP<const MV> Vprev;
	std::vector<int> curind(blockSize_);

	// z_ must be zeroed out in order to compute Givens rotations correctly
	Z_.putScalar(0.0);

	// Orthonormalize the new V_0
	for(int i = 0; i<blockSize_; i++){curind[i] = i;};
	
	Vnext = MVT::CloneViewNonConst(*V_,curind);
	//Orthonormalize Initial Columns
	//Store orthogonalization coefficients in Z0
	Teuchos::RCP<SDM > Z0 =
               Teuchos::rcp( new SDM(blockSize_,blockSize_) );
	int rank = ortho_->normalize(*Vnext,Z0);



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
	while(stest_->checkStatus(this) != Passed && curDim_+blockSize_-1 < numBlocks_*blockSize_){
		lclIter_++;
		iter_++;
//KMS
std::cout << "Iter=" << iter_ << std::endl << "lclIter=" << lclIter_ <<  std::endl;
std::cout << "Here are the current residuals" << std::endl;
{
	std::vector<MagnitudeType> norms;
	getNativeResiduals( &norms );
	for(int jj=0; jj<norms.size(); jj++){
		std::cout << "norms[" << jj << "]=" << norms[jj] << std::endl;
	}
}
		//int lclDim = curDim_ + 1;



		int HFirstCol = curDim_-blockSize_;//First column of H we need view of
		int HLastCol = HFirstCol + blockSize_-1 ;//last column of H we need a view of
		int HLastOrthRow = HLastCol;//Last row of H we will put orthog coefficients in
		int HFirstNormRow = HLastOrthRow + 1;//First row of H where normalization matrix goes
//KMS
std::cout << "curDim_ = " << curDim_ << ", HFirstCol = " << HFirstCol << ", HLastCol =  " << HLastCol <<", HLastOrthRow =  " << HLastOrthRow << ", HFirstNormRow =  " << HFirstNormRow << std::endl;		
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
		// Compute the next std::vector in the Krylov basis:  Vnext = Op*Vprev
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
//KMS

std::cout << "This is subR2 a "<< (lclIter_+1)*blockSize_ << " by " << blockSize_ << " submatrix starting at row "<< 0  << " and column " << HFirstCol << std::endl;
subR2.print(std::cout);
std::cout << "This is subH2 a "<< (lclIter_+1)*blockSize_ << " by " << blockSize_ << " submatrix starting at row "<< 0  << " and column " << HFirstCol << std::endl;
subH2.print(std::cout);
std::cout << "This is R!" << std::endl;
R_->print(std::cout);
std::cout << "This is H!" << std::endl;
H_->print(std::cout);

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
	if(curDim_==0) {
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


		//
		// Solve the least squares problem.
		// Observe that in calling TRSM, we use the value
		// curDim_ -blockSize_. This is because curDim_ has
		// already been incremented but the proper size of R is still 
		// based on the previous value.
		//
		blas.TRSM( Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
                Teuchos::NON_UNIT_DIAG, curDim_-blockSize_, blockSize_, one,
                R_->values(), R_->stride(), Y.values(), Y.stride() );
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
sprintf(filename,"Rtril.mat");
ofs.open(filename);
RblockView->matlab(ofs);
ofs.close();

sprintf(filename,"Houseitril.mat");
ofs.open(filename);
House_[i].matlab(ofs);
ofs.close();

			blas.GEMM(Teuchos::NO_TRANS,Teuchos::NO_TRANS, 2*blockSize_,blockSize_,2*blockSize_,one,House_[i].values(),House_[i].stride(), RblockCopy->values(),RblockCopy -> stride(), zero, RblockView->values(),RblockView -> stride());

sprintf(filename,"RtrilHousePrev.mat");
ofs.open(filename);
RblockView->matlab(ofs);
ofs.close();
		}


		//Get a view of last 2*blockSize entries of entire block to 
		//generate new reflections.
		Rblock = rcp(new SDM (Teuchos::View, *R_, 2*blockSize_,blockSize_, curDim_-blockSize_, curDim_-blockSize_));
//KMS
std::ofstream ofs;
char filename[30];



		//Calculate and apply the new reflections
		for(i=0; i<blockSize_; i++){
			//
			//Calculating Reflection
			//
			int curcol = (lclIter_ - 1)*blockSize_ + i;//current column of R_
			int lclCurcol = i;//current column of Rblock
std::cout << "Current Column = " << curcol << std::endl;
			ScalarType signDiag = (*R_)(curcol,curcol) / Teuchos::ScalarTraits<ScalarType>::magnitude((*R_)(curcol,curcol));
sprintf(filename,"Rtril.mat");
ofs.open(filename);
Rblock->matlab(ofs);
ofs.close();

sprintf(filename,"Vtril.mat");
ofs.open(filename);
MVT::MvPrint(*V_, ofs);
ofs.close();

                        // Norm of the vector to be reflected.
                        // BLAS returns a ScalarType, but it really should be a magnitude.
			ScalarType nvs = blas.NRM2(blockSize_+1,&((*R_)[curcol][curcol]),1);
			ScalarType alpha = -signDiag*nvs;

			//norm of reflection vector which is just the vector being reflected
sprintf(filename,"RBlockmodtril.mat");
ofs.open(filename);
Rblock->matlab(ofs);
ofs.close();		//with the first entry modified by alpha
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

sprintf(filename,"RBlockmodtril.mat");
ofs.open(filename);
Rblock->matlab(ofs);
ofs.close();

sprintf(filename,"vreflTril.mat");
ofs.open(filename);
v_refl->matlab(ofs);
ofs.close();
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



sprintf(filename,"workMatrixTril.mat");
ofs.open(filename);
workmatrix->matlab(ofs);
ofs.close();

std::cout << "numRows = " << workmatrix->numRows() << std::endl;
std::cout << "numCols = " << workmatrix->numCols() << std::endl;
std::cout << "stride = " << workmatrix->stride() << std::endl;
std::cout << "v_refl length = " << v_refl->length() << std::endl;
std::cout << "workvec length = " << workvec->length() << std::endl;

				blas.GEMV(Teuchos::TRANS, workmatrix->numRows(), workmatrix->numCols(), one, workmatrix->values(), workmatrix->stride(), v_refl->values(), 1, zero, workvec->values(), 1);

std::cout << "This is a workmatrix view of Rblock" << std::endl;
workmatrix->print(std::cout);
std::cout << "This is vrefl" << std::endl;
v_refl->print(std::cout);
std::cout << "This is workvec" << std::endl;
workvec->print(std::cout);
				blas.GER(workmatrix->numRows(),workmatrix->numCols(), -2.*one/nvs, v_refl->values(),1,workvec->values(),1,workmatrix->values(),workmatrix->stride());
workmatrix->print(std::cout);
			}

sprintf(filename,"v_reflTril.mat");
ofs.open(filename);
v_refl->matlab(ofs);
ofs.close();

sprintf(filename,"RtrilHouse.mat");
ofs.open(filename);
Rblock->matlab(ofs);
ofs.close();

			//
			//2.
			//
			workvec = Teuchos::rcp(new SDV(2*blockSize_));
			workmatrix = Teuchos::rcp(new SDM (Teuchos::View, House_[lclIter_ -1], blockSize_+1, 2*blockSize_, i, 0 ) );

sprintf(filename,"workvecTril.mat");
ofs.open(filename);
workvec->matlab(ofs);
ofs.close();

sprintf(filename,"v_reflTril.mat");
ofs.open(filename);
v_refl->matlab(ofs);
ofs.close();

sprintf(filename,"workMatrixTril.mat");
ofs.open(filename);
workmatrix->matlab(ofs);
ofs.close();
//exit(1);
			blas.GEMV(Teuchos::TRANS,workmatrix->numRows(),workmatrix->numCols(),one,workmatrix->values(),workmatrix->stride(), v_refl->values(), 1,zero,workvec->values(),1);
			blas.GER(workmatrix->numRows(),workmatrix->numCols(), -2*one/nvs, v_refl -> values(),1,workvec->values(),1,workmatrix->values(),(*workmatrix).stride());

sprintf(filename,"HouseTril.mat");
ofs.open(filename);
House_[lclIter_ - 1].matlab(ofs);
ofs.close();

sprintf(filename,"workvecTril.mat");
ofs.open(filename);
workvec->matlab(ofs);
ofs.close();
                        //
                        //3.
                        //
			workvec = Teuchos::rcp(new SDV(blockSize_));
			workmatrix = Teuchos::rcp(new SDM (Teuchos::View, Z_, blockSize_+1, blockSize_, curcol, 0 ) );
			blas.GEMV(Teuchos::TRANS, workmatrix->numRows(), workmatrix->numCols(), one, workmatrix-> values(), workmatrix->stride(), v_refl -> values(), 1, zero, workvec->values(), 1);
                        blas.GER((*workmatrix).numRows(),(*workmatrix).numCols(), -2*one/nvs,v_refl -> values(), 1,&((*workvec)[0]),1,(*workmatrix)[0],(*workmatrix).stride());

sprintf(filename,"ZTril.mat");
ofs.open(filename);
Z_.matlab(ofs);
ofs.close();


			//
			//4.
			//
			(*R_)[curcol][curcol] = alpha;
			for(int ii=1; ii<= blockSize_; ii++){
				(*R_)[curcol][curcol+ii] = 0;
			}
		}
sprintf(filename,"HouseTril.mat");
ofs.open(filename);
House_[lclIter_].matlab(ofs);
ofs.close();

sprintf(filename,"RtrilHouse.mat");
ofs.open(filename);
Rblock->matlab(ofs);
ofs.close();

sprintf(filename,"ZTril.mat");
ofs.open(filename);
Z_.matlab(ofs);
ofs.close();

std::cout << "******************This is Z after iteration " << lclIter_ << "**************************" << std::endl;
Z_.print(std::cout);
std::cout << "*******************************************************************" << std::endl;
/*
std::cout << "******************This is R after iteration " << lclIter_ << "**************************" << std::endl;
Rblock->print(std::cout);
std::cout << "*******************************************************************" << std::endl;
*/
//exit(1);
	}


    } // end updateLSQR()


}//End Belos Namespace
