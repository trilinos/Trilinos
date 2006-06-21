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

#ifndef BELOS_LINEAR_PROBLEM_HPP
#define BELOS_LINEAR_PROBLEM_HPP

// Define if you want to see a dump of information for detailed debugging
//#define BELOS_LINEAR_PROBLEM_DUMP_COUT

// Define if you want to see a dump multi-vectors as well
//#define BELOS_LINEAR_PROBLEM_DUMP_COUT_ALL

#include "Belos_LinearProblemSetup.hpp"
#include "Belos_deflate.hpp"
#include "TSFCoreLinearOpHandle.hpp"
#include "TSFCoreMultiVectorStdOps.hpp"
#include "TSFCoreMultiplicativeLinearOp.hpp"
#ifdef BELOS_LINEAR_PROBLEM_DUMP_COUT
#include "TSFCoreTestingTools.hpp"
#endif
#include "TSFCoreAssertOp.hpp"
#include "Teuchos_arrayArg.hpp"

namespace Belos {

///
/** Standard concrete implementation of a linear problem.
 *
 * ToDo: Finish documentation!
 */
template<class Scalar>
class LinearProblem : public LinearProblemSetup<Scalar> {
public:

	///
	/** Default constructor.
	 *
	 * Postconditions:<ul>
	 * <li><tt>this->getBlockSize()==1</tt>
	 * <li><tt>this->getTotalNumRhs()==0</tt>
	 * <li><tt>this->getCurrNumRhs()==0</tt>
	 * <li> ...
	 * </ul>
	 */
	LinearProblem();

	/** @name Overridden from LinearProblemState */
	//@{

	///
	int getCurrNumRhs() const;
	///
	int getCurrBlockSize() const;
	///
	void getCurrRhsIndexes(
		const int              currNumRhs
		,int                   currRhsIndexes[]
		) const;
	///
	const TSFCore::MultiVector<Scalar>& getCurrRhs() const;
	///
	const TSFCore::MultiVector<Scalar>& getCurrInitResidual() const;
	///
	bool isCurrInitResidualComputed() const;
	///
	const TSFCore::MultiVector<Scalar>& getCurrInitLhs() const;
	///
	const TSFCore::MultiVector<Scalar>& getCurrLhs() const;
	///
	bool isCurrLhsUpdated() const;
	///
	const TSFCore::MultiVector<Scalar>& getCurrResidual() const;
	///
	bool isCurrResidualComputed() const;

	//@}

	/** @name Overridden from LinearProblemIteration */
	//@{

	///
	int getTotalNumRhs() const;
	///
	TSFCore::LinearOpHandle<Scalar> getOperator() const;
	///
	EOpSymmetry getOperatorSymmetry() const;
	///
	TSFCore::LinearOpHandle<Scalar> getRightPrec() const;
	///
	EOpSymmetry getRightPrecSymmetry() const;
	///
	TSFCore::LinearOpHandle<Scalar> getLeftPrec() const;
	///
	EOpSymmetry getLeftPrecSymmetry() const;
	///
	TSFCore::LinearOpHandle<Scalar> getCompositeOperator() const;
	///
	const TSFCore::Vector<Scalar>* getRightScaling() const;
	///
	const TSFCore::Vector<Scalar>* getLeftScaling() const;
	///
	const TSFCore::MultiVector<Scalar>& getRhs() const;
	///
	TSFCore::MultiVector<Scalar>& getLhs() const;
	///
	void setBlockSize( const int blockSize );
	///
	int getBlockSize() const;
	///
	void setAugmentationAllowed( const bool augmentationAllowed );
	///
	bool getAugmentationAllowed() const;
	///
	StatusTest<Scalar>* getStatusTest();
	///
	bool setupCurrSystem();
	///
	RefCountPtr<TSFCore::MultiVector<Scalar> > getCurrLhs();
	///
	void setCurrLhsUpdated( const bool isCurrLhsUpdated );
	///
	void updateCurrLhs( const TSFCore::MultiVector<Scalar> &nativeLhsStep );
	///
	void setCurrLhs( const TSFCore::MultiVector<Scalar> &nativeLhs );
	///
	RefCountPtr<TSFCore::MultiVector<Scalar> > getCurrResidual();
	///
	void setCurrResidualComputed( const bool isCurrResidualComputed );
	///
	void computeCurrCompositeResidual(
		const TSFCore::MultiVector<Scalar>    *currLhs
		,const TSFCore::MultiVector<Scalar>   *currRhs
		,TSFCore::MultiVector<Scalar>         *currCompositeResidual
		);
	///
	void deflate(
		const BasicIterationState<Scalar>        &bis
		,const int                               numCurrRhsToRemove
		,const int                               currRhsIndexesToRemove[]
		);
	///
	void restart( const BasicIterationState<Scalar> &bis );
	///
	void finalizeCurrSystem( const BasicIterationState<Scalar> &bis );

	//@}

	/** @name Overridden from LinearProblemSetup */
	//@{

	///
	void initialize(
		const TSFCore::LinearOpHandle<Scalar>                    &op
		,const EOpSymmetry                                       symmetry
		,const RefCountPtr<const TSFCore::MultiVector<Scalar> >  &rhs
		,const RefCountPtr<TSFCore::MultiVector<Scalar> >        &lhs
		);
	///
	void setOperator( const TSFCore::LinearOpHandle<Scalar> &op, const EOpSymmetry symmetry );
	///
	TSFCore::LinearOpHandle<Scalar> getOperator();
	///
	void setRightPrec( const TSFCore::LinearOpHandle<Scalar> &rightPrec, const EOpSymmetry symmetry );
	///
	TSFCore::LinearOpHandle<Scalar> getRightPrec();
	///
	void setLeftPrec( const TSFCore::LinearOpHandle<Scalar> &leftPrec, const EOpSymmetry symmetry );
	///
	TSFCore::LinearOpHandle<Scalar> getLeftPrec();
	///
	void setRhs( const RefCountPtr<const TSFCore::MultiVector<Scalar> > &rhs );
	///
	RefCountPtr<const TSFCore::MultiVector<Scalar> > getRhs();
	///
	void setLhs( const RefCountPtr<TSFCore::MultiVector<Scalar> > &lhs );
	///
	RefCountPtr<TSFCore::MultiVector<Scalar> > getLhs();
	///
	void setStatusTest( const RefCountPtr<StatusTest<Scalar> > &statusTest );
	///
	void completeSetup();
	///
	void uninitialize();

	//@}

private:

	// /////////////////////////////////////
	// Private types

	typedef Teuchos::ScalarTraits<Scalar> ST;
	typedef typename ST::magnitudeType ScalarMag;

	// /////////////////////////////////////
	// Private data members

	// Members for overall set of systems
	int                                               totalNumRhs_;
	TSFCore::LinearOpHandle<Scalar>                   operator_;
	EOpSymmetry                                       operatorSymmetry_;
	TSFCore::LinearOpHandle<Scalar>                   rightPrec_;
	EOpSymmetry                                       rightPrecSymmetry_;
	TSFCore::LinearOpHandle<Scalar>                   leftPrec_;
	EOpSymmetry                                       leftPrecSymmetry_;
	RefCountPtr<const TSFCore::Vector<Scalar> >       rightScaling_;
	RefCountPtr<const TSFCore::Vector<Scalar> >       leftScaling_;
	TSFCore::LinearOpHandle<Scalar>                   compositeOperator_;
	RefCountPtr<const TSFCore::MultiVector<Scalar> >  rhs_;
	RefCountPtr<TSFCore::MultiVector<Scalar> >        lhs_;
	RefCountPtr<StatusTest<Scalar> >                  statusTest_;

	// Members for controlling setup of current block	
	int                                               blockSize_;
	bool                                              augmentationAllowed_;

	// Members for current block	
	int                                               currFirstRhsOffset_;
	int                                               currInitNumRhs_;
	int                                               currNumRhs_;
	int                                               currBlockSize_;
	std::vector<int>                                  currRhsIndexes_;
	RefCountPtr<TSFCore::MultiVector<Scalar> >        currRhs_store_;
	RefCountPtr<TSFCore::MultiVector<Scalar> >        currRhs_;
	RefCountPtr<TSFCore::MultiVector<Scalar> >        currInitResidual_store_;
	RefCountPtr<TSFCore::MultiVector<Scalar> >        currInitResidual_;
	mutable bool                                      isCurrInitResidualComputed_;
	RefCountPtr<TSFCore::MultiVector<Scalar> >        currInitLhs_store_;
	RefCountPtr<TSFCore::MultiVector<Scalar> >        currInitLhs_;
	RefCountPtr<TSFCore::MultiVector<Scalar> >        currLhs_store_;
	RefCountPtr<TSFCore::MultiVector<Scalar> >        currLhs_;
	mutable bool                                      currLhsIsCurrInitLhs_;
	mutable bool                                      isCurrLhsUpdated_;
	RefCountPtr<TSFCore::MultiVector<Scalar> >        currResidual_store_;
	RefCountPtr<TSFCore::MultiVector<Scalar> >        currResidual_;
	mutable bool                                      isCurrResidualComputed_;

	// RAB: 2004/09/17: Warning, do not access the values of any
	// xxx_store_ quantity directly because there is a view of these
	// objects that may not be in sync (see the behavior of
	// TSFCore::MultiVector::subView()).

	// //////////////////////////////////////////
	// Private member functions

	void assertTotalNumRhs() const;
	void assertCurrNumRhs() const;

};


// /////////////////////////////
// Inline member functions

template<class Scalar>
inline
void LinearProblem<Scalar>::assertTotalNumRhs() const
{
#ifdef TEUCHOS_DEBUG
	TEST_FOR_EXCEPTION( totalNumRhs_ <= 0, std::logic_error, "LinearProblem<Scalar>: Error, precondition not satisfied!" );
#endif
}

template<class Scalar>
inline
void LinearProblem<Scalar>::assertCurrNumRhs() const
{
#ifdef TEUCHOS_DEBUG
	TEST_FOR_EXCEPTION( currNumRhs_ <= 0, std::logic_error, "LinearProblem<Scalar>: Error, precondition not satisfied!" );
#endif
}

// /////////////////////////////
// Implementation

template<class Scalar>
LinearProblem<Scalar>::LinearProblem()
	:totalNumRhs_(0)
	,blockSize_(1)
	,augmentationAllowed_(true)
	,currNumRhs_(0)
{}

// Overridden from LinearProblemState

template<class Scalar>
int LinearProblem<Scalar>::getCurrNumRhs() const
{
	return currNumRhs_;
}

template<class Scalar>
int LinearProblem<Scalar>::getCurrBlockSize() const
{
	assertCurrNumRhs();
	return currBlockSize_;
}

template<class Scalar>
void LinearProblem<Scalar>::getCurrRhsIndexes(
	const int              currNumRhs
	,int                   currRhsIndexes[]
	) const
{
	assertCurrNumRhs();
#ifdef TEUCHOS_DEBUG
	TEST_FOR_EXCEPT( currNumRhs != currNumRhs_ );
#endif
	std::copy( &currRhsIndexes_[0], &currRhsIndexes_[0]+currNumRhs, currRhsIndexes );
}

template<class Scalar>
const TSFCore::MultiVector<Scalar>& LinearProblem<Scalar>::getCurrRhs() const
{
	assertCurrNumRhs();
	return *currRhs_;
}

template<class Scalar>
const TSFCore::MultiVector<Scalar>& LinearProblem<Scalar>::getCurrInitResidual() const
{
	assertCurrNumRhs();
	TEST_FOR_EXCEPT( leftScaling_.get() ); // Can't handle scaling yet!
	if(!isCurrInitResidualComputed_) {
		// currInitResidual = operator * curLhsInitGuess - currRhs
		assign( &*currInitResidual_, *currRhs_ );
		const ScalarMag lhs_norm_1 = norm_1(*lhs_);
		if( lhs_norm_1 != ST::zero() ) {
			operator_.apply(
				TSFCore::NOTRANS
				,*lhs_->subView(currNumRhs_,&currRhsIndexes_[0]) // Note, \bar{X} is 0 for augmented columns
				,&*currInitResidual_->subView(TSFCore::Range1D(1,currNumRhs_))
				,+ST::one()
				,-ST::one()
				);
		}
		else {
			scale( -ST::one(), &*currInitResidual_ );
		}
		isCurrInitResidualComputed_ = true;
	}
	return *currInitResidual_;
}

template<class Scalar>
bool LinearProblem<Scalar>::isCurrInitResidualComputed() const
{
	assertCurrNumRhs();
	return isCurrInitResidualComputed_;
}

template<class Scalar>
const TSFCore::MultiVector<Scalar>& LinearProblem<Scalar>::getCurrInitLhs() const
{
	assertCurrNumRhs();
	return *currInitLhs_;
}

template<class Scalar>
const TSFCore::MultiVector<Scalar>& LinearProblem<Scalar>::getCurrLhs() const
{
	assertCurrNumRhs();
	return *currLhs_;
}

template<class Scalar>
bool LinearProblem<Scalar>::isCurrLhsUpdated() const
{
	assertCurrNumRhs();
	return isCurrLhsUpdated_;
}

template<class Scalar>
const TSFCore::MultiVector<Scalar>& LinearProblem<Scalar>::getCurrResidual() const
{
#if defined(BELOS_LINEAR_PROBLEM_DUMP_COUT)
	std::cout << "\nLinearProblem<Scalar>::getCurrResidual(...):\n";
#endif
	assertCurrNumRhs();
	TEST_FOR_EXCEPT( leftScaling_.get() ); // Can't handle scaling yet!
#ifdef TEUCHOS_DEBUG
	TEST_FOR_EXCEPT( isCurrLhsUpdated_==false );
#endif
	if(!isCurrResidualComputed_) {
		// currResidual = operator * curLhs - currRhs
		assign( &*currResidual_, *currRhs_ );
#if defined(BELOS_LINEAR_PROBLEM_DUMP_COUT)
		std::cout << "\n||currResidual = currRhs||1 = " << norm_1(*currResidual_) << std::endl;
#endif
#if defined(BELOS_LINEAR_PROBLEM_DUMP_COUT) && (BELOS_LINEAR_PROBLEM_DUMP_COUT_ALL)
		std::cout << "\ncurrResidual = currRhs =\n" << *currResidual_;
#endif
		const ScalarMag currLhs_norm_1 = norm_1(*currLhs_);
		if( currLhs_norm_1 != ST::zero() ) {
#if defined(BELOS_LINEAR_PROBLEM_DUMP_COUT)
			std::cout << "\n||currLhs|| = " << currLhs_norm_1 << std::endl;
#endif
#if defined(BELOS_LINEAR_PROBLEM_DUMP_COUT) && (BELOS_LINEAR_PROBLEM_DUMP_COUT_ALL)
			std::cout << "\ncurrLhs = \n" << *currLhs_;
#endif
			operator_.apply(
				TSFCore::NOTRANS
				,*currLhs_
				,&*currResidual_
				,+ST::one()
				,-ST::one()
				);
#if defined(BELOS_LINEAR_PROBLEM_DUMP_COUT)
			std::cout << "\n||currResidual = operator*currLhs - currResidual|| = " << norm_1(*currResidual_) << std::endl;
#endif
#if defined(BELOS_LINEAR_PROBLEM_DUMP_COUT) && (BELOS_LINEAR_PROBLEM_DUMP_COUT_ALL)
			std::cout << "\ncurrResidual = operator*currLhs - currResidual = \n" << *currResidual_;
#endif
		}
		else {
			scale( -ST::one(), &*currResidual_ );
		}
		isCurrResidualComputed_ = true;
		if(currLhsIsCurrInitLhs_) {
			assign( &*currInitResidual_, *currResidual_ );
			isCurrInitResidualComputed_ = true;
		}
	}
	return *currResidual_;
}

template<class Scalar>
bool LinearProblem<Scalar>::isCurrResidualComputed() const
{
	assertCurrNumRhs();
	return isCurrResidualComputed_;
}

// Overridden from LinearProblemIteration

template<class Scalar>
int LinearProblem<Scalar>::getTotalNumRhs() const
{
	return totalNumRhs_;
}

template<class Scalar>
TSFCore::LinearOpHandle<Scalar> LinearProblem<Scalar>::getOperator() const
{
	assertTotalNumRhs();
	return operator_;
}

template<class Scalar>
EOpSymmetry LinearProblem<Scalar>::getOperatorSymmetry() const
{
	assertTotalNumRhs();
	return operatorSymmetry_;
}

template<class Scalar>
TSFCore::LinearOpHandle<Scalar> LinearProblem<Scalar>::getRightPrec() const
{
	assertTotalNumRhs();
	return rightPrec_;
}

template<class Scalar>
EOpSymmetry LinearProblem<Scalar>::getRightPrecSymmetry() const
{
	assertTotalNumRhs();
	return rightPrecSymmetry_;
}

template<class Scalar>
TSFCore::LinearOpHandle<Scalar> LinearProblem<Scalar>::getLeftPrec() const
{
	assertTotalNumRhs();
	return leftPrec_;
}

template<class Scalar>
EOpSymmetry LinearProblem<Scalar>::getLeftPrecSymmetry() const
{
	assertTotalNumRhs();
	return leftPrecSymmetry_;
}

template<class Scalar>
TSFCore::LinearOpHandle<Scalar> LinearProblem<Scalar>::getCompositeOperator() const
{
	assertTotalNumRhs();
	return compositeOperator_;
}

template<class Scalar>
const TSFCore::Vector<Scalar>* LinearProblem<Scalar>::getRightScaling() const
{
	assertTotalNumRhs();
	return rightScaling_.get();
}

template<class Scalar>
const TSFCore::Vector<Scalar>* LinearProblem<Scalar>::getLeftScaling() const
{
	assertTotalNumRhs();
	return leftScaling_.get();
}

template<class Scalar>
const TSFCore::MultiVector<Scalar>& LinearProblem<Scalar>::getRhs() const
{
	assertTotalNumRhs();
	return *rhs_;
}

template<class Scalar>
TSFCore::MultiVector<Scalar>& LinearProblem<Scalar>::getLhs() const
{
	assertTotalNumRhs();
	return *lhs_;
}

template<class Scalar>
void LinearProblem<Scalar>::setBlockSize( const int blockSize )
{
#ifdef TEUCHOS_DEBUG
	TEST_FOR_EXCEPT( blockSize <= 0 );
#endif
	blockSize_ = blockSize;
	currNumRhs_ = 0;
}

template<class Scalar>
int LinearProblem<Scalar>::getBlockSize() const
{
	return blockSize_;
}

template<class Scalar>
void LinearProblem<Scalar>::setAugmentationAllowed( const bool augmentationAllowed )
{
	augmentationAllowed_ = augmentationAllowed;
}

template<class Scalar>
bool LinearProblem<Scalar>::getAugmentationAllowed() const
{
	return augmentationAllowed_;
}

template<class Scalar>
StatusTest<Scalar>* LinearProblem<Scalar>::getStatusTest()
{
	return statusTest_.get();
}

template<class Scalar>
bool LinearProblem<Scalar>::setupCurrSystem()
{
#ifdef TEUCHOS_DEBUG
	TEST_FOR_EXCEPT( currFirstRhsOffset_  > totalNumRhs_ );
#endif
	// See if we are finished yet or not
	if( currFirstRhsOffset_  == totalNumRhs_ ) return false;
	// Determine the block size and the number of RHSs of the next set block
	const int numRhsRemaining = totalNumRhs_ - currFirstRhsOffset_;
	const int currNumRhs      = min( numRhsRemaining, blockSize_ );
	const int blockSize       = ( augmentationAllowed_ ? blockSize_ : currNumRhs );
	// Setup the current block system
	RefCountPtr<const TSFCore::VectorSpace<Scalar> > rhs_range = rhs_->range();
	RefCountPtr<const TSFCore::VectorSpace<Scalar> > lhs_range = lhs_->range();
	currRhs_store_ = rhs_range->createMembers(blockSize);
	currRhs_ = currRhs_store_;
	currInitLhs_store_ = lhs_range->createMembers(blockSize);
	currInitLhs_ = currInitLhs_store_;
	currLhs_store_ = lhs_range->createMembers(blockSize);
	currLhs_ = currLhs_store_;
	currInitResidual_store_ = rhs_range->createMembers(blockSize);
	currInitResidual_ = currInitResidual_store_;
	currResidual_store_ = rhs_range->createMembers(blockSize);
	currResidual_ = currResidual_store_;
	const TSFCore::Range1D currRhsRng(1,currNumRhs);
	const TSFCore::Range1D currOrigRhsRng = currRhsRng + currFirstRhsOffset_;
	assign( &*currRhs_->subView(currRhsRng), *rhs_->subView(currOrigRhsRng) );
	assign( &*currLhs_->subView(currRhsRng), *lhs_->subView(currOrigRhsRng) );
	isCurrLhsUpdated_ = true;
	currLhsIsCurrInitLhs_ = true;
	if( blockSize > currNumRhs ) {
		const TSFCore::Range1D currAugRng(currNumRhs+1,blockSize);
		randomize( Scalar(-ST::one()), Scalar(+ST::one()), &*currRhs_->subView(currAugRng) ); // Random RHS for augmented systems
		assign( &*currLhs_->subView(currAugRng), ST::zero() );                                // Zero initial guess for augmented systems
	}
	assign( &*currInitLhs_, *currLhs_ );
	// Initialize list of active RHSs (one-based)
	currRhsIndexes_.resize(currNumRhs);
	for( int k = 0; k < currNumRhs; ++k ) currRhsIndexes_[k] = currFirstRhsOffset_+k+1;
	// Mark derived quantities as uninitialized (unscaled unpreconditioned residuals)
	isCurrInitResidualComputed_ = false;
	isCurrResidualComputed_ = false;
	// Finish up initialization
	currInitNumRhs_  = currNumRhs;
	currNumRhs_      = currNumRhs;
	currBlockSize_   = blockSize;
	return true;
}

template<class Scalar>
RefCountPtr<TSFCore::MultiVector<Scalar> >
LinearProblem<Scalar>::getCurrLhs()
{
	assertCurrNumRhs();
	currLhsIsCurrInitLhs_ = false; // Just to be safe!
	return currLhs_;
}

template<class Scalar>
void LinearProblem<Scalar>::setCurrLhsUpdated( const bool isCurrLhsUpdated )
{
	assertCurrNumRhs();
	isCurrLhsUpdated_ = isCurrLhsUpdated;
	isCurrResidualComputed_ = false;
	if(isCurrLhsUpdated) currLhsIsCurrInitLhs_ = false;
}

template<class Scalar>
void LinearProblem<Scalar>::updateCurrLhs( const TSFCore::MultiVector<Scalar> &nativeLhsStep )
{
	assertCurrNumRhs();
	TEST_FOR_EXCEPT( rightScaling_.get() ); // Can't handle scaling yet!
	// Update from the current LHS!
	assign( &*currLhs_, *currInitLhs_ );
	if(rightPrec_.op().get())
		rightPrec_.apply(TSFCore::NOTRANS,nativeLhsStep,&*currLhs_,ST::one(),ST::one());
	else
		update( ST::one(), nativeLhsStep, &*currLhs_ );
	isCurrLhsUpdated_ = true;
	isCurrResidualComputed_ = false;
	currLhsIsCurrInitLhs_ = false;
}

template<class Scalar>
void LinearProblem<Scalar>::setCurrLhs( const TSFCore::MultiVector<Scalar> &nativeLhs )
{
#if defined(BELOS_LINEAR_PROBLEM_DUMP_COUT)
	std::cout << "\nLinearProblem<Scalar>::setCurrLhs(nativeLhs):\n";
	std::cout << "\n||nativeLhs|| = " << norm_1(nativeLhs) << std::endl;
#endif
	assertCurrNumRhs();
	TEST_FOR_EXCEPT( rightScaling_.get() ); // Can't handle scaling yet!
	if(rightPrec_.op().get())
		rightPrec_.apply(TSFCore::NOTRANS,nativeLhs,&*currLhs_);
	else
		assign( &*currLhs_, nativeLhs );
#if defined(BELOS_LINEAR_PROBLEM_DUMP_COUT)
	std::cout << "\n||currLhs|| = " << norm_1(*currLhs_) << std::endl;
#endif
	isCurrLhsUpdated_ = true;
	isCurrResidualComputed_ = false;
	currLhsIsCurrInitLhs_ = false;
}

template<class Scalar>
RefCountPtr<TSFCore::MultiVector<Scalar> >
LinearProblem<Scalar>::getCurrResidual()
{
	assertCurrNumRhs();
	return currResidual_;
}

template<class Scalar>
void LinearProblem<Scalar>::setCurrResidualComputed( const bool isCurrResidualComputed )
{
	assertCurrNumRhs();
	isCurrResidualComputed_ = isCurrResidualComputed;
}

template<class Scalar>
void LinearProblem<Scalar>::computeCurrCompositeResidual(
	const TSFCore::MultiVector<Scalar>    *currLhs
	,const TSFCore::MultiVector<Scalar>   *currRhs
	,TSFCore::MultiVector<Scalar>         *currCompositeResidual
	)
{
	assertCurrNumRhs();
	if( currLhs==NULL && currRhs==NULL ) {
#ifdef TEUCHOS_DEBUG
		TEST_FOR_EXCEPT( !isCurrLhsUpdated_ );
#endif
		const TSFCore::MultiVector<Scalar> &currResidual = const_cast<const LinearProblem<Scalar>*>(this)->getCurrResidual();
		if( leftPrec_.op().get() )
			leftPrec_.apply( TSFCore::NOTRANS, currResidual, currCompositeResidual );
		else
			assign( currCompositeResidual, currResidual );
	}
	else {
		TEST_FOR_EXCEPT(true); // Can't handle this case yet but it is faily easy to do so!
	}
}

template<class Scalar>
void LinearProblem<Scalar>::deflate(
	const BasicIterationState<Scalar>        &bis
	,const int                               numCurrRhsToRemove
	,const int                               currRhsIndexesToRemove[]
	)
{
	assertCurrNumRhs();
#ifdef BELOS_LINEAR_PROBLEM_DUMP_COUT
	std::cout
		<< "\nLinearProblem<Scalar>::deflate(...):\n";
	std::cout
		<< "  currRhsIndexes after ={";
	for(int k=0;k<currNumRhs_;++k) { if(k!=0) std::cout << ","; std::cout << currRhsIndexes_[k]; }
	std::cout
		<< "}\n"
		<< "  currRhsIndexesToRemove={";
	for(int k=0;k<numCurrRhsToRemove;++k) { if(k!=0) std::cout << ","; std::cout << currRhsIndexesToRemove[k]; }
	std::cout << "}\n";
#endif
	if(!numCurrRhsToRemove) return; // Nothing to deflate!
	// Get the indexes in the original RHS of the converged RHSs
	std::vector<int> origRhsIndexesToRemove(numCurrRhsToRemove);
	for( int k = 0; k < numCurrRhsToRemove; ++k ) origRhsIndexesToRemove[k] = currRhsIndexes_[currRhsIndexesToRemove[k]-1];
	// Copy converged LHS into full LHS
	if(!isCurrLhsUpdated_) bis.forceCurrLhsUpdate();
	assign( &*lhs_->subView(numCurrRhsToRemove,&origRhsIndexesToRemove[0]), *currLhs_->subView(numCurrRhsToRemove,&currRhsIndexesToRemove[0]) );
	// Release multi-vector views so that we can modify full storage multi-vectors
	currRhs_ = currInitResidual_ = currInitLhs_ = currResidual_ = null;
	// Deflate the multi-vectors
	::Belos::deflate(
		numCurrRhsToRemove, currRhsIndexesToRemove, currBlockSize_
		,5
		,Teuchos::arrayArg<TSFCore::MultiVector<Scalar>*>(
			&*currRhs_store_, &*currInitResidual_store_, &*currInitLhs_store_
			,&*currLhs_store_, &*currResidual_store_
			)()
		);
	// Deflate index set of current RHSs
	::Belos::deflate(
		numCurrRhsToRemove, currRhsIndexesToRemove, currBlockSize_
		,1,Teuchos::arrayArg<int*>(&currRhsIndexes_[0])()
		);
	// Update the dimensions
	currNumRhs_ -= numCurrRhsToRemove;
	currBlockSize_ -= numCurrRhsToRemove;
#ifdef BELOS_LINEAR_PROBLEM_DUMP_COUT
	std::cout
		<< "  currRhsIndexes after ={";
	for(int k=0;k<currNumRhs_;++k) { if(k!=0) std::cout << ","; std::cout << currRhsIndexes_[k]; }
	std::cout << "}\n";
#endif
	// Reform the views
	const TSFCore::Range1D currRng(1,currBlockSize_);
	currRhs_ = currRhs_store_->subView(currRng);
	currInitResidual_ = currInitResidual_store_->subView(currRng);
	currInitLhs_ = currInitLhs_store_->subView(currRng);
	currLhs_ = currLhs_store_->subView(currRng);
	currResidual_ = currResidual_store_->subView(currRng);
}

template<class Scalar>
void LinearProblem<Scalar>::restart( const BasicIterationState<Scalar> &bis )
{
	assertCurrNumRhs();
	if(!isCurrLhsUpdated_) bis.forceCurrLhsUpdate();
	assign( &*currInitLhs_, *currLhs_ );
}

template<class Scalar>
void LinearProblem<Scalar>::finalizeCurrSystem( const BasicIterationState<Scalar> &bis )
{
	if(!isCurrLhsUpdated_) bis.forceCurrLhsUpdate();
	assign( &*lhs_->subView(currNumRhs_,&currRhsIndexes_[0]), *currLhs_->subView(TSFCore::Range1D(1,currNumRhs_)) );
	currFirstRhsOffset_ += currInitNumRhs_;
	currNumRhs_ = 0; // Set as uninitialized
}

// Overridden from LinearProblemSetup

template<class Scalar>
void LinearProblem<Scalar>::initialize(
	const TSFCore::LinearOpHandle<Scalar>                    &op
	,const EOpSymmetry                                       symmetry
	,const RefCountPtr<const TSFCore::MultiVector<Scalar> >  &rhs
	,const RefCountPtr<TSFCore::MultiVector<Scalar> >        &lhs
	)
{
#ifdef TEUCHOS_DEBUG
	TEST_FOR_EXCEPT( op.op().get() == NULL );
	TEST_FOR_EXCEPT( rhs.get() == NULL );
	TEST_FOR_EXCEPT( lhs.get() == NULL );
#endif
	operator_           = op;
	operatorSymmetry_   = symmetry;
	rhs_                = rhs;
	lhs_                = lhs;
	completeSetup();
}

template<class Scalar>
void LinearProblem<Scalar>::setOperator( const TSFCore::LinearOpHandle<Scalar> &op, const EOpSymmetry symmetry )
{
#ifdef TEUCHOS_DEBUG
	TEST_FOR_EXCEPT( op.op().get() == NULL );
#endif
	operator_           = op;
	operatorSymmetry_   = symmetry;
	totalNumRhs_        = 0;
	currNumRhs_         = 0;
}

template<class Scalar>
TSFCore::LinearOpHandle<Scalar> LinearProblem<Scalar>::getOperator()
{
	return operator_;
}

template<class Scalar>
void LinearProblem<Scalar>::setRightPrec( const TSFCore::LinearOpHandle<Scalar> &rightPrec, const EOpSymmetry symmetry )
{
#ifdef TEUCHOS_DEBUG
	TEST_FOR_EXCEPT( rightPrec.op().get() == NULL );
#endif
	rightPrec_           = rightPrec;
	rightPrecSymmetry_   = symmetry;
	totalNumRhs_         = 0;
	currNumRhs_          = 0;
}

template<class Scalar>
TSFCore::LinearOpHandle<Scalar> LinearProblem<Scalar>::getRightPrec()
{
	return rightPrec_;
}

template<class Scalar>
void LinearProblem<Scalar>::setLeftPrec( const TSFCore::LinearOpHandle<Scalar> &leftPrec, const EOpSymmetry symmetry )
{
#ifdef TEUCHOS_DEBUG
	TEST_FOR_EXCEPT( leftPrec.op().get() == NULL );
#endif
	leftPrec_           = leftPrec;
	leftPrecSymmetry_   = symmetry;
	totalNumRhs_        = 0;
	currNumRhs_         = 0;
}

template<class Scalar>
TSFCore::LinearOpHandle<Scalar> LinearProblem<Scalar>::getLeftPrec()
{
	return leftPrec_;
}

template<class Scalar>
void LinearProblem<Scalar>::setRhs( const RefCountPtr<const TSFCore::MultiVector<Scalar> > &rhs )
{
#ifdef TEUCHOS_DEBUG
	TEST_FOR_EXCEPT( rhs.get() == NULL );
#endif
	rhs_ = rhs;
	totalNumRhs_ = 0;
	currNumRhs_ = 0;
}

template<class Scalar>
RefCountPtr<const TSFCore::MultiVector<Scalar> >
LinearProblem<Scalar>::getRhs()
{
	return rhs_;
}

template<class Scalar>
void LinearProblem<Scalar>::setLhs( const RefCountPtr<TSFCore::MultiVector<Scalar> > &lhs )
{
#ifdef TEUCHOS_DEBUG
	TEST_FOR_EXCEPT( lhs.get() == NULL );
#endif
	lhs_ = lhs;
	totalNumRhs_ = 0;
	currNumRhs_ = 0;
}

template<class Scalar>
RefCountPtr<TSFCore::MultiVector<Scalar> >
LinearProblem<Scalar>::getLhs()
{
	return lhs_;
}

template<class Scalar>
void LinearProblem<Scalar>::setStatusTest( const RefCountPtr<StatusTest<Scalar> > &statusTest )
{
	statusTest_ = statusTest;
	totalNumRhs_ = 0;
	currNumRhs_ = 0;
}

template<class Scalar>
void LinearProblem<Scalar>::completeSetup()
{
#ifdef TEUCHOS_DEBUG
	const char funcName[] = "LinearProblem<Scalar>::completeSetup()";
	TEST_FOR_EXCEPTION( operator_.op().get() == NULL, std::logic_error, funcName<<", Error!" );
	TEST_FOR_EXCEPTION( rhs_.get() == NULL, std::logic_error, funcName<<", Error!"  );
	TEST_FOR_EXCEPTION( lhs_.get() == NULL, std::logic_error, funcName<<", Error!"  );
	TSFCORE_ASSERT_VEC_SPACES( funcName, *rhs_->domain(), *lhs_->domain() );
	if( rightPrec_.op().get() ) {
		TSFCORE_ASSERT_VEC_SPACES( funcName, *rightPrec_.domain(), *lhs_->range() );
		TSFCORE_ASSERT_VEC_SPACES( funcName, *operator_.domain(), *rightPrec_.range() );
	}
	else {
		TSFCORE_ASSERT_VEC_SPACES( funcName, *operator_.domain(), *lhs_->range() );
	}
	if( leftPrec_.op().get() ) {
		TSFCORE_ASSERT_VEC_SPACES( funcName, *leftPrec_.range(), *rhs_->range() );
		TSFCORE_ASSERT_VEC_SPACES( funcName, *leftPrec_.domain(), *operator_.range() );
	}
	else {
		TSFCORE_ASSERT_VEC_SPACES( funcName, *operator_.range(), *rhs_->range() );
	}
	if( leftScaling_.get() ) {
		TSFCORE_ASSERT_VEC_SPACES( funcName, *rhs_->range(), *leftScaling_->space() );
	}
	if( rightScaling_.get() ) {
		TSFCORE_ASSERT_VEC_SPACES( funcName, *lhs_->range(), *rightScaling_->space() );
	}
#endif
	// Build the composite operator.
	if( leftPrec_.op().get() ) {
		if( rightPrec_.op().get() ) {
			compositeOperator_ = TSFCore::LinearOpHandle<Scalar>(
				rcp(
					new TSFCore::MultiplicativeLinearOp<Scalar>(
						3, Teuchos::arrayArg<TSFCore::LinearOpHandle<Scalar> >(leftPrec_,operator_,rightPrec_)()
						)
					)
				);
		}
		else {
			compositeOperator_ = TSFCore::LinearOpHandle<Scalar>(
				rcp(
					new TSFCore::MultiplicativeLinearOp<Scalar>(
						2, Teuchos::arrayArg<TSFCore::LinearOpHandle<Scalar> >(leftPrec_,operator_)()
						)
					)
				);
		}
	}
	else {
		if( rightPrec_.op().get() ) {
			compositeOperator_ = TSFCore::LinearOpHandle<Scalar>(
				rcp(
					new TSFCore::MultiplicativeLinearOp<Scalar>(
						2, Teuchos::arrayArg<TSFCore::LinearOpHandle<Scalar> >(operator_,rightPrec_)()
						)
					)
				);
		}
		else {
			compositeOperator_ = operator_;
		}
	}
	// Finish initialization
	totalNumRhs_ = rhs_->domain()->dim();
	currFirstRhsOffset_ = 0;
	currNumRhs_ = 0;
}

template<class Scalar>
void LinearProblem<Scalar>::uninitialize()
{
	totalNumRhs_ = 0;
	operator_ = TSFCore::LinearOpHandle<Scalar>();
	leftPrec_ = TSFCore::LinearOpHandle<Scalar>();
	rightPrec_ = TSFCore::LinearOpHandle<Scalar>();
	compositeOperator_ = TSFCore::LinearOpHandle<Scalar>();
	rhs_ = null;
	lhs_ = null;
	statusTest_ = null;
	currFirstRhsOffset_ = 0;
	currNumRhs_ = 0;
	currBlockSize_ = 0;
	currRhsIndexes_.resize(0);
	currRhs_store_ = null;
	currRhs_ = null;
	currInitResidual_ = null;
	currInitResidual_store_ = null;
	currInitLhs_ = null;
	currInitLhs_store_ = null;
	currLhs_ = null;
	currLhs_store_ = null;
	currResidual_ = null;
	currResidual_store_ = null;
	// Note: Must release the ***_ view objects before we release the
	// ***_storage_ objects as according to the specification of
	// TSFCore::MultiVector
}

} // end Belos namespace

#endif // BELOS_LINEAR_PROBLEM_HPP
