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

//#define BELOS_LINEAR_PROBLEM_DUMP_COUT

#include "Belos_LinearProblemSetup.hpp"
#include "Belos_deflate.hpp"
#include "TSFCoreLinOp.hpp"
#include "TSFCoreMultiVectorStdOps.hpp"
#include "TSFCoreMultiplicativeLinearOp.hpp"
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
	int getTotalNumRhs() const;
	///
	TSFCore::LinOpNonPersisting<Scalar> getOperator() const;
	///
	EOpSymmetry getOperatorSymmetry() const;
	///
	TSFCore::LinOpNonPersisting<Scalar> getRightPrec() const;
	///
	EOpSymmetry getRightPrecSymmetry() const;
	///
	TSFCore::LinOpNonPersisting<Scalar> getLeftPrec() const;
	///
	EOpSymmetry getLeftPrecSymmetry() const;
	///
	TSFCore::LinOpNonPersisting<Scalar> getCombinedOperator() const;
	///
	const TSFCore::Vector<Scalar>* getRightScaling() const;
	///
	const TSFCore::Vector<Scalar>* getLeftScaling() const;
	///
	const TSFCore::MultiVector<Scalar>& getRhs() const;
	///
	TSFCore::MultiVector<Scalar>& getLhs() const;
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
	int getBlockSize() const;
	///
	StatusTest<Scalar>* getStatusTest();
	///
	void setCurrSystem(
		const int                   firstRhsOffset
		,const int                  numRhs
		,const int                  blockSize
		);
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
	void computeCurrPrecResidual(
		const TSFCore::MultiVector<Scalar>    *currLhs
		,const TSFCore::MultiVector<Scalar>   *currRhs
		,TSFCore::MultiVector<Scalar>         *currPrecResidual
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
	void setCurrToFullLhs( const BasicIterationState<Scalar> &bis );

	//@}

	/** @name Overridden from LinearProblemSetup */
	//@{

	///
	void initialize(
		const TSFCore::LinOpPersisting<Scalar>                   &op
		,const EOpSymmetry                                       symmetry
		,const RefCountPtr<const TSFCore::MultiVector<Scalar> >  &rhs
		,const RefCountPtr<TSFCore::MultiVector<Scalar> >        &lhs
		);
	///
	void setOperator( const TSFCore::LinOpPersisting<Scalar> &op, const EOpSymmetry symmetry );
	///
	TSFCore::LinOpPersisting<Scalar> getOperator();
	///
	void setRightPrec( const TSFCore::LinOpPersisting<Scalar> &rightPrec, const EOpSymmetry symmetry );
	///
	TSFCore::LinOpPersisting<Scalar> getRightPrec();
	///
	void setLeftPrec( const TSFCore::LinOpPersisting<Scalar> &leftPrec, const EOpSymmetry symmetry );
	///
	TSFCore::LinOpPersisting<Scalar> getLeftPrec();
	///
	void setRhs( const RefCountPtr<const TSFCore::MultiVector<Scalar> > &rhs );
	///
	RefCountPtr<const TSFCore::MultiVector<Scalar> > getRhs();
	///
	void setLhs( const RefCountPtr<TSFCore::MultiVector<Scalar> > &lhs );
	///
	RefCountPtr<TSFCore::MultiVector<Scalar> > getLhs();
	///
	void setBlockSize( const int blockSize );
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
	int                                               blockSize_;
	int                                               totalNumRhs_;
	TSFCore::LinOpPersisting<Scalar>                  operator_;
	EOpSymmetry                                       operatorSymmetry_;
	TSFCore::LinOpPersisting<Scalar>                  rightPrec_;
	EOpSymmetry                                       rightPrecSymmetry_;
	TSFCore::LinOpPersisting<Scalar>                  leftPrec_;
	EOpSymmetry                                       leftPrecSymmetry_;
	RefCountPtr<const TSFCore::Vector<Scalar> >       rightScaling_;
	RefCountPtr<const TSFCore::Vector<Scalar> >       leftScaling_;
	TSFCore::LinOpPersisting<Scalar>                  combinedOperator_;
	RefCountPtr<const TSFCore::MultiVector<Scalar> >  rhs_;
	RefCountPtr<TSFCore::MultiVector<Scalar> >        lhs_;
	RefCountPtr<StatusTest<Scalar> >                  statusTest_;

	// Members for current block	
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

};

// /////////////////////////////
// Implementation

template<class Scalar>
LinearProblem<Scalar>::LinearProblem()
	:blockSize_(1)
	,totalNumRhs_(0)
	,currNumRhs_(0)
{}

// Overridden from LinearProblemState

template<class Scalar>
int LinearProblem<Scalar>::getTotalNumRhs() const
{
	return totalNumRhs_;
}

template<class Scalar>
TSFCore::LinOpNonPersisting<Scalar> LinearProblem<Scalar>::getOperator() const
{
	return operator_;
}

template<class Scalar>
EOpSymmetry LinearProblem<Scalar>::getOperatorSymmetry() const
{
	return operatorSymmetry_;
}

template<class Scalar>
TSFCore::LinOpNonPersisting<Scalar> LinearProblem<Scalar>::getRightPrec() const
{
	return rightPrec_;
}

template<class Scalar>
EOpSymmetry LinearProblem<Scalar>::getRightPrecSymmetry() const
{
	return rightPrecSymmetry_;
}

template<class Scalar>
TSFCore::LinOpNonPersisting<Scalar> LinearProblem<Scalar>::getLeftPrec() const
{
	return leftPrec_;
}

template<class Scalar>
EOpSymmetry LinearProblem<Scalar>::getLeftPrecSymmetry() const
{
	return leftPrecSymmetry_;
}

template<class Scalar>
TSFCore::LinOpNonPersisting<Scalar> LinearProblem<Scalar>::getCombinedOperator() const
{
	return combinedOperator_;
}

template<class Scalar>
const TSFCore::Vector<Scalar>* LinearProblem<Scalar>::getRightScaling() const
{
	return rightScaling_.get();
}

template<class Scalar>
const TSFCore::Vector<Scalar>* LinearProblem<Scalar>::getLeftScaling() const
{
	return leftScaling_.get();
}

template<class Scalar>
const TSFCore::MultiVector<Scalar>& LinearProblem<Scalar>::getRhs() const
{
	return *rhs_;
}

template<class Scalar>
TSFCore::MultiVector<Scalar>& LinearProblem<Scalar>::getLhs() const
{
	return *lhs_;
}

template<class Scalar>
int LinearProblem<Scalar>::getCurrNumRhs() const
{
	return currNumRhs_;
}

template<class Scalar>
int LinearProblem<Scalar>::getCurrBlockSize() const
{
	return currBlockSize_;
}

template<class Scalar>
void LinearProblem<Scalar>::getCurrRhsIndexes(
	const int              currNumRhs
	,int                   currRhsIndexes[]
	) const
{
#ifdef _DEBUG
	TEST_FOR_EXCEPT( currNumRhs != currNumRhs_ );
#endif
	std::copy( &currRhsIndexes_[0], &currRhsIndexes_[0]+currNumRhs, currRhsIndexes );
}

template<class Scalar>
const TSFCore::MultiVector<Scalar>& LinearProblem<Scalar>::getCurrRhs() const
{
	return *currRhs_;
}

template<class Scalar>
const TSFCore::MultiVector<Scalar>& LinearProblem<Scalar>::getCurrInitResidual() const
{
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
	return isCurrInitResidualComputed_;
}

template<class Scalar>
const TSFCore::MultiVector<Scalar>& LinearProblem<Scalar>::getCurrInitLhs() const
{
	return *currInitLhs_;
}

template<class Scalar>
const TSFCore::MultiVector<Scalar>& LinearProblem<Scalar>::getCurrLhs() const
{
	return *currLhs_;
}

template<class Scalar>
bool LinearProblem<Scalar>::isCurrLhsUpdated() const
{
	return isCurrLhsUpdated_;
}

template<class Scalar>
const TSFCore::MultiVector<Scalar>& LinearProblem<Scalar>::getCurrResidual() const
{
	TEST_FOR_EXCEPT( leftScaling_.get() ); // Can't handle scaling yet!
#ifdef _DEBUG
	TEST_FOR_EXCEPT( isCurrLhsUpdated_==false );
#endif
	if(!isCurrResidualComputed_) {
		// currResidual = operator * curLhs - currRhs
		assign( &*currResidual_, *currRhs_ );
		const ScalarMag currLhs_norm_1 = norm_1(*currLhs_);
		if( currLhs_norm_1 != ST::zero() ) {
			operator_.apply(
				TSFCore::NOTRANS
				,*currLhs_
				,&*currResidual_
				,+ST::one()
				,-ST::one()
				);
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
	return isCurrResidualComputed_;
}

// Overridden from LinearProblemIteration

template<class Scalar>
int LinearProblem<Scalar>::getBlockSize() const
{
	return blockSize_;
}

template<class Scalar>
StatusTest<Scalar>* LinearProblem<Scalar>::getStatusTest()
{
	return statusTest_.get();
}

template<class Scalar>
void LinearProblem<Scalar>::setCurrSystem(
	const int                   firstRhsOffset
	,const int                  numRhs
	,const int                  blockSize
	)
{
#ifdef _DEBUG
	TEST_FOR_EXCEPT( ! ( firstRhsOffset + numRhs <= totalNumRhs_ ) );
	TEST_FOR_EXCEPT( ! ( blockSize >= numRhs ) );
#endif
	RefCountPtr<const TSFCore::VectorSpace<Scalar> > rhs_range = rhs_->range();
	RefCountPtr<const TSFCore::VectorSpace<Scalar> > lhs_range = lhs_->range();
	// Initialize current RHS, LHS and initial LHS
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
	const TSFCore::Range1D currRhsRng(1,numRhs);
	const TSFCore::Range1D currOrigRhsRng = currRhsRng + firstRhsOffset;
	assign( &*currRhs_->subView(currRhsRng), *rhs_->subView(currOrigRhsRng) );
	assign( &*currLhs_->subView(currRhsRng), *lhs_->subView(currOrigRhsRng) );
	isCurrLhsUpdated_ = true;
	currLhsIsCurrInitLhs_ = true;
	if( blockSize > numRhs ) {
		const TSFCore::Range1D currAugRng(numRhs+1,blockSize);
		randomize( -ST::one(), +ST::one(), &*currRhs_->subView(currAugRng) ); // Random RHS for augmented systems
		assign( &*currLhs_->subView(currAugRng), ST::zero() );                // Zero initial guess for augmented systems
	}
	assign( &*currInitLhs_, *currLhs_ );
	// Initialize list of active RHSs (one-based)
	currRhsIndexes_.resize(numRhs);
	for( int k = 0; k < numRhs; ++k ) currRhsIndexes_[k] = firstRhsOffset+k+1;
	// Mark derived quantities as uninitialized (unscaled unpreconditioned residuals)
	isCurrInitResidualComputed_ = false;
	isCurrResidualComputed_ = false;
	// Finish up initialization
	currNumRhs_ = numRhs;
	currBlockSize_ = numRhs;
}

template<class Scalar>
RefCountPtr<TSFCore::MultiVector<Scalar> >
LinearProblem<Scalar>::getCurrLhs()
{
	currLhsIsCurrInitLhs_ = false; // Just to be safe!
	return currLhs_;
}

template<class Scalar>
void LinearProblem<Scalar>::setCurrLhsUpdated( const bool isCurrLhsUpdated )
{
	isCurrLhsUpdated_ = isCurrLhsUpdated;
	isCurrResidualComputed_ = false;
	if(isCurrLhsUpdated) currLhsIsCurrInitLhs_ = false;
}

template<class Scalar>
void LinearProblem<Scalar>::updateCurrLhs( const TSFCore::MultiVector<Scalar> &nativeLhsStep )
{
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
	TEST_FOR_EXCEPT( rightScaling_.get() ); // Can't handle scaling yet!
	if(rightPrec_.op().get())
		rightPrec_.apply(TSFCore::NOTRANS,nativeLhs,&*currLhs_,ST::one());
	else
		assign( &*currLhs_, nativeLhs );
	isCurrLhsUpdated_ = true;
	isCurrResidualComputed_ = false;
	currLhsIsCurrInitLhs_ = false;
}

template<class Scalar>
RefCountPtr<TSFCore::MultiVector<Scalar> >
LinearProblem<Scalar>::getCurrResidual()
{
	return currResidual_;
}

template<class Scalar>
void LinearProblem<Scalar>::setCurrResidualComputed( const bool isCurrResidualComputed )
{
	isCurrResidualComputed_ = isCurrResidualComputed;
}

template<class Scalar>
void LinearProblem<Scalar>::computeCurrPrecResidual(
	const TSFCore::MultiVector<Scalar>    *currLhs
	,const TSFCore::MultiVector<Scalar>   *currRhs
	,TSFCore::MultiVector<Scalar>         *currPrecResidual
	)
{
	if( currLhs==NULL && currRhs==NULL ) {
#ifdef _DEBUG
		TEST_FOR_EXCEPT( !isCurrLhsUpdated_ );
#endif
		const TSFCore::MultiVector<Scalar> &currResidual = const_cast<const LinearProblem<Scalar>*>(this)->getCurrResidual();
		if( leftPrec_.op().get() )
			leftPrec_.apply( TSFCore::NOTRANS, currResidual, currPrecResidual );
		else
			assign( currPrecResidual, currResidual );
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
	if(!isCurrLhsUpdated_) bis.forceCurrLhsUpdate();
	assign( &*currInitLhs_, *currLhs_ );
}

template<class Scalar>
void LinearProblem<Scalar>::setCurrToFullLhs( const BasicIterationState<Scalar> &bis )
{
	if(!isCurrLhsUpdated_) bis.forceCurrLhsUpdate();
	assign( &*lhs_->subView(currNumRhs_,&currRhsIndexes_[0]), *currLhs_->subView(TSFCore::Range1D(1,currNumRhs_)) );
	currNumRhs_ = 0; // Set as uninitialized
}

// Overridden from LinearProblemSetup

template<class Scalar>
void LinearProblem<Scalar>::initialize(
	const TSFCore::LinOpPersisting<Scalar>                   &op
	,const EOpSymmetry                                       symmetry
	,const RefCountPtr<const TSFCore::MultiVector<Scalar> >  &rhs
	,const RefCountPtr<TSFCore::MultiVector<Scalar> >        &lhs
	)
{
#ifdef _DEBUG
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
void LinearProblem<Scalar>::setOperator( const TSFCore::LinOpPersisting<Scalar> &op, const EOpSymmetry symmetry )
{
#ifdef _DEBUG
	TEST_FOR_EXCEPT( op.op().get() == NULL );
#endif
	operator_           = op;
	operatorSymmetry_   = symmetry;
	totalNumRhs_        = 0;
	currNumRhs_         = 0;
}

template<class Scalar>
TSFCore::LinOpPersisting<Scalar> LinearProblem<Scalar>::getOperator()
{
	return operator_;
}

template<class Scalar>
void LinearProblem<Scalar>::setRightPrec( const TSFCore::LinOpPersisting<Scalar> &rightPrec, const EOpSymmetry symmetry )
{
#ifdef _DEBUG
	TEST_FOR_EXCEPT( rightPrec.op().get() == NULL );
#endif
	rightPrec_           = rightPrec;
	rightPrecSymmetry_   = symmetry;
	totalNumRhs_         = 0;
	currNumRhs_          = 0;
}

template<class Scalar>
TSFCore::LinOpPersisting<Scalar> LinearProblem<Scalar>::getRightPrec()
{
	return rightPrec_;
}

template<class Scalar>
void LinearProblem<Scalar>::setLeftPrec( const TSFCore::LinOpPersisting<Scalar> &leftPrec, const EOpSymmetry symmetry )
{
#ifdef _DEBUG
	TEST_FOR_EXCEPT( leftPrec.op().get() == NULL );
#endif
	leftPrec_           = leftPrec;
	leftPrecSymmetry_   = symmetry;
	totalNumRhs_        = 0;
	currNumRhs_         = 0;
}

template<class Scalar>
TSFCore::LinOpPersisting<Scalar> LinearProblem<Scalar>::getLeftPrec()
{
	return leftPrec_;
}

template<class Scalar>
void LinearProblem<Scalar>::setRhs( const RefCountPtr<const TSFCore::MultiVector<Scalar> > &rhs )
{
#ifdef _DEBUG
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
#ifdef _DEBUG
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
void LinearProblem<Scalar>::setBlockSize( const int blockSize )
{
#ifdef _DEBUG
	TEST_FOR_EXCEPT( blockSize <= 0 );
#endif
	blockSize_ = blockSize;
	totalNumRhs_ = 0;
	currNumRhs_ = 0;
}

template<class Scalar>
void LinearProblem<Scalar>::setStatusTest( const RefCountPtr<StatusTest<Scalar> > &statusTest )
{
#ifdef _DEBUG
	TEST_FOR_EXCEPT( statusTest.get() == NULL );
#endif
	statusTest_ = statusTest;
	totalNumRhs_ = 0;
	currNumRhs_ = 0;
}

template<class Scalar>
void LinearProblem<Scalar>::completeSetup()
{
	// ToDo: Validate that a consistent linear problem has been specified!
	// Build the combined operator.
	if( leftPrec_.op().get() ) {
		if( rightPrec_.op().get() ) {
			combinedOperator_ = TSFCore::LinOpPersisting<Scalar>(
				rcp(
					new TSFCore::MultiplicativeLinearOp<Scalar>(
						3, Teuchos::arrayArg<TSFCore::LinOpPersisting<Scalar> >(leftPrec_,operator_,rightPrec_)()
						)
					)
				);
		}
		else {
			combinedOperator_ = TSFCore::LinOpPersisting<Scalar>(
				rcp(
					new TSFCore::MultiplicativeLinearOp<Scalar>(
						2, Teuchos::arrayArg<TSFCore::LinOpPersisting<Scalar> >(leftPrec_,operator_)()
						)
					)
				);
		}
	}
	else {
		if( rightPrec_.op().get() ) {
			combinedOperator_ = TSFCore::LinOpPersisting<Scalar>(
				rcp(
					new TSFCore::MultiplicativeLinearOp<Scalar>(
						2, Teuchos::arrayArg<TSFCore::LinOpPersisting<Scalar> >(operator_,rightPrec_)()
						)
					)
				);
		}
		else {
			combinedOperator_ = operator_;
		}
	}
	// Finish initialization
	totalNumRhs_ = rhs_->domain()->dim();
	currNumRhs_ = 0;
}

template<class Scalar>
void LinearProblem<Scalar>::uninitialize()
{
	totalNumRhs_ = 0;
	operator_ = TSFCore::LinOpPersisting<Scalar>();
	leftPrec_ = TSFCore::LinOpPersisting<Scalar>();
	rightPrec_ = TSFCore::LinOpPersisting<Scalar>();
	combinedOperator_ = TSFCore::LinOpPersisting<Scalar>();
	rhs_ = null;
	lhs_ = null;
	statusTest_ = null;
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
