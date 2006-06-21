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

#ifndef BELOS_NONBLOCK_GMRES_HPP
#define BELOS_NONBLOCK_GMRES_HPP

#include "Belos_BasicIteration.hpp"
#include "Belos_LinearProblemIteration.hpp"
#include "Belos_StatusTest.hpp"
#include "Belos_deflate.hpp"
#include "TSFCoreTestingTools.hpp"
#include "TSFCoreMultiVectorStdOps.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"
#include "Teuchos_arrayArg.hpp"
#include "StandardCompositionMacros.hpp"

namespace Belos {

///
/** Simple implementation of a non-block single RHS generized minimum
 * residual (GMRES) iterative linear solver for unsymmetric systems.
 *
 * ToDo: Finish documentation!
 */
template<class Scalar>
class NonblockGmres : public BasicIteration<Scalar> {
public:

	typedef Teuchos::ScalarTraits<Scalar>  ST;
	typedef typename ST::magnitudeType ScalarMag;

	/** @name Constructors and initializers */
	//@{

	///
	/** Determine the maximum Krylov subspace dimension before a restart
	 * is performed.
	 */
	STANDARD_MEMBER_COMPOSITION_MEMBERS(int,maxKrylovDim)

	///
	/** Set the tolerance for detecting solver breakdown
	 */
	STANDARD_MEMBER_COMPOSITION_MEMBERS(ScalarMag,breakdown_tol )

	///
	/** Set the output steam for iterative algorithm.
	 */
	STANDARD_COMPOSITION_MEMBERS(std::ostream,out)

	///
	/** Determine if we dump all info or not.
	 */
	STANDARD_MEMBER_COMPOSITION_MEMBERS(bool,dump_all)

	///
	/** Default initialization.
	 */
	NonblockGmres(
		const int                                                           maxKrylovDim  = 1000
		,const typename Teuchos::ScalarTraits<Scalar>::magnitudeType        breakdown_tol = 0.0
		,const Teuchos::RefCountPtr<std::ostream>                           &out          = Teuchos::null
		,const bool                                                         dump_all      = false
		);

	//@}

	/** @name Overridden from BasicIterationState */
	//@{

	///
	const LinearProblemState<Scalar>& getProblem() const;
	///
	int getCurrNumIters() const;
	///
	int getCurrNumRestarts() const;
	///
	void forceCurrLhsUpdate() const;
	///
	ENativeResidualType getCurrNativeResidualType() const;
	///
	void getCurrNativeResiduals(
		const int                              currBlockSize
		,Scalar                                norms[]
		,const TSFCore::MultiVector<Scalar>*   *residuals
		) const;

	//@}

	/** @name Overridden from BasicIteration */
	//@{

	/// Returns <tt>false</tt>
	bool adjointRequired() const;
	///
	void setProblem( const RefCountPtr<LinearProblemIteration<Scalar> > &lpi );
	///
	RefCountPtr<LinearProblemIteration<Scalar> > getProblem();
	///
	void initialize();
	///
	IterateReturn iterate( const int maxNumIter );
	///
	void finalize();

	//@}

private:

	// //////////////////////////////////
	// Private types

	typedef Teuchos::ScalarTraits<Scalar> ST;
	typedef typename ST::magnitudeType ScalarMag;

	// //////////////////////////////////
	// Private data members
	
	RefCountPtr<LinearProblemIteration<Scalar> >           lpi_;
	int                                                    currNumIters_;
	int                                                    currTotalNumIters_;
	int                                                    currNumRestarts_;
	Scalar                                                 r0_nrm_; 
	Scalar			                                           rel_r_nrm_;
	Teuchos::RefCountPtr< TSFCore::MultiVector<Scalar> >   V_;
	mutable std::valarray<Scalar>                          z_;
	mutable std::valarray<Scalar>                          z_tmp_;
	Teuchos::SerialDenseMatrix<int,Scalar>                 H_;
	std::valarray<Scalar>	                                 cs_, sn_;

	mutable Teuchos::BLAS<int,Scalar>                      blas_;

	// //////////////////////////////////
	// Private member functions

	///
	IterateReturn nextBlock( const int maxNumIter, bool *allFinished );
	/// Returns true of solver breakdown has occured
	bool computeIteration();
	///
	void cleanUp();

};

// //////////////////////////////
// Implementation

template<class Scalar>
NonblockGmres<Scalar>::NonblockGmres(
	const int                                                           maxKrylovDim
	,const typename Teuchos::ScalarTraits<Scalar>::magnitudeType        breakdown_tol
	,const Teuchos::RefCountPtr<std::ostream>                           &out
	,const bool                                                         dump_all
	)
	:maxKrylovDim_(maxKrylovDim)
	,breakdown_tol_(breakdown_tol)
	,out_(out)
	,dump_all_(dump_all)
{}

// Overridden from BasicIterationState

template<class Scalar>
const LinearProblemState<Scalar>& NonblockGmres<Scalar>::getProblem() const
{
	return *lpi_;
}

template<class Scalar>
int NonblockGmres<Scalar>::getCurrNumIters() const
{
	return currTotalNumIters_;
}

template<class Scalar>
int NonblockGmres<Scalar>::getCurrNumRestarts() const
{
	return currNumRestarts_;
}

template<class Scalar>
void NonblockGmres<Scalar>::forceCurrLhsUpdate() const
{
	const Scalar one = ST::one();
	// Solve least squares problem.
	if( !lpi_->isCurrLhsUpdated() && currNumIters_ > 0 ) {
		if(z_tmp_.size() < z_.size()) z_tmp_.resize(z_.size());
		std::copy( &z_[0], &z_[0] + currNumIters_, &z_tmp_[0] );
		blas_.TRSM(
			Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, Teuchos::NON_UNIT_DIAG
			,currNumIters_, 1, one, H_.values(), maxKrylovDim()+1, &z_tmp_[0], maxKrylovDim()+1
			);
		// Compute the new solution update
		RefCountPtr<const TSFCore::MultiVector<Scalar> >
			V = V_->subView(TSFCore::Range1D(1,currNumIters_));
		RefCountPtr<const TSFCore::Vector<Scalar> >
			z_vec	= V->domain()->createMemberView(RTOpPack::SubVectorT<Scalar>(0,currNumIters_,&z_tmp_[0],1));
		RefCountPtr<TSFCore::Vector<Scalar> >
			nativeLhsStep = V->range()->createMember();
		V->apply( TSFCore::NOTRANS, *z_vec, &*nativeLhsStep, -one );
		lpi_->updateCurrLhs( *nativeLhsStep );
	}
	else {
		lpi_->setCurrLhsUpdated(true);
	}
}

template<class Scalar>
ENativeResidualType NonblockGmres<Scalar>::getCurrNativeResidualType() const
{
	return NATIVE_RESIDUAL_PRECONDITIONED;
}

template<class Scalar>
void NonblockGmres<Scalar>::getCurrNativeResiduals(
	const int                              currBlockSize
	,Scalar                                norms[]
	,const TSFCore::MultiVector<Scalar>*   *residuals
	) const
{
	if(norms) {
#ifdef TEUCHOS_DEBUG
		TEST_FOR_EXCEPT( currBlockSize != 1 );
#endif
		norms[0] = rel_r_nrm_;
	}
	if(residuals) {
		*residuals = NULL;
	}
}

// Overridden from BasicIteration

template<class Scalar>
bool NonblockGmres<Scalar>::adjointRequired() const
{
	return false;
}

template<class Scalar>
void NonblockGmres<Scalar>::setProblem( const RefCountPtr<LinearProblemIteration<Scalar> > &lpi )
{
#ifdef TEUCHOS_DEBUG
	TEST_FOR_EXCEPT(lpi.get()==NULL);
#endif
	lpi_ = lpi;
}

template<class Scalar>
RefCountPtr<LinearProblemIteration<Scalar> > NonblockGmres<Scalar>::getProblem()
{
	return lpi_;
}

template<class Scalar>
void NonblockGmres<Scalar>::initialize()
{
	// I don't know what to do with this?
}

template<class Scalar>
IterateReturn NonblockGmres<Scalar>::iterate( const int maxNumIter )
{
	int numCumulativeIter = 0;
	bool allConverged = true;
	// Get the composite operator
	TSFCore::LinearOpHandle<Scalar> Op = lpi_->getCompositeOperator();
	if(get_out().get()) {
		*get_out() << "\n*** Entering NonblockGmres<Scalar>::iterate(...)\n" << std::setprecision(16);
		if(dump_all()) {
			*get_out() << "\nOp.op =\n" << *Op.op();
			*get_out() << "\nOp.defaultTrans = " << toString(Op.defaultTrans()) << std::endl;
			*get_out() << "\nrhs =\n" << lpi_->getRhs();
			*get_out() << "\nlhs =\n" << lpi_->getLhs();
		}
	}
	// Initialize workspace/data for GMRES
	const int maxKrylovDim = this->maxKrylovDim();
	V_ = Op.domain()->createMembers(maxKrylovDim+1);
	if( maxKrylovDim+1 != H_.numRows() ) {
		H_.shapeUninitialized(maxKrylovDim+1,maxKrylovDim);
		z_.resize(maxKrylovDim+1);
		cs_.resize(maxKrylovDim);
		sn_.resize(maxKrylovDim);
	}
	// Solve the systems as a set of single RHS solves
	lpi_->setBlockSize(1); // Only unit block sizes allowed!
	bool allFinished = false;
	while ( !allFinished  ) {
		IterateReturn nextBlockReturn = nextBlock( maxNumIter, &allFinished );
		numCumulativeIter += nextBlockReturn.numCumulativeIter;
		if( nextBlockReturn.iterateTermination != TERMINATION_STATUS_TEST ) allConverged = false;
	}
	cleanUp();
	return IterateReturn( allConverged ? TERMINATION_STATUS_TEST : TERMINATION_MAX_NUM_ITER, numCumulativeIter );
}

template<class Scalar>
void NonblockGmres<Scalar>::finalize()
{
	// I don't know what to do with this?
}

// private

template<class Scalar>
IterateReturn NonblockGmres<Scalar>::nextBlock( const int maxNumIter, bool *allFinished )
{
	*allFinished = !lpi_->setupCurrSystem();
	if(*allFinished) return IterateReturn();
	// Inform the status test that we are solving a new block of systems
	StatusTest<Scalar> *statusTest = lpi_->getStatusTest();
	if(statusTest) statusTest->reset();
	// Perform iterations on this block of RHSs
	bool breakdown = false;
	EIterateTermination doIterationReturn = TERMINATION_UNDEFINED;
	currNumRestarts_ = 0;
	currNumIters_ = 0;
	try {
		for( currTotalNumIters_ = 0 ; ; ++currNumIters_, ++currTotalNumIters_ ) {
			// Setup for first GMRES iteration for current restart if needed
			if( currNumIters_ == 0 || currNumIters_ == maxKrylovDim() - 1 ) {
				// Setup for the first iteration of GMRES for the current restart
				if( currNumIters_ == maxKrylovDim() - 1 ) {
					// Must do a restart!
					lpi_->restart(*this);
					currNumIters_ = 0;
					++currNumRestarts_;
				}
				Teuchos::RefCountPtr<TSFCore::Vector<Scalar> > v_1 = V_->col(1);	// Get a mutable view of the first column of V_
				RefCountPtr<TSFCore::Vector<Scalar> > r_curr_0 = v_1;             // Just a name change!
				lpi_->computeCurrCompositeResidual(NULL,NULL,&*r_curr_0);
				const ScalarMag r_curr_0_nrm = norm(*r_curr_0);
				if( currTotalNumIters_ == 0 ) r0_nrm_ = r_curr_0_nrm;             // Remember initial residual for convergence test!
				rel_r_nrm_ = r_curr_0_nrm / r0_nrm_;
				z_[0] = r_curr_0_nrm;
				// Note! v_1 *is* r_curr_0
				Vt_S( &*v_1, ST::one()/r_curr_0_nrm );                            // v_1 = r_curr_0 / ||r_curr_0||
			}
			// Ask the status test to check the status of the current RHS
			StatusTest<Scalar> *statusTest = lpi_->getStatusTest();
			if(statusTest) {
				const int currNumRhs = 1;
				std::vector<EStatusType> status(currNumRhs,STATUS_UNCHECKED);
				statusTest->checkStatus(*this,currNumRhs,currNumRhs,&status[0]);
				TEST_FOR_EXCEPTION( status[0]==STATUS_FAILED || status[0] == STATUS_NAN, Exceptions::SolverBreakdown, "Error!" );
				if(status[0]==STATUS_CONVERGED) {
					doIterationReturn = TERMINATION_STATUS_TEST;
					break;
				}
			}
			// If breakdown occurred but status test did not recognize convergence
			// then throw SolverBreakdown excpetion!
			TEST_FOR_EXCEPTION(
				breakdown, Exceptions::SolverBreakdown
				,"Belos::NonblockGmres<Scalar>::nextBlock(...): The solver breakdown tolerance"
				" has been reached indicating the system(s) is solved but status test did not indicate convergence!"
				);
			// Check if max iterations are exceeded
			if( currTotalNumIters_ == maxNumIter ) {
				doIterationReturn = TERMINATION_MAX_NUM_ITER;
				break;
			}
			// Perform linear algebra for current iteration
			breakdown = computeIteration();
			// Inform that the current solution is *not* up to date!
			lpi_->setCurrLhsUpdated(false);
		}
	}
	catch(...) {
		// Copy whatever is remaining in current LHS into full LHS
		lpi_->finalizeCurrSystem(*this);
		throw;
	}
	// Copy whatever is remaining in current LHS into full LHS
	lpi_->finalizeCurrSystem(*this);
	// Return status of this block
	return IterateReturn(doIterationReturn,currTotalNumIters_);
}

#define BELOS_NONBLOCK_CG_SOLVER_ERR_MSG "NonblockGmres<Scalar>::computeIteration(...): iteration = " << this->getCurrNumIters() << ": Error, "

template<class Scalar>
bool NonblockGmres<Scalar>::computeIteration()
{
	const Scalar one = ST::one();
	int i;
	TSFCore::LinearOpHandle<Scalar> Op = lpi_->getCompositeOperator();
	Teuchos::SerialDenseMatrix<int, Scalar> &H = H_;
	std::valarray<Scalar> &z = z_, &cs = cs_, &sn = sn_;
	const int jm1 = currNumIters_;
	// 
	Teuchos::RefCountPtr<TSFCore::Vector<Scalar> >
		w = V_->col(jm1+2);                                              // w = v_{j+1}
	Op.apply( TSFCore::NOTRANS, *V_->col(jm1+1), &*w );                // w = M * v_{j}	
	//
	// Perform MGS to orthogonalize new Krylov vector.
	//
	for( i = 0; i < jm1 + 1; i++ ) {
		H( i, jm1 ) = w->space()->scalarProd( *w, *V_->col(i+1) );       // h_{i,j} = ( w, v_{i} )
		Vp_StV( &*w, -H( i, jm1 ), *V_->col(i+1) );                      // w = w - h_{i,j} * v_{i}
	}
	const Scalar H_jp1_j = H( jm1+1, jm1 ) =
		norm( *w );                                                      // h_{j+1,j} = || w ||
	const bool breakdown = ST::magnitude(H_jp1_j) <= breakdown_tol();
	if(!breakdown) Vt_S( &*w, one / H_jp1_j );                         // v_{j+1} = w / h_{j+1,j}			
	//
	// Apply previous Givens rotations
	//
	for( i = 0; i < jm1; i++ ) {
		const Scalar temp = cs[i]*H( i, jm1 ) + sn[i]*H( i+1, jm1 );
		H( i+1, jm1 ) = -sn[i]*H( i, jm1 ) + cs[i]*H( i+1, jm1 );
		H( i, jm1 ) = temp;
	}
	//
	// Calculate new Givens rotation
	//
	blas_.ROTG(
		&H( jm1, jm1 ), &H( jm1+1, jm1 ), 
		&cs[jm1], &sn[jm1]
		);
	//
	// Update RHS and residual w/ new transform and compute residual norm.
	//
	z[jm1+1] = -sn[jm1]*z[jm1];
	z[jm1] *= cs[jm1];
	rel_r_nrm_ = ST::magnitude( z[jm1+1] ) / r0_nrm_;
	//
	// Return that the solver has not broken down
	//
	return breakdown;
}

#undef BELOS_NONBLOCK_CG_SOLVER_ERR_MSG

template<class Scalar>
void NonblockGmres<Scalar>::cleanUp()
{
	// ToDo: Resize arrays to zero?
}

} // end Belos namespace

#endif // BELOS_NONBLOCK_GMRES_HPP
