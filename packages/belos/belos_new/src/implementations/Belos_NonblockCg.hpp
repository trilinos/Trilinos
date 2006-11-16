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

#ifndef BELOS_NONBLOCK_CG_HPP
#define BELOS_NONBLOCK_CG_HPP

#include "Belos_BasicIteration.hpp"
#include "Belos_LinearProblemIteration.hpp"
#include "Belos_StatusTest.hpp"
#include "Belos_deflate.hpp"
#include "TSFCoreTestingTools.hpp"
#include "TSFCoreMultiVectorStdOps.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"
#include "Teuchos_arrayArg.hpp"
#include "Teuchos_getConst.hpp"
#include "StandardCompositionMacros.hpp"

namespace Belos {

///
/** Simple implementation of a non-block (but mulit-RHS) conjugate
 * gradient (CG) iterative linear solver for symmetric positive
 * definite systems.
 *
 * ToDo: Finish documentation!
 */
template<class Scalar>
class NonblockCg : public BasicIteration<Scalar> {
public:

	/** @name Constructors and initializers */
	//@{

	///
	/** Set the output steam for iterative algorithm.
	 */
	STANDARD_COMPOSITION_MEMBERS(std::ostream,out);

	///
	/** Determine if we dump all info or not.
	 */
	STANDARD_MEMBER_COMPOSITION_MEMBERS(bool,dump_all);

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
	enum EOpPrec { OP_ONLY, COMPOSITE_OP, LEFT_PREC, RIGHT_PREC };

	// //////////////////////////////////
	// Private data members
	
	RefCountPtr<LinearProblemIteration<Scalar> >  lpi_;
	int                                           currNumIters_;
	EOpPrec                                       opPrecType_;
	RefCountPtr<TSFCore::MultiVector<Scalar> >    X_; // Only used if we need it
	RefCountPtr<TSFCore::MultiVector<Scalar> >    R_; // Only used if we need it
	RefCountPtr<TSFCore::MultiVector<Scalar> >    Q_;
	RefCountPtr<TSFCore::MultiVector<Scalar> >    Z_;
	RefCountPtr<TSFCore::MultiVector<Scalar> >    P_;
	std::vector<Scalar>                           rho_;
	std::vector<Scalar>                           rho_old_;
	std::vector<Scalar>                           beta_;
	std::vector<Scalar>                           gamma_;
	std::vector<Scalar>                           alpha_;
	std::vector<ScalarMag>                        R_native_init_norms_;
	mutable std::vector<ScalarMag>                R_native_norms_;
	mutable bool                                  R_native_norms_updated_;

	// //////////////////////////////////
	// Private member functions

	///
	IterateReturn nextBlock( const int maxNumIter, bool *allFinished );
	///
	bool setupCurrSystem();
	///
	EIterateTermination doIteration();
	///
	void computeIteration();
	///
	TSFCore::LinearOpHandle<Scalar> getOperator() const;
	///
	TSFCore::LinearOpHandle<Scalar> getPrec() const;
	///
	RefCountPtr<TSFCore::MultiVector<Scalar> > getCurrLhs() const;
	///
	RefCountPtr<TSFCore::MultiVector<Scalar> > getCurrResidual() const;
	///
	void cleanUp();

};

// //////////////////////////////
// Implementation

// Overridden from BasicIterationState

template<class Scalar>
const LinearProblemState<Scalar>& NonblockCg<Scalar>::getProblem() const
{
	return *lpi_;
}

template<class Scalar>
int NonblockCg<Scalar>::getCurrNumIters() const
{
	return currNumIters_;
}

template<class Scalar>
int NonblockCg<Scalar>::getCurrNumRestarts() const
{
	return 0; // We never restart a CG method!
}

template<class Scalar>
void NonblockCg<Scalar>::forceCurrLhsUpdate() const
{
	if(opPrecType_==COMPOSITE_OP) {
		lpi_->setCurrLhs(*this->getCurrLhs());
	}
}

template<class Scalar>
ENativeResidualType NonblockCg<Scalar>::getCurrNativeResidualType() const
{
	switch(opPrecType_) {
		case OP_ONLY: return NATIVE_RESIDUAL_UNPRECONDITIONED;
		case COMPOSITE_OP: return  NATIVE_RESIDUAL_PRECONDITIONED;
		case RIGHT_PREC: case LEFT_PREC: return NATIVE_RESIDUAL_UNPRECONDITIONED;
		default:
			TEST_FOR_EXCEPT(true);
	}
	return NATIVE_RESIDUAL_UNPRECONDITIONED; // Should never be executed!
}

template<class Scalar>
void NonblockCg<Scalar>::getCurrNativeResiduals(
	const int                              currBlockSize
	,Scalar                                norms[]
	,const TSFCore::MultiVector<Scalar>*   *residuals
	) const
{
	if(norms) {
#ifdef TEUCHOS_DEBUG
		TEST_FOR_EXCEPT( currBlockSize != lpi_->getCurrBlockSize() );
#endif
		if(!R_native_norms_updated_) {
			TSFCore::norms( *this->getCurrResidual(), &R_native_norms_[0] );
			R_native_norms_updated_ = true;
			if(get_out().get()) {
				*get_out() << "\n||R|| =\n"; for(int j=0;j<currBlockSize;++j) *get_out() << " " << R_native_norms_[j]; *get_out() << std::endl;
				*get_out() << "\n||R0|| =\n"; for(int j=0;j<currBlockSize;++j) *get_out() << " " << R_native_init_norms_[j]; *get_out() << std::endl;
			}
		}
		for(int j=0;j<currBlockSize;++j) norms[j] = R_native_norms_[j] / R_native_init_norms_[j];
	}
	if(residuals) {
		*residuals = &*this->getCurrResidual();
	}
}

// Overridden from BasicIteration

template<class Scalar>
bool NonblockCg<Scalar>::adjointRequired() const
{
	return false;
}

template<class Scalar>
void NonblockCg<Scalar>::setProblem( const RefCountPtr<LinearProblemIteration<Scalar> > &lpi )
{
#ifdef TEUCHOS_DEBUG
	TEST_FOR_EXCEPT(lpi.get()==NULL);
#endif
	lpi_ = lpi;
}

template<class Scalar>
RefCountPtr<LinearProblemIteration<Scalar> > NonblockCg<Scalar>::getProblem()
{
	return lpi_;
}

template<class Scalar>
void NonblockCg<Scalar>::initialize()
{
	// I don't know what to do with this?
}

template<class Scalar>
IterateReturn NonblockCg<Scalar>::iterate( const int maxNumIter )
{
	int numCumulativeIter = 0;
	bool allConverged = true;
	try{
		// Determine the types of the operators
		const TSFCore::LinearOpHandle<Scalar>
			rightPrec = lpi_->getRightPrec(),
			leftPrec  = lpi_->getLeftPrec();
		if(rightPrec.op().get()) {
			if(leftPrec.op().get()) {
				// Here we can just hope that the overall operator is
				// symmetric since the P_L must equal P_R' but P_L and P_R
				// need not be symmetric
				opPrecType_ = COMPOSITE_OP;
			}
			else {
#ifdef TEUCHOS_DEBUG
				TEST_FOR_EXCEPT( lpi_->getRightPrecSymmetry() != OP_SYMMETRIC );
#endif
				opPrecType_ = RIGHT_PREC;
			}
		}
		else {
			if(leftPrec.op().get()) {
#ifdef TEUCHOS_DEBUG
				TEST_FOR_EXCEPT( lpi_->getLeftPrecSymmetry() != OP_SYMMETRIC );
#endif
				opPrecType_ = LEFT_PREC;
			}
			else {
				opPrecType_ = OP_ONLY;
			}
		}
		// Print linear system
		TSFCore::LinearOpHandle<Scalar> Op = lpi_->getOperator();
		if(get_out().get()) {
			*get_out() << "\n*** Entering NonblockCg<Scalar>::iterate(...)\n" << std::setprecision(16);
			if(dump_all()) {
				*get_out() << "\nOp.op =\n" << *Op.op();
				*get_out() << "\nOp.defaultTrans = " << toString(Op.defaultTrans()) << std::endl;
				*get_out() << "\nrhs =\n" << lpi_->getRhs();
				*get_out() << "\nlhs =\n" << lpi_->getLhs();
			}
		}
		// Solve the systems as a set of block solves
		lpi_->setAugmentationAllowed(false);
		bool allFinished = false; 
		while ( !allFinished ) {
			IterateReturn nextBlockReturn = nextBlock( maxNumIter, &allFinished );
			if(allFinished) break;
			numCumulativeIter += nextBlockReturn.numCumulativeIter;
			if( nextBlockReturn.iterateTermination != TERMINATION_STATUS_TEST ) allConverged = false;
		}
		cleanUp();
	}
	catch(...) {
		cleanUp();
		throw;
	}
	return IterateReturn( allConverged ? TERMINATION_STATUS_TEST : TERMINATION_MAX_NUM_ITER, numCumulativeIter );
}

template<class Scalar>
void NonblockCg<Scalar>::finalize()
{
	// I don't know what to do with this?
}

// private

template<class Scalar>
IterateReturn NonblockCg<Scalar>::nextBlock( const int maxNumIter, bool *allFinished )
{
	// Setup for next block of RHSs to solve
	*allFinished = !this->setupCurrSystem(); // We are all finished if there is no current block system
	if(*allFinished) return IterateReturn();
	// Inform the status test that we are solving a new block of systems
	StatusTest<Scalar> *statusTest = lpi_->getStatusTest();
	if(statusTest) statusTest->reset();
	// Perform iterations on this block of RHSs
	EIterateTermination doIterationReturn = TERMINATION_UNDEFINED;
	for( currNumIters_ = 0; currNumIters_ < maxNumIter; ++currNumIters_ ) {
		doIterationReturn = doIteration();
		if(doIterationReturn==TERMINATION_STATUS_TEST)
			break;
	}
	if( currNumIters_==maxNumIter && doIterationReturn!=TERMINATION_STATUS_TEST )
		doIterationReturn = TERMINATION_MAX_NUM_ITER;
	// Copy whatever is remaining in current LHS into full LHS
	lpi_->finalizeCurrSystem(*this);
	// Return status of this block
	return IterateReturn(doIterationReturn,currNumIters_);
}

template<class Scalar>
bool NonblockCg<Scalar>::setupCurrSystem()
{	
	if(!lpi_->setupCurrSystem()) return false; // All systems are solved!
	// Create private objects needed by the algorithm
	const int currInitNumRhs = lpi_->getCurrNumRhs();
	const RefCountPtr<const TSFCore::VectorSpace<Scalar> >
		space = this->getOperator().range(); // Range and domain should be the same for a symmetric operator!
	if(opPrecType_==COMPOSITE_OP) {
		X_ = space->createMembers(currInitNumRhs);
		R_ = space->createMembers(currInitNumRhs);
	}
	Q_ = space->createMembers(currInitNumRhs);
	Z_ = space->createMembers(currInitNumRhs);
	P_ = space->createMembers(currInitNumRhs);
	rho_.resize(currInitNumRhs);
	rho_old_.resize(currInitNumRhs);
	beta_.resize(currInitNumRhs);
	gamma_.resize(currInitNumRhs);
	alpha_.resize(currInitNumRhs);
	R_native_init_norms_.resize(currInitNumRhs);
	R_native_norms_.resize(currInitNumRhs);
	// Compute the initial residual and set it equal to the current
	// residual (for the computeIteration() function)
	if(opPrecType_==COMPOSITE_OP) {
		assign( &*X_, Teuchos::getConst(*lpi_).getCurrLhs() );
		lpi_->computeCurrCompositeResidual(NULL,NULL,&*R_);
		norms( *R_, &R_native_init_norms_[0] );
	}
	else {
		const TSFCore::MultiVector<Scalar> &R0 = lpi_->getCurrInitResidual();
		assign( &*lpi_->getCurrResidual(), R0 );
		lpi_->setCurrResidualComputed(true);
		// Remember the natural norms of this initial residual for the convergence tests
		norms( R0, &R_native_init_norms_[0] );
	}
	std::copy( &R_native_init_norms_[0], &R_native_init_norms_[0]+currInitNumRhs, &R_native_norms_[0] );
	R_native_norms_updated_ = true;
	return true;
}

template<class Scalar>
EIterateTermination NonblockCg<Scalar>::doIteration()
{
	// Ask the status test to check the status of the current systems
	const int currNumRhs = lpi_->getCurrNumRhs();
	StatusTest<Scalar> *statusTest = lpi_->getStatusTest();
	if(statusTest) {
		std::vector<EStatusType> status(currNumRhs,STATUS_UNCHECKED);
		statusTest->checkStatus(*this,currNumRhs,currNumRhs,&status[0]);
		//std::cout<<"\nNonblockCg::doIteration() : status = "; for(int k=0;k<currNumRhs;++k) std::cout<<","<<toString(status[0]); std::cout<<"\n";  
		// Get a list of the systems that have converged (according to status test)
		bool allConverged = true;
		int numToRemove = 0;
		std::vector<int> indexesToRemove(currNumRhs);
		for( int k = 0; k < currNumRhs; ++k ) {
			TEST_FOR_EXCEPTION( status[k]==STATUS_FAILED || status[k] == STATUS_NAN, Exceptions::SolverBreakdown, "Error!" );
			if(status[k]==STATUS_CONVERGED)
				indexesToRemove[numToRemove++] = k+1; // These are one-based
			else
				allConverged = false;
		}
		if(allConverged)
			return TERMINATION_STATUS_TEST;
		// Deflate systems that are converged (if not all converged)
		if(numToRemove) {
			lpi_->deflate(*this,numToRemove,&indexesToRemove[0]);
			// Deflate private workspace
			deflate(
				numToRemove, &indexesToRemove[0], currNumRhs
				,opPrecType_==COMPOSITE_OP ? 5 : 3
				,Teuchos::arrayArg<TSFCore::MultiVector<Scalar>*>(
					&*Q_, &*Z_, &*P_, X_.get(), R_.get()
					)()
				);
			deflate(
				numToRemove, &indexesToRemove[0], currNumRhs
				,5
				,Teuchos::arrayArg<Scalar*>(
					&rho_[0], &rho_old_[0], &beta_[0], &gamma_[0], &alpha_[0]
					)()
				);
			deflate(
				numToRemove, &indexesToRemove[0], currNumRhs
				,2
				,Teuchos::arrayArg<ScalarMag*>(
					&R_native_init_norms_[0], &R_native_norms_[0]
					)()
				);
		}
	}
	// Perform linear algebra for current iteration
	computeIteration();
	// Inform that the current solution and the residual are up to date
	lpi_->setCurrLhsUpdated(opPrecType_==COMPOSITE_OP?false:true);
	lpi_->setCurrResidualComputed(opPrecType_==COMPOSITE_OP?false:true);
	return TERMINATION_UNDEFINED;
}

#define BELOS_NONBLOCK_CG_SOLVER_ERR_MSG "NonblockCg<Scalar>::computeIteration(...): iteration = " << this->getCurrNumIters() << ": Error, "

template<class Scalar>
void NonblockCg<Scalar>::computeIteration()
{
	const int m = lpi_->getCurrNumRhs(); // Same as getCurrBlockSize()
	const TSFCore::Range1D currRng(1,m); // Range for active RHSs
	// Get the operator and preconditioner
	TSFCore::LinearOpHandle<Scalar> Op   = this->getOperator();
	TSFCore::LinearOpHandle<Scalar> Prec = this->getPrec();
	// Get current views of data
	RefCountPtr<TSFCore::MultiVector<Scalar> >
		X = this->getCurrLhs(),
		R = this->getCurrResidual(),
		Q = Q_->subView(currRng),
		Z = Z_->subView(currRng),
		P = P_->subView(currRng);
	// Print header
	if(get_out().get()) {
		*get_out() << "\n*** currNumIters = " << currNumIters_ << std::endl;
		if(dump_all()) {
			*get_out() << "\nX =\n" << *X;
			*get_out() << "\nR =\n" << *R;
		}
	}
	// Negate the residual to be currRhs - Operator * currLhs as according to standard CG
	scale( -ST::one(), &*R );
	RefCountPtr<const TSFCore::VectorSpace<Scalar> >
		space = R->range(); // Operator should be symmetric so any space will do!
	int j;
	if( Prec.op().get() ) {  // Preconditioner is available
		Prec.apply( TSFCore::NOTRANS, *R, &*Z );                   // Prec*R^{i-1}                    -> Z^{i-1}
	}
	else {                   // No preconditioner is available
		assign( &*Z, *R );                                         // R^{i-1}                          -> Z^{i-1}
	}
	if(get_out().get() && dump_all()) {
		*get_out() << "\nZ =\n" << *Z;
	}
	space->scalarProds( *Z, *R, &rho_[0] );                      // rho_{i-1}[j] = <Z^{i-1}[j],R^{i-1}[j]>
	for(j=0;j<m;++j) { 	// Check indefinite operator
		TEST_FOR_EXCEPTION(
			ST::isnaninf(rho_[j]), Exceptions::SolverBreakdown
			,BELOS_NONBLOCK_CG_SOLVER_ERR_MSG << "rho["<<j<<"] = " << rho_[j] << " is not valid number, the method has failed!"
			);
		TEST_FOR_EXCEPTION(
			rho_[j] <= 0.0, Exceptions::Indefinite
			,BELOS_NONBLOCK_CG_SOLVER_ERR_MSG << "rho["<<j<<"] = " << rho_[j] << " <= 0, the preconditioner is indefinite!"
			);
	}
	if(get_out().get() && dump_all()) {
		*get_out() << "\nrho =\n"; for(j=0;j<m;++j) *get_out() << " " << rho_[j]; *get_out() << std::endl;
		if(currNumIters_ > 0 )
			*get_out() << "\nrho_old =\n"; for(j=0;j<m;++j) *get_out() << " " << rho_old_[j]; *get_out() << std::endl;
	}
	for(j=0;j<m;++j) { // Check for failure: rho_{i-1} = 0
		TEST_FOR_EXCEPTION(
			rho_[j] == 0.0, Exceptions::SolverBreakdown
			,BELOS_NONBLOCK_CG_SOLVER_ERR_MSG << "rho["<<j<<"] = 0.0, the method has failed!"
			);
	}
	if( currNumIters_ == 0 ) {
		assign( &*P, *Z );                                 // Z^{i-1}                                 -> P^{i}
	}
	else {
		for(j=0;j<m;++j) beta_[j] = rho_[j]/rho_old_[j];   // rho_{i-1}[j]/rho_{i-2}[j]               -> beta_{i-1}[j]
		update( *Z, &beta_[0], 1.0, &*P );                 // Z^{i-1}[j] + beta_{i-1} * P^{i-1}[j]    -> P^{i}[j]
	}
	if(get_out().get() && dump_all()) {
		*get_out() << "\nP =\n" << *P;
	}
	Op.apply( TSFCore::NOTRANS, *P, &*Q );               // Op*P^{i}                                -> Q^{i}
	if(get_out().get() && dump_all()) {
		*get_out() << "\nQ =\n" << *Q;
	}
	space->scalarProds( *P, *Q, &gamma_[0] );            // <P^{i}[j],Q^{i}[j]>                     -> gamma_{i-1}[j]
	for(j=0;j<m;++j) { 	// Check indefinite operator
		TEST_FOR_EXCEPTION(
			ST::isnaninf(gamma_[j]), Exceptions::SolverBreakdown
			,BELOS_NONBLOCK_CG_SOLVER_ERR_MSG << "gamma["<<j<<"] = " << gamma_[j] << " is not valid number, the method has failed!"
			);
		TEST_FOR_EXCEPTION(
			gamma_[j] <= 0.0, Exceptions::Indefinite
			,BELOS_NONBLOCK_CG_SOLVER_ERR_MSG << "gamma["<<j<<"] = " << gamma_[j] << " <= 0, the operator is indefinite!"
			);
	}
	for(j=0;j<m;++j) alpha_[j] = rho_[j]/gamma_[j];      // rho_{i-1}[j] / gamma_{i-1}[j]           -> alpha_{i}[j]
	if(get_out().get() && dump_all()) {
		*get_out() << "\ngamma =\n"; for(j=0;j<m;++j) *get_out() << " " << gamma_[j]; *get_out() << std::endl;
		*get_out() << "\nalpha =\n"; for(j=0;j<m;++j) *get_out() << " " << alpha_[j]; *get_out() << std::endl;
	}
	update( &alpha_[0], +1.0, *P, &*X );                 // +alpha_{i}[j] * P^{i}[j] + X^{i-1}      -> X^{i}
	update( &alpha_[0], -1.0, *Q, &*R );                 // -alpha_{i}[j] * Q^{i}[j] + R^{i-1}      -> R^{i}
	R_native_norms_updated_ = false;
	if(get_out().get() && dump_all()) {
		*get_out() << "\nX =\n" << *X;
		*get_out() << "\nR =\n" << *R;
	}
	for(j=0;j<m;++j) rho_old_[j] = rho_[j];
	// Negate the residual to be Operator * currLhs - currRhs as according to LinearProblemState definition
	scale( -ST::one(), &*R );
}

#undef BELOS_NONBLOCK_CG_SOLVER_ERR_MSG

template<class Scalar>
TSFCore::LinearOpHandle<Scalar> NonblockCg<Scalar>::getOperator() const
{
	switch(opPrecType_) {
		case OP_ONLY: return lpi_->getOperator();
		case COMPOSITE_OP: return lpi_->getCompositeOperator();
		case RIGHT_PREC: case LEFT_PREC: return lpi_->getOperator();
		default:
			TEST_FOR_EXCEPT(true);
	}
	return TSFCore::LinearOpHandle<Scalar>(); // Should never be executed!
}

template<class Scalar>
TSFCore::LinearOpHandle<Scalar> NonblockCg<Scalar>::getPrec() const
{
	switch(opPrecType_) {
		case OP_ONLY: case COMPOSITE_OP: return TSFCore::LinearOpHandle<Scalar>();
		case RIGHT_PREC: return lpi_->getRightPrec();
		case LEFT_PREC: return lpi_->getLeftPrec();
		default:
			TEST_FOR_EXCEPT(true);
	}
	return TSFCore::LinearOpHandle<Scalar>(); // Should never be executed!
}

template<class Scalar>
RefCountPtr<TSFCore::MultiVector<Scalar> > NonblockCg<Scalar>::getCurrLhs() const
{
	if(opPrecType_ == COMPOSITE_OP)
		return X_->subView(TSFCore::Range1D(1,lpi_->getCurrBlockSize()));
	else
		return lpi_->getCurrLhs();
}

template<class Scalar>
RefCountPtr<TSFCore::MultiVector<Scalar> > NonblockCg<Scalar>::getCurrResidual() const
{
	if(opPrecType_ == COMPOSITE_OP)
		return R_->subView(TSFCore::Range1D(1,lpi_->getCurrBlockSize()));
	else
		return lpi_->getCurrResidual();
}

template<class Scalar>
void NonblockCg<Scalar>::cleanUp()
{
	X_ = null;
	R_ = null;
	Q_ = null;
	Z_ = null;
	P_ = null;
	// ToDo: Resize arrays to zero?
}

} // end Belos namespace

#endif // BELOS_NONBLOCK_CG_HPP
