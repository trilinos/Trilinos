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

#ifndef BELOS_RESIDUAL_NORM_STATUS_TEST_HPP
#define BELOS_RESIDUAL_NORM_STATUS_TEST_HPP

#include "Belos_ToleranceBasedStatusTestBase.hpp"

namespace Belos {

///
/** Defines a status test based on the relative norm of the unscaled
 * unpreconditioned residual.
 *
 * ToDo: Finish documentation!
 */
template<class Scalar>
class ResidualNormStatusTest : public ToleranceBasedStatusTestBase<Scalar> {
public:

	///
	typedef Teuchos::ScalarTraits<Scalar>  ST;

	/// Set the stream that output will be sent to
	STANDARD_COMPOSITION_MEMBERS( std::ostream, out );

	/// Set the leading string that will be printed at the beginning of each new line of output.
	STANDARD_MEMBER_COMPOSITION_MEMBERS( std::string, leadingOutputStr );

	/** @name Constructor/initializers/accessors */
	//@{

	///
	/** Construct to uninitialized.
	 */
	ResidualNormStatusTest();

	///
	/** Calls <tt>initialize(tol)</tt>
	 */
	ResidualNormStatusTest(
		const typename Teuchos::ScalarTraits<Scalar>::magnitudeType        tol
		);

	///
	/** Calls <tt>initialize(totalNumRhs,tols)</tt>
	 */
	ResidualNormStatusTest(
		const int                                                          totalNumRhs
		,const typename Teuchos::ScalarTraits<Scalar>::magnitudeType       tols[]
		);

	//@}

protected:

	/** @name Overridden from AttachStatusTestBase */
	//@{
	///
	void protectedCheckStatus(
		const BasicIterationState<Scalar>         &bis
		,const int                                currBlockSize
		,const int                                currNumRhs
		,EStatusType                              status[]
		);
	//@}

private:

	mutable std::vector<typename ST::magnitudeType>  R_bar_norms_;
	mutable std::vector<typename ST::magnitudeType>  R_bar_0_norms_;
	mutable std::vector<int>                         currRhsIndexes_;

};

// //////////////////////////////////
// Implementation

// Constructor/initializers/accessors

template<class Scalar>
ResidualNormStatusTest<Scalar>::ResidualNormStatusTest()
{}

template<class Scalar>
ResidualNormStatusTest<Scalar>::ResidualNormStatusTest(
	const typename Teuchos::ScalarTraits<Scalar>::magnitudeType        tol
	)
{
	this->initialize(tol);
}

template<class Scalar>
ResidualNormStatusTest<Scalar>::ResidualNormStatusTest(
	const int                                                          totalNumRhs
	,const typename Teuchos::ScalarTraits<Scalar>::magnitudeType       tols[]
	)
{
	this->initialize(totalNumRhs,tols);
}

// Overridden from AttachStatusTestBase

template<class Scalar>
void ResidualNormStatusTest<Scalar>::protectedCheckStatus(
	const BasicIterationState<Scalar>         &bis
	,const int                                currBlockSize
	,const int                                currNumRhs
	,EStatusType                              status[]
	)
{
	const std::vector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &tols = this->tols();
	const LinearProblemState<Scalar> &lps = bis.getProblem();
	// Force update of currLhs
	if(!lps.isCurrLhsUpdated()) bis.forceCurrLhsUpdate();	
	// Get the current and initial unscaled unpreconditioned residuals and compute their norms
	const TSFCore::MultiVector<Scalar>
		&R_bar_0 = lps.getCurrInitResidual(),
		&R_bar   = lps.getCurrResidual();
	if(static_cast<int>(R_bar_norms_.size()) < currBlockSize) {
		R_bar_norms_.resize(currBlockSize);
		R_bar_0_norms_.resize(currBlockSize);
	}
	TSFCore::norms( R_bar, &R_bar_norms_[0] );
	TSFCore::norms( R_bar_0, &R_bar_0_norms_[0] );
	if(tols.size()==1 && currNumRhs > 1) {
		// One tolerance for all RHSs.
		TEST_FOR_EXCEPT(true);
	}
	else {
		// Must match a tolerance to each specific RHS
		if( static_cast<int>(currRhsIndexes_.size()) < currNumRhs ) currRhsIndexes_.resize(currNumRhs);
		bis.getProblem().getCurrRhsIndexes( currNumRhs, &currRhsIndexes_[0] );
		for( int k = 0; k < currNumRhs; ++k ) {
			const int origRhsIndex = currRhsIndexes_[k];
#ifdef TEUCHOS_DEBUG
			TEST_FOR_EXCEPT( origRhsIndex > static_cast<int>(tols.size())  );
#endif
			const typename ST::magnitudeType R_rel_norm = R_bar_norms_[k] / R_bar_0_norms_[k];
			if( ST::isnaninf(R_rel_norm) )                 status[k] = STATUS_NAN;
			else if( R_rel_norm <= tols[origRhsIndex-1] )  status[k] = STATUS_CONVERGED;
			else                                           status[k] = STATUS_UNCONVERGED;
		}
	}
}

} // end Belos namespace

#endif // BELOS_RESIDUAL_NORM_STATUS_TEST_HPP
