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

#ifndef BELOS_TYPES_HPP
#define BELOS_TYPES_HPP

#include "BelosConfigDefs.hpp"
#include "TSFCoreTypes.hpp"
#include "Teuchos_RefCountPtr.hpp"

/*!	\file Belos_Types.hpp
*/

namespace Belos {

using Teuchos::RefCountPtr;
using Teuchos::rcp;
using Teuchos::rcp_const_cast;
const Teuchos::ENull null = Teuchos::null;

///
/** Determines the symmetry property of an operator.
 */
enum EOpSymmetry {
	OP_SYMMETRIC       ///< The operator is assumed to be symmetric
	,OP_UNSYMMETRIC    ///< The operator is assumed to be unsymmetric
};

///
enum EStatusType {
 	STATUS_UNCHECKED = 2,   ///< Initial state of status
	STATUS_UNCONVERGED = 1, ///< Convergence is not reached
	STATUS_CONVERGED = 0,   ///< Convergence is reached.
	STATUS_FAILED = -1,     ///< Some failure occured.  Should stop.
	STATUS_NAN = -2         ///< Result from test contains a NaN value.  Should stop.
};

///
inline
const char* toString( const EStatusType status )
{
 	if(status==STATUS_UNCHECKED) return "STATUS_UNCHECKED";
	if(status==STATUS_UNCONVERGED) return "STATUS_UNCONVERGED";
	if(status==STATUS_CONVERGED) return "STATUS_CONVERGED";
	if(status==STATUS_FAILED) return "STATUS_FAILED";
	if(status==STATUS_NAN) return "STATUS_NAN";
#ifdef TEUCHOS_DEBUG
	TEST_FOR_EXCEPT(true); // Should never get here!
#endif
	return "*** Invalid value for status ***";
}

///
enum EIterateTermination {
	TERMINATION_UNDEFINED       ///< Not defined (not a valid value)
	,TERMINATION_STATUS_TEST    ///< The status test terminated before maxNumIter was exceeded
	,TERMINATION_MAX_NUM_ITER   ///< The status test terminated because maxNumIter iterations was exceeded.
};

///
inline
const char* toString( const EIterateTermination iterateTermination )
{
	if( iterateTermination == TERMINATION_UNDEFINED )
		return "TERMINATION_UNDEFINED";
	if( iterateTermination == TERMINATION_STATUS_TEST )
		return "TERMINATION_STATUS_TEST";
	if( iterateTermination == TERMINATION_MAX_NUM_ITER )
		return "TERMINATION_MAX_NUM_ITER";
#ifdef TEUCHOS_DEBUG
	TEST_FOR_EXCEPT(true); // Should never get here!
#endif
	return "*** Invalid value for iterateTermination ***";
}

///
struct IterateReturn {
	IterateReturn()
		:iterateTermination(TERMINATION_UNDEFINED),numCumulativeIter(-1) {} 
	IterateReturn(EIterateTermination iterateTermination_in, int numCumulativeIter_in)
		:iterateTermination(iterateTermination_in),numCumulativeIter(numCumulativeIter_in) {} 
	EIterateTermination    iterateTermination;
	int                    numCumulativeIter;
};

///
enum ENativeResidualType {
	NATIVE_RESIDUAL_UNPRECONDITIONED  ///< The native residual is scaled and unpreconditioned \f$S_L \bar{R}\f$.
	,NATIVE_RESIDUAL_PRECONDITIONED   ///< The native residual is scaled and preconditioned \f$P_L S_L \bar{R}\f$.
};

///
inline
const char* toString( const ENativeResidualType nativeResidualType )
{
 	if(nativeResidualType==NATIVE_RESIDUAL_UNPRECONDITIONED) return "NATIVE_RESIDUAL_UNPRECONDITIONED";
	if(nativeResidualType==NATIVE_RESIDUAL_PRECONDITIONED) return "NATIVE_RESIDUAL_PRECONDITIONED";
#ifdef TEUCHOS_DEBUG
	TEST_FOR_EXCEPT(true); // Should never get here!
#endif
	return "*** Invalid value for nativeResidualType ***";
}

namespace Exceptions {

/** \defgroup BelosExceptions_grp Basic Belos exception types.
 */
//@{

/// Thrown if a linear solve failed
class FailureToConverge : public std::logic_error
{public: FailureToConverge(const std::string& what_arg) : std::logic_error(what_arg) {}};

/// Thrown if the solver just can not make any more progress for some reason
class SolverBreakdown : public FailureToConverge
{public: SolverBreakdown(const std::string& what_arg) : FailureToConverge(what_arg) {}};

/// Thrown if operator turns out to be indefinite
class Indefinite : public SolverBreakdown
{public: Indefinite(const std::string& what_arg) : SolverBreakdown(what_arg) {}};

//@}

} // namespace Exceptions

// Forward declarations for interface types
template<class Scalar> class StatusTest;
template<class Scalar> class LinearProblemState;
template<class Scalar> class LinearProblemIteration;
template<class Scalar> class LinearProblemSetup;
template<class Scalar> class BasicIterationState;
template<class Scalar> class BasicIteration;

} // end Belos namespace

#endif // BELOS_TYPES_HPP
