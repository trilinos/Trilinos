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

#ifndef BELOS_BASIC_ITERATION_HPP
#define BELOS_BASIC_ITERATION_HPP

#include "Belos_BasicIterationState.hpp"

namespace Belos {

///
/** Base interface for a linear solver for access by clients of linear solvers.
 *
 * ToDo: Finish documentation!
 */
template<class Scalar>
class BasicIteration : public BasicIterationState<Scalar> {
public:

	///
	/** Returns <tt>true</tt> of adjoints for the operator and preconditioners are required.
	 */
	virtual bool adjointRequired() const = 0;

	///
	/** Set the linear problem to be solved.
	 *
	 * @param  lpi  [in] The linear problem to be solved.  After this object is set,
	 *              the underlying definition of the linear problem should not be
	 *              distrubed until after <tt>this->finalize()</tt> is called
	 *              after <tt>this->iterate()</tt> is called.
	 *
	 * Preconditions:<ul>
	 * <li><tt>lpi.get()!=NULL</tt>
	 * <li><tt>*lpi</tt> is fully setup for the problem defintion.
	 * <li>[<tt>this->adjointRequired()==true]
	 *     <tt>op.op()->adjointSupported(TSFCore::not_trans(lpi.getOperator().defaultTrans()))==true</tt>
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li><tt>this->getProblem().get() == lpi.get()</tt>
	 * </ul>
	 */
	virtual void setProblem( const RefCountPtr<LinearProblemIteration<Scalar> > &lpi ) = 0;

	///
	/** Return smart pointer to the <tt>LinearProblemIteration</tt> object.
	 */
	virtual RefCountPtr<LinearProblemIteration<Scalar> > getProblem() = 0;

	///
	/** Initialize the linear solver ???
	 */
	virtual void initialize() = 0;

	///
	/** Iterate on the linear problem.
	 *
	 * \param  maxNumIter  [in] The maximum nubmer of iterations per block of RHSs
	 *                     allowed.
	 *
	 * Postconditions:<ul>
	 * <li>The linear problem should be left in valid state for any return
	 *     (i.e. standard return or if an exception is thrown).
	 * </ul>
	 *
	 * \exception Exceptions::SolverBreakdown
	 * This exception is thrown if the status test returns <tt>STATUS_FAILED</tt>
	 * or <tt>STATUS_NAN</tt> or if any other numerical problem is detected in
	 * the solution procedure.
	 *
	 * \exception Exceptions::Indefinite
	 * This exception is thrown if the underlying solver is only for (positive) definite
	 * systems but the operator and/or the preconditioner is found to be indefinite.
	 * 
	 * \return Returns
	 * <tt>return.iterationTermination==TERMINATION_STATUS_TEST</tt> if
	 * the status test caused termination.  Returns
	 * <tt>return.iterationTermination==TERMINATION_MAX_NUM_ITER</tt> if
	 * the maximum number of iterations <tt>maxNumIter</tt> was exceeded
	 * for any of the block solves.  In any case,
	 * <tt>return.numCumulativeIter</tt> gives the total cumulative
	 * number (over all right-hand-side blocks) of iterations performed
	 * by the iterative solver.  Note that therefore
	 * <tt>return.numCumulativeIter</tt> may be larger than
	 * <tt>maxNumIter</tt> if the block size is smaller than the total
	 * number of right-hand sides.  Also note that
	 * <tt>return.numCumulativeIter</tt> also really only makes since in
	 * the context of the maximum number of RHSs and the block size.
	 */
	virtual IterateReturn iterate( const int maxNumIter ) = 0;

	///
	/** Finalize the iterative solver after <tt>iterate()</tt> returns ???.
	 */
	virtual void finalize() = 0;
	
	// ToDo: Fill in the rest of the operations

};

} // end Belos namespace

#endif // BELOS_BASIC_ITERATION_HPP
