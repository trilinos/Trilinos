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

#ifndef BELOS_BASIC_ITERATION_STATE_HPP
#define BELOS_BASIC_ITERATION_STATE_HPP

#include "Belos_Types.hpp"

namespace Belos {

///
/** Base interface of a linear solver for access by status tests.
 *
 * ToDo: Finish documentation!
 */
template<class Scalar>
class BasicIterationState {
public:
	
	///
	virtual ~BasicIterationState() {}

	///
	/** Return the current state of the linear problem.
	 */
	virtual const LinearProblemState<Scalar>& getProblem() const = 0;

	///
	/** Return the current total number of iterations over all restarts
	 * for current block of systems.
	 *
	 * Preconditions:<ul>
	 * <li><tt>this->getProblem().getCurrNumRhs() > 0</tt>
	 * </ul>
	 */
	virtual int getCurrNumIters() const = 0;

	///
	/** Return the number of restarts that have been performed.
	 *
	 * Preconditions:<ul>
	 * <li><tt>this->getProblem().getCurrNumRhs() > 0</tt>
	 * </ul>
	 */
	virtual int getCurrNumRestarts() const = 0;

	///
	/** Force the computation of the current LHS \f$\bar{X}\f$ if not done so already.
	 *
	 * Preconditions:<ul>
	 * <li><tt>this->getProblem().getCurrNumRhs() > 0</tt>
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li><tt>this->getProblem().isCurrLhsUpdated() == true</tt>
	 * </ul>
	 */
	virtual void forceCurrLhsUpdate() const = 0;

	///
	/** Returns the type of the current native residual returned from <tt>getCurrNativeResiduals()</tt>.
	 */
	virtual ENativeResidualType getCurrNativeResidualType() const = 0;

	///
	/** Get the current native residual norms.
	 *
	 * \param  currBlockSize  [in] The current block size.
	 * \param  norms          [out] Array (length <tt>currBlockSize</tt>) of relative native norms
	 *                        for scaled preconditioned residual (see ???).  If <tt>norms==NULL</tt> then
	 *                        these norms will not be returned.
	 * \param  residuals      [out] Pointer to a pointer to the native scaled preconditioned
	 *                        residuals if it exits.  If <tt>residuals!=NULL</tt> is passed in then
	 *                        this is a request for the native residual \f$\breve{R}\f$ it is readily
	 *                        available.  If <tt>residuals==NULL</tt> is passed in then the native
	 *                        residual \f$\breve{R}\f$ is not requested.
	 *
	 * Preconditions:<ul>
	 * <li><tt>this->getProblem().getCurrNumRhs() > 0</tt>
	 * <li><tt>this->getProblem().getCurrBlockSize() == currBlockSize</tt>
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li>If <tt>norms!=NULL</tt> then <tt>norms[]</tt> will contain the relative norms
	 *     of the native residuals defined as:
	 *     <ul>
	 *     <li>If <tt>this->getCurrNativeResidualType()==NATIVE_RESIDUAL_UNPRECONDITIONED</tt> then
	 *       <ul><li><tt>norms[k]</tt> = \f$||S_L \bar{R}_{(:,k+1)}|| / ||S_L (\bar{R}_0)_{(:,k+1)}||\f$, for <tt>k = 0 ... currBlockSize</tt></ul>
	 *     <li>ElseIf <tt>this->getCurrNativeResidualType()==NATIVE_RESIDUAL_PRECONDITIONED</tt> then
	 *       <ul><li><tt>norms[k]</tt> = \f$||P_L S_L \bar{R}_{(:,k+1)}|| / ||P_L S_L (\bar{R}_0)_{(:,k+1)}||\f$, for <tt>k = 0 ... currBlockSize</tt></ul>
	 *     </ul>
	 * <li>If <tt>residuals!=NULL</tt> then
	 *     <ul>
	 *     <li>If <tt>*residuals==NULL</tt> then no native residual multi-vector exists (e.g. GMRES)
	 *     <li>If <tt>*residuals!=NULL</tt> the a native residual multi-vector is avaiable (e.g. CG).
   *         If this value is returned then this multi-vector should be queried quickly and then
	 *         forgoten!
	 *     </ul>
	 * </ul>
	 */
	virtual void getCurrNativeResiduals(
		const int                              currBlockSize
		,Scalar                                norms[]
		,const TSFCore::MultiVector<Scalar>*   *residuals  = NULL
		) const = 0;

};

} // end Belos namespace

#endif // BELOS_BASIC_ITERATION_STATE_HPP
