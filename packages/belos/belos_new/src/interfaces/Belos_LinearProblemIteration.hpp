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

#ifndef BELOS_LINEAR_PROBLEM_ITERATION_HPP
#define BELOS_LINEAR_PROBLEM_ITERATION_HPP

#include "Belos_LinearProblemState.hpp"

namespace Belos {

///
/** Base interface of a linear problem being solved for access by an
 * iterative solver implementation.
 *
 * ToDo: Finish documentation!
 */
template<class Scalar>
class LinearProblemIteration : public LinearProblemState<Scalar> {
public:

	/** @name Pure virtual methods */
	//@{

	///
	/** Get the block size as set by external client.
	 *
	 * Preconditions:<ul>
	 * <li><tt>this->getTotalNumRhs() > 0</tt>
	 * </ul>
	 */
	virtual int getBlockSize() const = 0;

	///
	/** Get the current status test if it exists.
	 *
	 * Preconditions:<ul>
	 * <li><tt>this->getTotalNumRhs() > 0</tt>
	 * </ul>
	 */
	virtual StatusTest<Scalar>* getStatusTest() = 0;

	///
	/** Setup the next set of linear systems to solve.
	 *
	 * @param  firstRhsOffset  [in] Offset to first RHS in next block of systems.
	 * @param  numRhs          [in] Number of RHS in next block of systems.
	 * @param  blockSize       [in] The block size to use.  If <tt>blockSize > numRhs</tt>
	 *                         then current block of systems will be augmented.
	 *
	 * Preconditions:<ul>
	 * <li><tt>firstRhsOffset+numRhs <= this->getTotalNumRhs()</tt>
	 * <li><tt>numRhs <= blockSize</tt>
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li><tt>this->getCurrNumRhs() == numRhs</tt>
	 * <li><tt>this->getCurrBlockSize() == blockSize</tt>
	 * <li><tt>this->isCurrInitResidualComputed() == false</tt>
	 * <li><tt>this->isCurrLhsUpdated() == true</tt>
	 * <li><tt>this->isCurrResidualComputed() == false</tt>
	 * <li><tt>currRhsIndexes[k] = k</tt>, for <tt>k=0...numRhs</tt> where <tt>currRhsIndexes</tt>
	 *     is returned from <tt>this->getCurrRhsIndexes((currRhsIndexes)</tt>.
	 * </ul>
	 */
	virtual void setCurrSystem(
		const int                   firstRhsOffset
		,const int                  numRhs
		,const int                  blockSize
		) = 0;

	///
	/** Get smart pointer to mutable current LHS \f$\bar{X}\f$.
	 *
	 * Preconditions:<ul>
	 * <li><tt>this->getCurrNumRhs() > 0</tt>
	 * </ul>
	 *
	 * This function gives an iterative solver implementation access to
	 * the internal current LHS \f$\bar{X}\f$ for direct updating.
	 * Directly updating \f$\bar{X}\f$ would be advantageous for any
	 * solver that directly updates the unscaled unpreconditioned
	 * solution \f$\bar{X}\f$ such as a CG, BiCG, BiCGStab or other
	 * solver where there is no right scaling or right preconditioning.
	 */
	virtual RefCountPtr<TSFCore::MultiVector<Scalar> > getCurrLhs() = 0;

	///
	/** Inform that current LHS is up to date or not.
	 *
	 * @param  isCurrLhsUpdated  [in] If <tt>true</tt> then <tt>this->getCurrLhs()</tt> is current.
	 *                         If <tt>false</tt> then <tt>this->getCurrLhs()</tt> is not current.
	 *
	 * Preconditions:<ul>
	 * <li><tt>this->getCurrNumRhs() > 0</tt>
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li><tt>this->isCurrLhsUpdated() == isCurrLhsUpdated</tt>
	 * <li><tt>this->isCurrResidualComputed() == false</tt>
	 * </ul>
	 */
	virtual void setCurrLhsUpdated( const bool isCurrLhsUpdated ) = 0;

	///
	/** Update the current solution given an update MultiVector.
	 *
	 * @param  nativeLhsStep  [in] Multivector \f$\bar{Z}\f$ that defines the
	 *                        update \f$\bar{X} = \bar{X}_{0,r} + S_R P_R \bar{Z}\f$.
	 *
	 * Preconditions:<ul>
	 * <li><tt>this->getCurrNumRhs() > 0</tt>
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li><tt>this->isCurrLhsUpdated() == true</tt>
	 * <li><tt>this->isCurrResidualComputed() == false</tt>
	 * </ul>
	 */
	virtual void updateCurrLhs( const TSFCore::MultiVector<Scalar> &nativeLhsStep ) = 0;

	///
	/** Set the current solution given a native LHS MultiVector.
	 *
	 * @param  nativeLhs  [in] Multivector \f$\breve{X}\f$ that defines the
	 *                    update \f$\bar{X} = S_R P_R \breve{X}\f$.
	 *
	 * Preconditions:<ul>
	 * <li><tt>this->getCurrNumRhs() > 0</tt>
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li><tt>this->isCurrLhsUpdated() == true</tt>
	 * <li><tt>this->isCurrResidualComputed() == false</tt>
	 * </ul>
	 */
	virtual void setCurrLhs( const TSFCore::MultiVector<Scalar> &nativeLhs ) = 0;

	///
	/** Get smart pointer to mutable current residual \f$\bar{R}\f$.
	 *
	 * Preconditions:<ul>
	 * <li><tt>this->getCurrNumRhs() > 0</tt>
	 * </ul>
	 * 
	 * This function allows an interative solver to directly update the
	 * current residual if it is convenient to do so.  For example, a
	 * standard CG algorithm that uses just one preconditioner and no
	 * scaling will directly update the unscaled, unpreconditioned
	 * residual \f$\bar{R}\f$ and therefore save on storage and
	 * computational cost.
	 */
	virtual RefCountPtr<TSFCore::MultiVector<Scalar> > getCurrResidual() = 0;

	///
	/** Inform that current Residual is up to date or not.
	 *
	 * @param  isCurrResidualComputed  [in] If <tt>true</tt> then <tt>this->getCurrResidual()</tt> is current.
	 *                                If <tt>false</tt> then <tt>this->getCurrResidual()</tt> is not current.
	 *
	 * Preconditions:<ul>
	 * <li><tt>this->getCurrNumRhs() > 0</tt>
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li><tt>this->isCurrResidualComputed() == isCurrResidualComputed</tt>
	 * </ul>
	 */
	virtual void setCurrResidualComputed( const bool isCurrResidualComputed ) = 0;

	///
	/** Compute the preconditioned residual \f$P_L S_L \bar{R}\f$
	 * consistent with the conbined operator returned from
	 * <tt>this->getCombinedOperator()>/tt>.
	 *
	 * @param  currLhs  [in] Pointer to current LHS to use.  If <tt>currLhs==NULL</tt> then
	 *                  <tt>this->getCurrLhs()</tt> will be used in its place.
	 * @param  currRhs  [in] Pointer to current RHS to use.  If <tt>currRhs==NULL</tt> then
	 *                  <tt>this->getCurrRhs()</tt> will be used in its place.
	 * @param  currPrecResidual
	 *                  [out] The preconditioned residual that must be been created
	 *                  using <tt>this->getCombinedOperator().range()->createMembers(this->getCurrBlockSize())</tt>.
	 *
	 * Preconditions:<ul>
	 * <li><tt>this->getCurrNumRhs() > 0</tt>
	 * <li>[<tt>currLhs==NULL</tt>] <tt>this->isCurrLhsUpdated() == true</tt>
	 * <li><tt>currPrecResidual != NULL</tt>
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li>[currLhs==NULL && currRhs==NULL</tt>] <tt>this->isCurrResidualComputed() == true</tt>
	 * </ul>
	 *
	 * Note, this function would be called by any iterative solver that
	 * wants to be oblivious as to what type of preconditioning or scaling
	 * is being used.
	 */
	virtual void computeCurrPrecResidual(
		const TSFCore::MultiVector<Scalar>    *currLhs
		,const TSFCore::MultiVector<Scalar>   *currRhs
		,TSFCore::MultiVector<Scalar>         *currPrecResidual
		) = 0;

	///
	/** Deflate out some of the current RHSs.
	 *
	 * @param  bis              [in] The iterative solver.
	 * @param  numCurrRhsToRemove
	 *                          [in] The number of current systems to deflate out.
	 * @param  currRhsIndexesToRemove
	 *                          [in] Array (length <tt>numCurrRhsToRemove</tt>) of one-based
	 *                          indexes in the current block to remove.
	 *                          The entries <tt>k=currRhsIndexesToRemove[i]</tt> give one-based
	 *                          columns in <tt>this->getCurrRhs()->col(k)</tt> (and other
	 *                          current quantities) to remove.
	 *
	 * Preconditions:<ul>
	 * <li><tt>this->getCurrNumRhs() > 0</tt>
	 * <li><tt>0 < currRhsIndexesToRemove[i] <= this->getCurrNumRhs()</tt>, for <tt>i=0...numCurrRhsToRemove-1</tt>
	 * <li><tt>currRhsIndexesToRemove[i] < currRhsIndexesToRemove[i+1], for <tt>i=0...numCurrRhsToRemove-2</tt>
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> ToDo: Fill these out!
	 * </ul>
	 */
	virtual void deflate(
		const BasicIterationState<Scalar>        &bis
		,const int                               numCurrRhsToRemove
		,const int                               currRhsIndexesToRemove[]
		) = 0;

	///
	/** Setup for a restart by copying \f$\bar{X}\f$ into \f$\bar{X}_{0,r}\f$.
	 *
	 * Preconditions:<ul>
	 * <li><tt>this->getCurrNumRhs() > 0</tt>
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li><tt>this->isCurrLhsUpdated() == true</tt>
	 * <li><tt>this->getCurrInitLhs()</tt> contains a copy of <tt>this->getCurrLhs()</tt>.
	 * </ul>
	 *
	 * If <tt>this->isCurrLhsUpdated()==false</tt> before this function
	 * is called then <tt>bis.forceCurrLhsUpdate()</tt> will be called
	 * before copying <tt>this->getCurrLhs()</tt> into
	 * <tt>this->getCurrInitLhs()</tt>.
	 */
	virtual void restart( const BasicIterationState<Scalar> &bis ) = 0;

	///
	/** Set what is left in <tt>this->getCurrLhs()</tt> into
	 * <tt>this->getLhs()</tt> and invalidate current block of systems.
	 *
	 * Preconditions:<ul>
	 * <li><tt>this->getCurrNumRhs() > 0</tt>
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li><tt>this->getCurrNumRhs() == 0</tt>, i.e. <tt>*this</tt> will
	 *     no longer have a current set of linear systems.
	 * </ul>
	 *
	 * If <tt>this->isCurrLhsUpdated()==false</tt> before this function
	 * is called then <tt>bis.forceCurrLhsUpdate()</tt> will be called
	 * before the solution is copied.
	 */
	virtual void setCurrToFullLhs( const BasicIterationState<Scalar> &bis ) = 0;

	// ToDo: Fill in the rest of the operations

	//@}

};

} // end Belos namespace

#endif // BELOS_LINEAR_PROBLEM_ITERATION_HPP
