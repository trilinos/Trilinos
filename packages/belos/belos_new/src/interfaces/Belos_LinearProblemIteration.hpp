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

	///
	using LinearProblemState<Scalar>::getCurrLhs;

	/** @name Overall linear system */
	//@{

	///
	/** Return the total number of RHSs to be solved.
	 *
	 * Postconditions:<ul>
	 * <li><tt>return == this->getRhs().domain()->dim()</tt>
	 * </ul>
	 *
	 * A return value of <tt>return==0</tt> is a flag that the linear
	 * problem is not sufficiently initialized.
	 */
	virtual int getTotalNumRhs() const = 0;

	///
	/** Return non-persisting reference to the operator \f$\tilde{A}\f$.
	 *
	 * Preconditions:<ul>
	 * <li><tt>this->getTotalNumRhs() > 0</tt>
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li><tt>return.op() != NULL</tt>
	 * </ul>
	 */
	virtual TSFCore::LinearOpHandle<Scalar> getOperator() const = 0;
	
	///
	/** Return if the operator \f$\tilde{A}\f$ is symmetric or not.
	 *
	 * Preconditions:<ul>
	 * <li><tt>this->getTotalNumRhs() > 0</tt>
	 * </ul>
	 */
	virtual EOpSymmetry getOperatorSymmetry() const = 0;

	///
	/** Return non-persisting reference to the right preconditioner \f$P_R\f$ if it exits.
	 *
	 * Preconditions:<ul>
	 * <li><tt>this->getTotalNumRhs() > 0</tt>
	 * </ul>
	 */
	virtual TSFCore::LinearOpHandle<Scalar> getRightPrec() const = 0;
	
	///
	/** Return if the right preconditioner \f$P_R\f$ is symmetric or not.
	 *
	 * Preconditions:<ul>
	 * <li><tt>this->getRightPrec() != NULL</tt>
	 * </ul>
	 */
	virtual EOpSymmetry getRightPrecSymmetry() const = 0;

	///
	/** Return non-persisting reference to the left preconditioner \f$P_L\f$ if it exits.
	 *
	 * Preconditions:<ul>
	 * <li><tt>this->getTotalNumRhs() > 0</tt>
	 * </ul>
	 */
	virtual TSFCore::LinearOpHandle<Scalar> getLeftPrec() const = 0;
	
	///
	/** Return if the left preconditioner \f$P_L\f$ is symmetric or not.
	 *
	 * Preconditions:<ul>
	 * <li><tt>this->getLeftPrec() != NULL</tt>
	 * </ul>
	 */
	virtual EOpSymmetry getLeftPrecSymmetry() const = 0;

	///
	/** Return non-persisting reference to the composite operator
	 * \f$\hat{A} = P_L \tilde{A} P_R\f$ as a composite
	 * <tt>TSFCore::LinearOp</tt> object.
	 *
	 * Preconditions:<ul>
	 * <li><tt>this->getCurrNumRhs() > 0</tt>
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li><tt>return.op() != NULL</tt>
	 * </ul>
	 */
	virtual TSFCore::LinearOpHandle<Scalar> getCompositeOperator() const = 0;

	///
	/** Return pointer to the right scaling vector \f$S_R\f$ if it exists.
	 *
	 * Preconditions:<ul>
	 * <li><tt>this->getTotalNumRhs() > 0</tt>
	 * </ul>
	 */
	virtual const TSFCore::Vector<Scalar>* getRightScaling() const = 0;

	///
	/** Return pointer to the left scaling vector \f$S_L\f$ if it exists.
	 *
	 * Preconditions:<ul>
	 * <li><tt>this->getTotalNumRhs() > 0</tt>
	 * </ul>
	 */
	virtual const TSFCore::Vector<Scalar>* getLeftScaling() const = 0;

	///
	/** Return const reference to the original RHS \f$B\f$.
	 *
	 * Preconditions:<ul>
	 * <li><tt>this->getTotalNumRhs() > 0</tt>
	 * </ul>
	 */
	virtual const TSFCore::MultiVector<Scalar>& getRhs() const = 0;

	///
	/** Return non-const reference to the original LHS \f$B\f$.
	 *
	 * Preconditions:<ul>
	 * <li><tt>this->getTotalNumRhs() > 0</tt>
	 * </ul>
	 *
	 * Warning! Manually modifying this while the iteration process is
	 * taking place may be dangerous and it not recommended in general.
	 */
	virtual TSFCore::MultiVector<Scalar>& getLhs() const = 0;

	//@}

	/** @name Properties for selecting current block system */
	//@{

	///
	/** Set the block size \f$b\f$.
	 *
	 * @param  blockSize  [in] The block size \f$b\f$.
	 *
	 * Preconditions:<ul>
	 * <li><tt>blockSize > 0</tt>
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li><tt>this->getBlockSize() == blockSize</tt>
	 * <li><tt>this->getCurrNumRhs() == 0</tt>
	 * </ul>
	 */
	virtual void setBlockSize( const int blockSize ) = 0;

	///
	/** Get the block size.
	 */
	virtual int getBlockSize() const = 0;

	///
	/** Set if augmentation is allowed or not.
	 *
	 * @param  augmentationAllowed  [in]
	 *
	 * Postconditions:<ul>
	 * <li><tt>this->augmentationAllowed() == augmentationAllowed</tt>
	 * <li><tt>this->getCurrNumRhs() == 0</tt>
	 * </ul>
	 */
	virtual void setAugmentationAllowed( const bool augmentationAllowed ) = 0;

	///
	/** Get if augmentation is allowed or not.
	 */
	virtual bool getAugmentationAllowed() const = 0;

	//@}

	/** @name Manipulating the current block system */
	//@{

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
	 * Preconditions:<ul>
	 * <li><tt>this->getTotalNumRhs() > 0</tt>
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li><tt>this->getCurrNumRhs() > 0</tt>
	 * <li>If <tt>this->getAugmentationAllowed() == true</tt> then
	 *     <ul><li><tt>this->getCurrBlockSize() == this->getBlockSize()</tt></ul>
	 *     else
	 *     <ul><li><tt>this->getCurrBlockSize() <= this->getBlockSize()</tt></ul>
	 * <li><tt>this->isCurrInitResidualComputed() == false</tt>
	 * <li><tt>this->isCurrLhsUpdated() == true</tt>
	 * <li><tt>this->isCurrResidualComputed() == false</tt>
	 * <li><tt>currRhsIndexes[k] = k</tt>, for <tt>k=0...this->getCurrNumRhs()-1</tt> where <tt>currRhsIndexes</tt>
	 *     is returned from <tt>this->getCurrRhsIndexes((currRhsIndexes)</tt>.
	 * </ul>
	 *
	 * @return Returns <tt>true</tt> if a new current block of RHSs has
	 * been setup and <tt>faluse</tt> if no more RHS remain.
	 */
	virtual bool setupCurrSystem() = 0;

	///
	/** Get smart pointer to mutable current unscaled unpreconditioned LHS \f$\bar{X}\f$.
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
	 *
	 * Note that if right scaling is present then direct access to this
	 * solution multi-vector probably should not be sought.  When in
	 * doubt, a concrete <tt>Belos::BasicIteration</tt> subclass should
	 * not call this function but should instead update the current
	 * unscaled unpreconditioned solution using
	 * <tt>this->setCurrLhs()</tt>.
	 */
	virtual RefCountPtr<TSFCore::MultiVector<Scalar> > getCurrLhs() = 0;

	///
	/** Inform that current unscaled unpreconditioned LHS \f$\bar{X}\f$ is up to date or not.
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
	 *
	 * This function should only be called if <tt>this->getCurrLhs()</tt> was
	 * called to get direct access to the unscaled unpreconditioned
	 * solution multi-vector for direct update by the iterative solver
	 * implementation.
	 */
	virtual void setCurrLhsUpdated( const bool isCurrLhsUpdated ) = 0;

	///
	/** Update the current unscaled unpreconditoned LHS \f$\bar{X}\f$ given an scaled preconditioned LHS update \f$\breve{Z}\f$.
	 *
	 * @param  nativeLhsStep  [in] Multivector \f$\breve{Z}\f$ that defines the
	 *                        update \f$\bar{X} = \bar{X}_{0,r} + S_R P_R \breve{Z}\f$.
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
	/** Set the current unscaled unpreconditioned LHS \f$\bar{X}\f$ given the scaled preconditioned LHS \f$\breve{X}\f$.
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
	/** Get smart pointer to mutable current unscaled unpreconditioned residual \f$\bar{R}\f$.
	 *
	 * Preconditions:<ul>
	 * <li><tt>this->getCurrNumRhs() > 0</tt>
	 * </ul>
	 * 
	 * This function allows an interative solver to directly update the
	 * current residual if it is convenient to do so.  For example, a
	 * standard CG algorithm that uses just one preconditioner and no
	 * scaling will directly update the unscaled unpreconditioned
	 * residual \f$\bar{R}\f$ and therefore save on storage and
	 * computational cost.
	 *
	 * Note that if left scaling is present then direct access to this
	 * residual probably should not be sought.  When in doubt, a
	 * concrete <tt>Belos::BasicIteration</tt> subclass should not call
	 * this function and instead just let <tt>*this</tt> compute it when
	 * requested by the status test.
	 */
	virtual RefCountPtr<TSFCore::MultiVector<Scalar> > getCurrResidual() = 0;

	///
	/** Inform that current unscaled unpreconditioned residual \f$\bar{R}\f$ is up to date or not.
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
	 *
	 * This function should only be called if <tt>this->getCurrResidual()</tt>
	 * was called to get direct access to the unscaled unpreconditioned
	 * residual multi-vector for direct update by the iterative solver
	 * implementation.
	 */
	virtual void setCurrResidualComputed( const bool isCurrResidualComputed ) = 0;

	///
	/** \brief Compute the scaled preconditioned residual \f$\breve{R}\f$ consistent
	 * with the composite operator returned from this->getCompositeOperator().
	 *
	 * @param  currLhs  [in] Pointer to current LHS to use.  If <tt>currLhs==NULL</tt> then
	 *                  <tt>this->getCurrLhs()</tt> will be used in its place.
	 * @param  currRhs  [in] Pointer to current RHS to use.  If <tt>currRhs==NULL</tt> then
	 *                  <tt>this->getCurrRhs()</tt> will be used in its place.
	 * @param  currCompositeResidual
	 *                  [out] The preconditioned residual that must be been created
	 *                  using <tt>this->getCompositeOperator().range()->createMembers(this->getCurrBlockSize())</tt>.
	 *
	 * Preconditions:<ul>
	 * <li><tt>this->getCurrNumRhs() > 0</tt>
	 * <li>[<tt>currLhs==NULL</tt>] <tt>this->isCurrLhsUpdated() == true</tt>
	 * <li><tt>currCompositeResidual != NULL</tt>
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
	virtual void computeCurrCompositeResidual(
		const TSFCore::MultiVector<Scalar>    *currLhs
		,const TSFCore::MultiVector<Scalar>   *currRhs
		,TSFCore::MultiVector<Scalar>         *currCompositeResidual
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
	virtual void finalizeCurrSystem( const BasicIterationState<Scalar> &bis ) = 0;

	//@}

};

} // end Belos namespace

#endif // BELOS_LINEAR_PROBLEM_ITERATION_HPP
