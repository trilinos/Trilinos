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

#ifndef BELOS_LINEAR_PROBLEM_STATE_HPP
#define BELOS_LINEAR_PROBLEM_STATE_HPP

#include "Belos_Types.hpp"
#include "TSFCoreLinOp.hpp"

namespace Belos {

///
/** Base interface of a linear problem being solved for access by status tests.
 *
 * ToDo: Finish documentation!
 */
template<class Scalar>
class LinearProblemState {
public:
	
	///
	virtual ~LinearProblemState() {}

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
	virtual TSFCore::LinOpNonPersisting<Scalar> getOperator() const = 0;
	
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
	virtual TSFCore::LinOpNonPersisting<Scalar> getRightPrec() const = 0;
	
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
	virtual TSFCore::LinOpNonPersisting<Scalar> getLeftPrec() const = 0;
	
	///
	/** Return if the left preconditioner \f$P_L\f$ is symmetric or not.
	 *
	 * Preconditions:<ul>
	 * <li><tt>this->getLeftPrec() != NULL</tt>
	 * </ul>
	 */
	virtual EOpSymmetry getLeftPrecSymmetry() const = 0;

	///
	/** Return non-persisting reference to the combined operator
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
	virtual TSFCore::LinOpNonPersisting<Scalar> getCombinedOperator() const = 0;

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

	/** @name Specific to current block of RHSs */
	//@{

	///
	/** Return the current number of RHSs being solved in the current block \f$\bar{m}\f$.
	 *
	 * Preconditions:<ul>
	 * <li><tt>this->getTotalNumRhs() > 0</tt>
	 * </ul>
	 *
	 * A return value of <tt>return==0</tt> is a flag that a current
	 * block of linear systems is not active.
	 */
	virtual int getCurrNumRhs() const = 0;

	///
	/** Return the current block size \f$\bar{b}\f$.
	 *
	 * Preconditions:<ul>
	 * <li><tt>this->getCurrNumRhs() > 0</tt>
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li><tt>return == getCurrRhs().domain()->dim()</tt>
	 * <li><tt>return >= this->getCurrNumRhs()</tt>
	 * </ul>
	 */
	virtual int getCurrBlockSize() const = 0;

	///
	/** Return the current set of RHS indexes.
	 *
	 * @param  currNumRhs      [in] The current number of RHSs in the current block of systems.
	 * @param  currRhsIndexes  [out] Array (length <tt>currNumRhs</tt>) that on output returns
	 *                         the indexes of the current RHS.  Specifically, the column vectors:
	 *                         <tt>*this->getRhs().col(currRhsIndexes[k])</tt> are equal to 
	 *                         <tt>this->getCurrRhs().col(k+1)</tt>, for <tt>k=0...currNumRhs-1</tt>.
	 *
	 * Preconditions:<ul>
	 * <li><tt>this->getCurrNumRhs() > 0</tt>
	 * <li><tt>currNumRhs == this->getCurrNumRhs()</tt>
	 * </ul>
	 *
	 * Note that the returned array <tt>currRhsIndexes</tt> does not
	 * contain entries for augmented RHSs since there is no relation to
	 * the original RHSs.
	 */
	virtual void getCurrRhsIndexes(
		const int              currNumRhs
		,int                   currRhsIndexes[]
		) const = 0;

	///
	/** Return <tt>const</tt> reference to the current RHS \f$\bar{B}\f$.
	 *
	 * Preconditions:<ul>
	 * <li><tt>this->getCurrNumRhs() > 0</tt>
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li><tt>return.range()->isCompatible(*this->getOperator().range())==true</tt>
	 * <li><tt>return.domain()->dim()==this->getCurrBlockSize()</tt>
	 * </ul>
	 */
	virtual const TSFCore::MultiVector<Scalar>& getCurrRhs() const = 0;

	///
	/** Return <tt>const</tt> reference to the current residual for the
	 * starting initial guess \f$\bar{R}_0 = A \bar{X}_0 - \bar{B}\f$.
	 *
	 * Preconditions:<ul>
	 * <li><tt>this->getCurrNumRhs() > 0</tt>
	 * </ul>
	 *
	 * If <tt>this->isCurrInitResidualComputed()==false</tt> before this
	 * function is called then this residual will be computed during
	 * this function call.
	 *
	 * Postconditions:<ul>
	 * <li><tt>this->isCurrInitResidualComputed()==true</tt>
	 * <li><tt>return.range()->isCompatible(*this->getOperator().range())==true</tt>
	 * <li><tt>return.domain()->dim()==this->getCurrBlockSize()</tt>
	 * </ul>
	 */
	virtual const TSFCore::MultiVector<Scalar>& getCurrInitResidual() const = 0;

	///
	/** Return if current residual for the starting initial guess \f$\bar{R}_0\f$ is computed yet or not.
	 *
	 * Preconditions:<ul>
	 * <li><tt>this->getCurrNumRhs() > 0</tt>
	 * </ul>
	 */
	virtual bool isCurrInitResidualComputed() const = 0;

	///
	/** Return <tt>const</tt> reference to the current initial guess for the current restart \f$\bar{X}_{0,r}\f$.
	 *
	 * Preconditions:<ul>
	 * <li><tt>this->getCurrNumRhs() > 0</tt>
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li><tt>return.range()->isCompatible(*this->getOperator().domain())==true</tt>
	 * <li><tt>return.domain()->dim()==this->getCurrBlockSize()</tt>
	 * </ul>
	 */
	virtual const TSFCore::MultiVector<Scalar>& getCurrInitLhs() const = 0;

	///
	/** Return <tt>const</tt> reference to the current LHS \f$\bar{X}\f$.
	 *
	 * Preconditons:<ul>
	 * <li><tt>this->isCurrLhsUpdated() == true</tt>.
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li><tt>return.range()->isCompatible(*this->getOperator().domain())==true</tt>
	 * <li><tt>return.domain()->dim()==this->getCurrBlockSize()</tt>
	 * </ul>
	 */
	virtual const TSFCore::MultiVector<Scalar>& getCurrLhs() const = 0;

	///
	/** Return if the current LHS \f$\bar{X}\f$ is updated or not.
	 *
	 * Preconditions:<ul>
	 * <li><tt>this->getCurrNumRhs() > 0</tt>
	 * </ul>
	 */
	virtual bool isCurrLhsUpdated() const = 0;

	///
	/** Return <tt>const</tt> reference to the current residual \f$\bar{R} = A \bar{X} - \bar{B}\f$.
	 *
	 * Preconditions:<ul>
	 * <li>[<tt>this->isCurrResidualComputed()==false</tt>] <tt>this->isCurrLhsUpdated()==true</tt>
	 * </ul>
	 *
	 * If <tt>this->isCurrResidualComputed()==false</tt> before this
	 * function is called then this residual will be computed during
	 * this function call.
	 *
	 * Postconditions:<ul>
	 * <li><tt>this->isCurrResidualComputed()==true</tt>
	 * <li><tt>return.range()->isCompatible(*this->getOperator().range())==true</tt>
	 * <li><tt>return.domain()->dim()==this->getCurrBlockSize()</tt>
	 * </ul>
	 */
	virtual const TSFCore::MultiVector<Scalar>& getCurrResidual() const = 0;

	///
	/** Return if current residual \f$\bar{R}\f$ is computed yet or not.
	 *
	 * Preconditions:<ul>
	 * <li><tt>this->getCurrNumRhs() > 0</tt>
	 * </ul>
	 */
	virtual bool isCurrResidualComputed() const = 0;

	// ToDo: Fill in the rest of the operations

	//@}

};

} // end Belos namespace

#endif // BELOS_LINEAR_PROBLEM_STATE_HPP
