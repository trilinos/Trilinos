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

#ifndef BELOS_LINEAR_PROBLEM_SETUP_HPP
#define BELOS_LINEAR_PROBLEM_SETUP_HPP

#include "Belos_LinearProblemIteration.hpp"

namespace Belos {

///
/** Interface of a linear problem being solved for access by a client
 * of an iterative solver.
 *
 * ToDo: Finish documentation!
 */
template<class Scalar>
class LinearProblemSetup : public LinearProblemIteration<Scalar> {
public:

	///
  using LinearProblemIteration<Scalar>::getOperator;

	/** @name Setup functions */
	//@{

	///
	/** Set operator, RHS and LHS.
	 *
	 * @param  operator  [in] Defines the operator \f$A\f$.
	 * @param  symmetry  [in] Defines the assumed symmetry of the operator \f$A\f$
	 * @param  rhs       [in] Defines the RHS multi-vector \f$B\f$.
	 * @param  lhs       [in] Defines the LHS multi-vector \f$X\f$.  On input, <tt>lhs</tt>
	 *                   must contain the initial guess of the solution (usually just zero).
	 *                   This multi-vector will contain be updated with the solution after
	 *                   linear solver is finished.
	 *
	 * Preconditions:<ul>
	 * <li><tt>operator.op().get() != NULL</tt>
	 * <li><tt>rhs.get() != NULL</tt>
	 * <li><tt>lhs.get() != NULL</tt>
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li><tt>this->getOperator().op().get() == operator.op().get()</tt>
	 * <li><tt>this->getOperatorSymmetry() == symmetry</tt>
	 * <li><tt>this->getRhs().get() == rhs.get()</tt>
	 * <li><tt>this->getLhs().get() == lhs.get()</tt>
	 * <li><tt>this->getTotalNumRhs() == this->getRhs().domain()->dim()</tt>
	 * <li><tt>this->getCurrNumRhs() == 0</tt>
	 * </ul>
	 */
	virtual void initialize(
		const TSFCore::LinearOpHandle<Scalar>                    &op
		,const EOpSymmetry                                       symmetry
		,const RefCountPtr<const TSFCore::MultiVector<Scalar> >  &rhs
		,const RefCountPtr<TSFCore::MultiVector<Scalar> >        &lhs
		) = 0;

	///
	/** Set the operator \f$A\f$ or \f$\tilde{A}\f$.
	 *
	 * @param  operator  [in] Defines the operator \f$A\f$ or \f$\tilde{A}\f$
	 * @param  symmetry  [in] Defines the assumed symmetry of the operator
	 *
	 * Preconditions:<ul>
	 * <li><tt>operator.op().get() != NULL</tt>
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li><tt>this->getOperator().op().get() == operator.op().get()</tt>
	 * <li><tt>this->getOperatorSymmetry() == symmetry</tt>
	 * <li><tt>this->getTotalNumRhs() == 0</tt>
	 * </ul>
	 */
	virtual void setOperator( const TSFCore::LinearOpHandle<Scalar> &op, const EOpSymmetry symmetry = OP_UNSYMMETRIC ) = 0;

	///
	/** Get persisting relationship with the operator \f$A\f$ or \f$\tilde{A}\f$.
	 */
	virtual TSFCore::LinearOpHandle<Scalar> getOperator() = 0;

	///
	/** Set the right preconditioner.
	 *
	 * @param  rightPrec  [in] Defines the right preconditioner \f$P_R\f$.
	 * @param  symmetry   [in] Defines the assumed symmetry of the operator.
	 *
	 * Preconditions:<ul>
	 * <li><tt>rightPrec.op().get() != NULL</tt>
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li><tt>const_cast<LinearProblemSetup<Scalar>*>(this)->getRightPrec().op().get() == rightPrec.op().get()</tt>
	 * <li><tt>const_cast<const LinearProblemSetup<Scalar>*>(this)->getRightPrec().op() == rightPrec.op().get()</tt>
	 * <li><tt>this->getRightPrecSymmetry() == symmetry</tt>
	 * <li><tt>this->getTotalNumRhs() == 0</tt>
	 * <li><tt>this->getCurrNumRhs() == 0</tt>
	 * </ul>
	 */
	virtual void setRightPrec( const TSFCore::LinearOpHandle<Scalar> &rightPrec, const EOpSymmetry symmetry = OP_UNSYMMETRIC ) = 0;

	///
	/** Get persisting relationship with the right preconditioner \f$P_R\f$.
	 */
	virtual TSFCore::LinearOpHandle<Scalar> getRightPrec() = 0;

	///
	/** Set the left preconditioner.
	 *
	 * @param  leftPrec  [in] Defines the left preconditioner \f$P_L\f$.
	 * @param  symmetry  [in] Defines the assumed symmetry of the operator
	 *
	 * Preconditions:<ul>
	 * <li><tt>leftPrec.op().get() != NULL</tt>
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li><tt>const_cast<LinearProblemSetup<Scalar>*>(this)->getLeftPrec().op().get() == rightPrec.op().get()</tt>
	 * <li><tt>const_cast<const LinearProblemSetup<Scalar>*>(this)->getLeftPrec().op() == rightPrec.op().get()</tt>
	 * <li><tt>this->getLeftPrecSymmetry() == symmetry</tt>
	 * <li><tt>this->getTotalNumRhs() == 0</tt>
	 * <li><tt>this->getCurrNumRhs() == 0</tt>
	 * </ul>
	 */
	virtual void setLeftPrec( const TSFCore::LinearOpHandle<Scalar> &leftPrec, const EOpSymmetry symmetry = OP_UNSYMMETRIC ) = 0;

	///
	/** Get persisting relationship with the left preconditioner \f$P_L\f$.
	 */
	virtual TSFCore::LinearOpHandle<Scalar> getLeftPrec() = 0;

	///
	/** Set the RHS multi-vector \f$B\f$.
	 *
	 * @param  rhs       [in] Defines the RHS multi-vector \f$B\f$.
	 *
	 * Preconditions:<ul>
	 * <li><tt>rhs.get() != NULL</tt>
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li><tt>this->getRhs().get() == rhs.get()</tt>
	 * <li><tt>this->getTotalNumRhs() == 0</tt>
	 * <li><tt>this->getCurrNumRhs() == 0</tt>
	 * </ul>
	 */
	virtual void setRhs( const RefCountPtr<const TSFCore::MultiVector<Scalar> > &rhs ) = 0;

	///
	/** Set the RHS multi-vector \f$B\f$.
	 */
	virtual RefCountPtr<const TSFCore::MultiVector<Scalar> > getRhs() = 0;

	///
	/** Set the LHS multi-vector \f$X\f$.
	 *
	 * @param  lhs       [in] Defines the LHS multi-vector \f$X\f$.  On input, <tt>lhs</tt>
	 *                   must contain the initial guess of the solution (usually just zero).
	 *                   This multi-vector will contain be updated with the solution after
	 *                   linear solver is finished.
	 *
	 * Preconditions:<ul>
	 * <li><tt>lhs.get() != NULL</tt>
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li><tt>this->getLhs().get() == lhs.get()</tt>
	 * <li><tt>this->getTotalNumRhs() == 0</tt>
	 * <li><tt>this->getCurrNumRhs() == 0</tt>
	 * </ul>
	 */
	virtual void setLhs( const RefCountPtr<TSFCore::MultiVector<Scalar> > &lhs ) = 0;

	///
	/** Set the LHS multi-vector \f$X\f$.
	 */
	virtual RefCountPtr<TSFCore::MultiVector<Scalar> > getLhs() = 0;

	///
	/** Set the status test.
	 *
	 * @param  statusTest  [in] Defines the status test.  It is allowed for <tt>statusTest.get()==NULL</tt>
	 *                     in which case the status test will be unset.
	 *
	 * Postconditions:<ul>
	 * <li><tt>getStatusTest.get() ==  statusTest.get()</tt>
	 * <li><tt>this->getTotalNumRhs() == 0</tt>
	 * <li><tt>this->getCurrNumRhs() == 0</tt>
	 * </ul>
	 */
	virtual void setStatusTest( const RefCountPtr<StatusTest<Scalar> > &statusTest ) = 0;

	///
	/** Complete setup of the linear problem
	 *
	 * Preconditions:<ul>
	 * <li><tt>this->getOperator().op().get() != NULL</tt>
	 * <li><tt>this->getRhs().get() != NULL</tt>
	 * <li><tt>this->getLhs().get() != NULL</tt>
	 * <li> Spaces for all opeators and multi-vectors must agree (ToDo: Make more explicit!)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li><tt>this->getTotalNumRhs() == this->getRhs().domain()->dim()</tt>
	 * <li><tt>this->getCurrNumRhs() == 0</tt>
	 * </ul>
	 */
	virtual void completeSetup() = 0;

	///
	/** Set to an uninitialized state and release smart pointers to all
	 * operators, multi-vectors, vectors and status tests that where
	 * set.
	 *
	 * Postconditions:<ul>
	 * <li><tt>this->getTotalNumRhs() == 0</tt>
	 * <li><tt>this->getCurrNumRhs() == 0</tt>
	 * </ul>
	 */
	virtual void uninitialize() = 0;

	//@}

};

} // end Belos namespace

#endif // BELOS_LINEAR_PROBLEM_SETUP_HPP
