// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_DIAGONAL_LINEAR_OP_WITH_SOLVE_DECL_HPP
#define THYRA_DIAGONAL_LINEAR_OP_WITH_SOLVE_DECL_HPP

#include "Thyra_DefaultDiagonalLinearOp.hpp"
#include "Thyra_LinearOpWithSolveBase.hpp"


namespace Thyra {


/** \brief Concrete <tt>LinearOpWithSolveBase</tt> subclass for diagonal linear
 * operators.
 *
 * This class represents a diagonal linear operator <tt>M</tt> of the form:
 \verbatim

 M = diag(diag)
 \endverbatim
 *
 * where <tt>diag</tt> is a <tt>VectorBase</tt> object.
 *
 * The defined operator implements <tt>this->apply()</tt> as follows:
 *
 \verbatim

 y = alpha*op(M)*x + beta*y
 
 =>

 y(i) = alpha*diag(i)*x(i) + beta*y(i), for i = 0 ... n-1
 \endverbatim
 *
 * where <tt>n = this->domain()->dim()</tt>.
 *
 * The defined operator implements <tt>this->solve()</tt> as follows:
 *
 \verbatim

 x = inv(op(M))*b
 
 =>

 x(i) = b(i)/diag(i), for i = 0 ... n-1
 \endverbatim
 *
 * where <tt>n = this->domain()->dim()</tt>.
 *
 * That is all there is to this subclass.
 */
template<class Scalar>
class DefaultDiagonalLinearOpWithSolve
  : virtual public DefaultDiagonalLinearOp<Scalar>,
    virtual public LinearOpWithSolveBase<Scalar>
{
public:

  /** @name Constructors/initializers/accessors */
  //@{

  /** \brief Constructs to uninitialized.
   *
   * Postconditions:<ul>
   * <li><tt>this->getDiag().get()==NULL</tt>
   * </ul>
   */
  DefaultDiagonalLinearOpWithSolve();

  /// Calls <tt>initialize()</tt>
  DefaultDiagonalLinearOpWithSolve(
    const RCP<const VectorBase<Scalar> >   &diag
    );

protected:

  /** @name Overridden from LinearOpWithSolveBase */
  //@{
  /** \brief . */
  bool solveSupportsImpl(EOpTransp M_trans) const;
  /** \brief . */
  bool solveSupportsSolveMeasureTypeImpl(
    EOpTransp M_trans, const SolveMeasureType& solveMeasureType) const;
  /** \brief . */
  SolveStatus<Scalar> solveImpl(
    const EOpTransp transp,
    const MultiVectorBase<Scalar> &B,
    const Ptr<MultiVectorBase<Scalar> > &X,
    const Ptr<const SolveCriteria<Scalar> > solveCriteria
    ) const;
  //@}

};


/** \brief Nonmember constructor.
 *
 * \relates DefaultDiagonalLinearOpWithSolve
 */
template<class Scalar>
RCP<DefaultDiagonalLinearOpWithSolve<Scalar> >
defaultDiagonalLinearOpWithSolve()
{
  return Teuchos::rcp(new DefaultDiagonalLinearOpWithSolve<Scalar>);
}


/** \brief Nonmember constructor.
 *
 * \relates DefaultDiagonalLinearOpWithSolve
 */
template<class Scalar>
RCP<DefaultDiagonalLinearOpWithSolve<Scalar> >
defaultDiagonalLinearOpWithSolve(
  const RCP<const VectorBase<Scalar> >   &diag
  )
{
  RCP<DefaultDiagonalLinearOpWithSolve<Scalar> > ddlows =
    defaultDiagonalLinearOpWithSolve<Scalar>();
  ddlows->initialize(diag);
  return ddlows;
}



}	// end namespace Thyra


#endif	// THYRA_DIAGONAL_LINEAR_OP_WITH_SOLVE_DECL_HPP
