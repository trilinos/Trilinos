// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
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

#ifndef THYRA_DIAGONAL_LINEAR_OP_DECL_HPP
#define THYRA_DIAGONAL_LINEAR_OP_DECL_HPP

#include "Thyra_LinearOpWithSolveBaseDecl.hpp"
#include "Thyra_SingleRhsLinearOpWithSolveBaseDecl.hpp"

namespace Thyra {

/** \brief Concrete <tt>LinearOpBase</tt> subclass for diagonal linear
 * operators.
 *
 * This class represents a diagonal linear operator <tt>M</tt> of the form:
 \verbatim

 M = gamma*diag(diag)
 \endverbatim
 *
 * where <tt>diag</tt> is a <tt>VectorBase</tt> object and <tt>gamma</tt>
 * is a <tt>Scalar</tt>.
 *
 * The defined operator implements <tt>this->apply()</tt> as follows:
 *
 \verbatim

 y = alpha*op(M)*x + beta*y
 
 =>

 y(i) = (alpha*gamma)*diag(i)*x(i) + beta*y(i), for i = 1 ... n
 \endverbatim
 *
 * where <tt>n = this->domain()->dim()</tt>.
 *
 * The defined operator implements <tt>this->solve()</tt> as follows:
 *
 \verbatim

 x = inv(op(M))*b
 
 =>

 x(i) = (1/gamma)*b(i)/diag(i), for i = 1 ... n
 \endverbatim
 *
 * where <tt>n = this->domain()->dim()</tt>.
 *
 * That is all there is to this subclass.
 *
 * \ingroup Thyra_Op_Vec_ANA_Development_grp
 */
template<class Scalar>
class DiagonalLinearOp
  : virtual public LinearOpWithSolveBase<Scalar>
  , virtual protected SingleRhsLinearOpWithSolveBase<Scalar>
{
public:

  /** \brief . */
  SingleRhsLinearOpBase<Scalar>::apply;

  /** @name Constructors/initializers/accessors */
  //@{

  /** \brief Constructs to uninitialized.
   *
   * Postconditions:<ul>
   * <li><tt>this->getDiag().get()==NULL</tt>
   * </ul>
   */
  DiagonalLinearOp();

  /// Calls <tt>initialize()</tt>
  DiagonalLinearOp(
    const Teuchos::RefCountPtr<const VectorBase<Scalar> >   &diag
    ,const Scalar                                           &gamma = Teuchos::ScalarTraits<Scalar>::one()
    );

  /** \brief Initialize given the diagonal.
   *
   * @param  diag   [in] Smart pointer to diagonal vector. 
   * @param  gamma  [in] Scalar multiplier.
   *
   * Preconditions:<ul>
   * <li><tt>diag.get()!=NULL</tt>
   * </ul>
   *
   * Postconditions:<ul>
   * <li><tt>this->diag().get()==diag.get()</tt>
   * <li><tt>this->gamma()==gamma</tt>
   * <li><tt>this->this->domain().get() == diag->space().get()</tt>
   * <li><tt>this->this->range().get() == diag->space().get()</tt>
   * </ul>
   */
  void initialize(
    const Teuchos::RefCountPtr<const VectorBase<Scalar> >   &diag
    ,const Scalar                                           &gamma = Teuchos::ScalarTraits<Scalar>::one()
    );

  /** \brief Returns the diagonal vector <tt>diag</tt>.
   *
   * A return value of <tt>return.get()==NULL</tt> indicates that
   * <tt>this</tt> is not fully initialized.
   */
  Teuchos::RefCountPtr<const VectorBase<Scalar> > diag() const;

  /** \brief Returns the scalar multiplier <tt>gamma</tt>.
   */
  Scalar gamma() const;

  /** \brief Set to uninitialized.
   *
   * @param  diag  [out] Optional pointer to smart pointer for diagonal.
   *               If <tt>diag!=NULL</tt> then on output <tt>*diag</tt>
   *               is set to <tt>this->diag()</tt> (before call).
   * @param  gamma [out] Optional pointer to scalar <tt>gamma</tt>.
   *               If <tt>gamma!=NULL</tt> then on output <tt>*gamma</tt>
   *               is set to <tt>this->gamma()</tt> (before call).
   *
   * Postconditions:<ul>
   * <li><tt>this->getDiag().get()==NULL</tt>
   * </ul>
   */
  void uninitialize(
    Teuchos::RefCountPtr<const VectorBase<Scalar> >  *diag   = NULL
    ,Scalar                                          *gamma  = NULL
    );

  /** @name Overridden from LinearOpBase */
  //@{
  /** \brief Returns <tt>this->getDiag()->space()</tt>.
   *
   * Preconditions:<ul>
   * <li><tt>this->getDiag().get()!=NULL</tt>
   * </ul>
   */
  Teuchos::RefCountPtr< const VectorSpaceBase<Scalar> > range() const;
  /** \brief Returns <tt>this->getDiag()->space()</tt>.
   *
   * Preconditions:<ul>
   * <li><tt>this->getDiag().get()!=NULL</tt>
   * </ul>
   */
  Teuchos::RefCountPtr< const VectorSpaceBase<Scalar> > domain() const;
  /** \brief . */
  Teuchos::RefCountPtr<const LinearOpBase<Scalar> > clone() const;
  //@}

protected:

  /** @name Overridden from SingleScalarLinearOpBase */
  //@{
  /** \brief . */
  bool opSupported(ETransp M_trans) const;
  //@}

  /** @name Overridden from SingleRhsLinearOpBase */
  //@{
  /** \brief . */
  void apply(
    const ETransp                M_trans
    ,const VectorBase<Scalar>    &x
    ,VectorBase<Scalar>          *y
    ,const Scalar                alpha
    ,const Scalar                beta
    ) const;
  //@}

  /** @name Overridden from SingleScalarLinearOpWithSolveBase */
  //@{
  /** \brief . */
  bool solveSupportsTrans(ETransp M_trans) const;
  /** \brief . */
  bool solveSupportsSolveTolType(ETransp M_trans, ESolveTolType solveTolType) const;
  //@}

  /** @name Overridden from SingleRhsLinearOpWithSolveBase */
  //@{
  /** \brief . */
  SolveStatus<Scalar> solve(
    const ETransp                         M_trans
    ,const VectorBase<Scalar>             &b
    ,VectorBase<Scalar>                   *x
    ,const SolveCriteria<Scalar>          *solveCriteria
    ) const;
  //@}

private:

  Teuchos::RefCountPtr< const VectorBase<Scalar> > diag_;
  Scalar                                           gamma_;

  void assertInitialized() const;

};

// /////////////////////////////////
// Inline members

template<class Scalar>
inline
Teuchos::RefCountPtr<const VectorBase<Scalar> >  
DiagonalLinearOp<Scalar>::diag() const
{
  return diag_;
}

template<class Scalar>
inline
Scalar DiagonalLinearOp<Scalar>::gamma() const
{
  return gamma_;
}

template<class Scalar>
inline
void DiagonalLinearOp<Scalar>::assertInitialized() const
{
#ifdef _DEBUG
  TEST_FOR_EXCEPT( diag_.get() == NULL );
#endif
}

}	// end namespace Thyra

#endif	// THYRA_DIAGONAL_LINEAR_OP_DECL_HPP
