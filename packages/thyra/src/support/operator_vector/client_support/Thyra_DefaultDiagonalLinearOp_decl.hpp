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

#include "Thyra_DiagonalLinearOpBase.hpp"
#include "Teuchos_ConstNonconstObjectContainer.hpp"


namespace Thyra {


/** \brief Default concrete <tt>LinearOpBase</tt> subclass for diagonal linear
 * operators.
 *
 * This class represents a diagonal linear operator <tt>M</tt> of the form:

 \verbatim

 M = diag(diag)

 \endverbatim
 
 * where <tt>diag</tt> is a <tt>VectorBase</tt> object.
 *
 * The defined operator implements <tt>this->apply()</tt> as follows:
 
 \verbatim

 y = alpha*op(M)*x + beta*y
 
 =>

 y(i) = alpha*diag(i)*x(i) + beta*y(i), for i = 0 ... n-1

 \endverbatim
 
 * where <tt>n = this->domain()->dim()</tt>.
 *
 * \ingroup Thyra_Op_Vec_ANA_Development_grp
 */
template<class Scalar>
class DefaultDiagonalLinearOp : virtual public DiagonalLinearOpBase<Scalar>
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
  DefaultDiagonalLinearOp();

  /** \brief Calls <tt>initialize()</tt> to construct given a vector space.
   */
  DefaultDiagonalLinearOp(
    const RCP<const VectorSpaceBase<Scalar> >  &space
    );

  /** \brief Calls <tt>initialize()</tt> to construct for a non-const diagonal
   * vector.
   */
  DefaultDiagonalLinearOp(
    const RCP<VectorBase<Scalar> >   &diag
    );

  /** \brief Calls <tt>initialize()</tt> to construct for a const diagonal
   * vector.
   */
  DefaultDiagonalLinearOp(
    const RCP<const VectorBase<Scalar> >   &diag
    );

  /** \brief Initialize given a vector space which allocates a vector internally.
   *
   * \param  space   [in] Smart pointer to vector space
   *
   * Preconditions:<ul>
   * <li><tt>space.get()!=NULL</tt>
   * </ul>
   *
   * Postconditions:<ul>
   * <li><tt>this->getNonconstDiag()->space()->isCompatible(*space)==true</tt>
   * <li><tt>this->getDiag()->space()->isCompatible(*space)==true</tt>
   * <li><tt>this->this->domain().get() == space.get()</tt>
   * <li><tt>this->this->range().get() == space.get()</tt>
   * </ul>
   */
  void initialize(
    const RCP<const VectorSpaceBase<Scalar> >  &space
    );

  /** \brief Initialize given a non-const diagonal vector.
   *
   * \param diag [in] Smart pointer to diagonal vector.
   *
   * Preconditions:<ul>
   * <li><tt>diag.get()!=NULL</tt>
   * </ul>
   *
   * Postconditions:<ul>
   * <li><tt>this->getNonconstDiag().get()==diag.get()</tt>
   * <li><tt>this->getDiag().get()==diag.get()</tt>
   * <li><tt>this->this->domain().get() == diag->space().get()</tt>
   * <li><tt>this->this->range().get() == diag->space().get()</tt>
   * </ul>
   */
  void initialize(const RCP<VectorBase<Scalar> > &diag);

  /** \brief Initialize given a const diagonal vector.
   *
   * \param  diag   [in] Smart pointer to diagonal vector. 
   *
   * Preconditions:<ul>
   * <li><tt>diag.get()!=NULL</tt>
   * </ul>
   *
   * Postconditions:<ul>
   * <li><tt>this->getNonconstDiag()</tt> with throw an exception if called
   * <li><tt>this->getDiag().get()==diag.get()</tt>
   * <li><tt>this->this->domain().get() == diag->space().get()</tt>
   * <li><tt>this->this->range().get() == diag->space().get()</tt>
   * </ul>
   */
  void initialize(
    const RCP<const VectorBase<Scalar> >   &diag
    );

  /** \brief Uninitialize.
   *
   * Postconditions:<ul>
   * <li><tt>this->getNonconstDiag().get()==NULL</tt>
   * <li><tt>this->getDiag().get()==NULL</tt>
   * <li><tt>this->this->domain().get()==NULL</tt>
   * <li><tt>this->this->range().get()==NULL</tt>
   * </ul>
   *
   * <b>Note:</b> If the client wants to hold on to the underlying wrapped
   * diagonal vector then they had better grab it using
   * <tt>this->getDiag()</tt> or <tt>this->getNonconstDiag()</tt> before
   * calling this function!
   */
  void uninitialize();

  //@}

  /** @name Overridden from DiagonalLinearOpBase */
  //@{

  /** \brief . */
  bool isDiagConst() const;
  /** \brief . */
  RCP<VectorBase<Scalar> > getNonconstDiag();
  /** \brief . */
  RCP<const VectorBase<Scalar> > getDiag() const;

  //@}

  /** @name Overridden from LinearOpBase */
  //@{
  /** \brief Returns <tt>this->getDiag()->space()</tt>.
   *
   * Preconditions:<ul>
   * <li><tt>this->getDiag().get()!=NULL</tt>
   * </ul>
   */
  RCP< const VectorSpaceBase<Scalar> > range() const;
  /** \brief Returns <tt>this->getDiag()->space()</tt>.
   *
   * Preconditions:<ul>
   * <li><tt>this->getDiag().get()!=NULL</tt>
   * </ul>
   */
  RCP< const VectorSpaceBase<Scalar> > domain() const;
  /** \brief . */
  RCP<const LinearOpBase<Scalar> > clone() const;
  //@}

protected:

  /** @name Protected functions overridden from LinearOpBase */
  //@{
  /** \brief . */
  bool opSupportedImpl(EOpTransp M_trans) const;
  /** \brief . */
  void applyImpl(
    const EOpTransp M_trans,
    const MultiVectorBase<Scalar> &X,
    const Ptr<MultiVectorBase<Scalar> > &Y,
    const Scalar alpha,
    const Scalar beta
    ) const;
  //@}

private:

  Teuchos::ConstNonconstObjectContainer<VectorBase<Scalar> > diag_;

};


/** \brief Nonmember constructor function.
 *
 * \relates DefaultDiagonalLinearOp
 */
template<class Scalar>
RCP<const LinearOpBase<Scalar> >
diagonal(
  const RCP<VectorBase<Scalar> > &diag,
  const std::string &label = ""
  )
{
  RCP<LinearOpBase<Scalar> > dlo =
    Teuchos::rcp(new DefaultDiagonalLinearOp<Scalar>(diag));
  if (label.length())
    dlo->setObjectLabel(label);
  return dlo;
}


}	// end namespace Thyra


#endif	// THYRA_DIAGONAL_LINEAR_OP_DECL_HPP
