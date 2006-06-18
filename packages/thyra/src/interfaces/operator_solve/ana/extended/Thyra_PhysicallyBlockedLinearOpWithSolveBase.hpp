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

#ifndef THYRA_PHYSICALLY_BLOCKED_LINEAR_OP_WITH_SOLVE_BASE_HPP
#define THYRA_PHYSICALLY_BLOCKED_LINEAR_OP_WITH_SOLVE_BASE_HPP

#include "Thyra_BlockedLinearOpWithSolveBase.hpp"

namespace Thyra {

/** \brief .
 *
 * ToDo: Finish documentation.
 *
 * \ingroup Thyra_Op_Vec_Interoperability_Extended_Interfaces_grp
 */
template<class RangeScalar, class DomainScalar=RangeScalar>
class PhysicallyBlockedLinearOpBase
  : virtual public BlockLinearOpBase<RangeScalar,DomainScalar>
{
public:

  /** \brief Begin a block fill where the product range and domain spaces will
   * be created on the fly.
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>this->blockFillIsActive()==true</tt>
   * <li><tt>this->productRange().get()==NULL</tt>
   * <li><tt>this->productDomain().get()==NULL</tt>
   * </ul>
   */
  virtual void beginBlockFill() = 0;

  /** \brief Begin a block fill where the product range and domain spaces
   * are set a priori.
   *
   * \param  productRange
   *           [in] The product space to use of the range space.
   * \param  productRange
   *           [in] The product space to use of the domain space.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>productRange.get()!=NULL && productRange->dim() > 0</tt>
   * <li><tt>productDomain.get()!=NULL && productDomain->dim() > 0</tt>
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>this->blockFillIsActive()==true</tt>
   * <li><tt>this->productRange().get()==productRange.get()</tt>
   * <li><tt>this->productDomain().get()==productDomain.get()</tt>
   * </ul>
   */
  virtual void beginBlockFill(
    const Teuchos::RefCountPtr<const ProductVectorSpaceBase<Scalar> >  &productRange
    ,const Teuchos::RefCountPtr<const ProductVectorSpaceBase<Scalar> > &productDomain
    ) = 0;
  
  /** \brief Determines if a block fill is active or not . */
  virtual bool blockFillIsActive() const = 0;

  /** \brief Determines if the block <tt>(i,j)</tt> can be filled or not.
   *
   * \param  i  [in] Zero-based index for the block row.
   * \param  j  [in] Zero-based index for the block column.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>i >= 0 && j >= 0</tt>
   * <li>[<tt>this->productRange().get()!=NULL</tt>]
   *       <tt>i < this->productRange()->numBlocks()</tt>
   * <li>[<tt>this->productDomain().get()!=NULL</tt>]
   *       <tt>j < this->productDomain()->numBlocks()</tt>
   * </ul>
   */
  virtual bool acceptsBlock(const int i, const int j) const = 0;
  
  /** \brief Set a non-const block linear operator.
   *
   * \param  i  [in] Zero-based index for the block row.
   * \param  j  [in] Zero-based index for the block column.
   * \param  block
   *            [in] The block operator being set.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>this->acceptsBlock(i,j)==true</tt>
   * </ul>
   */
  virtual void setBlock(
    const int i, const int j
    ,const Teuchos::RefCountPtr<LinearOpBase<Scalar> > &block
    ) = 0;
  
  /** \brief Set a const block linear operator.
   *
   * \param  i  [in] Zero-based index for the block row.
   * \param  j  [in] Zero-based index for the block column.
   * \param  block
   *            [in] The block operator being set.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>this->acceptsBlock(i,j)==true</tt>
   * </ul>
   */
  virtual void setBlock(
    const int i, const int j
    ,const Teuchos::RefCountPtr<const LinearOpBase<Scalar> > &block
    ) = 0;
  
  /** \brief End a block fill after which <tt>*this</tt> object can be used.
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>this->blockFillIsActive()==false</tt>
   * <li><tt>this->productRange().get()!=NULL</tt>
   * <li><tt>this->productDomain().get()!=NULL</tt>
   * </ul>
   */
  virtual void endBlockFill() = 0;

};

} // namespace Thyra

#endif // THYRA_PHYSICALLY_BLOCKED_LINEAR_OP_WITH_SOLVE_BASE_HPP
