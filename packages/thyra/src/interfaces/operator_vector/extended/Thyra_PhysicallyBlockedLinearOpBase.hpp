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

#ifndef THYRA_PHYSICALLY_BLOCKED_LINEAR_OP_BASE_HPP
#define THYRA_PHYSICALLY_BLOCKED_LINEAR_OP_BASE_HPP

#include "Thyra_BlockedLinearOpBase.hpp"


namespace Thyra {


/** \brief Base interface for physically blocked linear operators.
 *
 * This interface allows clients to fill a blocked linear operator and create
 * blocked objects with various numbers of row and column blocks.  There are
 * two modes to fill a blocked linear operator represented by the two forms of
 * the <tt>beginBlockFill()</tt> function.
 *
 * ToDo: Finish documentation.
 *
 * \ingroup Thyra_Op_Vec_extended_interfaces_code_grp
 */
template<class Scalar>
class PhysicallyBlockedLinearOpBase
  : virtual public BlockedLinearOpBase<Scalar>
{
public:

  /** \brief Begin a block fill where the product range and domain spaces will
   * be created on the fly and the number of block rows and columns is not
   * known in advance.
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>this->blockFillIsActive()==true</tt>
   * <li><tt>this->productRange().get()==NULL</tt>
   * <li><tt>this->productDomain().get()==NULL</tt>
   * </ul>
   */
  virtual void beginBlockFill() = 0;

  /** \brief Begin a block fill where the product range and domain spaces will
   * be created on the fly but the total number of block rows and block
   * columns is known in advance.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>numRowBlocks > 0</tt>
   * <li><tt>numColBlocks > 0</tt>
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>this->blockFillIsActive()==true</tt>
   * <li><tt>this->productRange().get()==NULL</tt>
   * <li><tt>this->productDomain().get()==NULL</tt>
   * </ul>
   */
  virtual void beginBlockFill(
    const int numRowBlocks, const int numColBlocks
    ) = 0;

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
    const Teuchos::RCP<const ProductVectorSpaceBase<Scalar> >    &productRange
    ,const Teuchos::RCP<const ProductVectorSpaceBase<Scalar> >  &productDomain
    ) = 0;
  
  /** \brief Determines if a block fill is active or not . */
  virtual bool blockFillIsActive() const = 0;

  /** \brief Determines if the block <tt>(i,j)</tt> can be filled or not.
   *
   * \param  i  [in] Zero-based index for the block row.
   * \param  j  [in] Zero-based index for the block column.
   *
   * <b>Preconditions:</b><ul>
   * <li>this->blockFillIsActive()==true</tt>
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
  virtual void setNonconstBlock(
    const int i, const int j
    ,const Teuchos::RCP<LinearOpBase<Scalar> > &block
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
    ,const Teuchos::RCP<const LinearOpBase<Scalar> > &block
    ) = 0;
  
  /** \brief End a block fill after which <tt>*this</tt> object can be used.
   *
   * If a valid linear operator object can not be formed from what was set
   * then a ??? exception will be thrown.
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>this->blockFillIsActive()==false</tt>
   * <li><tt>this->productRange().get()!=NULL</tt>
   * <li><tt>this->productDomain().get()!=NULL</tt>
   * </ul>
   */
  virtual void endBlockFill() = 0;
  
  /** \brief Set to uninitlaized.
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>this->blockFillIsActive()==false</tt>
   * <li><tt>this->productRange().get()==NULL</tt>
   * <li><tt>this->productDomain().get()==NULL</tt>
   * </ul>
   */
  virtual void uninitialize() = 0;

private:
  
  // Not defined and not to be called
  PhysicallyBlockedLinearOpBase<Scalar>&
  operator=(const PhysicallyBlockedLinearOpBase<Scalar>&);

};


} // namespace Thyra


#endif // THYRA_PHYSICALLY_BLOCKED_LINEAR_OP_BASE_HPP
