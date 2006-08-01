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
#include "Thyra_PhysicallyBlockedLinearOpBase.hpp"

namespace Thyra {

/** \brief Base interface for filling an implicit
 * <tt>LinearOpWithSolveBase</tt> object as a set of
 * <tt>LinearOpWithSolveBase<tt> and <tt>>LinearOpBase</tt> blocks.
 *
 * ToDo: Finish documentation.
 *
 * \ingroup Thyra_Op_Solve_Interoperability_Extended_Interfaces_grp
 */
template<class RangeScalar, class DomainScalar=RangeScalar>
class PhysicallyBlockedLinearOpWithSolveBase
  : virtual public BlockedLinearOpWithSolveBase<RangeScalar,DomainScalar>
  , virtual public PhysicallyBlockedLinearOpBase<RangeScalar,DomainScalar>
{
public:

  /** \brief Determines if the block <tt>(i,j)</tt> can be filled with a
   * <tt>LinearOpWithSolveBase</tt> object or not.
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
  virtual bool acceptsBlockLOWSB(const int i, const int j) const = 0;
  
  /** \brief Set a non-const block <tt>LinearOpWithSolveBase</tt> object.
   *
   * \param  i  [in] Zero-based index for the block row.
   * \param  j  [in] Zero-based index for the block column.
   * \param  block
   *            [in] The block operator being set.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>this->acceptsBlockLOWS(i,j)==true</tt>
   * </ul>
   */
  virtual void setNonconstBlockLOWS(
    const int i, const int j
    ,const Teuchos::RefCountPtr<LinearOpWithSolveBase<Scalar> > &block
    ) = 0;
  
  /** \brief Set a const block <tt>LinearOpWithSolveBase</tt> object.
   *
   * \param  i  [in] Zero-based index for the block row.
   * \param  j  [in] Zero-based index for the block column.
   * \param  block
   *            [in] The block operator being set.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>this->acceptsBlockLOWS(i,j)==true</tt>
   * </ul>
   */
  virtual void setBlockLOWS(
    const int i, const int j
    ,const Teuchos::RefCountPtr<const LinearOpWithSolveBase<Scalar> > &block
    ) = 0;

};

} // namespace Thyra

#endif // THYRA_PHYSICALLY_BLOCKED_LINEAR_OP_WITH_SOLVE_BASE_HPP
