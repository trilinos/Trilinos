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

#ifndef THYRA_BLOCKED_LINEAR_OP_WITH_SOLVE_BASE_HPP
#define THYRA_BLOCKED_LINEAR_OP_WITH_SOLVE_BASE_HPP

#include "Thyra_LinearOpWithSolveBase.hpp"
#include "Thyra_BlockedLinearOpBase.hpp"

namespace Thyra {

/** \brief Base interface for <tt>LinearOpWithSolveBase</tt> objects that can
 * be accessed in blocks.
 *
 * ToDo: Finish Documentation.
 *
 * \ingroup Thyra_Op_Vec_Interoperability_Extended_Interfaces_grp
 */
template<class RangeScalar, class DomainScalar=RangeScalar>
class BlockedLinearOpWithSolveBase
  : virtual public LinearOpWithSolveBase<RangeScalar,DomainScalar>
  , virtual public BlockedLinearOpBase<RangeScalar,DomainScalar>
{
public:

  /** \brief Return a non-const view of the block <tt>(i,j)</tt> if it exists.
   *
   * \param  i  [in] Zero-based index for the block row.
   * \param  j  [in] Zero-based index for the block column.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>this->blockIsConst(i,j)==false</tt>
   * <li><tt>0 <= i && i < this->productRange()->numBlocks()</tt>
   * <li><tt>0 <= j && j < this->productDomain()->numBlocks()</tt>
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   * <li>[<tt>this->blockExists(i,j)==true</tt>] <tt>return.get()!=NULL</tt>
   * <li>[<tt>this->blockExists(i,j)==false</tt>] <tt>return.get()==NULL</tt>
   * </ul>
   */
  virtual Teuchos::RefCountPtr<LinearOpWithSolveBase<RangeScalar,DomainScalar> >
  getBlockLOWS(const int i, const int j) = 0; 

  /** \brief Return a const view of the block <tt>(i,j)</tt> if it exists.
   *
   * \param  i  [in] Zero-based index for the block row.
   * \param  j  [in] Zero-based index for the block column.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>0 <= i && i < this->productRange()->numBlocks()</tt>
   * <li><tt>0 <= j && j < this->productDomain()->numBlocks()</tt>
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   * <li>[<tt>this->blockExists(i,j)==true</tt>] <tt>return.get()!=NULL</tt>
   * <li>[<tt>this->blockExists(i,j)==false</tt>] <tt>return.get()==NULL</tt>
   * </ul>
   */
  virtual Teuchos::RefCountPtr<const LinearOpWithSolveBase<RangeScalar,DomainScalar> >
  getBlockLOWS(const int i, const int j) const = 0; 

};

} // namespace Thyra

#endif // THYRA_BLOCKED_LINEAR_OP_WITH_SOLVE_BASE_HPP
