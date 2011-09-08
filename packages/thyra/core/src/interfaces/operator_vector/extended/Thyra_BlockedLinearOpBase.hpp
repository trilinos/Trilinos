// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_BLOCKED_LINEAR_OP_BASE_HPP
#define THYRA_BLOCKED_LINEAR_OP_BASE_HPP


#include "Thyra_LinearOpBase.hpp"


namespace Thyra {


/** \brief . */
template <class Scalar> class ProductVectorSpaceBase;

/** \brief Base interface for linear operators that can be accessed as
 * sub-blocks.
 *
 * ToDo: Finish Documentation.
 *
 * \ingroup Thyra_Op_Vec_extended_interfaces_code_grp
 */
template<class Scalar>
class BlockedLinearOpBase
  : virtual public LinearOpBase<Scalar>
{
public:

  /** \brief Return the product space for the range.
   *
   * A return value of <tt>return.get()==NULL</tt> is an indication that
   * <tt>*this</tt> is not fully initialized.
   */
  virtual Teuchos::RCP<const ProductVectorSpaceBase<Scalar> >
  productRange() const = 0;

  /** \brief Return the product space for the domain.
   *
   * A return value of <tt>return.get()==NULL</tt> is an indication that
   * <tt>*this</tt> is not fully initialized.
   */
  virtual Teuchos::RCP<const ProductVectorSpaceBase<Scalar> >
  productDomain() const = 0;

  /** \brief Return if the block <tt>(i,j)</tt> exists or not.
   *
   * \param  i  [in] Zero-based index for the block row.
   * \param  j  [in] Zero-based index for the block column.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>0 <= i && i < this->productRange()->numBlocks()</tt>
   * <li><tt>0 <= j && j < this->productDomain()->numBlocks()</tt>
   * </ul>
   */
  virtual bool blockExists(const int i, const int j) const = 0; 

  /** \brief Return if the block <tt>(i,j)</tt> is const only or not.
   *
   * \param  i  [in] Zero-based index for the block row.
   * \param  j  [in] Zero-based index for the block column.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>0 <= i && i < this->productRange()->numBlocks()</tt>
   * <li><tt>0 <= j && j < this->productDomain()->numBlocks()</tt>
   * </ul>
   */
  virtual bool blockIsConst(const int i, const int j) const = 0; 

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
  virtual Teuchos::RCP<LinearOpBase<Scalar> >
  getNonconstBlock(const int i, const int j) = 0; 

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
  virtual Teuchos::RCP<const LinearOpBase<Scalar> >
  getBlock(const int i, const int j) const = 0; 

};


} // namespace Thyra


#endif // THYRA_BLOCKED_LINEAR_OP_BASE_HPP
