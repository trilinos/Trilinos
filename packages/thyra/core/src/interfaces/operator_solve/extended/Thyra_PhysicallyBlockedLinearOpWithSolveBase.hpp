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

#ifndef THYRA_PHYSICALLY_BLOCKED_LINEAR_OP_WITH_SOLVE_BASE_HPP
#define THYRA_PHYSICALLY_BLOCKED_LINEAR_OP_WITH_SOLVE_BASE_HPP

#include "Thyra_BlockedLinearOpWithSolveBase.hpp"
#include "Thyra_PhysicallyBlockedLinearOpBase.hpp"


namespace Thyra {


/** \brief Base interface for linear operators with a solve that are composed
 * out of individual LOB and LOWSB objects.
 *
 * \ingroup Thyra_Op_Solve_extended_interfaces_code_grp
 *
 * ToDo: Finish Documentation.
 */
template<class Scalar>
class PhysicallyBlockedLinearOpWithSolveBase
  : virtual public BlockedLinearOpWithSolveBase<Scalar>,
    virtual public PhysicallyBlockedLinearOpBase<Scalar>
{
public:

  /** \brief Determines if the block <tt>(i,j)</tt> can be filled with a LOWDB
   * object or not.
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
  virtual bool acceptsLOWSBlock(const int i, const int j) const = 0;

  /** \brief . */
  virtual void setNonconstLOWSBlock(
    const int i, const int j,
    const Teuchos::RCP<LinearOpWithSolveBase<Scalar> > &block
    ) = 0;
  
  /** \brief . */
  virtual void setLOWSBlock(
    const int i, const int j,
    const Teuchos::RCP<const LinearOpWithSolveBase<Scalar> > &block
    ) = 0;

private:
  
  // Not defined and not to be called
  PhysicallyBlockedLinearOpWithSolveBase<Scalar>&
  operator=(const PhysicallyBlockedLinearOpWithSolveBase<Scalar>&);

};


} // namespace Thyra


#endif // THYRA_PHYSICALLY_BLOCKED_LINEAR_OP_WITH_SOLVE_BASE_HPP
