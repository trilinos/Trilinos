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

#ifndef THYRA_DIAGONAL_LINEAR_OP_BASE_HPP
#define THYRA_DIAGONAL_LINEAR_OP_BASE_HPP

#include "Thyra_LinearOpBase.hpp"

namespace Thyra {

/** \brief Interface class for for diagonal linear operators.
 *
 * This interface represents a diagonal linear operator <tt>M</tt> of the form:
 \verbatim

 M = diag(diag)
 \endverbatim
 *
 * where <tt>diag</tt> is a <tt>VectorBase</tt> object.
 *
 * The operator subclass must implement <tt>apply()</tt> as follows:
 *
 \verbatim

 y = alpha*op(M)*x + beta*y
 
 =>

 y(i) = alpha*diag(i)*x(i) + beta*y(i), for i = 0 ... n-1
 \endverbatim
 *
 * where <tt>n = this->domain()->dim()</tt>.
 *
 * \ingroup Thyra_Op_Vec_extended_interfaces_code_grp
 */
template<class Scalar>
class DiagonalLinearOpBase : virtual public LinearOpBase<Scalar> {
public:

  /** @name Pure virtual functions that must be overridden in subclass */
  //@{

  /** \brief Returns true if the diagonal vector is const-only. */
  virtual bool isDiagConst() const = 0;

  /** \brief Returns the non-const diagonal vector <tt>diag</tt>.
   *
   * <b>Preconditions:</b><ul>
   * <li>[<tt>getDiag().get()!=NULL</tt>] <tt>isDiagConst()==false</tt>
   * </ul>
   *
   * Note that <tt>*this</tt> is not guaranteed to be fully modified until the
   * RCP returned is deleted.
   *
   * A return value of <tt>return.get()==NULL</tt> indicates that
   * <tt>this</tt> is not fully initialized.
   */
  virtual Teuchos::RCP<VectorBase<Scalar> > getNonconstDiag() = 0;

  /** \brief Returns the const diagonal vector <tt>diag</tt>.
   *
   * A return value of <tt>return.get()==NULL</tt> indicates that
   * <tt>this</tt> is not fully initialized.
   */
  virtual Teuchos::RCP<const VectorBase<Scalar> > getDiag() const = 0;

  //@}

};

}	// end namespace Thyra

#endif	// THYRA_DIAGONAL_LINEAR_OP_BASE_HPP
