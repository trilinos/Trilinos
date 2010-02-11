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
