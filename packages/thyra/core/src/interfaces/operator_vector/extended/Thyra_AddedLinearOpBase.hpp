// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_ADDED_LINEAR_OP_BASE_HPP
#define THYRA_ADDED_LINEAR_OP_BASE_HPP

#include "Thyra_LinearOpBase.hpp"

namespace Thyra {

/** \brief Interface class for implicitly added linear operators.
 *
 * This interface represents a added linear operator <tt>M</tt> of the form:
 \verbatim
 
 M = Op[0] + Op[1] + ... + Op[numOps-1]
 \endverbatim
 *
 * where <tt>Op[]</tt> is an array of <tt>numOps</tt> <tt>LinearOpBase</tt>
 * objects.  Of course the operator <tt>M</tt> is not constructed explicitly
 * but instead just applies the constituent linear operators accordingly using
 * temporaries.
 *
 * In other words, subclasses define <tt>apply()</tt> as:
 *
 \verbatim

 y = alpha*M*x + beta*y
   = alpha * ( Op[0] * x + Op[1] * x + ...  + Op[numOps-1] * x ) + beta * y
 \endverbatim
 *
 * \ingroup Thyra_Op_Vec_extended_interfaces_code_grp
 */
template<class Scalar>
class AddedLinearOpBase : virtual public LinearOpBase<Scalar> {
public:

  /** @name Pure virtual functions that must be overridden by subclasses */
  //@{

  /** \brief Returns the number of constituent operators.
   *
   * A return value of <tt>0</tt> indicates that <tt>this</tt> is not fully
   * initialized.
   */
  virtual int numOps() const = 0;

  /** \brief Determine if the <tt>k</tt>th constituent operator is const-only or not.
   *
   * \param  k  [in] The zero-based index of the constituent operator to return.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt> 0 <= k < this->numOps()</tt>
   * </ul>
   */
  virtual bool opIsConst(const int k) const = 0;

  /** \brief Return the <tt>k</tt>th non-constant constituent operator.
   *
   * \param  k  [in] The zero-based index of the constituent operator to return.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt> 0 <= k < this->numOps()</tt>
   * <li><tt>this->opIsConst(k)==false</tt>
   * </ul>
   */
  virtual Teuchos::RCP<LinearOpBase<Scalar> > getNonconstOp(const int k) = 0;

  /** \brief Return the <tt>k</tt>th constant constituent operator.
   *
   * \param  k  [in] The zero-based index of the constituent operator to return.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt> 0 <= k < this->numOps()</tt>
   * </ul>
   */
  virtual Teuchos::RCP<const LinearOpBase<Scalar> > getOp(const int k) const = 0;

  //@}

};

} // namespace Thyra

#endif	// THYRA_ADDED_LINEAR_OP_BASE_HPP
